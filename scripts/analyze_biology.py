#!/usr/bin/env python3
"""Descriptive analysis AFTER sealed final scoring. Never fit, tune, or reselect.

Usage: python scripts/analyze_biology.py --root PROJECT --lock FINAL_PREDICTION_LOCK
Requires reports/results/final_metrics.json. Does not produce test predictions.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from sklearn.metrics import average_precision_score

from benchmark_common import threshold_metrics
from evaluate_locked import validate_prediction_rows
from selection_lock import MODEL_NAMES, safe_path, validate_selection_lock


def require(condition,message):
    if not condition: raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_csv(path):
    with Path(path).open(newline='') as handle: return list(csv.DictReader(handle))


def load_json(path):
    return json.loads(Path(path).read_text())


def similarity_band(value):
    value=float(value)
    require(math.isfinite(value) and 0 <= value <= 1,'Invalid pinned similarity')
    return '<0.5' if value < .5 else '0.5–<0.7' if value < .7 else '>=0.7'


def subgroup_metrics(y,score,threshold):
    y=np.asarray(y,dtype=int);score=np.asarray(score,dtype=float)
    require(y.shape==score.shape and y.ndim==1,'Invalid subgroup vectors')
    require(set(y.tolist()) <= {0,1} and np.isfinite(score).all(),'Invalid subgroup data')
    require(math.isfinite(threshold),'Invalid fixed threshold')
    n=len(y);positive=int(y.sum())
    counts=threshold_metrics(y,score,threshold) if n else {
        'threshold':float(threshold),'n':0,'tp':0,'fp':0,'tn':0,'fn':0,
        'sensitivity':None,'specificity':None,'balanced_accuracy':None}
    return {**counts,'positive':positive,'negative':n-positive,
            'average_precision':float(average_precision_score(y,score)) if 0 < positive < n else None,
            'average_precision_undefined_reason':None if 0 < positive < n else 'Both label classes required'}


def source_bins(evidence,label):
    values=set()
    records=json.loads(evidence)
    require(isinstance(records,list) and bool(records),'Missing source evidence')
    for record in records:
        require(isinstance(record.get('values'),dict) and bool(record['values']),'Missing original activity bins')
        values.update(record['values'].values())
    require(bool(values) and values <= {'+','++','+++'},'Unknown original activity bin')
    require(not ('+++' in values and len(values)>1),'Conflicting threshold bins in retained parent')
    require(int(label)==int(values=={'+++'}),'Source bins disagree with frozen label')
    return '|'.join(sorted(values))


def molecular_descriptors(smiles):
    mol=Chem.MolFromSmiles(smiles)
    require(mol is not None,'Invalid frozen structure')
    return {'mw':float(Descriptors.MolWt(mol)), 'logp':float(Crippen.MolLogP(mol)),
            'hbd':int(Lipinski.NumHDonors(mol)), 'hba':int(Lipinski.NumHAcceptors(mol)),
            'tpsa':float(rdMolDescriptors.CalcTPSA(mol))}


def collision_audit(rows):
    groups=defaultdict(list)
    for row in rows:
        mol=Chem.MolFromSmiles(row['canonical_smiles'])
        require(mol is not None,'Invalid collision input')
        groups[Chem.MolToSmiles(mol,canonical=True,isomericSmiles=False)].append(row)
    flags={};summary={'rows':len(rows),'distinct_inputs':len(groups),'collision_groups':0,
                     'collision_members':0,'mixed_label_groups':0,'mixed_label_members':0,
                     'empirical_minimum_errors':0}
    for members in groups.values():
        n=len(members);positive=sum(int(r['binary_y']) for r in members)
        mixed=0 < positive < n
        summary['collision_groups']+=int(n>1)
        summary['collision_members']+=n if n>1 else 0
        summary['mixed_label_groups']+=int(mixed)
        summary['mixed_label_members']+=n if mixed else 0
        summary['empirical_minimum_errors']+=min(positive,n-positive)
        for row in members:
            flags[row['compound_id']]={'nonisomeric_collision':n>1,
                'nonisomeric_group_size':n,'nonisomeric_mixed_label_collision':mixed}
    summary['empirical_minimum_error_fraction']=summary['empirical_minimum_errors']/len(rows) if rows else None
    return flags,summary


def training_neighbors(rows,test_ids,train_ids):
    result={}
    for row in rows:
        if (row['reference_split'],row['query_split']) != ('train','test'): continue
        cid=row['query_id']
        require(cid in test_ids and cid not in result,'Unexpected or duplicated train/test neighbor')
        require(row['nearest_id'] in train_ids,'Neighbor is not a training parent')
        value=float(row['max_tanimoto'])
        result[cid]={'nearest_train_id':row['nearest_id'],'nearest_train_tanimoto':value,
                     'similarity_band':similarity_band(value)}
    require(set(result)==test_ids,'Incomplete pinned training-neighbor coverage')
    return result


def seed_variation(encoder):
    require(encoder['validation_parents']==70,'Unexpected validation scope')
    require(len(encoder['trials'])==6,'Incomplete seed/LR budget')
    groups=[]
    for lr in (1e-5,3e-5):
        trials=sorted([t for t in encoder['trials'] if t['learning_rate']==lr],key=lambda t:t['seed'])
        require([t['seed'] for t in trials]==[17,29,43],'Missing or duplicated optimization seed')
        values=[float(t['best_validation_average_precision']) for t in trials]
        require(all(math.isfinite(x) and 0<=x<=1 for x in values),'Invalid validation AP')
        groups.append({'learning_rate':lr,'selected_learning_rate':lr==encoder['partial_fine_tuning']['selected_learning_rate'],
            'seeds':[{'seed':t['seed'],'selected_epoch':t['best_epoch'],'validation_average_precision':v}
                     for t,v in zip(trials,values)],
            'mean':float(np.mean(values)),'sample_standard_deviation':float(np.std(values,ddof=1)),
            'minimum':min(values),'maximum':max(values)})
    return {'metric':'Best-epoch validation average precision per optimization seed',
            'n':70,'positive':39,'negative':31,'learning_rates':groups,
            'interpretation':'Optimization variability on reused validation data; not test performance or a confidence interval. Distinct from scaffold bootstrap uncertainty.'}


def test_seed_variation(receipt,test_rows,ensemble_scores,selected_lr):
    ids=receipt['compound_ids'];expected=[r['compound_id'] for r in test_rows]
    require(len(ids)==len(set(ids)) and set(ids)==set(expected),'Seed receipt parent coverage differs')
    probabilities=receipt['seed_probabilities']
    require(set(probabilities)=={'17','29','43'},'Incomplete selected test seeds')
    indices={cid:i for i,cid in enumerate(ids)}
    y=np.array([int(r['binary_y']) for r in test_rows]);vectors=[];results=[]
    require(0 < y.sum() < len(y),'Both test classes required')
    for seed in (17,29,43):
        original=np.asarray(probabilities[str(seed)],dtype=float)
        require(original.shape==(len(ids),) and np.isfinite(original).all() and
                ((0<=original)&(original<=1)).all(),'Invalid per-seed probability vector')
        ordered=original[[indices[cid] for cid in expected]];vectors.append(ordered)
        results.append({'seed':seed,'test_average_precision':float(average_precision_score(y,ordered))})
    require(np.allclose(np.mean(vectors,axis=0),ensemble_scores,rtol=0,atol=1e-12),
            'Per-seed mean differs from sealed adapted-encoder probabilities')
    values=[r['test_average_precision'] for r in results]
    return {'metric':'Final-test average precision per selected optimization seed','learning_rate':selected_lr,
            'n':len(y),'positive':int(y.sum()),'negative':int(len(y)-y.sum()),'seeds':results,
            'mean_seed_ap':float(np.mean(values)),'sample_standard_deviation':float(np.std(values,ddof=1)),
            'minimum':min(values),'maximum':max(values),
            'ensemble_average_precision':float(average_precision_score(y,ensemble_scores)),
            'interpretation':'Descriptive variation across all three preselected seeds, without seed selection. Mean seed AP differs from AP of mean probabilities. This is not scaffold-bootstrap uncertainty.'}


def validate_final_bindings(final,lock,lock_sha):
    require(final.get('status')=='SCORED_LOCKED_PREDICTIONS','Final scoring must already be complete')
    require(final.get('selection_lock_sha256')==lock_sha,'Final result selection-lock hash differs')
    for key in ('cohort_sha256','protocol_sha256','predictions_sha256'):
        require(final.get(key)==lock.get(key) and bool(lock.get(key)),'Final result binding differs: '+key)
    require((final.get('n'),final.get('positive'),final.get('negative'),final.get('source_families'))==(246,69,177,1),
            'Unexpected frozen final-result scope')
    require(set(final.get('metrics',{}))==set(lock.get('models',{}))==MODEL_NAMES,'Incomplete final model set')
    for name in MODEL_NAMES:
        threshold=float(lock['models'][name]['validation_threshold'])
        require(math.isfinite(threshold),'Invalid locked threshold')
        metric=final['metrics'][name]
        require(metric['validation_selected_threshold']==threshold and metric['classification']['threshold']==threshold,
                'Final result threshold differs: '+name)


def describe(values):
    values=np.asarray(values,dtype=float)
    return {'n':len(values),'median':float(np.median(values)) if len(values) else None,
            'q25':float(np.quantile(values,.25)) if len(values) else None,
            'q75':float(np.quantile(values,.75)) if len(values) else None}


def analyze(root,lock_path,metrics_path):
    root=Path(root).resolve();lock_path=Path(lock_path);metrics_path=Path(metrics_path)
    # All gates precede cohort label parsing and any additional test metric calculation.
    lock_hash=sha(lock_path);metrics_hash=sha(metrics_path)
    lock=load_json(lock_path);final=load_json(metrics_path)
    validate_final_bindings(final,lock,lock_hash)
    validate_selection_lock(root,lock)
    encoder_path='data/models/encoder/encoder_development.json'
    preflight_path='reports/models/encoder_preflight.json'
    encoder=load_json(root/encoder_path);preflight=load_json(root/preflight_path)
    seed_receipt_path='data/predictions/encoder/receipt.json'
    require(seed_receipt_path in lock.get('prediction_receipts',{}),'Missing locked encoder prediction receipt')
    require(sha(root/seed_receipt_path)==lock['prediction_receipts'][seed_receipt_path],
            'Changed encoder prediction receipt')
    seed_receipt=load_json(root/seed_receipt_path)
    require(seed_receipt.get('status')=='LOCKED_TEST_PREDICTIONS_NO_METRICS' and seed_receipt.get('kind')=='encoder',
            'Incomplete encoder prediction receipt')
    for key in ('protocol_sha256','cohort_sha256','model_selection_lock_sha256'):
        require(seed_receipt.get(key)==lock.get(key),'Encoder prediction receipt binding differs: '+key)
    require(seed_receipt.get('models')=={name:spec for name,spec in lock['models'].items()
            if name in ('frozen_encoder','adapted_encoder')},'Encoder prediction receipt model choices differ')
    gate_path='reports/curation/cohort_gate.json';gate=load_json(root/gate_path)
    require(sha(root/gate_path)==encoder['bindings']['cohort_gate_sha256'],'Changed cohort gate')
    require(preflight['status']=='PASS' and preflight['input_gate']['status']=='PASS','Incomplete encoder preflight')
    require(rdBase.rdkitVersion==preflight['input_gate']['rdkit_version'],'RDKit differs from frozen input gate')
    require(gate['fingerprint']=={'radius':2,'bits':2048,'include_chirality':True},'Changed structural audit policy')
    for relative in ('data/processed/cohort.csv','data/processed/nearest_neighbors.csv'):
        require(sha(root/relative)==gate['outputs'][relative],'Changed pinned data: '+relative)
    require(gate['outputs']['data/processed/cohort.csv']==lock['cohort_sha256'],'Cohort gate/lock differs')
    paths={'selection_lock':lock_path,'final_metrics':metrics_path,
           'predictions':safe_path(root,lock['predictions_path']),
           'cohort':root/'data/processed/cohort.csv','neighbors':root/'data/processed/nearest_neighbors.csv',
           'cohort_gate':root/gate_path,'encoder_development':root/encoder_path,
           'encoder_preflight':root/preflight_path,'analysis_code':Path(__file__)}
    paths['encoder_prediction_receipt']=root/seed_receipt_path
    before={key:sha(path) for key,path in paths.items()}
    require(before['selection_lock']==lock_hash and before['final_metrics']==metrics_hash,
            'Final scoring inputs changed during gate validation')
    require(before['predictions']==lock['predictions_sha256'],'Prediction bytes changed during gate validation')
    expected_hashes={'encoder_development':lock['development_artifacts'][encoder_path],
        'encoder_preflight':lock['development_artifacts'][preflight_path],
        'encoder_prediction_receipt':lock['prediction_receipts'][seed_receipt_path],
        'cohort_gate':encoder['bindings']['cohort_gate_sha256'],
        'cohort':lock['cohort_sha256'],'neighbors':gate['outputs']['data/processed/nearest_neighbors.csv']}
    require(all(before[key]==digest for key,digest in expected_hashes.items()),
            'Pinned artifacts changed during gate validation')
    all_rows=load_csv(paths['cohort']);test=[r for r in all_rows if r['split']=='test']
    require(len(test)==246 and sum(int(r['binary_y']) for r in test)==69,'Wrong frozen test cohort')
    neighbors=training_neighbors(load_csv(paths['neighbors']),{r['compound_id'] for r in test},
                                 {r['compound_id'] for r in all_rows if r['split']=='train'})
    flags,collisions=collision_audit(test)
    for key,value in collisions.items():
        require(preflight['input_gate']['by_split']['test'][key]==value,'Collision audit differs from preflight: '+key)
    scores=validate_prediction_rows(load_csv(paths['predictions']),test,sorted(MODEL_NAMES))
    y=np.array([int(r['binary_y']) for r in test]);parents=[]
    for row in test:
        cid=row['compound_id']
        parents.append({'compound_id':cid,'canonical_smiles':row['canonical_smiles'],
            'scaffold_group':row['scaffold_group'],'source_family':row['source_family'],'binary_y':int(row['binary_y']),
            'original_activity_bins':source_bins(row['source_evidence'],row['binary_y']),
            **neighbors[cid],**flags[cid],**molecular_descriptors(row['canonical_smiles'])})
    model_rows=[];overall={}
    for name in sorted(MODEL_NAMES):
        threshold=float(lock['models'][name]['validation_threshold'])
        current=subgroup_metrics(y,scores[name],threshold);overall[name]=current
        recorded=final['metrics'][name]
        require(all(current[k]==recorded['classification'][k] for k in ('n','tp','fp','tn','fn','threshold')),
                'Final result confusion counts disagree with sealed bytes')
        require(math.isclose(current['average_precision'],recorded['ranking']['average_precision'],abs_tol=1e-12,rel_tol=0),
                'Final result AP disagrees with sealed bytes')
        for parent,score in zip(parents,scores[name]):
            predicted=int(score>=threshold);label=parent['binary_y']
            error=('TP' if predicted else 'FN') if label else ('FP' if predicted else 'TN')
            model_rows.append({**parent,'model':name,'score':float(score),'fixed_validation_threshold':threshold,
                               'predicted_y':predicted,'error_type':error})
    subgroups={}
    for field,categories in [('similarity_band',['<0.5','0.5–<0.7','>=0.7']),
                              ('original_activity_bins',sorted({p['original_activity_bins'] for p in parents})),
                              ('nonisomeric_collision',[False,True]),
                              ('nonisomeric_mixed_label_collision',[False,True])]:
        subgroups[field]=[]
        for category in categories:
            mask=np.array([p[field]==category for p in parents])
            subgroups[field].append({'category':category,'n':int(mask.sum()),'positive':int(y[mask].sum()),
                'negative':int(mask.sum()-y[mask].sum()),'models':{name:subgroup_metrics(y[mask],score[mask],
                    float(lock['models'][name]['validation_threshold'])) for name,score in scores.items()}})
    descriptors={}
    for field in ('mw','logp','hbd','hba','tpsa'):
        descriptors[field]={'by_frozen_label':{str(label):describe([p[field] for p in parents if p['binary_y']==label]) for label in (0,1)},
            'by_model_error':{name:{error:describe([r[field] for r in model_rows if r['model']==name and r['error_type']==error])
                                    for error in ('TP','FN','FP','TN')} for name in sorted(MODEL_NAMES)}}
    report={'status':'DESCRIPTIVE_POST_SCORE_ANALYSIS','n_test_parents':246,'positive':69,'negative':177,
        'n_model_parent_rows':len(model_rows),'models':sorted(MODEL_NAMES),'overall':overall,'subgroups':subgroups,
        'nonisomeric_collision_audit':collisions,'descriptors':descriptors,
        'descriptor_units':{'mw':'Da','logp':'RDKit Crippen estimate','hbd':'count','hba':'count','tpsa':'angstrom^2'},
        'source_activity_bin_intervals_nm':{'+++':'[1,100]','++':'(100,1000]','+':'(1000,10000]'},
        'source_bin_policy':'Preserve every original bin; joint weak bins remain joint, without averaging or within-bin potency ranking.',
        'similarity_policy':{'reference_split':'train','fingerprint':gate['fingerprint'],
            'interpretation':'Pinned structural audit uses chirality; primary model inputs are nonisomeric/achiral. Bands do not use recomputed alternative similarities.'},
        'validation_optimization_seed_variation':seed_variation(encoder),
        'test_optimization_seed_variation':test_seed_variation(seed_receipt,test,scores['adapted_encoder'],
            lock['models']['adapted_encoder']['learning_rate']),'input_sha256':before,
        'runtime_versions':{'rdkit':rdBase.rdkitVersion,'numpy':np.__version__},
        'limitations':['Descriptive post-score summaries only; no tuning, significance tests or causal/mechanistic attribution.',
            'One held-out patent family; source, chemistry, time and assay protocol remain confounded.',
            'Small or single-class strata do not support ranking comparisons; AP is null unless both classes occur.',
            'Nonisomeric collisions describe representation ambiguity; they do not establish stereochemical biology or exclude assay noise.',
            'Descriptor summaries repeat the same parents across models; model-parent rows are not independent compounds.']}
    require(before=={key:sha(path) for key,path in paths.items()},'Inputs changed during analysis')
    report['input_hashes_unchanged_after_analysis']=True
    json.dumps(report,allow_nan=False)
    return report,model_rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--lock',type=Path,required=True)
    parser.add_argument('--final-metrics',type=Path)
    parser.add_argument('--output-dir',type=Path,help='Separate empty directory for analysis replay')
    args=parser.parse_args()
    report_path=(args.output_dir/'biological_analysis.json') if args.output_dir else args.root/'reports/results/biological_analysis.json'
    rows_path=(args.output_dir/'rows.csv') if args.output_dir else args.root/'data/analysis/rows.csv'
    require(not report_path.exists() and not rows_path.exists(),'Analysis output exists; preserve it and use --output-dir for replay')
    report,rows=analyze(args.root,args.lock,args.final_metrics or args.root/'reports/results/final_metrics.json')
    buffer=io.StringIO(newline='');writer=csv.DictWriter(buffer,fieldnames=list(rows[0]))
    writer.writeheader();writer.writerows(rows)
    row_bytes=buffer.getvalue().encode()
    report['rows_csv_sha256']=hashlib.sha256(row_bytes).hexdigest()
    for path,payload in ((rows_path,row_bytes),(report_path,(json.dumps(report,indent=2,allow_nan=False)+'\n').encode())):
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('xb') as handle: handle.write(payload)
    print(json.dumps({'status':report['status'],'test_parents':246,'model_parent_rows':len(rows)}))


if __name__=='__main__': main()
