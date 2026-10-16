"""Derive locked model choices from complete, hash-verified development artifacts."""
import json
import csv
import argparse
from datetime import datetime, timezone
from pathlib import Path
from benchmark_common import COHORT_SHA256, PROTOCOL_SHA256, sha256

BASELINE_FAMILIES = {'constant_prevalence','five_neighbor','fingerprint_logistic',
                     'fingerprint_logistic_chiral','random_forest'}
MODEL_NAMES = BASELINE_FAMILIES | {'frozen_encoder','adapted_encoder'}
CODE_PATHS = ('scripts/benchmark_common.py','scripts/prepare_encoder.py',
              'scripts/train_encoder.py','scripts/selection_lock.py',
              'scripts/predict_locked.py','scripts/evaluate_locked.py')


def require(condition,message):
    if not condition:
        raise ValueError(message)


def load_json(path):
    return json.loads(Path(path).read_text())


def safe_path(root, relative):
    root=Path(root).resolve();path=(root/relative).resolve()
    require(path.is_relative_to(root) and path != root, 'Artifact must be inside the project')
    return path


def verify_file_map(root, file_map):
    require(bool(file_map), 'Empty artifact manifest')
    for relative,digest in file_map.items():
        require(sha256(safe_path(root,relative))==digest,'Artifact hash mismatch: '+relative)


def load_expected_selection(root):
    root=Path(root).resolve()
    require(sha256(root/'docs/PROTOCOL.md')==PROTOCOL_SHA256,'Changed frozen protocol')
    require(sha256(root/'data/processed/cohort.csv')==COHORT_SHA256,'Changed frozen cohort')
    base_root=root/'data/models/baselines'
    seal=load_json(base_root/'SEALED.json')
    require(seal.get('status')=='SEALED_VALIDATION_ONLY' and seal.get('test_scored') is False,
            'Baselines must be sealed before test scoring')
    require(seal['protocol_sha256']==PROTOCOL_SHA256 and seal['cohort_sha256']==COHORT_SHA256,
            'Baseline freeze mismatch')
    verify_file_map(base_root,seal['files'])
    base=load_json(base_root/'validation_baselines.json')
    require(set(base['selected'])==BASELINE_FAMILIES and base['candidate_count']==10,
            'Missing prescribed baseline choices')
    require(base['train_n']==604 and base['validation_n']==70 and base['test_scored'] is False,
            'Wrong baseline development scope')
    require(base['protocol_sha256']==PROTOCOL_SHA256 and base['cohort_sha256']==COHORT_SHA256,
            'Baseline report bindings differ')
    require(base['code_sha256']['benchmark_common.py']==sha256(root/'scripts/benchmark_common.py'),
            'Shared scientific implementation changed after baseline fitting')
    primary={k:v for k,v in base['selected'].items() if k!='fingerprint_logistic_chiral'}
    comparator=sorted(primary,key=lambda k:(-primary[k]['validation_average_precision'],k))[0]
    require(comparator==base['strongest_primary_classical_baseline'],'Invalid classical comparator')
    models={name:{'kind':'baseline','candidate_id':record['candidate_id'],
                  'validation_threshold':record['threshold'],'chirality':record['chirality'],
                  'model_relative_path':'data/models/baselines/'+record['model_file']}
            for name,record in base['selected'].items()}
    encoder_root=root/'data/models/encoder'
    enc=load_json(encoder_root/'encoder_development.json')
    require(enc.get('status')=='DEVELOPMENT_COMPLETE_NO_FINAL_TEST_SCORES','Incomplete encoder development')
    require(enc['bindings']['protocol_sha256']==PROTOCOL_SHA256 and enc['bindings']['cohort_sha256']==COHORT_SHA256,
            'Encoder freeze mismatch')
    require(enc['training_parents']==604 and enc['validation_parents']==70,'Wrong encoder development scope')
    require(enc['code_sha256']==sha256(root/'scripts/train_encoder.py'),'Training implementation changed after fitting')
    verify_file_map(encoder_root,enc['artifacts'])
    selection=enc['partial_fine_tuning']
    require(selection['selected_seeds']==[17,29,43] and len(enc['trials'])==6,
            'Incomplete predeclared seed/LR budget')
    expected_means={}
    for lr in (1e-5,3e-5):
        group=[t for t in enc['trials'] if t['learning_rate']==lr]
        require(sorted(t['seed'] for t in group)==[17,29,43],'Missing or duplicated seed')
        expected_means[lr]=sum(t['best_validation_average_precision'] for t in group)/3
    selected_lr=max(expected_means,key=lambda lr:(expected_means[lr],-lr))
    require(selected_lr==selection['selected_learning_rate'],'LR does not follow all-seed validation selection')
    head='data/models/encoder/frozen_head.npz'
    asset_sha=enc['model_receipt']['asset_manifest_sha256']
    models['frozen_encoder']={'kind':'frozen_encoder',
        'validation_threshold':enc['frozen_head']['selected_threshold'],
        'selected_C':enc['frozen_head']['selected_C'],'head_path':head,'model_asset_sha':asset_sha}
    trials=sorted((t for t in enc['trials'] if t['learning_rate']==selected_lr),key=lambda t:t['seed'])
    models['adapted_encoder']={'kind':'adapted_encoder',
        'validation_threshold':enc['selected_threshold'],'head_path':head,'model_asset_sha':asset_sha,
        'learning_rate':selected_lr,'seeds':[17,29,43],
        'checkpoint_paths':['data/models/encoder/'+t['artifact_directory']+'/best.pt' for t in trials],
        'prediction_rule':'mean_probability_all_three_selected_seeds'}
    development_paths=['data/models/baselines/SEALED.json','data/models/baselines/validation_baselines.json',
                       'data/models/encoder/encoder_development.json','reports/models/encoder_preflight.json']
    artifacts={path:sha256(root/path) for path in development_paths}
    for path,digest in seal['files'].items():
        artifacts['data/models/baselines/'+path]=digest
    for path,digest in enc['artifacts'].items():
        artifacts['data/models/encoder/'+path]=digest
    return {'models':models,'classical_comparator':comparator,'development_artifacts':artifacts}


def validate_lock_definition(lock, expected, current_code_hashes, expected_status="LOCKED_BEFORE_FINAL_TEST_SCORING"):
    require(lock.get('status')==expected_status,'Missing pre-test selection lock')
    require(lock.get('protocol_sha256')==PROTOCOL_SHA256 and lock.get('cohort_sha256')==COHORT_SHA256,
            'Lock does not use frozen protocol/cohort constants')
    require(set(lock.get('models',{}))==MODEL_NAMES,'Missing or unprescribed model in lock')
    require(lock['models']==expected['models'],'Lock choices or thresholds differ from development evidence')
    require(lock.get('classical_comparator')==expected['classical_comparator'],
            'Lock comparator differs from validation selection')
    require(bool(lock.get('development_artifacts')) and lock['development_artifacts']==expected['development_artifacts'],
            'Lock development evidence is empty or incomplete')
    require(set(lock.get('code_sha256',{}))==set(CODE_PATHS) and lock['code_sha256']==current_code_hashes,
            'Prediction/scoring code changed or is not bound')


def verified_arm_predictions(root, model_lock):
    """Bind scored bytes to both completed prediction arms, without label access."""
    root = Path(root)
    rows, receipts = [], {}
    lock_hash = sha256(root/'manifests/model_selection_lock.json')
    for arm, names in (('baselines', BASELINE_FAMILIES),
                       ('encoder', {'frozen_encoder','adapted_encoder'})):
        relative = 'data/predictions/'+arm+'/receipt.json'
        receipt = load_json(root/relative)
        require(receipt.get('status') == 'LOCKED_TEST_PREDICTIONS_NO_METRICS' and
                receipt.get('kind') == arm, 'Incomplete prediction arm')
        require(receipt.get('target_columns_accessed') is False and
                receipt.get('fitting_performed') is False and
                receipt.get('metrics_computed') is False, 'Prediction arm exceeded scope')
        require(receipt.get('model_selection_lock_sha256') == lock_hash,
                'Prediction arm used a different selection lock')
        for key in ('protocol_sha256','cohort_sha256','code_sha256','development_artifacts'):
            require(receipt.get(key) == model_lock[key], 'Prediction receipt binding differs: '+key)
        require(receipt.get('models') == {name:model_lock['models'][name] for name in names},
                'Prediction models differ from locked selection')
        ids = receipt.get('compound_ids', [])
        require(len(ids) == len(set(ids)) == 246 and receipt.get('test_input_count') == 246,
                'Wrong prediction input coverage')
        prediction_file = root/'data/predictions'/arm/'predictions.csv'
        require(sha256(prediction_file) == receipt.get('predictions_sha256'), 'Prediction bytes differ from receipt')
        with prediction_file.open(newline='') as handle:
            reader = csv.DictReader(handle)
            require(reader.fieldnames == ['model','compound_id','split','score'], 'Unexpected prediction schema')
            arm_rows = list(reader)
        require(len(arm_rows) == receipt.get('prediction_row_count') == len(names)*246,
                'Wrong arm prediction count')
        require({(r['model'], r['compound_id']) for r in arm_rows} ==
                {(name,cid) for name in names for cid in ids} and
                all(r['split'] == 'test' for r in arm_rows), 'Incomplete arm prediction coverage')
        rows.extend(arm_rows)
        receipts[relative] = sha256(root/relative)
    return sorted(rows, key=lambda r:(r['model'],r['compound_id'])), receipts


def require_merged_predictions(path, expected_rows):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ['model','compound_id','split','score'], 'Unexpected merged prediction schema')
        rows = list(reader)
    require(sorted(rows,key=lambda r:(r['model'],r['compound_id'])) == expected_rows,
            'Merged predictions differ from completed arm outputs')


def create_model_selection_lock(root):
    root=Path(root).resolve();expected=load_expected_selection(root)
    result={'status':'MODEL_SELECTION_LOCKED','locked_at_utc':datetime.now(timezone.utc).isoformat(),'protocol_sha256':PROTOCOL_SHA256,
            'cohort_sha256':COHORT_SHA256,**expected,
            'code_sha256':{name:sha256(root/name) for name in CODE_PATHS}}
    validate_lock_definition(result,expected,result['code_sha256'],expected_status='MODEL_SELECTION_LOCKED')
    return result


def load_locked_model_selection(root):
    root=Path(root).resolve()
    lock=load_json(root/'manifests/model_selection_lock.json')
    expected=load_expected_selection(root)
    current={name:sha256(root/name) for name in CODE_PATHS}
    validate_lock_definition(lock,expected,current,expected_status='MODEL_SELECTION_LOCKED')
    return lock


def create_selection_lock(root,predictions_path):
    root=Path(root).resolve();predictions_path=Path(predictions_path).resolve()
    relative=str(predictions_path.relative_to(root))
    model_lock=load_locked_model_selection(root)
    expected={key:model_lock[key] for key in ('models','classical_comparator','development_artifacts')}
    rows, receipts = verified_arm_predictions(root, model_lock)
    require_merged_predictions(predictions_path, rows)
    result={'status':'LOCKED_BEFORE_FINAL_TEST_SCORING','locked_at_utc':datetime.now(timezone.utc).isoformat(),'prediction_receipts':receipts,'protocol_sha256':PROTOCOL_SHA256,
            'cohort_sha256':COHORT_SHA256,**expected,
            'code_sha256':{name:sha256(root/name) for name in CODE_PATHS},
            'predictions_path':relative,'predictions_sha256':sha256(predictions_path),
            'model_selection_lock_sha256':sha256(root/'manifests/model_selection_lock.json')}
    validate_lock_definition(result,expected,result['code_sha256'])
    return result


def validate_selection_lock(root,lock):
    root=Path(root).resolve();expected=load_expected_selection(root)
    model_lock = load_locked_model_selection(root)
    rows, receipts = verified_arm_predictions(root, model_lock)
    require(lock.get('prediction_receipts') == receipts, 'Prediction receipts changed')
    require_merged_predictions(safe_path(root,lock['predictions_path']), rows)
    require(lock.get('model_selection_lock_sha256')==sha256(root/'manifests/model_selection_lock.json'),
            'Pre-prediction model choices changed')
    current={path:sha256(root/path) for path in CODE_PATHS}
    validate_lock_definition(lock,expected,current)
    verify_file_map(root,lock['development_artifacts'])
    require(sha256(safe_path(root,lock['predictions_path']))==lock['predictions_sha256'],
            'Changed sealed prediction bytes')
    return expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--stage', choices=('models','predictions'), required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.stage == 'models':
        path = root/'manifests/model_selection_lock.json'
        require(not path.exists(), 'Model selection lock already exists')
        result = create_model_selection_lock(root)
    else:
        path = root/'manifests/final_selection_lock.json'
        merged = root/'data/predictions/final_predictions.csv'
        require(not path.exists() and not merged.exists(), 'Final prediction lock/output already exists')
        lock = load_locked_model_selection(root)
        rows, _ = verified_arm_predictions(root,lock)
        with merged.open('x', newline='') as handle:
            writer = csv.DictWriter(handle,fieldnames=['model','compound_id','split','score'])
            writer.writeheader();writer.writerows(rows)
        result = create_selection_lock(root,merged)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as handle:
        handle.write(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'models':sorted(result['models'])}))


if __name__ == '__main__':
    main()
