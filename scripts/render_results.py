#!/usr/bin/env python3
"""Render the completed fixed benchmark; no fitting or model selection."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

NAMES = [
    ('constant_prevalence', 'Constant reference'),
    ('five_neighbor', 'Fingerprint 5-neighbor'),
    ('fingerprint_logistic', 'Fingerprint logistic'),
    ('random_forest', 'Random forest'),
    ('frozen_encoder', 'Frozen MoLFormer + head'),
    ('adapted_encoder', 'Partially adapted MoLFormer'),
    ('fingerprint_logistic_chiral', 'Chiral logistic (sensitivity)'),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    source = args.root/'reports/results/final_metrics.json'
    report = json.loads(source.read_text())
    if report['status'] != 'SCORED_LOCKED_PREDICTIONS' or report['n'] != 246:
        raise ValueError('Completed frozen final-test result required')
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False,
                         'svg.hashsalt':'prmt5-sourcebench-v1'})
    fig, axes = plt.subplots(1,2,figsize=(13.5,5.2),gridspec_kw={'width_ratios':[1.15,1]})
    y = np.arange(len(NAMES))
    colors = ['#8b949e','#7a9bb6','#7a9bb6','#245c85','#9467bd','#be5a72','#a6b5c2']
    for i, (key, label) in enumerate(NAMES):
        entry = report['metrics'][key]
        ap = entry['ranking']['average_precision']
        lo,hi = report['uncertainty']['models'][key]['average_precision_percentile_95']
        axes[0].plot([lo,hi],[i,i],color=colors[i],lw=3,alpha=.6)
        axes[0].scatter([ap],[i],color=colors[i],s=55,zorder=3)
        axes[0].text(1.015,i,f'{ap:.3f}',va='center',fontsize=9)
        screen = entry['ranking']['screening']
        expected = screen['expected_hits']
        axes[1].barh(i,expected,color=colors[i],height=.56)
        low, high = screen['minimum_hits'],screen['maximum_hits']
        if low != high:
            axes[1].plot([low,high],[i+.23,i+.23],color='#333333',lw=1.3)
            label_text=f'{expected:.1f} expected ({low}–{high})'
        else:
            label_text=f'{int(expected)}/25'
        axes[1].text(expected+.35,i,label_text,va='center',fontsize=9)
    prevalence = report['positive']/report['n']
    axes[0].axvline(prevalence,color='#5c6369',ls='--',lw=1)
    axes[0].set(yticks=y,yticklabels=[label for _,label in NAMES],xlim=(0,1),
                xlabel='Average precision (AP)',title='Ranking across 246 held-out compounds')
    axes[0].invert_yaxis()
    axes[1].set(yticks=y,yticklabels=[],xlim=(0,29),xlabel='Source-defined positives among top 25',
                title='Screening budget: 25 of 246 compounds')
    axes[1].invert_yaxis()
    axes[1].axvline(25*prevalence,color='#5c6369',ls='--',lw=1)
    fig.suptitle('PRMT5 SourceBench | 69 positives / 246 compounds | one held-out patent family',
                 x=.02,ha='left',fontweight='bold',fontsize=13)
    fig.text(.02,.075,'AP bars: 95% scaffold-bootstrap intervals (2,000 paired resamples). Dashed lines: source prevalence / random-selection expectation.',fontsize=9)
    fig.text(.02,.045,'Right-hand whiskers: attainable hits from boundary-tie ordering, not a confidence interval.',fontsize=9)
    fig.text(.02,.015,'Source-reported biochemical IC50 ≤100 nM; assay, chemistry and publication date change together. Not prospective or cellular validation.',fontsize=9)
    fig.tight_layout(rect=(0,.13,1,.94),w_pad=3)
    directory=args.root/'reports/figures';directory.mkdir(parents=True,exist_ok=True)
    outputs={}
    for extension in ('svg','png'):
        path=directory/('benchmark.'+extension)
        metadata={'Date':None} if extension=='svg' else {'Software':'PRMT5 SourceBench'}
        fig.savefig(path,dpi=180,metadata=metadata,bbox_inches='tight')
        outputs[str(path.relative_to(args.root))]=hashlib.sha256(path.read_bytes()).hexdigest()
    plt.close(fig)
    receipt={'input':'reports/results/final_metrics.json',
             'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
             'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'matplotlib':matplotlib.__version__,'outputs':outputs}
    (directory/'figure_manifest.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))


if __name__=='__main__':
    main()
