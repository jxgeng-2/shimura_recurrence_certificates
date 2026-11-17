"""Fresh, source-only reproduction of all twelve |D|<200 CM tables.

Start Sage in this directory, then:
    load('reproduce.sage')
    run_all()
Nothing runs merely by loading this file. No previous results are required.
"""
from sage.all import *
from pathlib import Path
import os

REPRO_ROOT=Path(globals().get('CM_REPRO_ROOT',Path.cwd())).resolve()
if not (REPRO_ROOT/'common/heegner_image_from_series.sage').exists():
    raise ValueError('Start Sage in the extracted source directory, or set CM_REPRO_ROOT before loading.')
CM_LABELS=('14a1','15a1','21a1','26a1','26b1','33a1','34a1','38a1','38b1','46a1','58a1','58b1')
CM_GENUS_TWO={'26a1','26b1','38a1','38b1','58a1','58b1'}
CM_BOUND=200
CM_PASS_SPECS={
    '14a1':(1000,1200,1000),'15a1':(1000,1200,1000),
    '21a1':(2400,2600,1000),'26a1':(1000,1200,1000),
    '26b1':(3000,3200,1600),'33a1':(4000,4400,1000),
    '34a1':(2400,2800,1000),'38a1':(1600,1800,1000),
    '38b1':(3000,3200,1000),'46a1':(4000,4400,1000),
    '58a1':(4000,4400,1000),'58b1':(4000,4400,1000)}
CM_NATIVE_TARGETS={
    '14a1':('14a2',2),'15a1':('15a1',1),'21a1':('21a2',2),
    '26a1':('26a1',1),'26b1':('26b2',7),'33a1':('33a1',1),
    '34a1':('34a3',3),'38a1':('38a1',1),'38b1':('38b2',5),
    '46a1':('46a2',2),'58a1':('58a1',1),'58b1':('58b2',5)}


def cm_output_dir(output_dir=None):
    return (REPRO_ROOT/'results' if output_dir is None else Path(output_dir)).resolve()


def cm_result_path(label,output_dir=None):
    if label not in CM_LABELS:raise ValueError('Unknown target: '+str(label))
    return cm_output_dir(output_dir)/('x%s_cm_images_%s_result.sobj'%(label[:-2],label))


def load_result(label,output_dir=None):
    return load(str(cm_result_path(label,output_dir)))


def run_all(labels=CM_LABELS,pass_specs=None,output_dir=None,write_table=True):
    """Compute from source, or resume this source package's generated caches.

    All classes, orientations and three nonidentity AL operators are used.
    pass_specs maps a target label to (low_Nper,high_Nper,prec_bits).
    A subset produces subset exports; it is not labeled a complete table.
    """
    labels=(labels,) if isinstance(labels,str) else tuple(labels)
    if not labels or len(labels)!=len(set(labels)) or not set(labels)<=set(CM_LABELS):
        raise ValueError('Choose distinct labels from CM_LABELS.')
    specs=dict(CM_PASS_SPECS);specs.update(pass_specs or {})
    destination=cm_output_dir(output_dir);destination.mkdir(parents=True,exist_ok=True)
    results={};previous=Path.cwd()
    try:
        for label in labels:
            if label in CM_GENUS_TWO:
                os.chdir(REPRO_ROOT/'genus_two');load('run_genus_two.sage')
                result=run_genus_two(labels=(label,),out_dir=str(destination),
                                     pass_specs={label:specs[label]})
                results[label]=result[label]
            else:
                N=int(label[:-2]);os.chdir(REPRO_ROOT/'six_curves');load('run_six_curves.sage')
                result=run_six_curves(levels=(N,),output_dir=str(destination),
                                      pass_specs={N:specs[label]})
                results[label]=result[N]
        if write_table:
            os.chdir(REPRO_ROOT);load('export_table.sage')
            export_table(labels=labels,output_dir=destination)
    finally:os.chdir(previous)
    return results


def cm_divisor_point(label,D,operator,orientation_index=1,class_index=1,output_dir=None):
    """Point over Q(x,y), coordinate fields, inclusion into H_D, and branch."""
    result=load_result(label,output_dir)
    key=(ZZ(D),ZZ(orientation_index),operator)
    item=next(r for r in result['numeric_rows'] if
              (r['D'],r['orientation_index'],r['operator'])==key and r['class_index']==class_index)
    row=next(r for r in result['rows'] if (r['D'],r['orientation_index'],r['operator'])==key)
    if row['status']=='identity':return dict(status='identity',point=EllipticCurve(label)(0),field=QQ)
    if row['status']!='recognized':raise ValueError('Unresolved entry: '+str(row['status']))
    exact=next(r for r in result['exact_rows'] if (r['D'],r['orientation_index'],r['operator'])==key)
    match=next(m for m in exact['numeric_matches'] if m['orbit_index']==item['orbit_index'])
    point=exact['points'][int(match['point_index'])-1]
    return dict(point,status='recognized',ambient_class_field=exact['ambient_class_field'],
                complex_embedding=exact['complex_embedding'],numerical_match=match)


def compact_point(entry):
    """Optionally simplify Q(x,y)'s presentation, preserving its inclusion."""
    if entry['status']=='identity':return dict(entry)
    P=entry['point'];L=entry['point_field']['field']
    if L is QQ:
        K=QQ;Q=P;inclusion=entry['point_field']['inclusion']
    else:
        K,fromK,toK=L.optimized_representation()
        Q=EllipticCurve(K,[QQ(a) for a in P.curve().ainvs()])([toK(P[0]),toK(P[1])])
        inclusion=entry['point_field']['inclusion']*fromK
    PH=entry['ambient_point']
    assert inclusion(Q[0])==PH[0] and inclusion(Q[1])==PH[1]
    return dict(point=Q,field=K,inclusion=inclusion,ambient_point=PH,
                complex_embedding=entry['complex_embedding']*inclusion)
