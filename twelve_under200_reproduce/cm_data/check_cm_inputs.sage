"""Small exact input sanity check; no Magma, periods, images or algdep.

Run from cm_data: load('check_cm_inputs.sage'); check_cm_inputs()
"""
from sage.all import *
from sage.repl.load import load as _cm_input_load
from pathlib import Path
import re

def check_cm_inputs(directory='.'):
    root=Path(directory).resolve()
    report=[]
    for N in (14,15,21,26,33,34,38,46,58):
        ns=dict(globals())
        _cm_input_load(str(root/('x%s_cm_data.sage'%N)),ns)
        assert ns['CM_DB']==N and ns['CM_ABS_D_MIN']==3 and ns['CM_ABS_D_MAX']==199
        primes=ZZ(N).prime_divisors()
        assert ns['CM_RAMIFIED_PRIMES']==primes
        Ds=[-d for d in range(3,200) if fundamental_discriminant(-d)==-d
            and all(kronecker(-d,p)!=1 for p in primes)]
        assert ns['CM_SELECTED_DS']==Ds
        B=QuaternionAlgebra(QQ,ns['CM_ALGEBRA_I2'],ns['CM_ALGEBRA_J2'])
        assert B.discriminant()==N
        basis=matrix(QQ,ns['CM_ORDER_BASIS']); inverse=basis.inverse()
        def inside(q):return all(c in ZZ for c in vector(QQ,list(q))*inverse)
        base=B(ns['CM_BASE_COORDS']);unit=B(ns['CM_ORIENTATION_UNIT'])
        assert inside(base) and inside(unit) and unit.reduced_norm()==-1
        AL={k:B(v) for k,v in ns['CM_AL_COORDS'].items()}
        assert set(AL)=={'w'+str(p) for p in list(primes)+[N]}
        for name,w in AL.items():
            assert w.reduced_norm()==ZZ(name[1:]) and inside(w*w/w.reduced_norm())
            assert all(inside(w*B(list(b))*(~w)) and inside((~w)*B(list(b))*w) for b in basis.rows())
        classes=ns['cm_points'];allpoints=ns['cm_points_all_orientations']
        assert classes==[dict(zip(ns['CM_CLASS_FIELDS'],r)) for r in ns['CM_CLASS_RECORDS']]
        assert allpoints==[dict(zip(ns['CM_ALL_FIELDS'],r)) for r in ns['CM_ALL_RECORDS']]
        assert len(classes)==ns['CM_EXPECTED_CLASS_REPS']
        assert len(allpoints)==ns['CM_EXPECTED_ALL_ORIENTATIONS']
        compact=0
        for D in Ds:
            crows=[r for r in classes if r['D']==D]
            arows=[r for r in allpoints if r['D']==D]
            h=crows[0]['h'];chi=[kronecker(D,p) for p in primes]
            assert h>0 and len(crows)==h
            assert sorted(r['orbit_index'] for r in crows)==list(range(1,h+1))
            assert all(r['h']==h and r['chi']==chi and r['beta_norm']>0 for r in crows)
            translators=[(1,B(1))]
            for p,c in zip(primes,chi):
                if c==-1:translators += [(p*d,AL['w'+str(p)]*w) for d,w in translators]
            assert len(translators)==prod(1-c for c in chi)
            assert len(arows)==h*len(translators)
            assert sorted(r['orbit_index'] for r in arows)==list(range(1,len(arows)+1))
            wanted={(o,c) for o in range(1,len(translators)+1) for c in range(1,h+1)}
            assert {(r['orientation_index'],r['class_index']) for r in arows}==wanted
            byclass={r['orbit_index']:B(r['mu_coords']) for r in crows}
            for r in crows+arows:
                q=B(r['mu_coords'])
                assert inside(q) and q.reduced_trace()^2-4*q.reduced_norm()==D
                assert r['h']==h and r['chi']==chi
            for r in arows:
                d,w=translators[r['orientation_index']-1]
                assert r['orientation_translate']==d
                assert r['orbit_index']==(r['orientation_index']-1)*h+r['class_index']
                assert B(r['mu_coords'])==w*byclass[r['class_index']]*(~w)
            compact += 3*len(translators)
        magma=(root/('x%s_cm_under200.m'%N)).read_text()
        assert magma.count('for absD := 3 to 199 do')==1
        assert 'assert 3 le -D and -D le 199' in magma and 'CM_ABS_D_MIN = 3' in magma
        assert 'CLASS_NUMBER_BOUND' not in magma
        assert not re.search(r'\b(Read|Write|PrintFile|Open|load)\s*\(',magma)
        report.append(dict(N=N,discriminants=len(Ds),class_representatives=len(classes),
                           oriented_CM_points=len(allpoints),divisors_per_factor=3*len(allpoints),
                           raw_orientation_operator_groups_per_factor=int(compact),passed=True))
    assert sum(r['divisors_per_factor']*(2 if r['N'] in (26,38,58) else 1) for r in report)==5616
    print(report)
    return report
