"""Light checks only: no full series, periods, divisor images or algdep.

From the bundle root: load('sanity_check.sage'); sanity_check()
"""
from sage.all import *
from pathlib import Path
from sage.repl.preparse import preparse_file
import os
if 'REPRO_ROOT' not in globals():load('reproduce.sage')


def sanity_check():
    previous=Path.cwd();checks=[]
    try:
        # Parse every Sage module, but do not execute its computation driver.
        for source in REPRO_ROOT.rglob('*.sage'):
            compile(preparse_file(source.read_text(encoding='utf-8')),str(source),'exec')
        os.chdir(REPRO_ROOT/'cm_data');load('check_cm_inputs.sage')
        inputs=check_cm_inputs()
        assert sum(r['divisors_per_factor']*(2 if r['N'] in (26,38,58) else 1) for r in inputs)==5616
        os.chdir(REPRO_ROOT/'six_curves');load('run_six_curves.sage')
        for N in (14,15,21,33,34,46):
            data=six_load_inputs(N);setup=data['builder'](data['side_path'],20,128)
            assert len(setup['b'])==20 and setup['b'][0]==1 and setup['tau0'].imag()>0
            record=data['records'][0];CC=ComplexField(128);RF=RealField(128)
            M=data['matrix'](record['mu_coords'],RF)
            tau=_cm_fixed_point_upper_half_plane(M,CC)
            assert tau.imag()>0 and abs(_six_mob(M,tau)-tau)<RF(2)**(-60)
            checks.append(dict(E=str(N)+'a1',terms=20,prec_bits=128,passed=True))
        load('x33_affine.sage')
        six33_exact_al_audit(six_load_inputs(33))
        E=EllipticCurve('33a1');T=E(0,0)
        assert six33_translate_point(E,{},source_is_identity=True)==T
        assert six33_translate_point(E,{'x':T[0],'y':T[1]},target_is_identity=True).is_zero()
        EC=E.change_ring(QQbar);P=EC.lift_x(QQbar(1))
        assert six33_translate_point(EC,{'x':P[0],'y':P[1]})==P+EC(T)
        os.chdir(REPRO_ROOT/'genus_two');load('run_genus_two.sage')
        for label in ('26a1','26b1','38a1','38b1','58a1','58b1'):
            checks.append(genus_two_sanity_check(label,Nper=20,prec_bits=128))
        assert len(checks)==12
        print('PASS: all source modules parse; exact CM inputs cover5616divisors; 20-term setups pass for all12forms.')
        print('No full coefficient run, period matching, divisor image integration or algebraic recognition was run.')
        return dict(passed=True,forms=checks,CM_inputs=inputs,full_computation_run=False)
    finally:os.chdir(previous)
