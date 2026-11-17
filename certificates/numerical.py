"""Assert numerical agreement; this does NOT certify infinite-sum error bounds."""
import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
import json
from pathlib import Path
import time
import nelson_cm_values as nelson

ROOT = Path(__file__).resolve().parents[1]

def exact_tau_check(E):
    """Check the quadratic fixed-point equation in an exact radical algebra."""
    R = nelson.Rad
    x,y,z,w,p,q = nelson.pq_coords(E.a, E.b, E.mu)
    a = R({1:x, p:y})
    b = R({1:z, p:w}) * R({q:1})
    c = R({1:-z, p:w}) * R({q:1})
    d = R({1:x, p:-y})
    t = E.tau_exact
    assert c*t*t + (d+a*(-1))*t + b*(-1) == R(), 'tau fixed-point equation'
    assert t == nelson.TABLE5_TAU[E.N], 'tau versus appendix'
    assert t.value()[1] > 0
    assert nelson.contains(E.R, (1,0,0,0))
    for v in E.R:
        assert (2*v[0]).denominator == 1 and E.A.nrd(v).denominator == 1
        for wv in E.R:
            assert nelson.contains(E.R, E.A.mul(v,wv))
    return dict(tau=str(t), fixed_point_equation_exact=True,
                order_closed_exact=True, reduced_order_discriminant=E.N,
                mu_discriminant=E.D, upper_half_plane=True)

def run_certificate(cutoffs=(144,196)):
    import mpmath as mp
    nelson.DIGITS = 90
    getcontext().prec = 90
    original_W = nelson.W_values
    nelson.W_values = lambda x,nus: original_W(x,nus,digits=80)
    PI = nelson.pi_dec()
    tolerance = Decimal('1e-40')
    report = dict(passed=False, numerical_only=True, interval_certified=False,
        scope='Observed discrepancy against conjectured exact values; no bound on infinite lattice tail or rounding error.',
        precision_digits=90, weight_target_digits=80, lattice_cutoffs=list(cutoffs),
        tolerance=str(tolerance), levels=[], forms=[])
    worst_c=worst_theta=worst_alpha=worst_stability=Decimal(0)
    start=time.time()
    for N in (14,15,21,26,33,34,38,46,58):
        E=nelson.CMPoint(N)
        level=exact_tau_check(E)
        level['N']=N
        P=nelson.pi_omega2(E.D,PI)
        with mp.workdps(100):
            p=-E.D
            logprod=sum((1 if pow(j,(p-1)//2,p)==1 else -1)*mp.loggamma(mp.mpf(j)/p) for j in range(1,p))
            gamma_P=mp.exp(logprod)/(2*p)
            gamma_error=abs(P-Decimal(mp.nstr(gamma_P,95)))
        assert gamma_error < Decimal('1e-75')
        level['period_eta_gamma_error']=str(gamma_error)
        report['levels'].append(level)
        for form in nelson.DATA[N]['forms']:
            row=dict(form=form,passes=[])
            previous=None
            for cutoff in cutoffs:
                c,theta=E.run(form,max(nelson.TABLE4[form]),Mcut=cutoff,PI=PI)
                comparisons=[]
                for k,expected in nelson.TABLE4[form].items():
                    ex=Decimal(expected.numerator)/Decimal(expected.denominator)
                    error=nelson.cabs((c[k][0]-ex,c[k][1]))
                    worst_c=max(worst_c,error)
                    assert error<tolerance, (form,k,error)
                    comparisons.append(dict(n=k,exact=str(expected),real=str(c[k][0]),imag=str(c[k][1]),observed_error=str(error)))
                alpha=nelson.TABLE5_THETA[form].value()
                theta_expected=nelson.cmul((P,Decimal(0)),alpha)
                theta_error=nelson.cabs((theta[0]-theta_expected[0],theta[1]-theta_expected[1]))
                alpha_error=nelson.cabs((theta[0]/P-alpha[0],theta[1]/P-alpha[1]))
                assert theta_error<tolerance and alpha_error<tolerance, (form,theta_error,alpha_error)
                worst_theta=max(worst_theta,theta_error)
                worst_alpha=max(worst_alpha,alpha_error)
                item=dict(cutoff=cutoff,coefficients=comparisons,theta_real=str(theta[0]),theta_imag=str(theta[1]),
                    theta_observed_error=str(theta_error),normalized_alpha_observed_error=str(alpha_error),
                    appendix_alpha_threshold_1e_50_observed=alpha_error<Decimal('1e-50'))
                if previous:
                    errors=[nelson.cabs((c[k][0]-previous[0][k][0],c[k][1]-previous[0][k][1])) for k in c]
                    errors.append(nelson.cabs((theta[0]-previous[1][0],theta[1]-previous[1][1])))
                    stability=max(errors)
                    assert stability<tolerance
                    item['successive_cutoff_difference']=str(stability)
                    worst_stability=max(worst_stability,stability)
                previous=(c,theta)
                row['passes'].append(item)
                print(f'{form}: cutoff={cutoff}; max c error={max(Decimal(x["observed_error"]) for x in comparisons):.3E}; alpha error={alpha_error:.3E}',flush=True)
            report['forms'].append(row)
            (ROOT/'results'/'numerical.partial.json').write_text(json.dumps(report,indent=2)+'\n')
    report.update(passed=True,max_coefficient_discrepancy=str(worst_c),max_theta_discrepancy=str(worst_theta),
        max_normalized_alpha_discrepancy=str(worst_alpha),max_cutoff_difference=str(worst_stability),seconds=round(time.time()-start,3))
    (ROOT/'results'/'numerical.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'results'/'numerical.partial.json').unlink(missing_ok=True)
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cutoffs',nargs='+',type=int,default=[144,196])
    args=parser.parse_args()
    result=run_certificate(tuple(args.cutoffs))
    print('PASS: numerical comparisons for all twelve forms. These are not interval-certified bounds.')
