"""Period normalization and the level33 affine correction."""
# six_curves/period_normalization.sage
"""Match all side periods, with an exact lattice-index calculation."""
from sage.all import ComplexField, RealField, EllipticCurve, ZZ, gcd


def six_lattice_coordinates(z, e1, e2):
    area = e1.real()*e2.imag() - e2.real()*e1.imag()
    return ((z.real()*e2.imag()-e2.real()*z.imag())/area,
            (e1.real()*z.imag()-z.real()*e1.imag())/area)


def six_match_periods(periods, E_label="33a1", prec_bits=1000,
                      search_bound=12, min_match_bits=60):
    """Find a scalar taking every side period into the target period lattice.

    The index is the gcd of ALL two-by-two minors of the integral lattice
    coordinates, not the determinant from an arbitrarily selected pair.
    Smallest index is selected before comparing numerical residuals.
    """
    CC, RR = ComplexField(prec_bits), RealField(prec_bits)
    periods = [CC(p) for p in periods]
    E = EllipticCurve(E_label)
    e1, e2 = map(CC, E.period_lattice().basis(prec=prec_bits))
    tol = RR(2)**(-min_match_bits)
    size = max([abs(p) for p in periods] + [RR(1)])
    nonzero = []
    for p in periods:
        if abs(p) < size*tol:
            continue
        if any(min(abs(p-q), abs(p+q)) < size*tol for q in nonzero):
            continue
        nonzero.append(p)
    nonzero.sort(key=abs)
    candidates, seen_mu = [], []
    for i, p1 in enumerate(nonzero):
        for p2 in nonzero[i+1:]:
            if abs((p2/p1).imag()) < RR("1e-10"):
                continue
            ratio = complex(p2/p1)
            c1, c2 = complex(e1), complex(e2)
            area = c1.real*c2.imag - c2.real*c1.imag
            for m in range(-int(search_bound), int(search_bound)+1):
                for n in range(-int(search_bound), int(search_bound)+1):
                    if m == 0 and n == 0:
                        continue
                    approx = ratio*(m*c1+n*c2)
                    a = (approx.real*c2.imag-c2.real*approx.imag)/area
                    b = (c1.real*approx.imag-approx.real*c1.imag)/area
                    if max(abs(a-round(a)), abs(b-round(b))) > 1e-7:
                        continue
                    mu = (m*e1+n*e2)/p1
                    if mu.real() < 0 or (mu.real() == 0 and mu.imag() < 0):
                        mu = -mu
                    if any(abs(mu-v) < tol*max(1, abs(mu)) for v in seen_mu):
                        continue
                    seen_mu.append(mu)
                    coords = [six_lattice_coordinates(mu*p, e1, e2)
                              for p in periods]
                    integers = [(ZZ(a.round()), ZZ(b.round())) for a,b in coords]
                    error = max([max(abs(a-u), abs(b-v))
                                 for (a,b),(u,v) in zip(coords,integers)] + [RR(0)])
                    if error > tol:
                        continue
                    index = gcd([abs(u*s-v*r) for u,v in integers for r,s in integers])
                    if index == 0:
                        continue
                    candidates.append(dict(mu=mu, index=index,
                        lattice_coordinates=integers, lattice_error=error,
                        lattice_error_bits=(-error.log(2)).floor() if error else prec_bits,
                        E_label=E_label, periods=periods, target_periods=(e1,e2)))
            primitive = [v for v in candidates if v["index"] == 1]
            if primitive:
                return min(primitive, key=lambda v: v["lattice_error"])
    if not candidates:
        raise ValueError("No all-period lattice match at %s bits for %s" %
                         (min_match_bits, E_label))
    return min(candidates, key=lambda v: (v["index"], v["lattice_error"]))


# genus_two/periods.sage
"""Balanced side periods and all-period lattice matching for the three levels."""
from sage.all import ComplexField, RealField, EllipticCurve, ZZ, gcd


def g2_lattice_coordinates(z, e1, e2):
    area = e1.real()*e2.imag() - e2.real()*e1.imag()
    return ((z.real()*e2.imag()-e2.real()*z.imag())/area,
            (e1.real()*z.imag()-z.real()*e1.imag())/area)


def g2_match_periods(periods, E_label, prec_bits=1000,
                      search_bound=12, min_match_bits=60):
    """Find a scalar taking every side period into the target period lattice.

    The index is the gcd of ALL two-by-two minors of the integral lattice
    coordinates, not the determinant from an arbitrarily selected pair.
    Smallest index is selected before comparing numerical residuals.
    """
    CC, RR = ComplexField(prec_bits), RealField(prec_bits)
    periods = [CC(p) for p in periods]
    E = EllipticCurve(E_label)
    e1, e2 = map(CC, E.period_lattice().basis(prec=prec_bits))
    tol = RR(2)**(-min_match_bits)
    size = max([abs(p) for p in periods] + [RR(1)])
    nonzero = []
    for p in periods:
        if abs(p) < size*tol:
            continue
        if any(min(abs(p-q), abs(p+q)) < size*tol for q in nonzero):
            continue
        nonzero.append(p)
    nonzero.sort(key=abs)
    candidates, seen_mu = [], []
    for i, p1 in enumerate(nonzero):
        for p2 in nonzero[i+1:]:
            if abs((p2/p1).imag()) < RR("1e-10"):
                continue
            ratio = complex(p2/p1)
            c1, c2 = complex(e1), complex(e2)
            area = c1.real*c2.imag - c2.real*c1.imag
            for m in range(-int(search_bound), int(search_bound)+1):
                for n in range(-int(search_bound), int(search_bound)+1):
                    if m == 0 and n == 0:
                        continue
                    approx = ratio*(m*c1+n*c2)
                    a = (approx.real*c2.imag-c2.real*approx.imag)/area
                    b = (c1.real*approx.imag-approx.real*c1.imag)/area
                    if max(abs(a-round(a)), abs(b-round(b))) > 1e-7:
                        continue
                    mu = (m*e1+n*e2)/p1
                    if mu.real() < 0 or (mu.real() == 0 and mu.imag() < 0):
                        mu = -mu
                    if any(abs(mu-v) < tol*max(1, abs(mu)) for v in seen_mu):
                        continue
                    seen_mu.append(mu)
                    coords = [g2_lattice_coordinates(mu*p, e1, e2)
                              for p in periods]
                    integers = [(ZZ(a.round()), ZZ(b.round())) for a,b in coords]
                    error = max([max(abs(a-u), abs(b-v))
                                 for (a,b),(u,v) in zip(coords,integers)] + [RR(0)])
                    if error > tol:
                        continue
                    index = gcd([abs(u*s-v*r) for u,v in integers for r,s in integers])
                    if index == 0:
                        continue
                    candidates.append(dict(mu=mu, index=index,
                        lattice_coordinates=integers, lattice_error=error,
                        lattice_error_bits=(-error.log(2)).floor() if error else prec_bits,
                        E_label=E_label, periods=periods, target_periods=(e1,e2)))
            primitive = [v for v in candidates if v["index"] == 1]
            if primitive:
                return min(primitive, key=lambda v: v["lattice_error"])
    if not candidates:
        raise ValueError("No all-period lattice match at %s bits for %s" %
                         (min_match_bits, E_label))
    return min(candidates, key=lambda v: (v["index"], v["lattice_error"]))

def g2_side_periods(setup,audit):
    """Evaluate the cocycle at each side's balanced hyperbolic midpoint."""
    assert min(r['residual_bits'] for r in audit)>=30
    tau,coeffs=setup['tau0'],setup['b'];bar=tau.conjugate()
    H=antiderivative_evaluator(coeffs,tau,setup['Nper'],setup['prec_bits'])
    def mob(g,z):return (g[0,0]*z+g[0,1])/(g[1,0]*z+g[1,1])
    def disk(z):return (z-tau)/(z-bar)
    rows=[]
    for i,g in enumerate(setup['S_matrices']):
        endpoint=disk(mob(g.inverse(),tau));w=endpoint/(1+(1-abs(endpoint)**2).sqrt())
        z=(tau-bar*w)/(1-w);gw=disk(mob(g,z))
        rows.append(dict(side_index=i,value=H(gw)-H(w),stats=dict(mode='balanced_period_cocycle',
            max_u=max(abs(w),abs(gw)),tail_abs_estimate=H.tail_estimate(w)+H.tail_estimate(gw),
            tail_is_rigorous=False,balanced_modularity_bits=audit[i]['residual_bits'])))
    return rows






# six_curves/x33_affine.sage
"""All-discriminant w11 preconditioning and P_w33=P_w3+(0,0).

The exact quaternion relations and rational2-torsion are checked afresh.
The +1 differential eigenvalue has numerical slash-operator evidence.
Both truncations retain their original coordinates and comparison evidence.
The origin and source-T exceptions are handled using both truncations.
"""
from sage.all import *
from pathlib import Path

def six33_exact_al_audit(inputs):
    B=QuaternionAlgebra(QQ,*inputs['algebra']);basis=matrix(QQ,inputs['order_basis']);inverse=basis.inverse()
    def coords(g):return vector(QQ,g)*inverse
    AL={k:B(v) for k,v in inputs['AL_coords'].items()};v=AL['w11'];square=v*v/v.reduced_norm()
    assert square.reduced_norm()==1 and all(c in ZZ for c in coords(square))
    for row in basis:assert all(c in ZZ for c in coords(v*B(list(row))*(~v)))
    commutators={}
    for op,w in AL.items():
        g=w*v*(~(v*w));assert g.reduced_norm()==1 and all(c in ZZ for c in coords(g))
        commutators[op]=dict(quaternion_coords=list(g),order_coords=list(coords(g)),norm=1)
    product=AL['w3']*AL['w11']*(~AL['w33'])
    assert product.reduced_norm()==1 and all(c in ZZ for c in coords(product))
    E=EllipticCurve('33a1');T=E(0,0);assert not T.is_zero() and (2*T).is_zero()
    return dict(N=33,E='33a1',kernel='w11',quaternion_invariants=inputs['algebra'],order_basis=inputs['order_basis'],
        kernel_normalizes_order=True,square_over_norm_coords=list(square),square_over_norm_order_coords=list(coords(square)),
        commutators=commutators,product_quaternion_coords=list(product),product_order_coords=list(coords(product)),
        product_norm=1,torsion_point=[0,0],exact_torsion_order=2,relation='P_w33=P_w3+(0,0)')

def _six33_eigenchecks(inputs,setup,reducer,bits):
    CC,RF=ComplexField(bits),RealField(bits);tau=setup['tau0']
    def value(z):
        w=(z-tau)/(z-tau.conjugate());v=CC(0)
        for bn in reversed(setup['b']):v=v*w+bn
        return (1-w)^2*v
    records=[]
    for op,coords in inputs['AL_coords'].items():
        W=inputs['matrix'](coords,RF);expected=1 if op=='w11' else -1;checks=[]
        for w0 in (CC(0),CC(RF(1)/10,RF(7)/100)):
            z=(tau-tau.conjugate()*w0)/(1-w0);wz=_six_mob(W,z);gamma,wr=reducer.reduce_upper(wz)
            zr=(tau-tau.conjugate()*wr)/(1-wr)
            eigen=W.det()*value(zr)/((W[1,0]*z+W[1,1])^2*(gamma[1,0]*wz+gamma[1,1])^2*value(z))
            err=abs(eigen-expected);accuracy=bits if err==0 else int(floor(-err.log2()))
            assert accuracy>=60, 'Insufficient AL eigenvalue accuracy.'
            checks.append(dict(eigenvalue=eigen,agreement_bits=accuracy,reduced_radius=abs(wr)))
        records.append(dict(E='33a1',operator=op,expected_eigenvalue=expected,checks=checks))
    return records

def six33_translate_point(curve,source,source_is_identity=False,target_is_identity=False):
    """Small reusable addition helper; identity flags require pairwise evidence."""
    T=curve(0,0)
    if source_is_identity:
        assert not target_is_identity
        return T
    if target_is_identity:return curve(0)
    return curve.point([source['x'],source['y'],curve.base_ring()(1)],check=False)+T

def six33_prepare_pair(low_nper=4000,high_nper=4400,bits=1000,output_dir=None):
    inputs=six_load_inputs(33);output=Path(output_dir).resolve() if output_dir is not None else SIX_SOURCE_DIR/'results'
    CC,RF=ComplexField(bits),RealField(bits);E=EllipticCurve('33a1');EC=E.change_ring(CC)
    reducer=DirichletReducer(inputs['side_path'],prec_bits=bits)
    AL={k:inputs['matrix'](v,RF) for k,v in inputs['AL_coords'].items()}
    paths={n:output/('x33_numeric_%s_%s.sobj'%(n,bits)) for n in (low_nper,high_nper)}
    audit=six33_exact_al_audit(inputs)
    audit['eigenchecks']=_six33_eigenchecks(inputs,load(str(output/('x33_setup_%s_%s.sobj'%(high_nper,bits)))),reducer,bits)
    audit['proof_boundary']='Exact quaternion and torsion arithmetic; numerical differential eigenvalue evidence, not interval certification.'
    e1,e2=map(CC,E.period_lattice().basis(prec=bits));halves=[e1/2,e2/2,(e1+e2)/2]
    zT=max(halves,key=lambda z:min(_six_agree(E.period_lattice().elliptic_exponential(z)[i],CC(0)) for i in (0,1)))
    assert min(_six_agree(E.period_lattice().elliptic_exponential(zT)[i],CC(0)) for i in (0,1))>=bits-30
    def radius(z):return abs(reducer.reduce_upper(z)[1])
    def original(r):return r.get('affine_AL_translation',{}).get('original',r.get('kernel_preconditioning',{}).get('original',r))
    def snapshot(r):return {k:r[k] for k in ('x','y','zE','zE_reduced','stats','identity_bits')}
    for nper,other_nper in ((low_nper,high_nper),(high_nper,low_nper)):
        payload=load(str(paths[nper]));other=load(str(paths[other_nper]));other_by={_six_key(r):r for r in other['rows']}
        setup=load(str(output/('x33_setup_%s_%s.sobj'%(nper,bits))))
        assert payload['geometry_signature']==inputs['signature']
        payload.update(AL_coords=inputs['AL_coords'],kernel_preconditioning_audit=audit,affine_AL_translation_audit=audit)
        for r in payload['rows']:
            if r['operator']!='w3' or r.get('kernel_preconditioning',{}).get('version')==1:continue
            z=r['tauD'];vz=_six_mob(AL['w11'],z);W=AL['w3']
            r0=max(radius(z),radius(_six_mob(W,z)));r1=max(radius(vz),radius(_six_mob(W,vz)))
            info=dict(version=1,kernel='w11',original_tau=z,evaluation_tau=z,original_radius=r0,alternate_radius=r1,selected='identity')
            if r1+RF('0.001')<r0 and r['identity_bits']<70:
                previous=snapshot(r);pair=original(other_by[_six_key(r)])
                known=min(_six_agree(previous[k],pair[k]) for k in ('x','y'))
                ans=heegner_image_from_series('w3',vz,b=setup['b'],tau0=setup['tau0'],AL=AL,
                    mu=payload['period_match']['mu'],E_label='33a1',Nper=nper,prec_bits=bits,reducer=reducer,
                    verbose=False,identity_min_bits=None,max_radius=.985)
                matches=[_six_agree(previous[k],ans['P'][j]) for j,k in enumerate(('x','y'))]
                assert min(matches)>=max(60,min(known,bits-40)-20)
                info.update(selected='w11',evaluation_tau=vz,original=previous,original_two_truncation_bits=known,original_xy_agreement_bits=matches)
                r.update(x=ans['P'][0],y=ans['P'][1],zE=ans['zE'],zE_reduced=ans['zE_reduced'],stats=ans['stats'],identity_bits=ans['identity_bits'])
            r['kernel_preconditioning']=info;_six_save(payload,paths[nper])
        own={_six_key(r):r for r in payload['rows']}
        for r in payload['rows']:
            if r['operator']!='w33' or r.get('affine_AL_translation',{}).get('version')==2:continue
            source_key=(int(r['D']),int(r['orbit_index']),'w3');source=own[source_key];source_other=other_by[source_key]
            previous=snapshot(r);pair=original(other_by[_six_key(r)])
            source_zero=min(source['identity_bits'],source_other['identity_bits'])>=70
            target_zero=min(previous['identity_bits'],pair['identity_bits'])>=70
            source_T_bits=min(_six_agree(s[k],CC(0)) for s in (source,source_other) for k in ('x','y'))
            if target_zero:assert not source_zero and source_T_bits>=70, 'Origin requires both source-T and original-period evidence.'
            P=six33_translate_point(EC,source,source_zero,target_zero)
            info=dict(version=2,source_operator='w3',torsion_point=[0,0],original=previous,
                translation_formula='P_w33=P_w3+(0,0)',source_is_identity=source_zero,target_is_identity=target_zero)
            stats=dict(source['stats']);stats['mode']='affine_AL_translation_from_w3'
            if target_zero:
                known=min(previous['identity_bits'],pair['identity_bits'],source_T_bits)
                info.update(identity_case='target_identity',original_two_truncation_bits=known,original_xy_agreement_bits=None,source_T_match_bits=source_T_bits)
                r.update(x=CC(0),y=CC(1),zE=CC(0),zE_reduced=CC(0),identity_bits=known,is_infinity=True,
                    identity_method='two_truncations_and_exact_affine_relation')
            else:
                assert not P.is_zero(), 'Unexpected identity without two-truncation confirmation.'
                known=min(_six_agree(previous[k],pair[k]) for k in ('x','y'))
                matches=[_six_agree(previous[k],P[j]) for j,k in enumerate(('x','y'))]
                assert min(matches)>=max(60,min(known,bits-40)-20)
                z=zT if source_zero else source['zE']+zT;c1,c2=six_lattice_coordinates(z,e1,e2)
                err=max(abs(c1-c1.round()),abs(c2-c2.round()));zr=z-ZZ(c1.round())*e1-ZZ(c2.round())*e2
                r.update(x=P[0],y=P[1],zE=z,zE_reduced=zr,identity_bits=bits if err==0 else max(0,int(floor(-err.log2()))),is_infinity=False)
                info.update(identity_case='source_identity' if source_zero else 'finite',original_two_truncation_bits=known,original_xy_agreement_bits=matches)
                if source_zero:r['accuracy_cap_bits']=min(known,min(matches),source['identity_bits'],source_other['identity_bits'])-15
            r['stats']=stats;r['affine_AL_translation']=info;_six_save(payload,paths[nper])
        _six_save(payload,paths[nper])
    return [load(str(paths[n])) for n in (low_nper,high_nper)]
