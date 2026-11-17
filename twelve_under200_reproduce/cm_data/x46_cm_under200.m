/* All admissible negative fundamental D with 3 <= |D| <= 199.
   No class-number cutoff. Frozen algebra, maximal order, base point and AL lifts.
   Paste into https://magma.maths.usyd.edu.au/calc/; print-only, no file I/O.
   Fresh class representatives may differ from the bundled exact Sage input. */
SetColumns(0); SetSeed(46);
DB := 46; RAMIFIED_PRIMES := [ 2, 23 ];
A := QuaternionAlgebra<Rationals()|46,-3>;
order_basis_coords := [
    [ 1, 0, 0, 0 ],
    [ 1/2, 0, 1/6, 1/3 ],
    [ 0, 0, 0, 1 ],
    [ 0, 1/2, 0, 1/2 ]
];
O := QuaternionOrder([A!c:c in order_basis_coords]);
assert Discriminant(O) eq DB and IsMaximal(O);
mu_base := A![ 1/2, 0, -35/6, -2/3 ];
orientation_unit := A![ 6, 1, 7, 1 ];
assert mu_base in O and orientation_unit in O and Norm(orientation_unit) eq -1;
W2 := A![ 1, 0, 7/3, -1/3 ]; assert Norm(W2) eq 2;
W23 := A![ -345, 89/2, 230/3, -109/6 ]; assert Norm(W23) eq 23;
W46 := A![ -46, -6, -46, -7 ]; assert Norm(W46) eq 46;
for W in [W2,W23,W46] do
 assert &and[W*(A!o)*W^(-1) in O and W^(-1)*(A!o)*W in O:o in Basis(O)];
 assert W^2/Norm(W) in O;
end for;

function PositiveGenerator(I,U)
 ok,w := IsPrincipal(I); assert ok;
 if Norm(w) lt 0 then w *:= U; end if;
 assert Norm(w) gt 0; return w;
end function;
selected := [* *];
for absD := 3 to 199 do
 D := -absD;
 if not IsFundamentalDiscriminant(D) then continue; end if;
 chis := [KroneckerSymbol(D,p):p in RAMIFIED_PRIMES];
 if 1 in chis then continue; end if;
 Append(~selected,<D,ClassNumber(D),chis>);
end for;
print "CM_DB = 46";
print "CM_RAMIFIED_PRIMES = [ 2, 23 ]";
print "CM_ABS_D_MIN = 3";
print "CM_ABS_D_MAX = 199";
print "CM_ALGEBRA_I2 = 46";
print "CM_ALGEBRA_J2 = -3";
print "CM_ORDER_BASIS = [";
print "    [ 1, 0, 0, 0 ],";
print "    [ 1/2, 0, 1/6, 1/3 ],";
print "    [ 0, 0, 0, 1 ],";
print "    [ 0, 1/2, 0, 1/2 ],";
print "]";
print "CM_BASE_COORDS = [ 1/2, 0, -35/6, -2/3 ]";
print "CM_ORIENTATION_UNIT = [ 6, 1, 7, 1 ]";
print "CM_AL_COORDS = {'w2': [ 1, 0, 7/3, -1/3 ], 'w23': [ -345, 89/2, 230/3, -109/6 ], 'w46': [ -46, -6, -46, -7 ]}";
printf "CM_SELECTED_DS = %o\n",[s[1]:s in selected];
printf "CM_EXPECTED_CLASS_REPS = %o\n",&+([0] cat [s[2]:s in selected]);
printf "CM_EXPECTED_ALL_ORIENTATIONS = %o\n",&+([0] cat [s[2]*&*[1-c:c in s[3]]:s in selected]);
records := [* *]; allrecords := [* *];
for entry in selected do
 D,h,chis := Explode(entry); assert 3 le -D and -D le 199;
 K<s> := QuadraticField(D); OK := MaximalOrder(K);
 assert Discriminant(OK) eq D and ClassNumber(OK) eq h;
 assert HasEmbedding(K,A:ComputeEmbedding:=false);
 mu0,phi0 := Embed(OK,O); f0 := MinimalPolynomial(A!mu0);
 class_map := ClassGroupPrimeRepresentatives(OK,ideal<OK|DB>);
 C := Domain(class_map); assert #C eq h;
 translators := [<1,A!1>];
 if chis[1] eq -1 then translators cat:= [<2*t[1],W2*t[2]>:t in translators]; end if;
 if chis[2] eq -1 then translators cat:= [<23*t[1],W23*t[2]>:t in translators]; end if;
 idx := 0;
 for c in C do
  idx +:= 1; aa := class_map(c); assert GCD(Integers()!Norm(aa),DB) eq 1;
  Iaa := rideal<O|[A!phi0(OK!x):x in Basis(aa)]>;
  beta := A!PositiveGenerator(Iaa,orientation_unit);
  muc := beta^(-1)*(A!mu0)*beta;
  assert muc in O and MinimalPolynomial(muc) eq f0;
  Append(~records,<D,h,idx,chis,Eltseq(muc),Norm(beta)>);
  for ti:=1 to #translators do
   t:=translators[ti]; mut:=t[2]*muc*t[2]^(-1);
   assert mut in O and MinimalPolynomial(mut) eq f0;
   Append(~allrecords,<D,h,(ti-1)*h+idx,chis,Eltseq(mut),t[1],ti,idx>);
  end for;
 end for;
 assert idx eq h;
end for;
assert #records eq &+([0] cat [s[2]:s in selected]);
assert #allrecords eq &+([0] cat [s[2]*&*[1-c:c in s[3]]:s in selected]);
print "CM_CLASS_FIELDS = ('D','h','orbit_index','chi','mu_coords','beta_norm')";
print "CM_CLASS_RECORDS = [";
for r in records do printf " [%o,%o,%o,%o,%o,%o],\n",r[1],r[2],r[3],r[4],r[5],r[6];end for;
print "]";
print "cm_points = [dict(zip(CM_CLASS_FIELDS,r)) for r in CM_CLASS_RECORDS]";
print "CM_ALL_FIELDS = ('D','h','orbit_index','chi','mu_coords','orientation_translate','orientation_index','class_index')";
print "CM_ALL_RECORDS = [";
for r in allrecords do printf " [%o,%o,%o,%o,%o,%o,%o,%o],\n",r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8];end for;
print "]";
print "cm_points_all_orientations = [dict(zip(CM_ALL_FIELDS,r)) for r in CM_ALL_RECORDS]";
printf "# All assertions passed; full range 3..199; %o class representatives; %o CM points.\n",#records,#allrecords;

