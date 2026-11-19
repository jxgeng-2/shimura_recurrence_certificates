/* Optional regeneration of exact CM inputs, not numerical Phi values.
   Load this file in Magma, then call GenerateCMInputs(N).
   It prints Sage assignments for one level. Use the corresponding block
   in CM Quaternion and Embedding Data.sage when comparing/replacing class representatives.
   This optional generator was not rerun during the sanity-only check. */

procedure GenerateCMInputs14()
  /* All admissible negative fundamental D with 3 <= |D| <= 199.
     No class-number cutoff. Frozen algebra, maximal order, base point and AL lifts.
     Paste into https://magma.maths.usyd.edu.au/calc/; print-only, no file I/O.
     Fresh class representatives may differ from the bundled exact Sage input. */
  SetColumns(0); SetSeed(14);
  DB := 14; RAMIFIED_PRIMES := [ 2, 7 ];
  A := QuaternionAlgebra<Rationals()|7,-2>;
  order_basis_coords := [
      [ 1, 0, 0, 0 ],
      [ 1/2, 1/2, 0, 1/2 ],
      [ 1/2, -1/2, 0, 1/2 ],
      [ 0, 0, 1/2, 1/2 ]
  ];
  O := QuaternionOrder([A!c:c in order_basis_coords]);
  assert Discriminant(O) eq DB and IsMaximal(O);
  mu_base := A![ 1/2, -1/2, -3/2, 0 ];
  orientation_unit := A![ 2, 1, 1, 0 ];
  assert mu_base in O and orientation_unit in O and Norm(orientation_unit) eq -1;
  W2 := A![ 0, 0, 1, 0 ]; assert Norm(W2) eq 2;
  W7 := A![ 7, 2, 0, 1 ]; assert Norm(W7) eq 7;
  W14 := A![ 0, 2, 7, -2 ]; assert Norm(W14) eq 14;
  for W in [W2,W7,W14] do
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
  print "CM_DB = 14";
  print "CM_RAMIFIED_PRIMES = [ 2, 7 ]";
  print "CM_ABS_D_MIN = 3";
  print "CM_ABS_D_MAX = 199";
  print "CM_ALGEBRA_I2 = 7";
  print "CM_ALGEBRA_J2 = -2";
  print "CM_ORDER_BASIS = [";
  print "    [ 1, 0, 0, 0 ],";
  print "    [ 1/2, 1/2, 0, 1/2 ],";
  print "    [ 1/2, -1/2, 0, 1/2 ],";
  print "    [ 0, 0, 1/2, 1/2 ],";
  print "]";
  print "CM_BASE_COORDS = [ 1/2, -1/2, -3/2, 0 ]";
  print "CM_ORIENTATION_UNIT = [ 2, 1, 1, 0 ]";
  print "CM_AL_COORDS = {'w2': [ 0, 0, 1, 0 ], 'w7': [ 7, 2, 0, 1 ], 'w14': [ 0, 2, 7, -2 ]}";
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
   if chis[2] eq -1 then translators cat:= [<7*t[1],W7*t[2]>:t in translators]; end if;
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


end procedure;

procedure GenerateCMInputs15()
  /* All admissible negative fundamental D with 3 <= |D| <= 199.
     No class-number cutoff. Frozen algebra, maximal order, base point and AL lifts.
     Paste into https://magma.maths.usyd.edu.au/calc/; print-only, no file I/O.
     Fresh class representatives may differ from the bundled exact Sage input. */
  SetColumns(0); SetSeed(15);
  DB := 15; RAMIFIED_PRIMES := [ 3, 5 ];
  A := QuaternionAlgebra<Rationals()|15,-30>;
  order_basis_coords := [
      [ 1, 0, 0, 0 ],
      [ 1/2, 0, 1/4, 1/20 ],
      [ 1/2, 0, 1/4, -1/60 ],
      [ 0, 1/2, 1/4, 1/60 ]
  ];
  O := QuaternionOrder([A!c:c in order_basis_coords]);
  assert Discriminant(O) eq DB and IsMaximal(O);
  mu_base := A![ 1/2, 0, 1/4, -1/60 ];
  orientation_unit := A![ 1, 0, 0, 1/15 ];
  assert mu_base in O and orientation_unit in O and Norm(orientation_unit) eq -1;
  W3 := A![ 3/2, 0, 1/4, 1/20 ]; assert Norm(W3) eq 3;
  W5 := A![ -5/2, 0, 1/4, -1/12 ]; assert Norm(W5) eq 5;
  W15 := A![ -15/2, -1, -1/4, -1/4 ]; assert Norm(W15) eq 15;
  for W in [W3,W5,W15] do
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
  print "CM_DB = 15";
  print "CM_RAMIFIED_PRIMES = [ 3, 5 ]";
  print "CM_ABS_D_MIN = 3";
  print "CM_ABS_D_MAX = 199";
  print "CM_ALGEBRA_I2 = 15";
  print "CM_ALGEBRA_J2 = -30";
  print "CM_ORDER_BASIS = [";
  print "    [ 1, 0, 0, 0 ],";
  print "    [ 1/2, 0, 1/4, 1/20 ],";
  print "    [ 1/2, 0, 1/4, -1/60 ],";
  print "    [ 0, 1/2, 1/4, 1/60 ],";
  print "]";
  print "CM_BASE_COORDS = [ 1/2, 0, 1/4, -1/60 ]";
  print "CM_ORIENTATION_UNIT = [ 1, 0, 0, 1/15 ]";
  print "CM_AL_COORDS = {'w3': [ 3/2, 0, 1/4, 1/20 ], 'w5': [ -5/2, 0, 1/4, -1/12 ], 'w15': [ -15/2, -1, -1/4, -1/4 ]}";
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
   if chis[1] eq -1 then translators cat:= [<3*t[1],W3*t[2]>:t in translators]; end if;
   if chis[2] eq -1 then translators cat:= [<5*t[1],W5*t[2]>:t in translators]; end if;
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


end procedure;

procedure GenerateCMInputs21()
  /* All admissible negative fundamental D with 3 <= |D| <= 199.
     No class-number cutoff. Frozen algebra, maximal order, base point and AL lifts.
     Paste into https://magma.maths.usyd.edu.au/calc/; print-only, no file I/O.
     Fresh class representatives may differ from the bundled exact Sage input. */
  SetColumns(0); SetSeed(21);
  DB := 21; RAMIFIED_PRIMES := [ 3, 7 ];
  A := QuaternionAlgebra<Rationals()|21,-1>;
  order_basis_coords := [
      [ 1, 0, 0, 0 ],
      [ 1/2, 1/2, 0, 0 ],
      [ 0, 0, 1/2, 1/2 ],
      [ 0, 0, 0, 1 ]
  ];
  O := QuaternionOrder([A!c:c in order_basis_coords]);
  assert Discriminant(O) eq DB and IsMaximal(O);
  mu_base := A![ 1/2, 1/2, 10, 2 ];
  orientation_unit := A![4,1,2,0];
  assert mu_base in O and orientation_unit in O and Norm(orientation_unit) eq -1;
  W3 := A![ -3, 1, 33/2, -7/2 ]; assert Norm(W3) eq 3;
  W7 := A![ 0, 0, 7/2, -1/2 ]; assert Norm(W7) eq 7;
  W21 := A![ -21, 4, -21, 5 ]; assert Norm(W21) eq 21;
  for W in [W3,W7,W21] do
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
  print "CM_DB = 21";
  print "CM_RAMIFIED_PRIMES = [ 3, 7 ]";
  print "CM_ABS_D_MIN = 3";
  print "CM_ABS_D_MAX = 199";
  print "CM_ALGEBRA_I2 = 21";
  print "CM_ALGEBRA_J2 = -1";
  print "CM_ORDER_BASIS = [";
  print "    [ 1, 0, 0, 0 ],";
  print "    [ 1/2, 1/2, 0, 0 ],";
  print "    [ 0, 0, 1/2, 1/2 ],";
  print "    [ 0, 0, 0, 1 ],";
  print "]";
  print "CM_BASE_COORDS = [ 1/2, 1/2, 10, 2 ]";
  print "CM_ORIENTATION_UNIT = [4, 1, 2, 0]";
  print "CM_AL_COORDS = {'w3': [ -3, 1, 33/2, -7/2 ], 'w7': [ 0, 0, 7/2, -1/2 ], 'w21': [ -21, 4, -21, 5 ]}";
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
   if chis[1] eq -1 then translators cat:= [<3*t[1],W3*t[2]>:t in translators]; end if;
   if chis[2] eq -1 then translators cat:= [<7*t[1],W7*t[2]>:t in translators]; end if;
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


end procedure;

procedure GenerateCMInputs26()
  /* All admissible negative fundamental D with 3 <= |D| <= 199.
     No class-number cutoff. Frozen algebra, maximal order, base point and AL lifts.
     Paste into https://magma.maths.usyd.edu.au/calc/; print-only, no file I/O.
     Fresh class representatives may differ from the bundled exact Sage input. */
  SetColumns(0); SetSeed(26);
  DB := 26; RAMIFIED_PRIMES := [ 2, 13 ];
  A := QuaternionAlgebra<Rationals()|26,-5>;
  order_basis_coords := [
      [ 1, 0, 0, 0 ],
      [ 1/2, 0, 1/10, 1/10 ],
      [ 0, 0, 0, 1 ],
      [ 0, 1/2, 0, 1/2 ]
  ];
  O := QuaternionOrder([A!c:c in order_basis_coords]);
  assert Discriminant(O) eq DB and IsMaximal(O);
  mu_base := A![ 1/2, 1, -7/2, -1/2 ];
  orientation_unit := A![ 2, 0, -1/5, -1/5 ];
  assert mu_base in O and orientation_unit in O and Norm(orientation_unit) eq -1;
  W2 := A![ 6, -1, 8/5, -2/5 ]; assert Norm(W2) eq 2;
  W13 := A![ 0, 0, -13/5, 2/5 ]; assert Norm(W13) eq 13;
  W26 := A![ 0, -2, -26, 5 ]; assert Norm(W26) eq 26;
  for W in [W2,W13,W26] do
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
  print "CM_DB = 26";
  print "CM_RAMIFIED_PRIMES = [ 2, 13 ]";
  print "CM_ABS_D_MIN = 3";
  print "CM_ABS_D_MAX = 199";
  print "CM_ALGEBRA_I2 = 26";
  print "CM_ALGEBRA_J2 = -5";
  print "CM_ORDER_BASIS = [";
  print "    [ 1, 0, 0, 0 ],";
  print "    [ 1/2, 0, 1/10, 1/10 ],";
  print "    [ 0, 0, 0, 1 ],";
  print "    [ 0, 1/2, 0, 1/2 ],";
  print "]";
  print "CM_BASE_COORDS = [ 1/2, 1, -7/2, -1/2 ]";
  print "CM_ORIENTATION_UNIT = [ 2, 0, -1/5, -1/5 ]";
  print "CM_AL_COORDS = {'w2': [ 6, -1, 8/5, -2/5 ], 'w13': [ 0, 0, -13/5, 2/5 ], 'w26': [ 0, -2, -26, 5 ]}";
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
   if chis[2] eq -1 then translators cat:= [<13*t[1],W13*t[2]>:t in translators]; end if;
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


end procedure;

procedure GenerateCMInputs33()
  /* All admissible negative fundamental D with 3 <= |D| <= 199.
     No class-number cutoff. Frozen algebra, maximal order, base point and AL lifts.
     Paste into https://magma.maths.usyd.edu.au/calc/; print-only, no file I/O.
     Fresh class representatives may differ from the bundled exact Sage input. */
  SetColumns(0); SetSeed(33);
  DB := 33; RAMIFIED_PRIMES := [ 3, 11 ];
  A := QuaternionAlgebra<Rationals()|33,-1>;
  order_basis_coords := [
      [ 1, 0, 0, 0 ],
      [ 1/2, 1/2, 0, 0 ],
      [ 0, 0, 1/2, 1/2 ],
      [ 0, 0, 0, 1 ]
  ];
  O := QuaternionOrder([A!c:c in order_basis_coords]);
  assert Discriminant(O) eq DB and IsMaximal(O);
  mu_base := A![ 1/2, -1/2, -5, 0 ];
  orientation_unit := A![ 4, 1, 4, 0 ];
  assert mu_base in O and orientation_unit in O and Norm(orientation_unit) eq -1;
  W3 := A![ 3, 0, -3/2, 1/2 ]; assert Norm(W3) eq 3;
  W11 := A![ -11, -6, -77/2, 7/2 ]; assert Norm(W11) eq 11;
  W33 := A![ -33, -4, 0, -4 ]; assert Norm(W33) eq 33;
  for W in [W3,W11,W33] do
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
  print "CM_DB = 33";
  print "CM_RAMIFIED_PRIMES = [ 3, 11 ]";
  print "CM_ABS_D_MIN = 3";
  print "CM_ABS_D_MAX = 199";
  print "CM_ALGEBRA_I2 = 33";
  print "CM_ALGEBRA_J2 = -1";
  print "CM_ORDER_BASIS = [";
  print "    [ 1, 0, 0, 0 ],";
  print "    [ 1/2, 1/2, 0, 0 ],";
  print "    [ 0, 0, 1/2, 1/2 ],";
  print "    [ 0, 0, 0, 1 ],";
  print "]";
  print "CM_BASE_COORDS = [ 1/2, -1/2, -5, 0 ]";
  print "CM_ORIENTATION_UNIT = [ 4, 1, 4, 0 ]";
  print "CM_AL_COORDS = {'w3': [ 3, 0, -3/2, 1/2 ], 'w11': [ -11, -6, -77/2, 7/2 ], 'w33': [ -33, -4, 0, -4 ]}";
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
   if chis[1] eq -1 then translators cat:= [<3*t[1],W3*t[2]>:t in translators]; end if;
   if chis[2] eq -1 then translators cat:= [<11*t[1],W11*t[2]>:t in translators]; end if;
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


end procedure;

procedure GenerateCMInputs34()
  /* All admissible negative fundamental D with 3 <= |D| <= 199.
     No class-number cutoff. Frozen algebra, maximal order, base point and AL lifts.
     Paste into https://magma.maths.usyd.edu.au/calc/; print-only, no file I/O.
     Fresh class representatives may differ from the bundled exact Sage input. */
  SetColumns(0); SetSeed(34);
  DB := 34; RAMIFIED_PRIMES := [ 2, 17 ];
  A := QuaternionAlgebra<Rationals()|34,-3>;
  order_basis_coords := [
      [ 1, 0, 0, 0 ],
      [ 1/2, 0, 1/6, 1/3 ],
      [ 0, 0, 0, 1 ],
      [ 0, 1/2, 0, 1/2 ]
  ];
  O := QuaternionOrder([A!c:c in order_basis_coords]);
  assert Discriminant(O) eq DB and IsMaximal(O);
  mu_base := A![ 1/2, 1/2, -7/2, 1/2 ];
  orientation_unit := A![ 21, 1, 0, 2 ];
  assert mu_base in O and orientation_unit in O and Norm(orientation_unit) eq -1;
  W2 := A![ -26, 10, 118/3, -13/3 ]; assert Norm(W2) eq 2;
  W17 := A![ 0, -1, 17/3, -2/3 ]; assert Norm(W17) eq 17;
  W34 := A![ -714, 21, -1564/3, 340/3 ]; assert Norm(W34) eq 34;
  for W in [W2,W17,W34] do
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
  print "CM_DB = 34";
  print "CM_RAMIFIED_PRIMES = [ 2, 17 ]";
  print "CM_ABS_D_MIN = 3";
  print "CM_ABS_D_MAX = 199";
  print "CM_ALGEBRA_I2 = 34";
  print "CM_ALGEBRA_J2 = -3";
  print "CM_ORDER_BASIS = [";
  print "    [ 1, 0, 0, 0 ],";
  print "    [ 1/2, 0, 1/6, 1/3 ],";
  print "    [ 0, 0, 0, 1 ],";
  print "    [ 0, 1/2, 0, 1/2 ],";
  print "]";
  print "CM_BASE_COORDS = [ 1/2, 1/2, -7/2, 1/2 ]";
  print "CM_ORIENTATION_UNIT = [ 21, 1, 0, 2 ]";
  print "CM_AL_COORDS = {'w2': [ -26, 10, 118/3, -13/3 ], 'w17': [ 0, -1, 17/3, -2/3 ], 'w34': [ -714, 21, -1564/3, 340/3 ]}";
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
   if chis[2] eq -1 then translators cat:= [<17*t[1],W17*t[2]>:t in translators]; end if;
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


end procedure;

procedure GenerateCMInputs38()
  /* All admissible negative fundamental D with 3 <= |D| <= 199.
     No class-number cutoff. Frozen algebra, maximal order, base point and AL lifts.
     Paste into https://magma.maths.usyd.edu.au/calc/; print-only, no file I/O.
     Fresh class representatives may differ from the bundled exact Sage input. */
  SetColumns(0); SetSeed(38);
  DB := 38; RAMIFIED_PRIMES := [2,19];
  A := QuaternionAlgebra<Rationals()|-19,2>;
  order_basis_coords := [
    [ 1, 0, 0, 0 ],
    [ 1/2, 1/2, 0, 0 ],
    [ 0, 0, 1/2, 1/2 ],
    [ 0, 0, 0, 1 ]
  ];
  O := QuaternionOrder([A!c:c in order_basis_coords]);
  assert Discriminant(O) eq DB and IsMaximal(O);
  mu_base := A![ 1/2, -1/2, -1, 0 ];
  orientation_unit := A![ 1, 0, 1, 0 ];
  assert mu_base in O and orientation_unit in O and Norm(orientation_unit) eq -1;
  W2 := A![ 2, 0, 1, 0 ]; assert Norm(W2) eq 2;
  W19 := A![ 0, -1, 0, 0 ]; assert Norm(W19) eq 19;
  W38 := A![ 0, -2, 0, 1 ]; assert Norm(W38) eq 38;
  for W in [W2,W19,W38] do
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
  print "CM_DB = 38";
  print "CM_RAMIFIED_PRIMES = [2,19]";
  print "CM_ABS_D_MIN = 3";
  print "CM_ABS_D_MAX = 199";
  print "CM_ALGEBRA_I2 = -19";
  print "CM_ALGEBRA_J2 = 2";
  print "CM_ORDER_BASIS = [";
  print "  [ 1, 0, 0, 0 ],";
  print "  [ 1/2, 1/2, 0, 0 ],";
  print "  [ 0, 0, 1/2, 1/2 ],";
  print "  [ 0, 0, 0, 1 ],";
  print "]";
  print "CM_BASE_COORDS = [ 1/2, -1/2, -1, 0 ]";
  print "CM_ORIENTATION_UNIT = [ 1, 0, 1, 0 ]";
  print "CM_AL_COORDS = {'w2': [ 2, 0, 1, 0 ], 'w19': [ 0, -1, 0, 0 ], 'w38': [ 0, -2, 0, 1 ]}";
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
   if chis[2] eq -1 then translators cat:= [<19*t[1],W19*t[2]>:t in translators]; end if;
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


end procedure;

procedure GenerateCMInputs46()
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


end procedure;

procedure GenerateCMInputs58()
  /* All admissible negative fundamental D with 3 <= |D| <= 199.
     No class-number cutoff. Frozen algebra, maximal order, base point and AL lifts.
     Paste into https://magma.maths.usyd.edu.au/calc/; print-only, no file I/O.
     Fresh class representatives may differ from the bundled exact Sage input. */
  SetColumns(0); SetSeed(58);
  DB := 58; RAMIFIED_PRIMES := [2,29];
  A := QuaternionAlgebra<Rationals()|2,-29>;
  order_basis_coords := [
    [ 1, 0, 0, 0 ],
    [ 1/2, 0, 1/2, 1/2 ],
    [ 1/2, 0, 1/2, -1/2 ],
    [ 0, 1/2, 0, 1/2 ]
  ];
  O := QuaternionOrder([A!c:c in order_basis_coords]);
  assert Discriminant(O) eq DB and IsMaximal(O);
  mu_base := A![ 1/2, -3/2, -1/2, 0 ];
  orientation_unit := A![ 1, 1, 0, 0 ];
  assert mu_base in O and orientation_unit in O and Norm(orientation_unit) eq -1;
  W2 := A![ 2, 1, 0, 0 ]; assert Norm(W2) eq 2;
  W29 := A![ 0, 0, -1, 0 ]; assert Norm(W29) eq 29;
  W58 := A![ 0, 0, -2, -1 ]; assert Norm(W58) eq 58;
  for W in [W2,W29,W58] do
   assert &and[W*(A!o)*W^(-1) in O and W^(-1)*(A!o)*W in O:o in Basis(O)];
   assert W^2/Norm(W) in O;
  end for;
  side_coords := [
    [ 3, -2, 0, 0 ],
    [ 10, -13/2, 0, 1/2 ],
    [ 10, -8, -1, 0 ],
    [ 5/2, -5/2, -1/2, 0 ],
    [ 2, -4, -1, 0 ],
    [ 4, -1/2, 0, 1/2 ],
    [ 1/2, -19/2, -5/2, 0 ],
    [ 1/2, 19/2, 5/2, 0 ],
    [ 4, 1/2, 0, 1/2 ],
    [ 2, 4, 1, 0 ],
    [ 5/2, 5/2, 1/2, 0 ],
    [ 10, 8, 1, 0 ],
    [ 10, 13/2, 0, 1/2 ],
    [ 3, 2, 0, 0 ],
    [ 10, 13/2, 0, -1/2 ],
    [ 10, 8, -1, 0 ],
    [ 5/2, 5/2, -1/2, 0 ],
    [ 4, 1/2, 0, -1/2 ],
    [ 4, -1/2, 0, -1/2 ],
    [ 5/2, -5/2, 1/2, 0 ],
    [ 10, -8, 1, 0 ],
    [ 10, -13/2, 0, -1/2 ]
  ];
  assert &and[A!c in O and Norm(A!c) eq 1:c in side_coords];
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
  print "CM_DB = 58";
  print "CM_RAMIFIED_PRIMES = [2,29]";
  print "CM_ABS_D_MIN = 3";
  print "CM_ABS_D_MAX = 199";
  print "CM_ALGEBRA_I2 = 2";
  print "CM_ALGEBRA_J2 = -29";
  print "CM_ORDER_BASIS = [";
  print "  [ 1, 0, 0, 0 ],";
  print "  [ 1/2, 0, 1/2, 1/2 ],";
  print "  [ 1/2, 0, 1/2, -1/2 ],";
  print "  [ 0, 1/2, 0, 1/2 ],";
  print "]";
  print "CM_BASE_COORDS = [ 1/2, -3/2, -1/2, 0 ]";
  print "CM_ORIENTATION_UNIT = [ 1, 1, 0, 0 ]";
  print "CM_SIDE_COORDS = [";
  print "  [ 3, -2, 0, 0 ],";
  print "  [ 10, -13/2, 0, 1/2 ],";
  print "  [ 10, -8, -1, 0 ],";
  print "  [ 5/2, -5/2, -1/2, 0 ],";
  print "  [ 2, -4, -1, 0 ],";
  print "  [ 4, -1/2, 0, 1/2 ],";
  print "  [ 1/2, -19/2, -5/2, 0 ],";
  print "  [ 1/2, 19/2, 5/2, 0 ],";
  print "  [ 4, 1/2, 0, 1/2 ],";
  print "  [ 2, 4, 1, 0 ],";
  print "  [ 5/2, 5/2, 1/2, 0 ],";
  print "  [ 10, 8, 1, 0 ],";
  print "  [ 10, 13/2, 0, 1/2 ],";
  print "  [ 3, 2, 0, 0 ],";
  print "  [ 10, 13/2, 0, -1/2 ],";
  print "  [ 10, 8, -1, 0 ],";
  print "  [ 5/2, 5/2, -1/2, 0 ],";
  print "  [ 4, 1/2, 0, -1/2 ],";
  print "  [ 4, -1/2, 0, -1/2 ],";
  print "  [ 5/2, -5/2, 1/2, 0 ],";
  print "  [ 10, -8, 1, 0 ],";
  print "  [ 10, -13/2, 0, -1/2 ],";
  print "]";
  print "CM_AL_COORDS = {'w2': [ 2, 1, 0, 0 ], 'w29': [ 0, 0, -1, 0 ], 'w58': [ 0, 0, -2, -1 ]}";
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
   if chis[2] eq -1 then translators cat:= [<29*t[1],W29*t[2]>:t in translators]; end if;
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


end procedure;

procedure GenerateCMInputs(N)
  if N eq 14 then GenerateCMInputs14();
  elif N eq 15 then GenerateCMInputs15();
  elif N eq 21 then GenerateCMInputs21();
  elif N eq 26 then GenerateCMInputs26();
  elif N eq 33 then GenerateCMInputs33();
  elif N eq 34 then GenerateCMInputs34();
  elif N eq 38 then GenerateCMInputs38();
  elif N eq 46 then GenerateCMInputs46();
  elif N eq 58 then GenerateCMInputs58();
  else error "Unsupported level";
  end if;
end procedure;
