(* Claim: for N=48 cell-centred points, even reflected +/-2 stencils
   have half-width 2 and no repeated independent column modulo 5.
   Exact witness: exhaust all rows/offsets in Z, including both poles.
   Fields and Cartesian derivatives in Theta_plus are point-local (source audit).
   This verifies the support/color claim, not the physical expansion formula. *)
Module[{n = 48, mirror, supports, bandwidth, unique, d1, d2},
  mirror[j_] := Which[j < 0, -j-1, j >= n, 2*n-1-j, True, j];
  supports = Table[Union[mirror /@ Range[i-2, i+2]], {i, 0, n-1}];
  bandwidth = Max[Flatten[Table[Abs[supports[[i+1]]-i], {i, 0, n-1}]]];
  unique = And @@ (DuplicateFreeQ[Mod[#, 5]] & /@ supports);
  d1 = {1, -8, 0, 8, -1}/12;
  d2 = {-1, 16, -30, 16, -1}/12;
  Print["N=48; integer support witness; even poles; coefficients over Q"];
  Print[{"half_width", bandwidth, "five_color_independence", unique,
    "first_moments", Table[Sum[d1[[k+3]] If[p == 0, 1, k^p], {k,-2,2}], {p,0,4}],
    "second_moments", Table[Sum[d2[[k+3]] If[p == 0, 1, k^p], {k,-2,2}], {p,0,4}]}];
  If[bandwidth != 2 || !TrueQ[unique] ||
    Table[Sum[d1[[k+3]] If[p == 0, 1, k^p], {k,-2,2}], {p,0,4}] != {0,1,0,0,0} ||
    Table[Sum[d2[[k+3]] If[p == 0, 1, k^p], {k,-2,2}], {p,0,4}] != {0,0,2,0,0}, Exit[1]];
  Print["PROVED: N48 reflected stencil/color witness (CLI-fallback)"];
];
