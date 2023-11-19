ClearAll["Global`*"];
Off[General::"spell",
  General::"spell1", InterpolatingFunction::"dmval"];

(*XPmatrices[dim_, parity_Integer:-1] :=
  Module[
    {x2array, p2array, x4array, x2, p2, x4, x1, x1array, p1, p1array},
*)

dim = 5;
parity = 0;

    Clear[i,n];
    
    singlepower = False;
    If[parity == 0 || parity == 1,
      Print["This is True"];
      n = 2*i - 2 + parity;
      shift = 1;,
      Print["This is False"];
      n = i - 1;
      shift = 2;
      singlepower = True
    ];

Print["singlepower = ", singlepower];
Print["n = ", n];
Print["shift = ", shift];


    x2 = SparseArray[
      { {i_, i_} -> (2*n + 1)/4, 
        {i_, k_} /; k - i == shift ->
          Sqrt[(n+1)*(n+2)]/2}, {dim,dim}
    ];

Print["x2 = "];
Print[MatrixForm[x2]];

    x2 = (x2 + Transpose[x2]);

Print["x2 = "];
Print[MatrixForm[x2]];


p2 = SparseArray[ {
  {i_, i_} -> (2*n + 1)/4,
  {i_, k_} /; k - i == shift -> -Sqrt[(n+1)*(n+1)]/2
}, {dim, dim} ];

Print["\np2"]
Print[MatrixForm[p2]];

p2 = p2 + Transpose[p2];

Print["\np2"]
Print[MatrixForm[p2]];

x4 = SparseArray[ {
 {i_, i_} -> 3/8*(2*n^2 + 2*n + 1),
 {i_, k_} /; k - i == shift -> (4*n + 6)/4 * Sqrt[(n+1)*(n+2)],
 {i_, k_} /; k - i == 2*shift -> 1/4 * Sqrt[(n+1)*(n+2)*(n+3)*(n+4)]
}, {dim,dim} ];

Print["\nx4"]
Print[MatrixForm[x4]];


(* ]; *)
(* End of module *)



Exit[]
