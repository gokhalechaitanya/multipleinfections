(* ::Package:: *)

BeginPackage["Tools`"];


MakeSystem::usage = "MakeSystem[var, vart, derivative, sys] return a dynamical system to be used for NDSolve 
var: variables in an ode system (form: {var1, var2});
vart: variables in the ode system under the form of functions of time (form: {var1[t], var2[t]};
derivative: derivative of variables with respect to t (form: {var1'[t], var2'[t]};
sys: the ode system (form: {a var1 + b var2,  x var2 - y var2 var3});
"


FollowRoot::usage = "FollowRoot[system, commonpars, followPar, range, variables, initialEq] return a list of equilibrium with respect to values of a parameter (using function FindRoot)"


ListStableMark::usage = "ListStableMark[jacobmatrix, parcommon, parfollow, range, equiList] return a list of symbols that correspond to the stability of the equilibrium, symbol '*' if the equilibrium is unstable, symbol '.' if the equilibrium is stable(form: {*,*,.,.}) 
jacobmatrix: jacobian matrix of the system (form: {{a, b}, {c, d}});
parcommon: values of fixed parameters (form: {p1 -> 2, p2-> 3};
parfollow: bifurcation parameter (form: symbol of the bifurcation parameter);
range: values of the bifurcation parameter (form: {0.1, 0.2});
equiList: list of values of equilibrium corresponds with different values of bifurcation parameter (form: {{var1 -> 4.5, var2 -> 44}, {var1 -> 4.6, var2-> 45}};
"


NSolveCodim2Positive::usage = "NSolveCodim2Positive[system, commonpars, bifurpar1, bifurpar2, bfparsName, equisymbol, variables] Find only positive solution varying two parameters, using NSolve (To be used in Parallele Table to follow the positive equilibrium)
system: the odes to be solved (form: dXdt = A * X, where X is the vector of variables and A is the characteristic matrix);\[IndentingNewLine]commonpars: list of the value of fixed parameters (form: {p1 -> 4, p2 -> 2, ..., pn -> 43});\[IndentingNewLine]bifurpar1: the value of the bifurcation parameter 1 (form: {pb1 -> 3.4});\[IndentingNewLine]bifurpar2: the value of the bifurcation parameter 2 (form: {pb2 -> 1.2});\[IndentingNewLine]bfparsName: symbol of the two bifurcation parameter (form: pb1pb2);\[IndentingNewLine]variables: list of the variables to be solved (form: {v1, v2, v3});\[IndentingNewLine]return a list of the values of the two bifurcation parameters and its corresponding positive equilibrium\[IndentingNewLine](form: {pb1pb2 -> {2, 2}, v1-> 4, v2-> 3, v3-> 4.5}}
"


ListStableMarkCodim2::usage = "ListStableMarkTwoParameters[jacobmatrix_, parcommon_, bifurParNames_,listbifurParsValues_, equiList_, useColor_] Mark a list of equilibrium corresponding to two bifurcation values, depending on its stability\[IndentingNewLine]jacobmatrix: jacobian matrix of the system (form: {{a, b}, {c, d}});\[IndentingNewLine]parcommon: values of fixed parameters (form: {p1 -> 2, p2-> 3};\[IndentingNewLine]bifurParNames: name of the two bifurcation parameters (form: {bifurpar1 , bifurpar2});\[IndentingNewLine]listbifurParsValues: list of corresponding pair of bifurcation parameters with the list of equilibrium (form: {{3.4, 3}, {2.4, 4}})\[IndentingNewLine]equiList: list of values of equilibrium corresponds with different values of bifurcation parameters (form: {{var1 -> 4.5, var2 -> 44}, {var1 -> 4.6, var2-> 45}};\[IndentingNewLine]useColor: True if the stability is marked as different color, False if the stability is marked as different markers\[IndentingNewLine]a list of symbols that correspond to the stability of the equilibrium, symbol '*' if the equilibrium is unstable, symbol ' . ' if the equilibrium is stable\[IndentingNewLine](form: {*,*,.,.})\[IndentingNewLine]If colorUse == True, the stable equilibrium will have different color than the unstable equilibrium
"


Begin["`Private`"]


MakeSystem[var_, t_, sys_]:=
Module[{varOft, sysOft},
varOft = Through[var[t]];
sysOft = sys/.Thread[var-> varOft]; 
Thread[D[varOft, t] == sysOft]]


FollowRoot[system_, commonpars_, followPar_, range_, variables_, initialEq_]:=
Module[{eq, fullpars, results, initValues, valuesIfFalse, i, initEq},
initEq = initialEq;
results={};

Do[eq=initEq;
	fullpars = Join[commonpars, {followPar -> i}];
	initValues = Thread[{variables, variables/.initEq, 0, Infinity}];
	valuesIfFalse = Thread[variables -> -1];

	initEq = Check[
					  FindRoot[Thread[system == 0]/.fullpars, initValues, MaxIterations->1000], 
					  valuesIfFalse, 
					  FindRoot::reged
					 ];

	If[AllTrue[variables/.initEq, #==-1&], Break[]];
	
	results = AppendTo[results, initEq], {i, range}
	];

results]


SingleStableMark[jacobmatrix_, parcommon_, parsfollow_, equilibrium_]:=
(* Mark an equilibrium depending on its stability 
Input:
jacobmatrix: jacobian matrix of the system (form: {{a, b}, {c, d}});
parcommon: values of fixed parameters (form: {p1 -> 2, p2-> 3};
parsfollow: bifurcation parameter which could be one or two parameters(form: {p4 -> 3.4} for one parameter, {p2 -> 1, p4 -> 4} for two parameters);
equilibrium: values of equilibrium (form: {var1 -> 4.5, var2 -> 44};
Output:
 symbol "*" if the equilibrium is unstable, symbol "." if the equilibrium is stable
*)
Module[{eiv, anyZero, anyPos, allNeg},

On[Assert];
Assert[Head[jacobmatrix] == List, "the jacobian matrix has not been defined"];

eiv = Eigenvalues[jacobmatrix/.parcommon/.parsfollow/.equilibrium];
anyZero = AnyTrue[Thread[-10^-10<=Re[eiv]<=10^-10], TrueQ];
anyPos = AnyTrue[Thread[Re[eiv]>10^-10], TrueQ];
allNeg = AllTrue[Thread[Re[eiv]<-10^-10], TrueQ];
Which[anyZero, "*",anyPos, "*", allNeg, "."]]


SingleStableColor[jacobmatrix_, parcommon_, parsfollow_, equilibrium_, colorlist_]:=
Module[{eiv, anyZero, anyPos, allNeg},

On[Assert];
Assert[Head[jacobmatrix] == List, "the jacobian matrix has not been defined"];

eiv = Eigenvalues[jacobmatrix/.parcommon/.parsfollow/.equilibrium];
anyZero = AnyTrue[Thread[-10^-10<=Re[eiv]<=10^-10], TrueQ];
anyPos = AnyTrue[Thread[Re[eiv]>10^-10], TrueQ];
allNeg = AllTrue[Thread[Re[eiv]<-10^-10], TrueQ];
Which[anyZero, colorlist[[2]],anyPos, colorlist[[1]], allNeg, colorlist[[1]]]]


ListStableMark[jacobmatrix_, parcommon_, parfollow_, range_, equiList_]:= 
MapThread[
		SingleStableMark,
		{ConstantArray[jacobmatrix, Length[equiList]], 
		ConstantArray[parcommon, Length[equiList]],
		Thread[parfollow->range], 
		equiList}
		]


NSolveCodim2Positive[system_, commonpars_, bifurpar1_, bifurpar2_,bfparsName_, equisymbol_,variables_]:=
Module[
{sysForNSolve, parspairVal, eqAll, eqpos, numberEq, i},

sysForNSolve = Thread[(system/.commonpars/. bifurpar1/. bifurpar2)==0];
eqAll = NSolve[sysForNSolve, variables, Reals];
eqpos = Select[eqAll, AllTrue[variables/.#, Positive]&];
parspairVal = {bifurpar1[[1]][[2]], bifurpar2[[1]][[2]]};
numberEq = Length[eqpos];
Print[{parspairVal, eqpos}];
If[
	numberEq == 1, 
	Join[{bfparsName -> parspairVal}, {equisymbol -> eqpos[[1]]}],
	If[numberEq > 1,
		For[i = 1, 
			i <= Length[eqpos], 
			i++, 
			Join[{bfparsName -> parspairVal}, {equisymbol -> eqpos[[i]]}]
			],
			{}
		]
	]
]


ListStableMarkCodim2[jacobmatrix_, parcommon_, bifurParNames_,listbifurParsValues_, equiList_, useColor_:False, colorlist_:Null]:=
Module[{eiv, marklist, listval},

listval = Thread[bifurParNames->#]&/@listbifurParsValues;

If[
	useColor == True, 
	On[Assert];
	colorlist != Null //Assert
	];
	
If[
	useColor,
	marklist = MapThread[
						SingleStableColor,
						{ConstantArray[jacobmatrix, Length[equiList]], 
						 ConstantArray[parcommon, Length[equiList]],
						 listval, 
						 equiList,
						 colorlist}
						],
	marklist = MapThread[
						SingleStableMark,
						{ConstantArray[jacobmatrix, Length[equiList]], 
						 ConstantArray[parcommon, Length[equiList]],
						 listval, 
						 equiList}
						 ]
	]
]


End[]


EndPackage[]


(* ::Input:: *)
(**)
