(* ::Package:: *)

n 


BeginPackage["Tools`"];


MakeSystem::usage = "MakeSystem[var, t, sys] return a dynamical system to be used for NDSolve 
var: variables in an ode system (form: {var1, var2});
t: time;
sys: the ode system (form: {a var1 + b var2,  x var2 - y var2 var3});
"


FollowRoot::usage = "FollowRoot[system, commonpars, followPar, range, variables, initialEq] 
return a list of equilibrium with respect to values of a parameter (using function FindRoot)"


NSolvePositive::usage = "NSolvePositive[system_, commonpars_, followPar_, variables_, equisymbol_]
Find all positive solutions using NSolve (can be use in parallel computation)
system: the dynamical system
commonpars: common parameters
followPar: bifurcation parameter and its value (form: {a -> 0.1}
variables: names of the variables of the dynamical system (form: {N1, N2})
equisymbol: the name of the solution to make rule (used for latter extraction of the solutions)

return a list of bifurcation parameter values and their corresponding solutions
"


ListStableMark::usage = "ListStableMark[jacobmatrix, parcommon, parfollow, range, equiList] return a list of symbols that correspond to the stability of the equilibrium, symbol '*' if the equilibrium is unstable, symbol '.' if the equilibrium is stable(form: {*,*,.,.}) 
jacobmatrix: jacobian matrix of the system (form: {{a, b}, {c, d}});
parcommon: values of fixed parameters (form: {p1 -> 2, p2-> 3};
parfollow: bifurcation parameter (form: symbol of the bifurcation parameter);
range: values of the bifurcation parameter (form: {0.1, 0.2});
equiList: list of values of equilibrium corresponds with different values of bifurcation parameter (form: {{var1 -> 4.5, var2 -> 44}, {var1 -> 4.6, var2-> 45}};
"


NSolveCodim2Positive::usage = "NSolveCodim2Positive[system, commonpars, bifurpar1, bifurpar2, bfparsName, equisymbol, variables]
 Find only positive solution varying two parameters, using NSolve (To be used in Parallele Table to follow the positive equilibrium)
system: the odes to be solved (form: dXdt = A * X, where X is the vector of variables and A is the characteristic matrix);\[IndentingNewLine]commonpars: list of the value of fixed parameters (form: {p1 -> 4, p2 -> 2, ..., pn -> 43});\[IndentingNewLine]bifurpar1: the value of the bifurcation parameter 1 (form: {pb1 -> 3.4});\[IndentingNewLine]bifurpar2: the value of the bifurcation parameter 2 (form: {pb2 -> 1.2});\[IndentingNewLine]bfparsName: symbol of the two bifurcation parameter (form: pb1pb2);\[IndentingNewLine]variables: list of the variables to be solved (form: {v1, v2, v3});\[IndentingNewLine]return a list of the values of the two bifurcation parameters and its corresponding positive equilibrium\[IndentingNewLine](form: {pb1pb2 -> {2, 2}, v1-> 4, v2-> 3, v3-> 4.5}}
"


ListStableMarkTwoParameters::usage = "ListStableMarkTwoParameters[jacobmatrix, parcommon, bifurParNames, listbifurParsValues, equiList, useColor]
Mark a list of equilibrium corresponding to two bifurcation values, depending on its stability
jacobmatrix: jacobian matrix of the system (form: {{a, b}, {c, d}});
parcommon: values of fixed parameters (form: {p1 -> 2, p2-> 3};
bifurParNames: name of the two bifurcation parameters (form: {bifurpar1 , bifurpar2});
listbifurParsValues: list of corresponding pair of bifurcation parameters with the list of equilibrium (form: {{3.4, 3}, {2.4, 4}})
equiList: list of values of equilibrium corresponds with different values of bifurcation parameters (form: {{var1 -> 4.5, var2 -> 44}, {var1 -> 4.6, var2-> 45}};
useColor: True if the stability is marked as different color, False if the stability is marked as different markers
Output: 
a list of symbols that correspond to the stability of the equilibrium, symbol '*' if the equilibrium is unstable, symbol '.' if the equilibrium is stable
(form: {*,*,.,.})
If colorUse == True, the stable equilibrium will have different color than the unstable equilibrium
"


Begin["`Private`"]


MakeSystem[var_, t_, sys_]:=
Module[{varOft, sysOft},
varOft = Through[var[t]];
sysOft = sys/.Thread[var-> varOft]; 
Thread[D[varOft, t] == sysOft]]


FollowRoot[system_,commonpars_,followPar_,range_,variables_, initialEq_]:=
Module[{iEq, eq, pars, results, initValues, isResultsWithinRange, checkresults},
iEq = initialEq; 
results={};

Do[
	eq=iEq;
	pars = Join[commonpars, {followPar -> i}];
	initValues = Thread[{variables, variables/.iEq, 0, Infinity}];
	isResultsWithinRange = Thread[variables->-1];
	iEq = Check[
				FindRoot[Thread[system == 0]/.pars, initValues, MaxIterations->2000, PrecisionGoal->Infinity], 
				isResultsWithinRange, 
				FindRoot::reged
				];
	checkresults = system/.pars/.iEq;
	If[
		AllTrue[variables/.iEq, #==-1&] || AnyTrue[checkresults, # > 10^-10&], 
		Break[]
		];
	results = AppendTo[results, iEq], {i, range}];
	results
	]


NSolvePositive[system_, commonpars_, followPar_, variables_, equisymbol_]:=
Module[
{eqAll, eqpos, nbpos, fparlist},
On[Assert];
Assert[Head[followPar]=== Rule];
eqAll = NSolve[Thread[(system/.commonpars/.followPar) == 0], variables, Reals];
eqpos = Select[eqAll, AllTrue[variables/.# , Positive]&];
nbpos = Length[eqpos];
If[
	eqpos == {}, 
	Unevaluated@Sequence[],
	If[
		nbpos == 1,
		{Join[{followPar}, {equisymbol -> eqpos[[1]]}]},
		fparlist = ConstantArray[followPar, nbpos];
		Thread[{fparlist, Thread[equisymbol -> eqpos]}]
		]
	]
]


SingleStableMark[jacobmatrix_, parcommon_, parsfollow_, equilibrium_, markerlist_:{"@", "*", "."}]:=
(* Mark an equilibrium depending on its stability 
Input:
jacobmatrix: jacobian matrix of the system (form: {{a, b}, {c, d}});
parcommon: values of fixed parameters (form: {p1 -> 2, p2-> 3};
parsfollow: bifurcation parameter which could be one or two parameters(form: {p4 -> 3.4} for one parameter, {p2 -> 1, p4 -> 4} for two parameters);
equilibrium: values of equilibrium (form: {var1 -> 4.5, var2 -> 44};
Output:
 symbol "*" if the equilibrium is unstable, symbol "." if the equilibrium is stable
*)
Module[
{eiv, anyZero, anyPos, allNeg, ps},
eiv = Eigenvalues[jacobmatrix/.parcommon/.parsfollow/.equilibrium];
anyZero = AnyTrue[Thread[-10^-10<=Re[eiv]<=10^-10], TrueQ];
anyPos = AnyTrue[Thread[Re[eiv]>10^-10], TrueQ];
allNeg = AllTrue[Thread[Re[eiv]<-10^-10], TrueQ];
Which[anyZero, markerlist[[1]], anyPos, markerlist[[2]], allNeg, markerlist[[3]]]
]


SingleStableColor[jacobmatrix_, parcommon_, parsfollow_, equilibrium_, colorlist_, opacity_:0.3, pointsize_:0.03]:=
Module[{eiv, anyZero, anyPos, allNeg},
eiv = Eigenvalues[jacobmatrix/.parcommon/.parsfollow/.equilibrium];
anyZero = AnyTrue[Thread[-10^-10<=Re[eiv]<=10^-10], TrueQ];
anyPos = AnyTrue[Thread[Re[eiv]>10^-10], TrueQ];
allNeg = AllTrue[Thread[Re[eiv]<-10^-10], TrueQ];
Which[
	anyZero, 
	Directive[colorlist[[3]], Opacity[opacity], PointSize[pointsize]],
	anyPos, 
	Directive[colorlist[[2]], Opacity[opacity], PointSize[pointsize]], 
	allNeg, 
	Directive[colorlist[[1]], Opacity[opacity], PointSize[pointsize]]
	]
]


ListStableMark[jacobmatrix_, parcommon_, parfollow_, range_,equiList_, markerCode_:{"@", "*", "."}]:=
 MapThread[SingleStableMark,
			{ConstantArray[jacobmatrix, Length[equiList]], 
			 ConstantArray[parcommon, Length[equiList]], 
			 Thread[parfollow -> range], 
			 equiList}
]


NSolveCodim2Positive[system_, commonpars_, bifurpar1_, bifurpar2_,bfparsName_, equisymbol_,variables_]:=
Module[
{eqAll, parspairVal, eqpos, nbpos, fparlist},
eqAll = NSolve[Thread[(system/.commonpars/.bifurpar1/.bifurpar2)==0], variables, Reals];
eqpos = Select[eqAll,And@@Thread[(variables/.# )> 0]&];
nbpos = Length[eqpos];
parspairVal = {bifurpar1[[1]][[2]], bifurpar2[[1]][[2]]};
If[
	nbpos == 0, 
	Unevaluated@Sequence[],
	If[
		nbpos == 1,
		{Join[{bfparsName -> parspairVal}, {equisymbol -> eqpos[[1]]}]},
		fparlist = ConstantArray[bfparsName -> parspairVal, nbpos];
		Thread[{fparlist, Thread[equisymbol -> eqpos]}]
		]		
	]
]


ListStableMarkTwoParameters[jacobmatrix_, parcommon_, bifurParNames_,listbifurParsValues_, 
							equiList_, markerCode_:{"@", "*", "."}, useColor_:False]:=
Module[{marklist, listval, markers, lenEqList},
listval = Thread[bifurParNames->#]&/@ listbifurParsValues;
lenEqList = Length[equiList];
markers = ConstantArray[markerCode, lenEqList];
If[useColor,
	Head[markerCode[[1]]] === RGBColor;
	marklist = MapThread[
						SingleStableColor,
						{ConstantArray[jacobmatrix, lenEqList], 
						 ConstantArray[parcommon, lenEqList],
						 listval, 
						 equiList,
						 markers}],
	marklist = MapThread[SingleStableMark,
						{ConstantArray[jacobmatrix, lenEqList], 
						 ConstantArray[parcommon, lenEqList],
						 listval, 
						 equiList,
						 markers}]]
]


End[]


EndPackage[]


(* ::Input:: *)
(**)
