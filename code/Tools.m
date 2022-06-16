(* ::Package:: *)

BeginPackage["Tools`"];


MakeSystem::usage = "MakeSystem[var, vart, derivative, sys] return a dynamical system to be used for NDSolve 
var: variables in an ode system (form: {var1, var2});
vart: variables in the ode system under the form of functions of time (form: {var1[t], var2[t]};
derivative: derivative of variables with respect to t (form: {var1'[t], var2'[t]};
sys: the ode system (form: {a var1 + b var2,  x var2 - y var2 var3});
"


FollowRoot::usage = "FollowRoot[system, commonpars, followPar, range, variables, initialEq] return a list of equilibrium with respect to values of a parameter (using function FindRoot)"


FollowRootTwoParameters::usage = "FollowRootTwoParameters[system, commonpars, followPar1, range1, followPar2, range2, variables, initialEq, minEqlimit, maxiteration:1000, precisionGoal:MachinePrecision] return a list of equilibrium with respect to values of two parameters"


ListStableMark::usage = "ListStableMark[jacobmatrix, parcommon, parfollow, range, equiList] return a list of symbols that correspond to the stability of the equilibrium, symbol '*' if the equilibrium is unstable, symbol '.' if the equilibrium is stable(form: {*,*,.,.}) 
jacobmatrix: jacobian matrix of the system (form: {{a, b}, {c, d}});
parcommon: values of fixed parameters (form: {p1 -> 2, p2-> 3};
parfollow: bifurcation parameter (form: symbol of the bifurcation parameter);
range: values of the bifurcation parameter (form: {0.1, 0.2});
equiList: list of values of equilibrium corresponds with different values of bifurcation parameter (form: {{var1 -> 4.5, var2 -> 44}, {var1 -> 4.6, var2-> 45}};
"


NSolveCodim2PositiveTest::usage = "NSolveCodim2Positive[system, commonpars, bifurpar1, bifurpar2, bfparsName, equisymbol, variables] Find only positive solution varying two parameters, using NSolve (To be used in Parallele Table to follow the positive equilibrium)
system: the odes to be solved (form: dXdt = A * X, where X is the vector of variables and A is the characteristic matrix);\[IndentingNewLine]commonpars: list of the value of fixed parameters (form: {p1 -> 4, p2 -> 2, ..., pn -> 43});\[IndentingNewLine]bifurpar1: the value of the bifurcation parameter 1 (form: {pb1 -> 3.4});\[IndentingNewLine]bifurpar2: the value of the bifurcation parameter 2 (form: {pb2 -> 1.2});\[IndentingNewLine]bfparsName: symbol of the two bifurcation parameter (form: pb1pb2);\[IndentingNewLine]variables: list of the variables to be solved (form: {v1, v2, v3});\[IndentingNewLine]return a list of the values of the two bifurcation parameters and its corresponding positive equilibrium\[IndentingNewLine](form: {pb1pb2 -> {2, 2}, v1-> 4, v2-> 3, v3-> 4.5}}
"


Begin["`Private`"]


MakeSystem[var_, t_, sys_]:=
Module[{varOft, sysOft},
varOft = Through[var[t]];
sysOft = sys/.Thread[var-> varOft]; 
Thread[D[varOft, t] == sysOft]]


FollowRoot[system_,commonpars_,followPar_,range_,variables_, initialEq_]:=
Module[{sys,cpar,fp, r,v,en, ieq, eq, pars, results, initValues, valuesIfFalse},
sys=system;
cpar=commonpars;
fp=followPar; 
r = range;
v=variables;
ieq = initialEq;
results={};

Do[eq=ieq;
pars = Join[cpar, {fp-> i}];
initValues = Thread[{v, v/.ieq, 0, Infinity}];
valuesIfFalse = Thread[v->-1];
ieq = Check[FindRoot[Thread[sys == 0]/.pars, initValues, MaxIterations->1000], valuesIfFalse, FindRoot::reged];
If[AllTrue[v/.ieq, #==-1&], Break[]];
results = AppendTo[results, ieq], {i, range}];
results]


FollowRootTwoParameters[system_,commonpars_,followPar1_,range1_,followPar2_, range2_,variables_, initialEq_, minEqlimit_, maxiteration_:1000, precisionGoal_:MachinePrecision]:=Module[{sys,cpar,fp1, r1, fp2, r2,v,en, ieq,mineqlim,maxiter,precgoal, pars, results},sys=system;
cpar=commonpars;
fp1=followPar1;
r1 = range1;
fp2 = followPar2;
r2 = range2;
v=variables;
ieq = initialEq;
mineqlim = minEqlimit;
maxiter = maxiteration;
precgoal = precisionGoal;
results={};
Do[pars = Join[cpar, {fp1-> i, fp2-> j}];
ieq = Check[FindRoot[Thread[sys == ConstantArray[0, Length[v]]]/.pars,Thread[{v, v/.ieq, minEqlimit, Infinity}], MaxIterations->maxiter, PrecisionGoal->precgoal], Thread[v-> -1]];
If[AllTrue[v/.ieq, # == -1&], Break[]];
results = AppendTo[results, ieq], {i, range1}, {j, range2}];
results];


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
Module[{jmat, pc, pf, r,eq, eiv, anyZero, anyPos, allNeg, ps},
jmat = jacobmatrix;
pc = parcommon;
pf = parsfollow;
eq = equilibrium;
eiv = Eigenvalues[jmat/.pc/.pf/.eq];
anyZero = AnyTrue[Thread[-10^-10<=Re[eiv]<=10^-10], TrueQ];
anyPos = AnyTrue[Thread[Re[eiv]>10^-10], TrueQ];
allNeg = AllTrue[Thread[Re[eiv]<-10^-10], TrueQ];
Which[anyZero, "*",anyPos, "*", allNeg, "."]]


ListStableMark[jacobmatrix_, parcommon_, parfollow_, range_,equiList_]:=
Module[{jmat, pc, pf, r,eql, eiv, ps},
On[Assert];
Assert[Head[jacobmatrix] == List, "the jacobian matrix has not been defined"];
jmat = jacobmatrix;
pc = parcommon;
pf = parfollow;
r=range;
eql = equiList;
MapThread[SingleStableMark,{ConstantArray[jmat, Length[eql]], ConstantArray[pc, Length[eql]],Thread[pf->r], eql}]]


NSolveCodim2PositiveTest[system_, commonpars_, bifurpar1_, bifurpar2_,bfparsName_, equisymbol_,variables_]:=
Module[
{sys, cpar,bfpar1, bfpar2,  vars, bfpName,eqsym,parspairVal, eqAll, eqpos},
sys = system;
cpar = commonpars;
bfpar1 = bifurpar1;
bfpar2 = bifurpar2;
bfpName = bfparsName;
eqsym = equisymbol;
vars = variables;
eqAll = NSolve[Thread[(sys/.cpar/.bfpar1/. bfpar2)==0], vars, Reals];
eqpos = Select[eqAll,And@@Thread[(vars/.# )> 0]&];
parspairVal = {bfpar1[[1]][[2]], bfpar2[[1]][[2]]};
If[eqpos != {}, Join[{bfpName -> parspairVal}, {eqsym -> eqpos[[1]]}] , Unevaluated@Sequence[]]
]


End[]


EndPackage[]


(* ::Input:: *)
(**)
