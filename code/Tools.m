(* ::Package:: *)

BeginPackage["Tools`"];


MakeSystem::usage = "MakeSystem[var, t, sys] return a dynamical system to be used for NDSolve 
var (List): variables in an ode system (form: {var1, var2});
t (Symbol): time;
sys (List): the ode system (form: {a var1 + b var2,  x var2 - y var2 var3});
"


FollowRoot::usage = "FollowRoot[system,commonpars,followPar,range,variables, initialEq] 
system (List): the dynamical system
commonpars (List): fixed parameters
followPar (Symbol): name of the bifurcation parameter
range (List): values of the bifurcation parameter
variables (List): list of Symbol of the variables
initialEq (List): initial values of the equilibrium (from: {x1 -> value1, x2 -> value2})
return a list of equilibrium with respect to values of a parameter (using function FindRoot)"


NSolvePositive::usage = "NSolvePositive[system, commonpars, followPar, variables, equisymbol]
Find all positive solutions using NSolve (can be use in parallel computation)

system (List): the dynamical system
commonpars (List of Rule): common parameters
followPar (Rule): bifurcation parameter value (form: a -> 0.1)
variables (List of Symbol): names of the variables of the dynamical system (form: {N1, N2})
equisymbol (Symbol): the name of the solution to make rule (used for latter extraction of the solutions)

return a list of bifurcation parameter values and their corresponding solutions
"


ListStableMark::usage = "ListStableMark[jacobmatrix, parcommon, parfollow, range, equiList, markerCode:{'@', '*', '.'}
return a list of symbols that correspond to the stability of the equilibrium 

jacobmatrix (Nested List): jacobian matrix of the system (form: {{a, b}, {c, d}});
parcommon (List of Rule): values of fixed parameters (form: {p1 -> 2, p2-> 3};
parfollow (Symbol): bifurcation parameter (form: symbol of the bifurcation parameter);
range (List): values of the bifurcation parameter (form: {0.1, 0.2});
equiList (Nested List of Rule): list of values of equilibrium corresponds with different values of bifurcation parameter (form: {{var1 -> 4.5, var2 -> 44}, {var1 -> 4.6, var2-> 45}};
markerCode (List): by default cycle: '@', unstable equilibrium '*', stable equilibrium: '.'
"


NSolveCodim2Positive::usage = "NSolveCodim2Positive[system, commonpars, bifurpar1, bifurpar2, bfparsName, equisymbol, variables]
 Find only positive solution varying two parameters, using NSolve (To be used in Parallele Table to follow the positive equilibrium)

system (List): the odes to be solved (form: dXdt = A * X, where X is the vector of variables and A is the characteristic matrix);\[IndentingNewLine]commonpars (List of Rule): list of the value of fixed parameters (form: {p1 -> 4, p2 -> 2, ..., pn -> 43});\[IndentingNewLine]bifurpar1 (Rule): the value of the bifurcation parameter 1 (form: {pb1 -> 3.4});\[IndentingNewLine]bifurpar2 (Rule): the value of the bifurcation parameter 2 (form: {pb2 -> 1.2});\[IndentingNewLine]bfparsName: symbol of the two bifurcation parameter (form: pb1pb2);\[IndentingNewLine]variables (List of Symbol): list of the variables to be solved (form: {v1, v2, v3});\[IndentingNewLine]return a list of the values of the two bifurcation parameters and its corresponding positive equilibrium\[IndentingNewLine](form: {pb1pb2 -> {2, 2}, v1-> 4, v2-> 3, v3-> 4.5}}
"


ListStableMarkTwoParameters::usage = "ListStableMarkTwoParameters[jacobmatrix, parcommon, bifurParNames, listbifurParsValues, equiList, markerCode:{'@', '*', '.'}, useColor:False]
Mark a list of equilibrium corresponding to two bifurcation values, depending on its stability

jacobmatrix (Nested List): jacobian matrix of the system (form: {{a, b}, {c, d}});
parcommon (List of Rule): values of fixed parameters (form: {p1 -> 2, p2-> 3};
bifurParNames (List of Symbol): name of the two bifurcation parameters (form: {bifurpar1 , bifurpar2});
listbifurParsValues (Nested List): list of corresponding pair of bifurcation parameters with the list of equilibrium (form: {{3.4, 3}, {2.4, 4}})
equiList (Nested List of Rule): list of values of equilibrium corresponds with different values of bifurcation parameters (form: {{var1 -> 4.5, var2 -> 44}, {var1 -> 4.6, var2-> 45}};
useColor (Boolean) (default = False): True if colors are used to mark stability, False if the markers are used to mark the stability
"


MakeListPlotData::usage = "MakeListPlotData[xcoord, ycoord]
Create data for list plot where each point can be a different color using PlotStyle option"


GetBoundaryLineBiStable::usage = "GetBoundaryLineBiStable[rawdat]
Get the boundary line for area with bistability
rawdat: (Nested List) List of bifurcation parameters (form: {{0, 1}, {1, 3}}
each element of the List is a List of the value of first parameter and second parameter

return a Nested List, first List is the upper boundary and second List is the lower boundary
"


GetBoundaryLineSingle::usage = "GetBoundaryLineSingle[rawdat]
Get the boundary line for area with bistability
rawdat: (Nested List) List of bifurcation parameters (form: {{0, 1}, {1, 3}}
each element of the List is a List of the value of first parameter and second parameter

return a Nested List, first List is the upper boundary and second List is the lower boundary
"


DirectionParametricPlot::usage="DirectionParametricPlot[{\!\(\*SubscriptBox[
StyleBox[\"f\", \"TI\"], 
StyleBox[\"x\", \"TI\"]]\),\!\(\*SubscriptBox[
StyleBox[\"f\", \"TI\"], 
StyleBox[\"y\", \"TI\"]]\)},{\!\(\*
StyleBox[\"t\", \"TI\"]\),\!\(\*SubscriptBox[
StyleBox[\"t\", \"TI\"], 
StyleBox[\"min\", \"TI\"]]\),\!\(\*SubscriptBox[
StyleBox[\"t\", \"TI\"], 
StyleBox[\"max\", \"TI\"]]\)}]
develop by Dennis M Schneider https://resources.wolframcloud.com/FunctionRepository/resources/DirectionParametricPlot/
"


Begin["`Private`"]


MakeSystem[var_, t_, sys_]:=
Module[{varOft, sysOft},

On[Assert];
Assert[Head@t===Symbol];
Thread[(Head/@var)===ConstantArray[Symbol, Length@var]];

varOft = Through[var[t]];
sysOft = sys/.Thread[var-> varOft]; 
Thread[D[varOft, t] == sysOft]]


FollowRoot[system_,commonpars_,followPar_,range_,variables_, initialEq_]:=
Module[{nbvar, iEq, eq, pars, results, initValues, isResultsWithinRange, checkresults},

nbvar = Length@initialEq;
On[Assert];
Thread[(Head/@variables) === ConstantArray[Rule, nbvar]];
Thread[(Head/@initialEq)=== ConstantArray[Rule, nbvar]];

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


NSolvePositive[system_, commonpars_, followPar_, variables_, equisymbol_, precision_:MachinePrecision]:=
Module[
{eqAll, eqpos, nbpos, fparlist},

On[Assert];
Assert[Head[followPar]=== Rule];

eqAll = NSolve[Thread[(system/.commonpars/.followPar) == 0], variables, Reals, WorkingPrecision->precision];
eqpos = Select[eqAll, AllTrue[variables/.# , Positive]&]//Sort;
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
markerlist: list of marker that one wish to set, markerlist[[1]]: cycle, markerlist[[2]]: unstable, markerlist[[3]]: stable
*)
Module[
{eiv, anyZero, anyPos, allNeg, ps},
eiv = Eigenvalues[jacobmatrix/.parcommon/.parsfollow/.equilibrium];
anyZero = AnyTrue[Thread[-10^-10<=Re[eiv]<=10^-10], TrueQ];
anyPos = AnyTrue[Thread[Re[eiv]>10^-10], TrueQ];
allNeg = AllTrue[Thread[Re[eiv]<-10^-10], TrueQ];
Which[anyZero, markerlist[[1]], anyPos, markerlist[[2]], allNeg, markerlist[[3]]]
]


SingleStableColor[jacobmatrix_, parcommon_, parsfollow_, equilibrium_, colorlist_, opacity_:0.3, pointsize_:0.02]:=
(*Similar to the SingleStableMark function, but instead of returning the mark, this function returns the colors*)
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


ListStableMark[jacobmatrix_, parcommon_, parfollow_, range_,equiList_, markerCode_:{"@", "*", "."}, useColor_:False]:=
(* MapThread[SingleStableMark,
			{ConstantArray[jacobmatrix, Length[equiList]], 
			 ConstantArray[parcommon, Length[equiList]], 
			 Thread[parfollow -> range], 
			 equiList}
]*)
Module[{lenEqList, markers},
lenEqList = Length[equiList];
markers = ConstantArray[markerCode, lenEqList];
If[useColor,
	On[Assert];
	Head[markerCode[[1]]] === RGBColor//Assert;
	marklist = MapThread[SingleStableColor,
						{ConstantArray[jacobmatrix, lenEqList], 
						 ConstantArray[parcommon, lenEqList],
						 Thread[parfollow -> range], 
						 equiList,
						 markers}],
	marklist = MapThread[SingleStableMark,
						{ConstantArray[jacobmatrix, lenEqList], 
						 ConstantArray[parcommon, lenEqList],
						 Thread[parfollow -> range], 
						 equiList,
						 markers}]
	]
]


NSolveCodim2Positive[system_, commonpars_, bifurpar1_, bifurpar2_,bfparsName_, equisymbol_,variables_, precision_:MachinePrecision]:=
Module[
{eqAll, parspairVal, eqpos, nbpos, fparlist},
eqAll = NSolve[Thread[(system/.commonpars/.bifurpar1/.bifurpar2)==0], variables, Reals, WorkingPrecision->precision];
eqpos = Select[eqAll,And@@Thread[(variables/.# )> 0]&]//Sort;
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
	On[Assert];
	Head[markerCode[[1]]] === RGBColor//Assert;
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


MakeListPlotData[xcoord_, ycoord_]:={#}& /@Transpose[{xcoord, ycoord}]


GetBoundaryLineBiStable[rawdat_ ]:=Module[
{tallyDat, duplicateDat,upperBound, lowerBound},
tallyDat = Tally[rawdat];
duplicateDat = Select[tallyDat, #[[2]]==2&][[All, 1]];
lowerBound = Sequence@@{#[[1]]} &/@GatherBy[duplicateDat, First];
upperBound = Sequence@@{#[[-1]]} &/@GatherBy[duplicateDat, First];
{lowerBound, upperBound}
]


GetBoundaryLineSingle[rawdat_ ]:=Module[
{tallyDat, singleDat, duplicateDat,lowerBoundSingleDat, lowerBoundDuplicateDat, joinDat, tallyByXaxisParam},
tallyDat = Tally[rawdat];
singleDat = Select[tallyDat, #[[2]]==1&][[All, 1]];
duplicateDat = Select[tallyDat, #[[2]]==2&][[All, 1]];
lowerBoundSingleDat = Sequence@@{#[[1]]} &/@GatherBy[singleDat, First];
lowerBoundDuplicateDat = Sequence@@{#[[1]]} &/@GatherBy[duplicateDat, First];
joinDat = Join[lowerBoundSingleDat, lowerBoundDuplicateDat]//Sort;
tallyByXaxisParam = Tally[joinDat, #1[[1]]==#2[[1]]&];
Sort[tallyByXaxisParam][[All, 1]]
]


(* ::Input::Initialization:: *)
DirectionParametricPlot[fun_, {t_, t0_, t1_}, opts___?OptionQ] := 
	Module[{x,y,z,funl, fp, aheads, asize, anumber, color, pp, gr, vars, arrow, data, dt = N[t1 - t0], ps, psdata, dpts, plotfun, plotrange}, 
	funl = ReleaseHold[Hold[fun] /. {Integrate -> NIntegrate, Tooltip[a_]->a}]//Quiet;
	If[VectorQ[funl],funl = ReleaseHold[{Hold[funl]}]];
	If[Head[funl] === Piecewise, funl = {funl}];
	If[MatrixQ[First[funl]], funl = Flatten[funl,1]];
	If[Not[MatrixQ[funl/.t->t0]] & , Return[Message[dpp::improperinput]]];
        {aheads, anumber, color} = 
	    {"DrawArrowheads", "ArrowNumber", ColorFunction} /. {opts} /.  
               Options[DirectionParametricPlot];		
     If[anumber <= 0, aheads = False];
     {ps,asize}=AppendOptions[DirectionParametricPlot,{PlotStyle,"ArrowSize"}, Length[funl],opts]; 
asize = If[Head[#] === Symbol,{#},{.05#}] & /@ (Last/@asize);
     vars = {x,y};
     arrow = (Arrow[#]&); 

     If[Not[FreeQ[ps, RGBColor[_,_,_] | Hue[_] | GrayLevel[_]]], color =None];
     If[MatchQ[color, Automatic],  
                 color = Function[Flatten[{vars,t}]//Evaluate, Hue[.9 t]]];		
     aheads = 
	If[aheads,fp = D[funl, t];
          If[FreeQ[fp, Derivative, Infinity], 			
	If[color === None, psdata = 
 {RGBColor[0.368417, 0.506779, 0.709798], RGBColor[0.880722, 0.611041, 0.142051],
RGBColor[0.560181, 0.691569, 0.194885], RGBColor[0.922526, 0.385626, 0.209179], RGBColor[0.528488, 0.470624, 0.701351], RGBColor[0.772079, 0.431554, 0.102387], 
		RGBColor[0.363898, 0.618501, 0.782349], RGBColor[1, 0.75, 0], RGBColor[0.647624, 0.37816, 0.614037], RGBColor[0.571589, 0.586483, 0.], 
		RGBColor[0.915, 0.3325, 0.2125], RGBColor[0.40082222609, 0.52200666434, 0.85], 
		RGBColor[0.9728288904374106, 0.621644452187053, 0.07336199581899142], RGBColor[0.736782672705901, 0.358, 0.5030266573755369], 
		RGBColor[0.28026441037696703`, 0.715, 0.4292089322474965]
		};
	     ps = MapThread[Flatten[{#1,#2}]&,{psdata[[1;;Length[funl]]],ps}];
	     ps1 = DeleteCases[ ps, Thickness[_]|Tube[_], Infinity]; 
data = Table[arrow /@ Transpose[{funl, funl + .001 Map[Normalize[#] &, fp]}//Chop], {t, t0 + dt/anumber, t1, dt/(anumber+.0001)}];
data = Map[MapThread[Flatten[{Arrowheads[Last[Flatten[{#2}]]],#1,#3}]&,{ps1,asize,#}]&,data],
			
acolor= Table[color[Sequence@@Flatten[{vars,(t - t0)/dt}]],{t, t0 + dt/anumber, t1, dt/(anumber + .0001)}];
asizes=Map[Arrowheads,asize,{2}]; 
arrows=Table[arrow /@ Transpose[{funl, funl + .001 Map[Normalize[#] &, fp]}//Chop], 
										{t, t0 + dt/anumber, t1, dt/(anumber + .0001)}];
data=MapThread[{#1,#2}&,{acolor,(Transpose[Flatten/@{asizes,#}]&/@arrows)}]];
              Select[data, FreeQ[#, Indeterminate, Infinity] &, Infinity],
              Message[dpp::noderiv];{}],{}];

Show[	
	ParametricPlot[funl//Evaluate, {t, t0, t1}, ColorFunction -> color, PlotStyle -> ps,
				Sequence @@ FilterRules[{opts}, Options[ParametricPlot]] // Evaluate], 
	Graphics[aheads],  Sequence@@FilterRules[{opts},Options[Graphics]]//Evaluate,
   PlotRangePadding->Scaled[.04]]
	]

Options[DirectionParametricPlot]=
{"ArrowNumber"->15,"ArrowSize"->Large,ColorFunction->Automatic,"DrawArrowheads"->True,PlotStyle->Automatic};

dpp::noderiv="Mathematica was unable to find the derivative which is used to produce the arrowheads. Plotting will proceed without arrowheads. ";

dpp::improperinput="The input function is not a vector or a list of vectors. ";
AppendOptions[command_Symbol,options_,rest__]:=AppendOptions[command,{options},rest]
AppendOptions[command_Symbol,{options__},number_?(IntegerQ[#]&&NonNegative[#]&),opts___Rule]:=AppendOptions[command,{options},Table[number,{Length[{options}]}],opts]
AppendOptions[command_Symbol,{options__},numbers_?(VectorQ[#,(IntegerQ[#]&&NonNegative[#])&]&),opts___Rule]:=Module[{defaultstyles,userstyles,lengths,newstyles},defaultstyles={options}/.Options[command];
userstyles={options}/.{opts}/.Options[command];
userstyles=(If[ListQ[#],#,{#}]&/@userstyles)/.{}->{{}};
lengths=Length/@userstyles;
newstyles=MapThread[MapThread[Flatten[{#1,#2}]&,{Table[#1,{#2}],#3}]&,{defaultstyles,lengths,userstyles}];
newstyles=MapThread[Join[Flatten[Table[#1,{Floor[#2/#3]}],1],Take[#1,Mod[#2,#3]]]&,{newstyles,numbers,lengths}];
With[{temp=DeleteCases[newstyles,Automatic,Infinity]},If[Length[{options}]===1,First[temp],temp]]]



End[]


EndPackage[]


(* ::Input:: *)
(**)
