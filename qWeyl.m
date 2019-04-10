BeginPackage["qWeyl`",{"qmatrix`"}]

(*Finite Abelian Groups*)

randomEndo::usage = "";

endomorphismQ::usage = "";

automorphismQ::usage = "";

inverseAut::usage = "";

inverseAut2::usage = "";

reduceToBlocks::usage = "A handy tool to create an automorphism by hand; shows the reduced matrices so you can fine-tune one that is invertible.";

adjugate::usage = "";

adjointEndo::usage = "";

symp::usage = "The symplectic bracket, usually done via vec1 ~symp~ vec2";

(*qWeyl*)

qInit::usage = "This function initializes your Hilbert space. numSystems must be a natural and listOfFactors must be a list (even if it's just {3})";

testfcn::usage = "This is a simple function designed for testing purposes";

ranlist::usage = "a helper function useful for testing. ranlist[sys,num] generates num numbers from system sys";

\[Chi]::usage = "";

X::usage = "The 'shift' operator";

Z::usage = "The 'boost' operator";

pi::usage = "The pi operator";

piSparse::usage = "";

WSparse::usage = "This is Gross' w (Weyl) operator, implemented sparsely using piSparse";

ASparse::usage = "Phase Space Point Operator";

AngleBracket::usage = "";

lim::usage = "";

numSystems::usage = "this global variable stores the number of systems we are currently considering";

nlist::usage = "this is a list of the factors that make up the finite abelian group.
Each factor must divide the next";

Nmat::usage = "this is a global variable; it's the factors in the finite abelian group listed on the diagonal of an numSystems x numSystems matrix";

Nsymp::usage = "this is a global variable; it's the symplectic transformation matrix except with Nmat instead of Identity"

primelist::usage = "this is a global variable; lists the prime numbers that are factors in nlist";

factorlist::usage = "this is a global variable; lists the factors that go into nlist."

FourierU::usage= "";

SFourier::usage = "";

SFourierInv::usage = "";

phiFourier::usage = "";

SDilation::usage = "";

phiDilation::usage = "";

DilationU::usage = "";

DilationUSparse::usage = "";

phiChirp::usage = "";

vecMod::usage = "This computes the modulo of a list where each component is reduced mod it's corresponding element in nlist. 
	It also happens to reduce a matrix's rows mod nlist, which is convenient.";

sub::usage = "This is a shortcut to create a new symbol; sub[var,i] creates symbol vari for use in equations";

(*Stone-von Neumann Operators*)

wJ::usage = "";

thetaJ::usage = "goes inside momentum; Arg(wj^n-1)=thetaJ";

phiJ::usage = "goes inside position; Arg(w^j)=phiJ";

position::usage = "If given a number, returns the position observable on that system. If given a binary list, returns an observable on system that 
  is identity on systems corresponding to 0's and position on systems corresponding to 1's";

momentum::usage = "";

shiftU::usage = "";

boostV::usage = "";

subtractMean::usage = "X^0";

generateR::usage = "";

meanVector::usage = "m";

covarianceMatrix::usage = "B_S(X)";

commutatorMatrix::usage = "C_S(X)";

variance::usage = "D_S(X)";

ASparse::usage = "";

indicatorFcn::usage = "";

gaussianKet::usage = "";

charfcn::usage = "";

Harper::usage = "";

minUncerVal::usage = "";

minUncerVec::usage = "";

ABstate::usage = "";


(*Phase Space*)

lineInPhaseSpace::usage = "Gives the coordinates of a line through point in the lattice GxG, where G is 
the group corresponding to the subsystem syst. ";

q::usage = "A protected character designed to make it easy to get access to a subsystem; q[i] is the ith subsystem.";
qWeyl`helpMe := "This package implements the Weyl representation 
of the Heisenberg group H(G) for a finite abelian group G";

Begin["`Private`"]

testfcn[c_] := dimension[Subscript[q,c]]

ranlist[sys_, num_] := RandomInteger[{0, nlist[[sys]] - 1}, num]

qInit[numberOfSystems_?IntegerQ, listOfFactors_?ListQ] := ( 
  numSystems = numberOfSystems; nlist = listOfFactors; 
  Nmat = DiagonalMatrix[listOfFactors]; Nsymp = ArrayFlatten[{{ConstantArray[0,{numSystems,numSystems}],Inverse[Nmat]},{-Inverse[Nmat],ConstantArray[0,{numSystems,numSystems}]}}];
  primelist = FactorInteger[Last[nlist]][[;; , 1]];
  factorlist = 
 Table[{primelist[[i]], 
   Table[IntegerExponent[nlist[[j]], primelist[[i]]], {j, 
     numSystems}]}, {i, Length[primelist]}];
     Unprotect[q];
  Map[(Subscript[q, #] := Symbol["q" <> ToString[#]]) &, 
   Table[i, {i, numSystems}]]; 
   Map[(q[#] := Symbol["q" <> ToString[#]]) &, 
   Table[i, {i, numSystems}]]; 
   Protect[q];
  Map[setModeType[Subscript[q, #], {discrete, listOfFactors[[#]]}] &, 
   Table[i, {i, numSystems}]]; 
  setSystem[Table[Subscript[q, i], {i, numSystems}]])

Protect[q]

(*Finite Abelian Group Stuff*)

(*still to come: randomAuto,*)
randomEndo[] := 
 Table[If[r <= s, RandomInteger[{0, nlist[[r]]-1}], 
   nlist[[r]] /nlist[[s]] RandomInteger[{0, nlist[[s]]-1}]], {r, 
   numSystems}, {s, numSystems}]

reduceToBlocks[A_] := 
 Table[{Mod[
     Det[A[[(val1 = (Count[factorlist[[i, 2]], x_ /; x == 0] + 1)) ;;,
         val1 ;;]]], primelist[[i]]](* !=0 *), primelist[[i]], 
   vecMod[A, primelist[[i]]^factorlist[[i, 2]]] // MatrixForm}, {i, 
   Length[primelist]}]

endomorphismQ[A_] := 
 AllTrue[Flatten[
   Table[Mod[A[[i, j]], nlist[[i]]/nlist[[j]]] == 0, {j, 1, 
     numSystems - 1}, {i, j + 1, numSystems}]], TrueQ]

automorphismQ[A_] := 
 endomorphismQ[A] && 
  AllTrue[Table[
    Mod[Det[A[[(val1 = (Count[factorlist[[i, 2]], x_ /; x == 0] + 
              1)) ;;, val1 ;;]]], primelist[[i]]] != 0, {i, Length[primelist]}], 
   TrueQ]

inverseAut[A_] := 
 If[! automorphismQ[A], 
  Message["I don't think A is invertible, you should double check"], 
  Module[{vals, blocks}, 
   Do[vals[i] = Count[factorlist[[i, 2]], x_ /; x == 0] + 1, {i, 
     Length[primelist]}];
   blocks = 
    Table[vecMod[A, primelist[[i]]^factorlist[[i, 2]]], {i, 
      Length[primelist]}];
   Do[blocks[[i]][[vals[i] ;;, vals[i] ;;]] = 
     inverseHp[blocks[[i]][[vals[i] ;;, vals[i] ;;]], i], {i, 
     Length[primelist]}];
   Table[ChineseRemainder[blocks[[;; , i, j]], 
     primelist[[ ;; ]]^factorlist[[;; , 2, i]]], {i, numSystems}, {j, 
     numSystems}]]]

adjugate[m_] := 
 Map[Reverse, Minors[Transpose[m], Length[m] - 1], {0, 1}]*
  Table[(-1)^(i + j), {i, Length[m]}, {j, Length[m]}]

inverseHp[A_, primePos_] := 
 vecMod[PowerMod[Det[A], -1, primelist[[primePos]]^
    Last[factorlist[[primePos, 2]]]] adjugate[A], 
  primelist[[primePos]]^Select[factorlist[[primePos, 2]],# !=0 &]]

adjointEndo[A_] := Nmat.Transpose[A].Inverse[Nmat]

symp[vec1_, vec2_] := vec1.SFourier[Length[vec1]/2].vec2

(*Weyl Stuff*)

\[Chi][syst_,x_] := Exp[2 Pi I x /dimension[syst]]

\[Chi][syst_?IntegerQ,x_] := \[Chi][Subscript[q,syst],x]

(* \[Chi][x_] := Product[\[Chi][i,x[[i]]],{i,numSystems}]*)

X[syst_, shift_?IntegerQ] := 
 matrix[op[permutation, syst, 
   Table[Mod[i + shift, dimension[syst], 1], {i, 
     dimension[syst]}]]] 

X[system_?IntegerQ, shift_?IntegerQ] := 
 X[q[system], shift] 

Z[syst_, boost_?IntegerQ] := 
 matrix[DiagonalMatrix[
   Table[\[Chi][syst,Mod[i*boost, dimension[syst]]], {i, 0, 
     dimension[syst] - 1}]], {ket[syst], bra[syst]}]

Z[syst_?IntegerQ, boost_?IntegerQ] := Z[Subscript[q, syst], boost]

WSparse[m_, l_] := piSparse[{l}, {m}] \[Chi][1, -1/2 m l]

pi[lambda_] := Fold[diracMatrixProduct, 
 Table[Z[i, lambda[[numSystems + i]]] ** X[i, lambda[[i]]], {i, 
   numSystems}]]

pi[l_,m_] := pi[Join[l,m]]

piHelp[glist_, llist_] := Riffle[vecMod[glist + llist] + 1, glist + 1]

piSparse[l_, m_] := 
 matrix[SparseArray[
   Flatten[Table[
       piHelp[Table[g[i], {i, numSystems}], 
         l] -> \[LeftAngleBracket]Table[g[i], {i, numSystems}], 
          m\[RightAngleBracket] \[LeftAngleBracket]l, 
          m\[RightAngleBracket], ##] & @@ 
     Evaluate[Table[{g[i], 0, nlist[[i]] - 1}, {i, numSystems}]]], 
   Riffle[nlist, nlist]], 
  Riffle[Table[ket[Symbol["q" <> ToString[i]]], {i, numSystems}], 
   Table[bra[Symbol["q" <> ToString[i]]], {i, numSystems}]]]

piSparse[lambda_] := 
 piSparse[lambda[[;; numSystems]], lambda[[numSystems + 1 ;;]]]

AHelp[glist_, llist_] := 
 Riffle[vecMod[glist + 2 llist] + 1, vecMod[-glist] + 1]

ASparse[m_, l_] := 
 matrix[SparseArray[
   Flatten[Table[
       AHelp[Table[g[i], {i, numSystems}], 
         l] -> \[LeftAngleBracket]Table[g[i], {i, numSystems}], 
          2 m\[RightAngleBracket] \[LeftAngleBracket]2 l, 
          m\[RightAngleBracket], ##] & @@ 
     Evaluate[Table[{g[i], 0, nlist[[i]] - 1}, {i, numSystems}]]], 
   Riffle[nlist, nlist]], 
  Riffle[Table[ket[q[i]], {i, numSystems}], 
   Table[bra[q[i]], {i, numSystems}]]]

ASparse[lambda_] := 
 ASparse[lambda[[;; numSystems]], lambda[[numSystems + 1 ;;]]]

ASparse[m_?IntegerQ,l_?IntegerQ] := ASparse[{m},{l}]

AngleBracket[m_, k_] := Exp[2 Pi I m.Inverse[Nmat].k]

FourierUBlock[sys_?IntegerQ] := 
  matrix[1/Sqrt[nlist[[sys]]] Table[
     Conjugate[\[Chi][sys,i j]], {i, 0, nlist[[sys]] - 1}, {j, 0, nlist[[sys]] - 1}], {ket[
     Subscript[q, sys]], bra[Subscript[q, sys]]}];

(* FourierU2 = 
  matrix[1/Sqrt[p] Table[
     Conjugate[\[Chi][sys,i j]], {i, 0, p - 1}, {j, 0, p - 1}], {ket[
     Subscript[q, 2]], bra[Subscript[q, 2]]}]; *)

FourierU[systems___] := If[Length[{systems}]==0,
	Fold[diracMatrixProduct, Table[FourierUBlock[i], {i, numSystems}]], 
	Fold[diracMatrixProduct, Table[FourierUBlock[i], {i, {systems}}]] ]
 
SFourier[numSys_] := 
 ArrayFlatten[{{ConstantArray[0, {numSys, numSys}], 
    IdentityMatrix[numSys]}, {-IdentityMatrix[numSys], 
    ConstantArray[0, {numSys, numSys}]}}]

SFourierInv[numSys_] := Inverse[SFourier[numSys]]

phiFourier[l_, m_] := \[LeftAngleBracket]m, l\[RightAngleBracket]

phiFourier[lambda_] := 
 phiFourier[lambda[[;;numSystems]],lambda[[numSystems+1;;]]]

(* DilationU[r_, system_] := 
 matrix[op[permutation, Subscript[q, system], 
   Table[Mod[PowerMod[r, -1, p] k, p] + 1, {k, 0, p - 1}]]] *)

(* DilationU[A_] := 
 matrix[Table[  
   If[{k1, k2} == Mod[A.{l1, l2}, 4], 1, 0], {l1, 0, 3}, {k1, 0, 
    3}, {l2, 0, 3}, {k2, 0, 3}], {ket[q1], bra[q1], ket[q2], bra[q2]}] *)

sub[symbol_, index_] := Symbol[ToString[symbol] <> ToString[index]]

vecMod[list_] := Table[Mod[list[[i]], nlist[[i]]], {i, 1, numSystems}]

vecMod[list_,mods_] := Table[Mod[list[[i]], mods[[i]]], {i, 1, Length[mods]}]


DilationU[A_] := 
 matrix[Table[
     If[Table[sub[vark, i], {i, numSystems}] == 
       vecMod[A.Table[sub[varl, i], {i, numSystems}]], 1, 0], ##] & @@
    Flatten[
    Table[{{sub[varl, i], 0, nlist[[i]] - 1}, {sub[vark, i], 0, 
       nlist[[i]] - 1}}
     , {i, 1, numSystems}], 1], 
  Flatten[Table[{ket[Symbol["q" <> ToString[i]]], 
     bra[Symbol["q" <> ToString[i]]]}, {i, numSystems}]]]

sendsTo[autInv_, jlist_] := 
 Riffle[vecMod[autInv.jlist] + 1, jlist + 1]

DilationUSparse[A_] := 
 matrix[SparseArray@
   Flatten[Table[
       sendsTo[A, Table[j[i], {i, 1, numSystems}]] -> 
        1, ##] & @@ 
     Evaluate[Table[{j[i], 0, nlist[[i]] - 1}, {i, 1, numSystems}]]], 
  Riffle[Table[ket[Symbol["q" <> ToString[i]]], {i, numSystems}], 
   Table[bra[Symbol["q" <> ToString[i]]], {i, numSystems}]]]

SDilation[A_] := 
 ArrayFlatten[{{A, 
    ConstantArray[0, Dimensions[A]]}, {ConstantArray[0, 
     Dimensions[A]], inverseAut@adjointEndo[A]}}]

phiDilation[l_, m_] := 1

phiDilation[lambda_] := phiDilation[lambda[[;;numSystems]], lambda[[numSystems+1;;]]]

phiChirp[c_, k_] := 
 Exp[Pi I k.(IdentityMatrix[numSystems] + Inverse[Nmat]).c.(IdentityMatrix[numSystems] + Nmat).k]

ChirpU12[c_] := 
 matrix[Table[
   If[l1 == k1 && l2 == k2, phiChirp[c, {l1, l2}], 0], {l1, 0, 
    p - 1}, {k1, 0, p - 1}, {l2, 0, p - 1}, {k2, 0, p - 1}], {ket[q1],
    bra[q1], ket[q2], bra[q2]}]

ChirpU[c_] := 
 matrix[SparseArray[
   Flatten[Table[
       Riffle[Table[j[i], {i, numSystems}] + 1, 
         Table[j[i], {i, numSystems}] + 1] -> 
        phiChirp[c, Table[j[i], {i, numSystems}]], ##] & @@ 
     Evaluate[Table[{j[i], 0, nlist[[i]] - 1}, {i, 1, numSystems}]]]],
   Riffle[Table[ket[Symbol["q" <> ToString[i]]], {i, numSystems}], 
   Table[bra[Symbol["q" <> ToString[i]]], {i, numSystems}]]]

SChirp[C_?ListQ] := 
 ArrayFlatten[{{IdentityMatrix[numSystems], 
    ConstantArray[0, {numSystems, numSystems}]}, {C, 
    IdentityMatrix[numSystems]}}]

InverseMod[val_, p_] := If[GCD[val, p] == 1, PowerMod[val, -1, p], 0]
SetAttributes[InverseMod, Listable];

(*Stone-von Neumann Operators*)

SmallCircle[X_, Y_] := 1/2 (X ** Y + Y ** X)

subtractMean[state_, opt_] := 
 opt - trace[state ** opt] identityMatrix[opt[[2]]]

(*goes inside momentum, Arg(wj^n-1)=thetaJ*)
thetaJ[syst_,j_] := Arg[Exp[I 2 Pi j (nlist[[syst]] - 1)/nlist[[syst]]]]

(*goes inside position, Arg(w^j)=phiJ*)
phiJ[syst_, j_] := Arg[Exp[I 2 Pi j/nlist[[syst]]]]

wJ[syst_,j_] := Exp[I 2 Pi j / nlist[[syst]]]

momentum[syst_?IntegerQ] := 
 matrix[Chop@N@Table[
   Sum[thetaJ[syst,j] wJ[syst,j]^(i - k)/nlist[[syst]], {j, 0, 
     nlist[[syst]] - 1}], {i, 0, nlist[[syst]] - 1}, {k, 0, 
    nlist[[syst]] - 1}], {ket[Subscript[q, syst]], bra[Subscript[q, syst]]}]

position[syst_?IntegerQ] := 
 matrix[DiagonalMatrix[
   Table[phiJ[syst,j], {j, 0, nlist[[syst]] - 1}]], {ket[Subscript[q, syst]], 
   bra[Subscript[q, syst]]}]

momentum[systList_?ListQ] := 
 Fold[diracMatrixProduct, 
  Table[If[systList[[j]] == 1, momentum[j], 
    identityMatrix[Subscript[q, j]]], {j, numSystems}]]

position[systList_?ListQ] := 
 Fold[diracMatrixProduct, 
  Table[If[systList[[i]] == 1, position[i], 
    identityMatrix[Subscript[q, i]]], {i, numSystems}]]

shiftU[shifts_] := 
 Fold[diracMatrixProduct, Table[X[i, shifts[[i]]], {i, numSystems}]]

boostV[boosts_] := 
 Fold[diracMatrixProduct, Table[Z[i, boosts[[i]]], {i, numSystems}]]

generateR[] := 
 Join[Table[position[UnitVector[numSystems, i]], {i, numSystems}], 
  Table[momentum[UnitVector[numSystems, i]], {i, numSystems}]]

meanVector[state_] := 
 Module[{R = generateR[]}, 
  Table[trace[state ** R[[i]]], {i, 2 numSystems}]]

covarianceMatrix[state_, obs_] := 
 Table[trace[
   state ** (subtractMean[state,obs[[j]]]\[SmallCircle]subtractMean[
       state,obs[[k]]])], {j, 2numSystems}, {k, 2numSystems}]

variance[state_, obs_] := trace[state ** subtractMean[state, obs]^2]

Unprotect[Power]
Power[state_?matrixQ, number_?IntegerQ] := 
 If[number!=0,Fold[diracMatrixProduct, Table[state, {i, number}]],identityMatrix[state[[2]]]]
Protect[Power]

commutatorMatrix[state_, obs_] := 
 I Table[trace[state ** commutator[obs[[j]], obs[[k]]]], {j, 
    2numSystems}, {k, 2numSystems}]

(*Phase Space*)

lineInPhaseSpace[syst_, slope_, point_: {0, 0}] := 
 Table[{Mod[slope q + point[[1]], nlist[[syst]]], 
   Mod[q + point[[2]], nlist[[syst]]]}, {q, 0, nlist[[syst]] - 1}]

indicatorFcn[set_, element_] := If[MemberQ[set, element], 1, 0]

gaussianKet[syst_?IntegerQ, covar_, mean_] := 
 matrix[1/Sqrt[
    nlist[[syst]]] Table[\[Chi][syst, j covar  j + mean j], {j, 0, 
     nlist[[syst]] - 1}], {ket[sub[q, syst]]}]

charfcn[state_] := 
 1/nlist[[1]] Table[
   trace[state ** hc[WSparse[i, j]]], {i, 0, nlist[[1]] - 1}, {j, 0, 
    nlist[[1]] - 1}]

Harper[] := 
 Module[{harp = (Z[1, 1] + hc[Z[1, 1]] + X[1, 1] + hc[X[1, 1]])/4}, 
  matrix[Chop@N[harp[[1]]], harp[[2]]]]

minUncerVal[] := normalizedEigensystem[Harper[]][[1, 2]]

minUncerVec[check_: 1] := 
 Module[{gam = normalizedEigensystem[Harper[]][[1, 1]]}, (If[
    check == 1, Print[gam],]; matrix[Abs[gam[[1]]], gam[[2]]])]

 ABstate[alpha_, beta_, minUncerVect_: minUncerVec[0]] := 
 Z[1, alpha] ** X[1, beta] ** minUncerVect


(* Useful Code Snippets that I want to keep *)


(* Here is the calculations for characteristic functions for G-Gaussian kets. *)
charfcn2[theta_, m_] := 
 1/nlist[[1]] Table[(-1)^(i j)
     indicatorFcn[lineInPhaseSpace[1, 2 theta], {i, j}] \[Chi][1, 
     m j], {i, 0, nlist[[1]] - 1}, {j, 0, nlist[[1]] - 1}]

charfcn3[theta_, x_] := 
 1/nlist[[1]] Table[(-1)^(i j)
     indicatorFcn[lineInPhaseSpace[1, 2 theta], {i, j}] \[Chi][
     1, ( x1 i + x2 j)], {i, 0, nlist[[1]] - 1}, {j, 0, 
    nlist[[1]] - 1}]


End[]

EndPackage[]