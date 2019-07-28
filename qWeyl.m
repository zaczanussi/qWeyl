BeginPackage["qWeyl`",{"qmatrix`"}]

(*Finite Abelian Groups*)

randomEndo::usage = "";

endomorphismQ::usage = "";

automorphismQ::usage = "";

symplecticQ::usage = "";

generateSymplectic::usage = "Computes a symplectic matrix like in Kaiblinger and Neuhauser.";

selfAdjointEndoQ::usage = "Returns True if endo is a self-adjoint endomorphism of G";

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

DilationU::usage = "Don't use this! use DilationUSparse instead. This should probably be deleted soon.";

DilationUSparse::usage = "This is the proper Dilation. A must be an automorphism of the group.";

phiChirp::usage = "";

ChirpU::usage = "C must be a self-adjoint endomorphism of the group. ";

SChirp::usage = "";

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

annihilation::usage= "";

creation::usage= "";

subtractMean::usage = "X^0";

generateR::usage = "";

meanVector::usage = "m";

covarianceMatrix::usage = "B_S(X)";

commutatorMatrix::usage = "C_S(X)";

variance::usage = "D_S(X)";

ASparse::usage = "";

indicatorFcn::usage = "";

sympWithSet::usage = "";

phaseSpace::usage = "The parameter matrix determines whether the phase space is returned as a plane (if true) or as a list of coordinates (if false).";

sympComplement::usage = "";

isotropicQ::usage = "";

coisotropicQ::usage = "";

lagrangianQ::usage = "";

subsetState::usage = "";

im::usage = "Use to get the imaginary part of a matrix object";

re::usage = "Use to get the real part of a matrix object";

gaussianKet::usage = "";

charfcn::usage = "";

invCharfcn::usage = "";

wigner::usage = "";

Harper::usage = "Generates the Harper operator H_theta as given in Hashimoto2019. Theta is Pi / 4 by default, but can be
  any value in the range [0,Pi/2].";

minUncerVal::usage = "";

minUncerVec::usage = "";

ABstate::usage = "";

scalarMultipleQ::usage = "If the two vectors are scalar multiples of each other, return the
  scalar. If not, prints a message informing you of this. Put the larger vector first to get a value greater than 1.";

vectorSolve::usage = "Wraps up LinearSolve to be more convenient for kets; put b in the first component, and the vector
  that makes up m in a list in the second component.";

(*Phase Space*)

lineInPhaseSpace::usage = "Gives the coordinates of a line through point in the lattice GxG, where G is 
the group corresponding to the subsystem syst. Right now, the line is of the form (theta q , q), q in Z_d. 
Can set theta = Infinity to get a horizontal line.";

q::usage = "A protected character designed to make it easy to get access to a subsystem; q[i] is the ith subsystem.";
qWeyl`helpMe := "This package implements the Weyl representation 
of the Heisenberg group H(G) for a finite abelian group G";

(*Hermite-Gaussians*)

generateG::usage = "";

generateHarperFunctions::usage = "";

generateHarperValues::usage = "";

unitKet::usage = "Gives the standard basis, with 0 begin the first element and nlist[[1]] - 1 being the last. Be careful;
  currently unitKet and vjKet have a different ordering, so that F|ej> = |v{n-j}> but F|vj> = |ej>.";

vjKet::usage = "Gives the Fourier basis, with 0 being the first element. This might seem like it is ordered backward, but it is like this to 
  be consistent with the X and Z matrices.";

(*Plotting*)

reorderKet::usage = "";

ketPlot::usage = "";

complexKetPlot::usage = "";

ketLinePlot::usage = "";

reorderCharfcn::usage = "";

blockify::usage = "";

charfcnPlot::usage = "";

randomProbDistro::usage = "";

hermiteGaussian::usage = "";

kravchukPoly::usage = "Here n is the dimension, s is a natural number in {0,...,n-1}, and z is between 0 and 2l, with n = 2l + 1";

kravchukFunc::usage = "Here n is the dimension, s is a natural number in {0,...,n-1}, and z is between -l and l, with n = 2l + 1";

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

symplecticQ[mat_] := 
 Block[{dim = Dimensions[mat][[1]], A, B, C, D}, 
  A = mat[[;; dim/2, ;; dim/2]]; B = mat[[;; dim/2, dim/2 + 1 ;;]]; 
  C = mat[[dim/2 + 1 ;;, ;; dim/2]]; 
  D = mat[[dim/2 + 1 ;;, dim/2 + 1 ;;]];
  If[vecMod[A.adjointEndo[B]] == vecMod[B.adjointEndo[A]] && 
    vecMod[C.adjointEndo[D]] == vecMod[D.adjointEndo[C]] && 
    vecMod[A.adjointEndo[D] - B.adjointEndo[C]] == 
     IdentityMatrix[dim/2], True, False]]

generateSymplectic[sa_, aut_, B_, theta_] := 
 Block[{mat = 
    SChirp[sa.inverseAut[aut]].SDilation[aut].SFourierInv[
      2].SChirp[-inverseAut[aut].B].SFourier[2].SChirp[-theta]}, 
  ArrayFlatten[{{vecMod[mat[[;; 2, ;; 2]]], 
     vecMod[mat[[;; 2, 3 ;;]]]}, {vecMod[mat[[3 ;;, ;; 2]]], 
     vecMod[mat[[3 ;;, 3 ;;]]]}}]]

selfAdjointEndoQ[endo_] := endo == vecMod[adjointEndo[endo]]

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

(*I have added a factor of -1 to this!*)
WSparse[m_, l_] := piSparse[{l}, {m}] \[Chi][1, -1/2 m l] (-1)^(l m)

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
 1/nlist[[1]] matrix[SparseArray[
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

unitKet[num_] := matrix[UnitVector[nlist[[1]], Mod[num,nlist[[1]]]+1], {ket[q[1]]}]

vjKet[num_] := 
 1/Sqrt[nlist[[1]]] matrix[
   Table[wJ[1, num n], {n, 0, nlist[[1]] - 1}], {ket[q[1]]}]

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

annihilation[syst_: 1] := 1/Sqrt[2] (position[syst] + I momentum[syst])

creation[syst_: 1] := 1/Sqrt[2] (position[syst] - I momentum[syst])

number[syst_: 1] := creation[syst] ** annihilation[syst]

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
 If[slope == Infinity, 
  Table[{Mod[ q + point[[1]], nlist[[syst]]], 
    Mod[point[[2]], nlist[[syst]]]}, {q, 0, nlist[[syst]] - 1}], 
  Table[{Mod[slope q + point[[1]], nlist[[syst]]], 
    Mod[q + point[[2]], nlist[[syst]]]}, {q, 0, nlist[[syst]] - 1}]]

indicatorFcn[set_, element_] := If[MemberQ[set, element], 1, 0]

sympWithSet[subset1_, subset2_] := 
 AllTrue[Flatten[Outer[symp, subset1, subset2, 1]], Mod[#,nlist[[1]]] == 0 &]

phaseSpace[matrix_:False] := 
 Flatten[Outer[List, Range[0, nlist[[1]] - 1], 
   Range[0, nlist[[1]] - 1]], If[matrix==True,0,1]]

sympComplement[set_] := 
 Select[phaseSpace[], sympWithSet[set, {#}] == True &]

isotropicQ[set_] := SubsetQ[sympComplement[set], set]

coisotropicQ[set_] := SubsetQ[set, sympComplement[set]]

lagrangianQ[set_] := isotropicQ[set] && coisotropicQ[set]

gaussianKet[syst_?IntegerQ, covar_, mean_] := 
 matrix[1/Sqrt[
    nlist[[syst]]] Table[\[Chi][syst, j covar  j + mean j], {j, 0, 
     nlist[[syst]] - 1}], {ket[sub[q, syst]]}]

(* charfcn[state_] := 
 1/nlist[[1]] Table[
   trace[state ** WSparse[i, j]], {i, 0, nlist[[1]] - 1}, {j, 0, 
    nlist[[1]] - 1}]
 *)
charfcn[state_, i_, j_] := 1/nlist[[1]] trace[state ** WSparse[i, j]]

charfcn[state_] := 
 Table[charfcn[state, i, j], {i, 0, nlist[[1]] - 1}, {j, 0, 
   nlist[[1]] - 1}]

charfcnGross[state_] := 
 1/nlist[[1]] Table[
   trace[state ** hc[WSparse[i, j]]], {i, 0, nlist[[1]] - 1}, {j, 0, 
    nlist[[1]] - 1}]

invCharfcn[charfcn_?MatrixQ] := 
 Sum[hc[WSparse[z[[1]], z[[2]]]] charfcn[[z[[1]] + 1, 
     z[[2]] + 1]], {z, phaseSpace[]}]

wigner[state_] := 
 1/nlist[[1]] Table[
   trace[state ** ASparse[i, j]], {i, 0, nlist[[1]] - 1}, {j, 0, 
    nlist[[1]] - 1}]

Harper[sys_Integer: 1, theta_: Pi/4] := 
 Block[{harp = (Sin[theta] (Z[sys, 1] + hc[Z[sys, 1]]) + 
       Cos[theta] (X[sys, 1] + hc[X[sys, 1]]))/2},
  matrix[Chop@N[harp[[1]]], harp[[2]]]]

minUncerVal[] := 
 Module[{sys}, sys := normalizedEigensystem[Harper[]]; 
  sys[[Ordering[sys[[;; , 2]], nlist[[1]], Greater][[1]], 2]]]

minUncerVec[theta_: Pi/4(* , check_: 0 *)] := 
 Module[{gam, sys}, (sys := normalizedEigensystem[Harper[theta]]; 
   gam := sys[[Ordering[sys[[;; , 2]], nlist[[1]], Greater][[1]], 1]];
    (* If[check == 1, Print[gam],]; *) matrix[Abs[gam[[1]]], gam[[2]]])]

ABstate[alpha_, beta_, theta_: Pi/4] := X[1, alpha] ** Z[1, beta] ** minUncerVec[theta] (*I switched Z and X, July 15 19*)

subsetState[line_] := 
 1/Length[line] Sum[
   ABstate[line[[i, 1]], line[[i, 2]]] ** 
    hc[ABstate[line[[i, 1]], line[[i, 2]]]], {i, 1, Length[line]}]

subsetState[line_, coeffs_] := 
 Block[{}, 
  If[Total[coeffs] != 1, 
   Print["Sum of coefficients is not equal to 1;"],]; 
  If[Length[line] != Length[coeffs], 
   Print["The number of elements in line is different than the number \
of coefficients"],]; 
  Sum[coeffs[[i]] ABstate[line[[i, 1]], line[[i, 2]]] ** 
     hc[ABstate[line[[i, 1]], line[[i, 2]]]], {i, 1, Length[line]}]]

re[ket_] := matrix[Re[ket[[1]]], ket[[2]]]

im[ket_] := matrix[Im[ket[[1]]], ket[[2]]]

scalarMultipleQ[vec1_?ListQ, vec2_?ListQ] := 
 Which[VectorAngle[vec1, vec2] == 0, Norm[vec1]/Norm[vec2], 
 VectorAngle[vec1, -vec2] == 0, -Norm[vec1]/Norm[vec2], 
 True,Print["The vectors are not multiples of each other"]]

scalarMultipleQ[vec1_?matrixQ, vec2_?matrixQ] := 
 scalarMultipleQ[vec1[[1]], vec2[[1]]]

vectorSolve[b_, vecs_?ListQ] := 
 LinearSolve[Table[vecs[[num, 1]], {num, 1, Length[vecs]}], b]
vectorSolve[b_?matrixQ, vecs_?ListQ] := 
 LinearSolve[Table[vecs[[num, 1]], {num, 1, Length[vecs]}], b[[1]]]

(*Hermite-Gaussians*)

generateG[] := 
 matrix[SparseArray[
    Join[Table[{x + 1, Mod[x + 1, nlist[[1]]] + 1} -> 1, {x, 0, 
       nlist[[1]] - 1}], 
     Table[{x + 1, Mod[x - 1, nlist[[1]]] + 1} -> 1, {x, 0, 
       nlist[[1]] - 1}], 
     Table[{x + 1, x + 1} -> 2 Cos[2 Pi/nlist[[1]] x], {x, 0, 
       nlist[[1]] - 1}]]] // N, {ket[q[1]], bra[q[1]]}]

generateHarperFunctions[sys_Integer:1,theta_: Pi/4] := 
 Block[{esyst = normalizedEigensystem[Harper[sys,theta]] // Chop, pos, 
   neg}, pos = 
   Sort[Select[esyst, #[[1, 1, 1]] != 0 &], #1[[2]] > #2[[2]] &];
  neg = Sort[
    Select[esyst, #[[1, 1, 1]] == 0 &], #1[[2]] > #2[[2]] &];
  DeleteDuplicates[
   Riffle[ReplacePart[pos, 
      Table[{x + 1, 1} -> 
        If[pos[[x + 1, 1, 1, Ceiling[(Length[esyst])/2]]] < 
          0, -pos[[x + 1, 1]], pos[[x + 1, 1]]], {x, 0, 
        Length[pos] - 1}]], 
     ReplacePart[neg, 
      Table[{x + 1, 1} -> 
        If[neg[[x + 1, 1, 1, Ceiling[(Length[esyst])/2]]] < 
          0, -neg[[x + 1, 1]], neg[[x + 1, 1]]], {x, 0, 
        Length[neg] - 1}]]][[;; , 1]]]]

generateHarperValues[sys_Integer:1,theta_:Pi / 4] := 
 Block[{esyst = normalizedEigensystem[Harper[sys,theta]] // Chop}, 
  pos = Sort[
    Select[esyst, #[[1, 1, 1]] != 0 &], #1[[2]] > #2[[2]] &];
  neg = Sort[Select[esyst, #[[1, 1, 1]] == 0 &], #1[[2]] > #2[[2]] &];
   DeleteDuplicates[
   Riffle[ReplacePart[pos, 
      Table[{x + 1, 1} -> 
        If[pos[[x + 1, 1, 1, 1]] < 0, -pos[[x + 1, 1]], 
         pos[[x + 1, 1]]], {x, 0, Length[pos] - 1}]], 
     ReplacePart[neg, 
      Table[{x + 1, 1} -> 
        If[neg[[x + 1, 1, 1, 2]] < 0, -neg[[x + 1, 1]], 
         neg[[x + 1, 1]]], {x, 0, Length[neg] - 1}]]][[;; , 2]]]]

hermiteGaussian[k_, x_] := 
 1/(Surd[Pi, 4] Sqrt[2^k Factorial[k] Exp[x^2]]) HermiteH[k, x]

kravchukPoly[s_, n_, z_?IntegerQ] := 
 Sum[(-2)^k Binomial[z, k] Binomial[s, k] Binomial[n - 1, k]^-1, {k, 
   0, Min[s, z]}]

kravchukFunc[s_, n_, z_] := (-1)^s/2^
  Floor[n/2] Sqrt[
   Binomial[n - 1, s] Binomial[n - 1, z + Floor[n/2]]] kravchukPoly[s,
    n, z + Floor[n/2]]

kravKet[s_, n_] := 
 matrix[kravchukFunc[s, n, #] & /@ 
   Join[Range[0, Ceiling[n/2] - 1], Range[-Floor[n/2], -1]], {ket[
    q[1]]}]

xBasis[s_, n_, z_] := (-1)^s kravchukFunc[Floor[n/2] - s, n, z]

xKet[s_, n_] := 
 matrix[xBasis[s, n, #] & /@ 
   Join[Range[0, Ceiling[n/2] - 1], Range[-Floor[n/2], -1]], {ket[
    q[1]]}]

(*Plotting*)

reorderKet[ket_] := 
Join[ket[[1, Ceiling[Length[ket[[1]]]/2] + 1 ;;]], 
 ket[[1, ;; Ceiling[Length[ket[[1]]]/2]]]]

reorderCharfcn[charfcn_] := 
 ArrayFlatten[{{charfcn[[Ceiling[nlist[[1]]/2] + 1 ;;, 
      Ceiling[nlist[[1]]/2] + 1 ;;]], 
    charfcn[[Ceiling[nlist[[1]]/2] + 1 ;;, ;; 
       Ceiling[nlist[[1]]/2]]]}, {charfcn[[;; Ceiling[nlist[[1]]/2], 
      Ceiling[nlist[[1]]/2] + 1 ;;]], 
    charfcn[[;; Ceiling[nlist[[1]]/2], ;; Ceiling[nlist[[1]]/2]]]}}]

blockify[mat_] := 
 ArrayFlatten[{{mat[[;; Floor[nlist[[1]]/2] + 1, ;; 
       Floor[nlist[[1]]/2] + 1]], 
    Reverse[mat[[;; Floor[nlist[[1]]/2] + 1, 
       2 ;; Floor[nlist[[1]]/2] + 1]], 2]}, {Reverse[
     mat[[2 ;; Floor[nlist[[1]]/2] + 1, ;; Floor[nlist[[1]]/2] + 1]]],
     Reverse[
     mat[[2 ;; Floor[nlist[[1]]/2] + 1, 
       2 ;; Floor[nlist[[1]]/2] + 1]], {1, 2}]}}]

ketPlot[ket_, offset_: False, squared_: False] := 
 Block[{ketList = ket[[1]], exponent = If[squared == True, 2, 1], 
   off = If[offset == True, 1, 0]}, 
  ListPlot[Join[ketList[[Ceiling[Length[ket[[1]]]/2] + 1 ;;]], 
     ketList[[;; Ceiling[Length[ket[[1]]]/2]]]]^exponent, 
   DataRange -> {-Floor[Length[ket[[1]]]/2] + 0.5 off, 
     Ceiling[Length[ket[[1]]]/2] - 1 + 0.5 off}, Filling -> Axis,PlotRange -> Full]]

ketPlot[kets_?ListQ, squared_: False,legends_:False] := 
 Block[{ketList = kets[[;; , 1]], 
   exponent = If[squared == True, 2, 1],leg = If[legends === False, Range[1,Length[kets]],legends]}, 
  ListPlot[Table[
     Join[ketList[[i]][[Ceiling[nlist[[1]]/2] + 1 ;;]], 
      ketList[[i]][[;; Ceiling[nlist[[1]]/2]]]], {i, 1, 
      Length[kets]}]^exponent, 
   DataRange -> {-Floor[nlist[[1]]/2], Ceiling[nlist[[1]]/2] - 1}, 
   Filling -> Axis, PlotLegends -> {leg},PlotRange -> Full]]

complexKetPlot[ket_, offset_: False, squared_: False] := 
 Block[{ketList = reorderKet[ket], 
   exponent = If[squared == True, 2, 1], 
   off = If[offset == True, 1, 0]}, 
  ListPlot[{Re[ketList], Im[ketList]}^exponent, 
   DataRange -> {-Floor[Length[ket[[1]]]/2] + 0.5 off, 
     Ceiling[Length[ket[[1]]]/2] - 1 + 0.5 off}, Filling -> Axis, 
   PlotRange -> Full, PlotLegends -> {Re, Im}]]

ketLinePlot[ket_, squared_: False] := 
 Block[{ketList = ket[[1]], exponent = If[squared == True, 2, 1]}, 
  ListLinePlot[
   Join[ketList[[Ceiling[nlist[[1]]/2] + 1 ;;]], 
    ketList[[;; Ceiling[nlist[[1]]/2]]]]^exponent, 
   DataRange -> {-Floor[nlist[[1]]/2], Ceiling[nlist[[1]]/2] - 1}, 
   Filling -> Axis]]

ketLinePlot[kets_?ListQ, squared_: False] := 
 Block[{ketList = kets[[;; , 1]], 
   exponent = If[squared == True, 2, 1]}, 
  ListLinePlot[Table[
     Join[ketList[[i]][[Ceiling[nlist[[1]]/2] + 1 ;;]], 
      ketList[[i]][[;; Ceiling[nlist[[1]]/2]]]], {i, 1, 
      Length[kets]}]^exponent, 
   DataRange -> {-Floor[nlist[[1]]/2], Ceiling[nlist[[1]]/2] - 1}, 
   Filling -> Axis, PlotLegends -> {Range[1, Length[kets]]}]]

charfcnPlot[charfcn_] := 
 Block[{dim = Dimensions[charfcn][[1]], 
   off = If[Mod[Dimensions[charfcn][[1]], 2] == 0, .5, 0]}, 
   ListPlot3D[reorderCharfcn[charfcn], PlotRange -> Full, 
    DataRange -> {{-dim/2 - off, dim/2 - off}, {-dim/2 - off, 
       dim/2 - off}}]]

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

randomProbDistro[length_] := 
 Normalize[Table[RandomReal[], {i, length}], Total]


End[]

EndPackage[]