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

(*qWeyl*)

qInit::usage = "This function initializes your Hilbert space. numSystems must be a natural and listOfFactors must be a list (even if it's just {3})";

testfcn::usage = "This is a simple function designed for testing purposes";

ranlist::usage = "a helper function useful for testing. ranlist[sys,num] generates num numbers from system sys";

\[Chi]::usage = "";

X::usage = "The 'shift' operator";

Z::usage = "The 'boost' operator";

pi::usage = "The pi operator"

AngleBracket::usage = "";

lim::usage = "";

numSystems::usage = "this global variable stores the number of systems we are currently considering";

nlist::usage = "this is a list of the factors that make up the finite abelian group.
Each factor must divide the next";

Nmat::usage = "this is a global variable; it's the factors in the finite abelian group listed on the diagonal of an numSystems x numSystems matrix";

Nsymp::usage = "this is a global variable; it's the symplectic transformation matrix except with Nmat instead of Identity"

primelist::usage = "this is a global variable; lists the prime numbers that are factors in nlist";

factorlist::usage = "this is a global variable; lists the factors that go into nlist."

FourierUBlock::usage = "this should only be temporary";

FourierU::usage= "";

SFourier::usage = "";

SFourierInv::usage = "";

phiFourier::usage = "";

SDilation::usage = "";

phiDilation::usage = "";

DilationU::usage = "";

vecMod::usage = "This computes the modulo of a list where each component is reduced mod it's corresponding element in nlist. 
	It also happens to reduce a matrix's rows mod nlist, which is convenient.";

sub::usage = "This is a shortcut to create a new symbol; sub[var,i] creates symbol vari for use in equations";

qWeyl`helpMe := "This package implements the Weyl representation 
of the Heisenberg group H(G) for a finite abelian group G";

Begin["`Private`"]

testfcn[c_] := dimension[Subscript[q,c]]
(* p = 4;
 *)

ranlist[sys_, num_] := RandomInteger[{0, nlist[[sys]] - 1}, num]

qInit[numberOfSystems_?IntegerQ, listOfFactors_?ListQ] := ( 
  numSystems = numberOfSystems; nlist = listOfFactors; 
  Nmat = DiagonalMatrix[listOfFactors]; Nsymp = ArrayFlatten[{{ConstantArray[0,{numSystems,numSystems}],Inverse[Nmat]},{-Inverse[Nmat],ConstantArray[0,{numSystems,numSystems}]}}];
  primelist = FactorInteger[Last[nlist]][[;; , 1]];
  factorlist = 
 Table[{primelist[[i]], 
   Table[IntegerExponent[nlist[[j]], primelist[[i]]], {j, 
     numSystems}]}, {i, Length[primelist]}];
  Map[(Subscript[q, #] := Symbol["q" <> ToString[#]]) &, 
   Table[i, {i, numSystems}]]; 
  Map[setModeType[Subscript[q, #], {discrete, listOfFactors[[#]]}] &, 
   Table[i, {i, numSystems}]]; 
  setSystem[Table[Subscript[q, i], {i, numSystems}]])

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


(*Weyl Stuff*)

\[Chi][syst_,x_] := N[Exp[2 Pi I x /dimension[syst]]]

\[Chi][syst_?IntegerQ,x_] := N[Exp[2 Pi I x /dimension[Subscript[q, syst]]]]

\[Chi][x_] := Product[\[Chi][i,x[[i]]],{i,numSystems}]

X[syst_, shift_?IntegerQ] := 
 matrix[op[permutation, syst, 
   Table[Mod[i + shift, dimension[syst], 1], {i, 
     dimension[syst]}]]] 

X[system_?IntegerQ, shift_?IntegerQ] := 
 X[Subscript[q, system], shift] 

Z[syst_, boost_?IntegerQ] := 
 matrix[DiagonalMatrix[
   Table[\[Chi][syst,Mod[i*boost, dimension[syst]]], {i, 0, 
     dimension[syst] - 1}]], {ket[syst], bra[syst]}]

Z[syst_?IntegerQ, boost_?IntegerQ] := Z[Subscript[q, syst], boost]

pi[lambda_] := Fold[diracMatrixProduct, 
 Table[Z[i, lambda[[numSystems + i]]] ** X[i, lambda[[i]]], {i, 
   numSystems}]]

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

SDilation[A_] := 
 ArrayFlatten[{{A, 
    ConstantArray[0, Dimensions[A]]}, {ConstantArray[0, 
     Dimensions[A]], inverseAut@adjointEndo[A]}}]

phiDilation[l_, m_] := 1

phiDilation[lambda_] := phiDilation[lambda[[;;numSystems]], lambda[[numSystems+1;;]]]

Id = IdentityMatrix[2];

phiChirp[c_, k_] := 
 Exp[Pi I k.(Id + Inverse[Nmat]).c.(Id + Nmat).k]

ChirpU12[c_] := 
 matrix[Table[
   If[l1 == k1 && l2 == k2, phiChirp[c, {l1, l2}], 0], {l1, 0, 
    p - 1}, {k1, 0, p - 1}, {l2, 0, p - 1}, {k2, 0, p - 1}], {ket[q1],
    bra[q1], ket[q2], bra[q2]}]

SChirp[c_] := 
 ArrayFlatten[{{IdentityMatrix[2], 
    ConstantArray[0, {2, 2}]}, {{{c, 0}, {0, c}}, 
    IdentityMatrix[2]}}]

SChirp[C_?ListQ] := 
 ArrayFlatten[{{IdentityMatrix[2], ConstantArray[0, {2, 2}]}, {C, 
    IdentityMatrix[2]}}]

InverseMod[val_, p_] := If[GCD[val, p] == 1, PowerMod[val, -1, p], 0]
SetAttributes[InverseMod, Listable];


End[]

EndPackage[]