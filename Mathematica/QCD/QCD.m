(* ::Package:: *)

BeginPackage["QCD`"]

beta::usage =
"beta[nf, i] returns the i-th beta-function coefficients with flavour number nf"
betaFunction::usage =
"betaFunction[a, i, nf] returns the beta-function to the ith order with flavour number nf"
gamma::usage =
"gamma[nf, i] return the i-th anomalous dimension coefficient with flavour nf"
gammaFunction::usage =
"gammaFunction[a, i, nf] returns the anoamlous dimension to the ith order with flavour number nf"

Begin["`Private`"]
beta[nf_Integer, i_Integer]:= {
		11/2-1/3*nf, 
		51/4-19/12*nf, 
		2857/64-5033/576*nf+325/1728*nf^2,
		149753/768+891/32*Zeta[3]-(1078361/20736+1627/864*Zeta[3])*nf+(50065/20736+809/1296*Zeta[3])*nf^2+1093/93312*nf^3,
		2/4^5 (8157455/16 + 621885/2* Zeta[3] - 88209/2* Zeta[4] - 288090* Zeta[5]
			-(336460813/1944 + 4811164/81 Zeta[3] - 33935/6Zeta[4] - 1358995/27 Zeta[5]) nf 
			+(25960913/1944 + 698531/81 Zeta[3] - 10526/9 Zeta[4] - 381760/81 Zeta[5]) nf^2 
			-(630559/5832 + 48722/243 Zeta[3] - 1618/27 Zeta[4] - 460/9 Zeta[5]) nf^3 
			+(1205/2916 - 152/81 Zeta[3]) nf^4)
}[[i]];
betaFunction[a_, i_Integer, nf_Integer]:=Sum[beta[nf, j+1]*a^(j+2),{j,0,i-1}];

gamma[nf_Integer, i_Integer] := {
		2, 
		101/12 - 5/18*nf, 
		1249/32 - (277/108 + 5/3*Zeta[3])*nf - 35/648*nf*nf,
		4603055/20736+1060/27*Zeta[3]-275/4*Zeta[5]+(11/144*Pi^4-91723/3456-2137/72*Zeta[3]+575/36*Zeta[5])*nf + (2621/15552-Pi^4/216+25/36*Zeta[3])*nf^2 - (83/7776-Zeta[3]/54)*nf^3
}[[i]];
gammaFunction[a_, i_Integer, nf_Integer]:=Sum[gamma[nf, j+1]*a^(j+1),{j,0,i-1}];		
End[]
EndPackage[]




