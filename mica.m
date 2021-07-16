%% calculates mica structural formula 
%Uses 11 Oxygen normalization scheme
%Assumes all Fe is FeO

clear all;
close all;

D=load('Mica.txt'); %EPMA data is loaded from a text file here

%Input data:
%Column1: SiO2
%Column2: TiO2
%Column3: Al2O3
%Column4: Cr2O3
%Column5: FeO
%Column6: MnO
%Column7: MgO
%Column8: CaO
%Column9: Na2O
%Column10: K2O
%Column11: BaO
%Column12: F
%Column13: Cl

%Output: StrctFrm 
%Column1: Si (T)
%Column2: Al (T)
%Column3: T total
%Column4: Al (M)
%Column5: Ti (M)
%Column6: Cr (M)
%Column7: Fe (M)
%Column8: Mn (M)
%Column9: Mg (M)
%Column10: M total
%Column11: Ca (I)
%Column12: Na (I)
%Column13: K (I)
%Column14: Ba (I)
%Column15: Total Cations
%Column16: F (A)
%Column17: Cl (A)
% %Column18: OH (A)

%Output: DiOct_Endmembers
%XMuscovite
%XFe-Celadonite (ferroaluminoceladonite)
%XMg-Celadonite (Magnesioaluminoceladonite)
%XParagonite 
%XMargarite
%XTrioctahedral

%Output: TriOct_Endmembers 
%Xphlogopite
%Xannite
%Xeastonite 
%Xsiderophyllite
%Xdioctohedral component

%% Molecular weights

SiO2=60.084;
TiO2=79.866;
Al2O3=101.961;
Cr2O3=151.9904;
FeO=71.8444;
MnO=70.937;
MgO=40.304;
CaO=56.077;
Na2O=61.979;
K2O=94.196;
BaO=153.3264;
F=18.998403; 
Cl=35.453; 

W=[SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O,BaO,F,Cl];

%% Moles of cations

[m,n]=size(D); %finds the x and y size of the input data matrix

MC=D./W;

%% Moles of O2 units
O2(:,1)=MC(:,1).*2; %SiO2
O2(:,2)=MC(:,2).*2; %TiO2
O2(:,3)=MC(:,3).*3; %Al2O3
O2(:,4)=MC(:,4).*3; %Cr2O3
O2(:,5)=MC(:,5); %FeO
O2(:,6)=MC(:,6); %MnO
O2(:,7)=MC(:,7); %MgO
O2(:,8)=MC(:,8); %CaO
O2(:,9)=MC(:,9); %Na2O
O2(:,10)=MC(:,10); %K2O
O2(:,11)=MC(:,11); %BaO
O2(:,12)=MC(:,12); %F
O2(:,13)=MC(:,13); %Cl

O2_Sum=sum(O2(:,1:11),2); %sum of O2, not including F and Cl

O2_N=11./O2_Sum; %normalization factor

%moles of anions
N_Ox=O2.*O2_N;

%% Structural formula

APFU(:,1)=N_Ox(:,1)./2; %Si
APFU(:,2)=N_Ox(:,2)./2; %Ti
APFU(:,3)=N_Ox(:,3).*(2/3); %Al
APFU(:,4)=N_Ox(:,4).*(2/3); %Cr
APFU(:,5)=N_Ox(:,5); %Fe
APFU(:,6)=N_Ox(:,6); %Mn
APFU(:,7)=N_Ox(:,7); %Mg
APFU(:,8)=N_Ox(:,8); %Ca
APFU(:,9)=N_Ox(:,9).*2; %Na
APFU(:,10)=N_Ox(:,10).*2; %K
APFU(:,11)=N_Ox(:,11); %Ba
APFU(:,12)=N_Ox(:,12); %F
APFU(:,13)=N_Ox(:,13); %Cl

%T site
StrctFrm(:,1)=APFU(:,1); %Si (T)

%Al (T)
for c=1:m
    if 4-StrctFrm(c,1) > APFU(c,3)
        StrctFrm(c,2)=APFU(c,3); 
    else
        StrctFrm(c,2)=4-StrctFrm(c,1);
    end
end

StrctFrm(:,3)=StrctFrm(:,1)+StrctFrm(:,2); %Sum of T

StrctFrm(:,4)=APFU(:,3)-StrctFrm(:,2); %Al(M)
StrctFrm(:,5)=APFU(:,2); %Ti (M)
StrctFrm(:,6)=APFU(:,4); %Cr (M)
StrctFrm(:,7)=APFU(:,5); %Fe2+ (M)
StrctFrm(:,8)=APFU(:,6); %Mn (M)
StrctFrm(:,9)=APFU(:,7); %Mg (M)
StrctFrm(:,10)=sum(StrctFrm(:,4:1:9),2); %sum (M)

StrctFrm(:,11)=APFU(:,8); %Ca (I)
StrctFrm(:,12)=APFU(:,9); %Na (I)
StrctFrm(:,13)=APFU(:,10); %K (I)
StrctFrm(:,14)=APFU(:,11); %Ba (I)

StrctFrm(:,15)=sum(StrctFrm(:,11:1:14),2)+StrctFrm(:,10)+StrctFrm(:,3); %cation sum

StrctFrm(:,16)=APFU(:,12); % F (A)
StrctFrm(:,17)=APFU(:,13); % Cl (A)
StrctFrm(:,18)=2-(APFU(:,12)+APFU(:,13)); % OH (A)

%% Endmembers

%Dioctahedral Micas
for c=1:m
   if StrctFrm(c,10)>2
       XTriOct(c,:)=StrctFrm(c,10)-2;
   else
       XTriOct(c,:)=0;
   end
end


A1(:,1)=APFU(:,10); %K
A1(:,2)=APFU(:,9); %Na
A1(:,3)=APFU(:,8); %Ca
A1(:,4)=APFU(:,5); %Fe
A1(:,5)=APFU(:,7); %Mg

M1=[1 1 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 1 0 0 0; 0 0 1 0 0];

AT1=transpose(A1); %transpose of A1

X1=zeros(5,m);
for c=1:m
    X1(:,c)=inv(M1)*AT1(:,c); %calculates endmembers
end

XT1=transpose(X1./sum(X1)); %normalizes and transposes back 
Xfact1=sum(XT1,2)+XTriOct; %normalization factor including trioct mica endmembers
Xnorm1=XT1;
Xnorm1(:,6)=XTriOct(:);
DiOct_Endmembers=Xnorm1./(Xfact1);

%Trioctahedral Mica
%this only works if Si (APFU) is between 3 and 2

XDiOct=1-XTriOct(:); %fraction of dioctahedral mica


%pholgopite-annite solid solution
XPhlAnn=APFU(:,1)-2; %fraction of the Phl+Ann endmembers
                     %Phl and Ann has 3 Si, whereas Sid and East has 2
XMg=APFU(:,7)./(APFU(:,7)+APFU(:,5)); %calculates XMg
XPhl=XPhlAnn.*XMg; %unnormalized fraction of Phl
XAnn=XPhlAnn-XPhl; %unnormalized fraction of Ann

%siderophyllite-eastonite solid solution
XSidEast=1-XPhlAnn; %fraction of the Sid+Eastendmembers
XEast=XSidEast.*XMg;
XSid=XSidEast-XEast;

xfact2=1./(XPhl+XAnn+XEast+XSid+XDiOct);

TriOct_Endmembers(:,1)=XPhl.*xfact2; %phlogopite
TriOct_Endmembers(:,2)=XAnn.*xfact2; %annite
TriOct_Endmembers(:,3)=XEast.*xfact2; %eastonite 
TriOct_Endmembers(:,4)=XSid.*xfact2; %siderophyllite
TriOct_Endmembers(:,5)=XDiOct.*xfact2; %dioctohedral component
