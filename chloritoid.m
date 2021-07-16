%%calculates chloritoid formulae
%Fe3+ calculated following method 2 of Schumacher (1991)
% Copyright (c) 2021 Jesse Walters

clear all;
close all;

D=load('chloritoid.txt'); %EPMA data is loaded from a text file here

%Input wt % oxides in the following order: 
%column1: SiO2
%column2: TiO2
%column3: Al2O3
%column4: FeO
%column5: MnO
%column6: MgO
%column7: CaO
%column8: Na2O

%OUTPUT: cations (APFU)
%column1: Si
%column2: Ti
%column3: Al
%column4: Fe3+
%column5: Fe2+
%column6: Mn
%column7: Mg
%column8: Ca
%column9: Na
%column12: total
%column13: O2 deficiency 

cat=8.0; %cations per formula unit
Opfu=12.0; %oxygens per formula unit

%% Molecular weights

SiO2=60.084;
TiO2=79.866;
Al2O3=101.961;
FeO=71.844;
MnO=70.937;
MgO=40.304;
CaO=56.077;
Na2O=61.979;

W=[SiO2,TiO2,Al2O3,FeO,MnO,MgO,CaO,Na2O];

%% Calculate cations units

[m,n]=size(D); %finds the x and y size of the input data matrix
MC=zeros(size(D)); %creates a matrix of zeroes the size of the input data
MC(:,1)=D(:,1)./W(:,1); %for SiO2
MC(:,2)=D(:,2)./W(:,2); %for TiO2
MC(:,3)=(D(:,3)./W(:,3)).*2; %for Al2O3
MC(:,4)=D(:,4)./W(:,4); %for FeO
MC(:,5)=D(:,5)./W(:,5); %for MnO
MC(:,6)=D(:,6)./W(:,6); %for MgO
MC(:,7)=D(:,7)./W(:,7); %for CaO
MC(:,8)=(D(:,8)./W(:,8)).*2; %for Na2O

MCnormfact=zeros(m,1); %creates a zeromatrix for the totals of the mole cation matrix 
MCnormfact=cat./sum(MC,2); %normalization factor

%% Calculate normalized cations units

MCnorm=MCnormfact.*MC; %creates a matrix of normalized cations

%% Calculate Oxygen Units

O2=zeros(size(D));
O2(:,1)=MCnorm(:,1).*2; %for SiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4); %for FeO
O2(:,5)=MCnorm(:,5); %for MnO
O2(:,6)=MCnorm(:,6); %for MgO
O2(:,7)=MCnorm(:,7); %for CaO
O2(:,8)=MCnorm(:,8)./2; %for Na2O

O2total=sum(O2,2); %O2 totals

%% Atoms pfu

APFU=zeros(m,n+3); %matrix of zeroes to be filled, n+2 for Fe3+ and total

APFU(:,1)=MCnorm(:,1); %for Si
APFU(:,2)=MCnorm(:,2); %for Ti
APFU(:,3)=MCnorm(:,3); %for Al
APFU(:,6)=MCnorm(:,5); %for Mn
APFU(:,7)=MCnorm(:,6); %for Mg
APFU(:,8)=MCnorm(:,7); %for Ca
APFU(:,9)=MCnorm(:,8); %for Na


%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 12
%if so, then there is no Fe3+
%if totalO2 < 12, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(12-totalO2) then the amount
%of Fe3+ = 2*(12-totalO2), if false then, all Fe is Fe3+

for c=1:m
    if (Opfu-O2total(c,1)) > 0;
        if MCnorm(c,4) > 2*(Opfu-O2total(c,1));
            APFU(c,4)=2*(Opfu-O2total(c,1)); 
        else
            APFU(c,4)=MCnorm(c,4);
        end
    else
        APFU(c,4)=MCnorm(c,4);
    end
end

APFU(:,5)=MCnorm(:,4)-APFU(:,4); %the APFU of Fe2+ equals totalFe-Fe3+

APFU(:,12)=sum(APFU,2); %calculations the total, which should be 8

% Oxygen deficiency 
APFU(:,13)=Opfu-O2total; %must be greater than zero

%% Structural Formula

%T SITE
StrctFrm(:,1)=APFU(:,1); %Si

%Al2O3 Layer (L1)
StrctFrm(:,2)=3; %Al

%(Al, Ti, Fe3+) + (Fe, Mg, Mn) layer (L2)
StrctFrm(:,3)=APFU(:,3)-3; %Al
StrctFrm(:,4)=APFU(:,2); %Ti
StrctFrm(:,5)=APFU(:,4); %Fe3+
StrctFrm(:,6)=StrctFrm(:,5)+StrctFrm(:,4)+StrctFrm(:,3); %Al + Ti + Fe3+
StrctFrm(:,7)=APFU(:,5); %Fe2+
StrctFrm(:,8)=APFU(:,6); %Mn
StrctFrm(:,9)=APFU(:,7); %Mg
StrctFrm(:,10)=APFU(:,8); %Ca
StrctFrm(:,11)=APFU(:,9); %Na
StrctFrm(:,12)=StrctFrm(:,11)+StrctFrm(:,10)+StrctFrm(:,9)+StrctFrm(:,8)+StrctFrm(:,7); %Sum of Fe2+, Mn, Mg, Ca, and Na 
