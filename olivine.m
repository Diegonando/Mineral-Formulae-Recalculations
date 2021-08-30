%%Calculates Olivine structural formula, 
%Schumacher (1991) for Fe3+ by stoichiometry 

%Copyright (c) 2021 Jesse Walters
%The DOI for this release is: 10.5281/zenodo.5110201

clear all;
close all;

D=load('Ol_EPMA.txt'); %EPMA data is loaded from a text file here

%input wt % oxide in the following order
%column1: SiO2
%column2: TiO2
%column3: Al2O3
%column4: Cr2O3
%column5: NiO
%column6: FeO
%column7: MnO
%column8: MgO
%column9: CaO


%OUTPUT: cations (APFU)
%column1: Si
%column2: Ti
%column3: Al
%column4: Fe3+
%column5: Cr
%column6: Ni
%column7: Fe2+
%column8: MnO
%column9: MgO
%column10: CaO
%column11: total
%column12: O2 deficiency 

%OUTPUT: Structural formula (APFU)
%column1: Si (T)
%column2: Al (T)
%column3: Fe3+ (T)
%column4: T sum
%column5: Ti (M)
%column6: Al (M)
%column7: Fe3+ (M)
%column8: Cr (M)
%column9: Ni (M)
%column10: Fe (M)
%column11: Mn (M)
%column12: Mg (M)
%column13: Ca (M)
%column14: M site total

%Output: endmembers
%column1: XFo
%column2: XFa
%column3: XTe
%column4: XCa-Ol


cat=3.0; %cations per formula unit
Opfu=4.0; %oxygens per formula unit

%% Molecular weights

SiO2=60.084;
TiO2=79.866;
Al2O3=101.961;
Fe2O3=159.688;
Cr2O3=151.99;
NiO=74.6928;
FeO=71.844;
MnO=70.937;
MgO=40.304;
CaO=56.077;


W=[SiO2,TiO2,Al2O3,Cr2O3,NiO,FeO,MnO,MgO,CaO];

%% Calculate cations units

[m,n]=size(D); %finds the x and y size of the input data matrix
MC=zeros(size(D)); %creates a matrix of zeroes the size of the input data
MC(:,1)=D(:,1)./W(:,1); %for SiO2
MC(:,2)=D(:,2)./W(:,2); %for TiO2
MC(:,3)=(D(:,3)./W(:,3)).*2; %for Al2O3
MC(:,4)=(D(:,4)./W(:,4)).*2; %for Cr2O3
MC(:,5)=D(:,5)./W(:,5); %for NiO
MC(:,6)=D(:,6)./W(:,6); %for FeO
MC(:,7)=D(:,7)./W(:,7); %for MnO
MC(:,8)=D(:,8)./W(:,8); %for MgO
MC(:,9)=D(:,9)./W(:,9); %for CaO

MCnormfact=zeros(m,1); %creates a zeromatrix for the totals of the mole cation matrix 
MCnormfact=cat./sum(MC,2); %normalization factor

%% Calculate normalized cations units

MCnorm=MCnormfact.*MC; %creates a matrix of normalized cations

%% Calculate Oxygen Units

O2=zeros(size(D));
O2(:,1)=MCnorm(:,1).*2; %for SiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4).*(3/2); %for Cr2O3
O2(:,5)=MCnorm(:,5); %for NiO
O2(:,6)=MCnorm(:,6); %for FeO
O2(:,7)=MCnorm(:,7); %for MnO
O2(:,8)=MCnorm(:,8); %for MgO
O2(:,9)=MCnorm(:,9); %for CaO

O2total=sum(O2,2); %O2 totals

%% Atoms pfu

APFU=zeros(m,n+2); %matrix of zeroes to be filled, n+2 for Fe3+ and total

APFU(:,1)=MCnorm(:,1); %for Si
APFU(:,2)=MCnorm(:,2); %for Ti
APFU(:,3)=MCnorm(:,3); %for Al
APFU(:,5)=MCnorm(:,4); %for Cr2O3
APFU(:,6)=MCnorm(:,5); %for NiO
APFU(:,8)=MCnorm(:,7); %for MnO
APFU(:,9)=MCnorm(:,8); %for MgO
APFU(:,10)=MCnorm(:,9); %for CaO


%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 4
%if so, then there is no Fe3+
%if totalO2 < 4, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(4-totalO2) then the amount
%of Fe3+ = 2*(4-totalO2), if false then, all Fe is Fe3+
for c=1:m
    if (Opfu-O2total(c,1)) >= 0;
        if MCnorm(c,6) > 2*(Opfu-O2total(c,1));
            APFU(c,4)=2*(Opfu-O2total(c,1)); 
        else
            APFU(c,4)=MCnorm(c,6);
        end
    else
        APFU(c,4)=0;
    end
end

APFU(:,7)=MCnorm(:,6)-APFU(:,4); %the APFU of Fe2+ equals totalFe-Fe3+

APFU(:,11)=sum(APFU,2); %calculations the total, which should be 4

% Oxygen deficiency 
APFU(:,12)=Opfu-O2total; %must be greater than zero

%% structural formula calculation

StrctFrm=zeros(m,n+8);%creates a matrix to be filled

%T SITE
%Si 
for c=1:m
    if APFU(c,1)<1.000
    StrctFrm(c,1)=APFU(c,1);
    else
    StrctFrm(c,1)=1;
    end
end

%Al(T)
for c=1:m
    if 1-StrctFrm(c,1)>APFU(c,3) 
        StrctFrm(c,2)=APFU(c,3);
    else
        if 1-StrctFrm(c,1)>0 
           StrctFrm(c,2)=2-StrctFrm(c,1); 
        else
            StrctFrm(c,2)=0; %if Si is 1 or more, no Al in T
        end
    end
end

%Fe3+(T)
for c=1:m
    if 1-StrctFrm(c,1)-StrctFrm(c,2)>0 && 1-StrctFrm(c,1)-StrctFrm(c,2)<=APFU(c,4) 
         StrctFrm(c,3)=1-StrctFrm(c,1)-StrctFrm(c,2);
    else
        if 1-StrctFrm(c,1)-StrctFrm(c,2)>0
            StrctFrm(c,3)=APFU(c,4);
        else
            StrctFrm(c,3)=0;
        end
    end
end

%Sum of T site
StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3);

%M SITE
%Ti
StrctFrm(:,5)=APFU(:,2);

%Al
StrctFrm(:,6)=APFU(:,3)-StrctFrm(:,2);

%Fe3+
StrctFrm(:,7)=APFU(:,4)-StrctFrm(:,3);

%Cr
StrctFrm(:,8)=APFU(:,5);

%Ni
StrctFrm(:,9)=APFU(:,6);

%Fe
StrctFrm(:,10)=APFU(:,7);

%Mn
StrctFrm(:,11)=APFU(:,8);

%Mg
StrctFrm(:,12)=APFU(:,9);

%Ca
StrctFrm(:,13)=APFU(:,10);

%M total
StrctFrm(:,14)=sum(StrctFrm(:,5:13),2);

%% Endmember calculation

Endmembers(:,1)=APFU(:,9)./(APFU(:,4)+APFU(:,7)+APFU(:,8)+APFU(:,9)+APFU(:,10)); % XFo
Endmembers(:,2)=(APFU(:,4)+APFU(:,7))./(APFU(:,4)+APFU(:,7)+APFU(:,8)+APFU(:,9)+APFU(:,10)); % XFa
Endmembers(:,3)=(APFU(:,8))./(APFU(:,4)+APFU(:,7)+APFU(:,8)+APFU(:,9)+APFU(:,10)); % XTe
Endmembers(:,4)=(APFU(:,10))./(APFU(:,4)+APFU(:,7)+APFU(:,8)+APFU(:,9)+APFU(:,10)); % XCa-Ol

