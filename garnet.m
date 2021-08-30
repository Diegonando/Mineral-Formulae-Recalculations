%%Garnet Structural Formula
%Copyright (c) 2021 Jesse Walters
%The DOI for this release is: 10.5281/zenodo.5110201

clear all;
close all;

D=load('grt_EPMA.txt'); %EPMA data is loaded from a text file here

%input wt % oxide in the following order
%column1: SiO2
%column2: TiO2
%column3: Al2O3
%column4: Cr2O3
%column5: Y2O3
%column6: FeO
%column7: MnO
%column8: MgO
%column9: CaO
%column10: Na2O

%OUTPUT: cations (APFU)
%column1: Si
%column2: Ti
%column3: Al
%column4: Fe3+
%column5: Cr
%column6: Y
%column7: Fe2+
%column8: Mn
%column9: Mg
%column10: Ca
%column11: Na
%column12: total
%column13: O2 deficiency 

%OUTPUT: Structural formula (APFU)
%column1: Si (T)
%column2: Al (T)
%column3: Fe3+ (T)
%column4: sum of T site
%column5: Al (Y)
%column6: Ti (Y)
%column7: Cr (Y)
%column8: Fe3+ (Y)
%column9: Sum Y
%column10: Y (X)
%column11: Fe2+ (X)
%column12: Mn (X)
%column13: Mg (X)
%column14: Ca (X)
%column15: Na (X)
%column16: Sum X

%OUTPUT: Endmembers (fractions
%column1: XAlm
%column2: XPrp
%column3: XSps
%column4: XGrs
%column5: XAdr
%column6: XUv

cat=8.0; %cations per formula unit
Opfu=12.0; %oxygens per formula unit


%% Molecular weights

SiO2=60.084;
TiO2=79.866;
Al2O3=101.961;
Fe2O3=159.688;
Cr2O3=151.99;
Y2O3=225.81;
FeO=71.844;
MnO=70.937;
MgO=40.304;
CaO=56.077;
Na2O=61.979;
K2O=94.196;

W=[SiO2,TiO2,Al2O3,Cr2O3,Y2O3,FeO,MnO,MgO,CaO,Na2O];

%% Calculate cations units

[m,n]=size(D); %finds the x and y size of the input data matrix
MC=zeros(size(D)); %creates a matrix of zeroes the size of the input data
MC(:,1)=D(:,1)./W(:,1); %for SiO2
MC(:,2)=D(:,2)./W(:,2); %for TiO2
MC(:,3)=(D(:,3)./W(:,3)).*2; %for Al2O3
MC(:,4)=(D(:,4)./W(:,4)).*2; %for Cr2O3
MC(:,5)=(D(:,5)./W(:,5)).*2; %for Y2O3
MC(:,6)=D(:,6)./W(:,6); %for FeO
MC(:,7)=D(:,7)./W(:,7); %for MnO
MC(:,8)=D(:,8)./W(:,8); %for MgO
MC(:,9)=D(:,9)./W(:,9); %for CaO
MC(:,10)=(D(:,10)./W(:,10)).*2; %for Na2O


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
O2(:,5)=MCnorm(:,5).*(3/2); %for Y2O3
O2(:,5)=MCnorm(:,6); %for FeO
O2(:,6)=MCnorm(:,7); %for MnO
O2(:,7)=MCnorm(:,8); %for MgO
O2(:,8)=MCnorm(:,9); %for CaO
O2(:,9)=MCnorm(:,10)./2; %for Na2O


O2total=sum(O2,2); %O2 totals

%% Atoms pfu

APFU=zeros(m,n+3); %matrix of zeroes to be filled, n+2 for Fe3+ and total

APFU(:,1)=MCnorm(:,1); %for Si
APFU(:,2)=MCnorm(:,2); %for Ti
APFU(:,3)=MCnorm(:,3); %for Al
APFU(:,5)=MCnorm(:,4); %for Cr
APFU(:,6)=MCnorm(:,5); %for Y
APFU(:,8)=MCnorm(:,7); %for MnO
APFU(:,9)=MCnorm(:,8); %for MgO
APFU(:,10)=MCnorm(:,9); %for CaO
APFU(:,11)=MCnorm(:,10); %for Na2O


%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 12
%if so, then there is no Fe3+
%if totalO2 < 12, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(12-totalO2) then the amount
%of Fe3+ = 2*(12-totalO2), if false then, all Fe is Fe3+

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

APFU(:,12)=sum(APFU,2); %calculations the total, which should be 8

% Oxygen deficiency 
APFU(:,13)=Opfu-O2total; %must be greater than zero

%% structural formula calculation

StrctFrm=zeros(m,n+6);%creates a matrix to be filled
%T SITE
%Si
StrctFrm=APFU(:,1);

%Al(T)
for c=1:m
    if 3-StrctFrm(c,1)>APFU(c,3) %if low Al Grt, all Al may be in T
        StrctFrm(c,2)=APFU(c,3);
    else
        if 3-StrctFrm(c,1)>0 
           StrctFrm(c,2)=3-StrctFrm(c,1); %for most Grt's, only some Al will be in T
        else
            StrctFrm(c,2)=0; %if Si is 3 or more, no Al in T
        end
    end
end

%Fe3+(T), for Al poor Grt
for c=1:m
    if 3-(StrctFrm(c,1)+StrctFrm(c,1))>0
        StrctFrm(c,3)=3-(StrctFrm(c,1)+StrctFrm(c,1));
    else
        StrctFrm(c,3)=0;
    end
end

%Sum of T site
StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3);

%Y SITE

%Al (Y)
StrctFrm(:,5)=APFU(:,3)-StrctFrm(:,2);

%Ti (Y)
StrctFrm(:,6)=APFU(:,2);

%Cr (Y) 
StrctFrm(:,7)=APFU(:,5);

%Fe3+ (Y)
StrctFrm(:,8)=APFU(:,4)-StrctFrm(:,3);

%Y sum
StrctFrm(:,9)=StrctFrm(:,8)+StrctFrm(:,7)+StrctFrm(:,6)+StrctFrm(:,5);

%X SITE

%Y (X)
StrctFrm(:,10)=APFU(:,6);

%Fe2+ (X)
StrctFrm(:,11)=APFU(:,7);

%Mn (X)
StrctFrm(:,12)=APFU(:,8);

%Mg (X)
StrctFrm(:,13)=APFU(:,9);

%Ca (X)
StrctFrm(:,14)=APFU(:,10);

%Na (X)
StrctFrm(:,15)=APFU(:,11);

% X Sum
StrctFrm(:,16)=StrctFrm(:,15)+StrctFrm(:,14)+StrctFrm(:,13)+StrctFrm(:,12)+StrctFrm(:,11)+StrctFrm(:,10);

%% end member calculations

Endmembers=zeros(m,6); %matrix to be filled with endmember data
A=zeros(m,6); %A vector to be filled
A(:,1)=APFU(:,10); %Ca
A(:,2)=APFU(:,9); %Mg
A(:,3)=APFU(:,4)+APFU(:,7); %Fetotal
A(:,4)=APFU(:,5); %Cr
A(:,5)=APFU(:,8); %Mn
A(:,6)=APFU(:,3); %Al
AT=transpose(A); %transpose of A

M=[0 0 0 3 3 3; 0 3 0 0 0 0; 3 0 0 0 2 0; 0 0 0 0 0 2; 0 0 3 0 0 0; 2 2 2 2 0 0];

X=zeros(6,m);
for c=1:m
    X(:,c)=inv(M)*AT(:,c); %calculates endmembers
end

Xtot=sum(X);%sum of endmembers
Xnorm=X./sum(X); %normalizes the endmembers to 1
XnormT=transpose(Xnorm); %transposes back 

Endmembers(:,1)=XnormT(:,1); % XAlm
Endmembers(:,2)=XnormT(:,2); % XPrp
Endmembers(:,3)=XnormT(:,3); % XSps
Endmembers(:,4)=XnormT(:,4); % XGrs
Endmembers(:,5)=XnormT(:,5); % XAdr 
Endmembers(:,6)=XnormT(:,6); % XUv

