%%calculate chlorite structural formulas and endmembers 
% Copyright (c) 2021 Jesse Walters
%assumes all Fe is FeO

clear all;
close all;

D=load('Chl_EPMA.txt');
%input wt % oxide in the following order
%column1: SiO2
%column2: Al2O3
%column3: NiO
%column4: FeO
%column5: MnO
%column6: MgO

%OUTPUT: cations (APFU)
%column1: Si
%column2: Al
%column3: Ni
%column4: Fe total
%column5: Mn
%column6: Mg
%column7: total

%OUTPUT: Structural formula (APFU)
%column1: Si (T)
%column2: Al (T)
%column3: Sum of T site
%column4: Al (M)
%column5: Ni (M)
%column6: Fe total (M)
%column7: Mn (M)
%column8: Mg (M)
%column9: Sum M

%OUTPUT: Endmembers (fractions
%column1: XClc (clinochlore)
%column2: XSud (Sudoite)
%column3: XAme (Amesite)
%column4: XChm (Chamosite)
%column5: XPen (Pennantite)

Opfu=14.0; %oxygens per formula unit


%% Molecular weights

SiO2=60.0843;
Al2O3=101.9612;
NiO=74.6928;
FeO=71.8444;
MnO=70.9374;
MgO=40.3044;


%% Calculate cations units

[m,n]=size(D); %finds the x and y size of the input data matrix

MC(:,1)=D(:,1)./SiO2; %for SiO2
MC(:,2)=D(:,2)./Al2O3; %for Al2O3
MC(:,3)=D(:,3)./NiO; %for NiO
MC(:,4)=D(:,4)./FeO; %for FeO
MC(:,5)=D(:,5)./MnO; %for MnO
MC(:,6)=D(:,6)./MgO; %for MgO


%% Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*3; %for Al2O3
O2(:,3)=MC(:,3); %for NiO
O2(:,4)=MC(:,4); %for FeO
O2(:,5)=MC(:,5); %for MnO
O2(:,6)=MC(:,6); %for MgO

O2total=sum(O2,2); %O2 totals


%% Normalized Oxygen Units

O2norm=O2.*(Opfu./O2total);

%% Atoms per formula unit



APFU(:,1)=O2norm(:,1)./2; %for Si
APFU(:,2)=O2norm(:,2).*(2/3); %for Al
APFU(:,3)=O2norm(:,3); %for Ni
APFU(:,4)=O2norm(:,4); %for Fe
APFU(:,5)=O2norm(:,5); %for Mn
APFU(:,6)=O2norm(:,6); %for Mg
APFU(:,7)=sum(APFU,2); %cation sum (should be ~ 10)

%% Structural formula

% Si (T)
StrctFrm(:,1)=APFU(:,1); 

% Al (T)
for c=1:m
    if 4-StrctFrm(c,1)>APFU(c,2)
        StrctFrm(c,2)=APFU(c,2);
    else
        StrctFrm(c,2)=4-StrctFrm(c,1);
    end
end

%T site sum
StrctFrm(:,3)=StrctFrm(:,1)+StrctFrm(:,2);

% Al (M)
for c=1:m
    if 4-StrctFrm(c,1)>APFU(c,2)
        StrctFrm(c,4)=0;
    else
        StrctFrm(c,4)=APFU(c,2)-StrctFrm(c,2);
    end
end

StrctFrm(:,5)=APFU(:,3); %Ni (M)
StrctFrm(:,6)=APFU(:,4); %Fe total (M)
StrctFrm(:,7)=APFU(:,5); %Mn (M)
StrctFrm(:,8)=APFU(:,6); %Mg (M)
StrctFrm(:,9)=StrctFrm(:,8)+StrctFrm(:,7)+StrctFrm(:,6)+StrctFrm(:,5)+StrctFrm(c,4); %Total of M site

%% end member calculations

Endmembers=zeros(m,5); %matrix to be filled with endmember data

A(:,1)=APFU(:,6); %Mg
A(:,2)=APFU(:,5); %Mn
A(:,3)=APFU(:,4); %Fe
A(:,4)=APFU(:,3); %Ni
A(:,5)=APFU(:,2); %Al

M=[2 5 0 0 0; 0 0 0 0 5; 0 0 5 0 0; 0 0 0 5 0; 4 2 2 2 2]; 

AT=transpose(A); %transpose of A

X=zeros(5,m);
for c=1:m
    X(:,c)=inv(M)*AT(:,c); %calculates endmembers
end

Xtot=sum(X);%sum of endmembers
Xnorm=X./sum(X); %normalizes the endmembers to 1
XnormT=transpose(Xnorm); %transposes back 

Endmembers(:,1)=XnormT(:,1); % XSudonite
Endmembers(:,2)=XnormT(:,2); % XClinochlore
Endmembers(:,3)=XnormT(:,3); % XChamosite
Endmembers(:,4)=XnormT(:,4); % XNimite
Endmembers(:,5)=XnormT(:,5); % XPennantite 
