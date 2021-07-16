%%Epidote structural formula
% Copyright (c) 2021 Jesse Walters

clear all;
close all;

Ep=load('Ep_EPMA.txt'); %EPMA data is loaded from a text file here

%input wt % oxide in the following order
%column1: SiO2
%column2: TiO2
%column3: Al2O3
%column4: Cr2O3
%column5: Fe2O3
%column6: MnO
%column7: MgO
%column8: CaO
%column9: Na2O
%column10: K2O

%OUTPUTS

%StrctFrm: Structural formula for individual analyses
%column1: Si (T)
%column2: Ti (T)
%column3: Al (T)
%column4: T sum
%column5: Ti (M)
%column6: Al (M)
%column7: Cr (M)
%column8: Mn (M)
%column9: Fe3+ (M)
%column10: M sum
%column11: Fe2+ (A)
%column12: Mg (A)
%column13: Ca (A)
%column14: Na (A)
%column15: K (A)
%column16: A sum

Opfu=12.5; %oxygens per formula unit


%% Molecular weights

SiO2=60.084;
TiO2=79.866;
Al2O3=101.961;
Cr2O3=151.99;
Fe2O3=159.688;
FeO=71.844;
MnO=70.937;
MgO=40.304;
CaO=56.077;
Na2O=61.979;
K2O=94.196;

W=[SiO2,TiO2,Al2O3,Cr2O3,Fe2O3,MnO,MgO,CaO,Na2O,K2O];


%% Calculate cations units

[m,n]=size(Ep); %finds the x and y size of the input data matrix
MC=zeros(size(Ep)); %creates a matrix of zeroes the size of the input data
MC(:,1)=Ep(:,1)./SiO2; %for SiO2
MC(:,2)=Ep(:,2)./TiO2; %for TiO2
MC(:,3)=Ep(:,3)./Al2O3; %for Al2O3
MC(:,4)=Ep(:,4)./Cr2O3; %for Cr2O3
MC(:,5)=Ep(:,5)./(2*FeO); %for Fe2O3
MC(:,6)=Ep(:,6)./MnO; %for MnO
MC(:,7)=Ep(:,7)./MgO; %for MgO
MC(:,8)=Ep(:,8)./CaO; %for CaO
MC(:,9)=Ep(:,9)./Na2O; %for Na2O
MC(:,10)=Ep(:,10)./K2O; %for K2O

%% Calculate Oxygen Units

O2=zeros(size(Ep));
O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3); %for Al2O3
O2(:,4)=MC(:,4).*(3); %for Cr2O3
O2(:,5)=MC(:,5).*(3); %for Fe2O3
O2(:,6)=MC(:,6); %for MnO
O2(:,7)=MC(:,7); %for MgO
O2(:,8)=MC(:,8); %for CaO
O2(:,9)=MC(:,9); %for Na2O
O2(:,10)=MC(:,10); %for K2O

O2total=sum(O2,2); %O2 totals

%% Normalized Oxygen Units

O2norm=O2.*(Opfu./O2total);

%% Atoms pfu

APFU=zeros(m,n+1); %matrix of zeroes to be filled, n+2 for Fe3+ and total

for c=1:m
APFU(c,1)=O2norm(c,1)./2; %for Si
end

for c=1:m
APFU(c,2)=O2norm(c,2)./2; %for Ti
end

for c=1:m
APFU(c,3)=O2norm(c,3).*(2/3); %for Al
end

for c=1:m
APFU(c,4)=O2norm(c,4).*(2/3); %for Cr3+
end

for c=1:m
APFU(c,5)=O2norm(c,5).*(2/3); %for Fe3+
end

for c=1:m
APFU(c,6)=O2norm(c,6); %for Mn
end

for c=1:m
APFU(c,7)=O2norm(c,7); %for Mg
end

for c=1:m
APFU(c,8)=O2norm(c,8); %for Ca
end

for c=1:m
APFU(c,9)=O2norm(c,9).*2; %for Na
end

for c=1:m
APFU(c,10)=O2norm(c,10).*2; %for K
end

APFU(:,11)=sum(APFU,2);

%% Structural Formula
%T site - Si, Ti, Al (sums to 3)
%M sites - Ti, Al, Cr3+, Fe3+, Mn3+
%A site - Ca, Na, K, Mg, Fe2+ (Sums to 2)

% if Al + Fe3+ + Cr + Ti + Mn in the M site is > 3 then the remaining Fe3+ is
% converted to Fe2+ in the A site

StrctFrm=zeros(m,16);%creates a matrix to be filled

%T SITE
%Si 
for c=1:m
    StrctFrm(c,1)=APFU(c,1);
end

%Ti(T)
for c=1:m
    if StrctFrm(c,1)<3
        if (StrctFrm(c,1)+APFU(c,2))>3
            StrctFrm(c,2)=3-APFU(c,1);
        else
            StrctFrm(c,2)=APFU(c,2);
        end
    else
        StrctFrm(c,2)=0;
    end
end

%Al(T)
for c=1:m
    if (StrctFrm(c,1)+StrctFrm(c,2))<3
        StrctFrm(c,3)=3-(StrctFrm(c,1)+StrctFrm(c,2));
    else
        StrctFrm(c,3)=0;
    end
end

%T-site sum
for c=1:m
    StrctFrm(c,4)=StrctFrm(c,3)+StrctFrm(c,2)+StrctFrm(c,1);
end

%M site

%Ti (M)
for c=1:m
    StrctFrm(c,5)=APFU(c,2)-StrctFrm(c,2);
end

%Al (M)
for c=1:m
    StrctFrm(c,6)=APFU(c,3)-StrctFrm(c,3);
end

%Cr (M)
for c=1:m
    StrctFrm(c,7)=APFU(c,4);
end

%Mn (M)
for c=1:m
    StrctFrm(c,8)=APFU(c,6);
end

%Fe3+ (M)
for c=1:m
    if (StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+APFU(c,5))>3
        StrctFrm(c,9)=3-(StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8));
    else
        StrctFrm(c,9)=APFU(c,5);
    end
end

%M site sum
for c=1:m
    StrctFrm(c,10)=StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9);
end

%A sites

%Fe2+ (A)
for c=1:m
    StrctFrm(c,11)=APFU(c,5)-StrctFrm(c,9);
end

%Mg (A)
for c=1:m
    StrctFrm(c,12)=APFU(c,7);
end

%Ca (A)
for c=1:m
    StrctFrm(c,13)=APFU(c,8);
end

%Na (A)
for c=1:m
    StrctFrm(c,14)=APFU(c,9);
end

%K (A)
for c=1:m
    StrctFrm(c,15)=APFU(c,10);
end

%A site sum
for c=1:m
    StrctFrm(c,16)=StrctFrm(c,11)+StrctFrm(c,12)+StrctFrm(c,13)+StrctFrm(c,14)+StrctFrm(c,15);
end
