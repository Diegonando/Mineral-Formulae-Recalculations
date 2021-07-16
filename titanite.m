%% calculates titanite structural formula 
%Uses cation normalization scheme

clear all;
close all;

D=load('Ttn_EPMA.txt'); %EPMA data is loaded from a text file here

%input wt % oxide in the following order
%column1: SiO2
%column2: TiO2
%column3: Al2O3
%column4: Y2O3
%column5: FeO
%column6: MnO
%column7: MgO
%column8: CaO
%column9: Na2O
%column10: K2O
%column11: F

%StrctFrm output (APFU): 
%column1: Si4+ (T)
%column2: Ti4+ (T)
%column3: Ti4+ (Oct) 
%column4: Al3+ (Oct)
%column5: Fe3+ (Oct)
%column6: Mn (Oct)
%column7: Mg (Oct)
%column8: Ca (Dec)
%column9: Y3+ (Dec)
%column10: Na+ (Dec)
%column11: K+ (Dec)
%column12: cation sum
%column13: F
%column14: OH
%column15: O
%column16: anion sum

%XTtn - fraction of CaTiSiO5 

%% Molecular weights

SiO2=60.084;
TiO2=79.866;
Al2O3=101.961;
Fe2O3=159.688;
Y2O3=225.81;
MnO=70.937;
MgO=40.304;
CaO=56.077;
Na2O=61.979;
K2O=94.196;
F=19; 

%% moles of cations

[m,n]=size(D); %finds the x and y size of the input data matrix
MC=zeros(size(D)); %creates a matrix of zeroes the size of the input data

MC(:,1)=D(:,1)./SiO2; %for SiO2
MC(:,2)=D(:,2)./TiO2; %for TiO2
MC(:,3)=(D(:,3)./Al2O3)*2; %for Al2O3
MC(:,4)=(D(:,4)./Y2O3)*2; %for Y2O3
MC(:,5)=((D(:,5).*1.111378)./Fe2O3)*2; %for Fe2O3
MC(:,6)=D(:,6)./MnO; %for MnO
MC(:,7)=D(:,7)./MgO; %for MgO
MC(:,8)=D(:,8)./CaO; %for CaO
MC(:,9)=(D(:,9)./Na2O)*2; %for Na2O
MC(:,10)=(D(:,10)./K2O)*2; %for K2O
MC(:,11)=D(:,11)./F; %for F

f=2./(MC(:,1)+MC(:,2)+MC(:,3)+MC(:,5)+MC(:,6)+MC(:,7)); %normalization factor

%% normalize moles of cations

[m,n]=size(D); %finds the x and y size of the input data matrix
NMC=zeros(size(D)); %creates a matrix of zeroes the size of the input data
NMC(:,1)=MC(:,1).*f(:); %SiO2
NMC(:,2)=MC(:,2).*f(:); %TiO2
NMC(:,3)=MC(:,3).*f(:); %Al2O3
NMC(:,4)=MC(:,4).*f(:); %Y2O3
NMC(:,5)=MC(:,5).*f(:); %Fe2O3
NMC(:,6)=MC(:,6).*f(:); %MnO
NMC(:,7)=MC(:,7).*f(:); %MgO
NMC(:,8)=MC(:,8).*f(:); %CaO
NMC(:,9)=MC(:,9).*f(:); %Na2O
NMC(:,10)=MC(:,10).*f(:); %K2O
NMC(:,11)=MC(:,11).*f(:); %F

%% structural formula

%cations
StrctFrm(:,1)=NMC(:,1); %Si4+ (T)
StrctFrm(:,2)=1-NMC(:,1); % Ti4+ (T)

StrctFrm(:,3)=NMC(:,2)-StrctFrm(:,2); %Ti4+ (Oct) 
StrctFrm(:,4)=NMC(:,3); %Al3+ (Oct)
StrctFrm(:,5)=NMC(:,5); %Fe3+ (Oct)
StrctFrm(:,6)=NMC(:,6); %Mn (Oct)
StrctFrm(:,7)=NMC(:,7); %Mg (Oct)

StrctFrm(:,8)=NMC(:,8); %Ca (Dec)
StrctFrm(:,9)=NMC(:,4); %Y3+ (Dec)
StrctFrm(:,10)=NMC(:,9); %Na+ (Dec)
StrctFrm(:,11)=NMC(:,10); %K+ (Dec)

%cation sum
StrctFrm(:,12)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3)+StrctFrm(:,4)+StrctFrm(:,5)+StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)+StrctFrm(:,9)+StrctFrm(:,10)+StrctFrm(:,11);

%anions
StrctFrm(:,13)=NMC(:,11); %F 
StrctFrm(:,14)=StrctFrm(:,5)+StrctFrm(:,4)-StrctFrm(:,13); %OH
StrctFrm(:,15)=5-(10-(2.*StrctFrm(:,8)+StrctFrm(:,11)+StrctFrm(:,10)+2.*StrctFrm(:,7)+2.*StrctFrm(:,6)+3.*StrctFrm(:,5)+3.*StrctFrm(:,4)+4.*StrctFrm(:,3)+4.*StrctFrm(:,2)+4.*StrctFrm(:,1)));

%anion sum
StrctFrm(:,16)=StrctFrm(:,13)+StrctFrm(:,14)+StrctFrm(:,15);

XTtn(:,1)=StrctFrm(:,3)./(StrctFrm(:,3)+StrctFrm(:,4)+StrctFrm(:,5)+StrctFrm(:,6)+StrctFrm(:,7));





