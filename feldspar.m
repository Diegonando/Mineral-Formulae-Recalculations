%% feldspar 
%using normalization to 5 cations

clear all;
close all;

D=load('SSP18-1B_fsp.txt'); %EPMA data is loaded from a text file here

%input wt % oxide in the following order
%column1: SiO2
%column2: Al2O3
%column3: FeO
%column4: MnO
%column5: MgO
%column6: CaO
%column7: Na2O
%column8: K2O
%column9: BaO

%% Molecular weights

SiO2=60.084;
Al2O3=101.9612;
FeO=71.8444;
MnO=70.9374;
MgO=40.3044;
CaO=56.0774;
Na2O=61.979;
K2O=94.196;
BaO=153.3264; 

%% moles of cations

[m,n]=size(D); %finds the x and y size of the input data matrix
MC=zeros(size(D)); %creates a matrix of zeroes the size of the input data

MC(:,1)=D(:,1)./SiO2; %for SiO2
MC(:,2)=(D(:,2)./Al2O3)*2; %for Al2O3
MC(:,3)=D(:,3)./FeO; %for FeO
MC(:,4)=D(:,4)./MnO; %for MnO
MC(:,5)=D(:,5)./MgO; %for MgO
MC(:,6)=D(:,6)./CaO; %for CaO
MC(:,7)=(D(:,7)./Na2O).*2; %for Na2O
MC(:,8)=(D(:,8)./K2O).*2; %for K2O
MC(:,9)=D(:,9)./BaO; %for BaO

f=5./(MC(:,1)+MC(:,2)+MC(:,3)+MC(:,5)+MC(:,6)+MC(:,7)+MC(:,8)+MC(:,9)); %normalization factor

%%

APFU=zeros(size(D)); %creates a matrix of zeroes the size of the input data
APFU(:,1)=MC(:,1).*f(:); %SiO2
APFU(:,2)=MC(:,2).*f(:); %Al2O3
APFU(:,3)=MC(:,3).*f(:); %FeO
APFU(:,4)=MC(:,4).*f(:); %MnO
APFU(:,5)=MC(:,5).*f(:); %MgO
APFU(:,6)=MC(:,6).*f(:); %CaO
APFU(:,7)=MC(:,7).*f(:); %Na2O
APFU(:,8)=MC(:,8).*f(:); %K2O
APFU(:,9)=MC(:,9).*f(:); %BaO

%% Oxygen Units
O2=zeros(size(D));
O2(:,1)=APFU(:,1).*2; %for SiO2
O2(:,2)=APFU(:,2).*(3/2); %for Al2O3
O2(:,3)=APFU(:,3); %for FeO
O2(:,4)=APFU(:,4); %for MnO
O2(:,5)=APFU(:,5); %for MgO
O2(:,6)=APFU(:,6); %for CaO
O2(:,7)=APFU(:,7)*(1/2); %for Na2O
O2(:,8)=APFU(:,8); %for K2O
O2(:,9)=APFU(:,9); %for BaO

O2total=sum(O2,2); %O2 totals

%% endmember calculation
Endmembers=zeros(m,3);
Endmembers(:,1)=APFU(:,6)./(APFU(:,6)+APFU(:,7)+APFU(:,8)); %anorthite 
Endmembers(:,2)=APFU(:,7)./(APFU(:,6)+APFU(:,7)+APFU(:,8)); %albite
Endmembers(:,3)=APFU(:,8)./(APFU(:,6)+APFU(:,7)+APFU(:,8)); %orthoclase

%% Ternary plot 
figure('Name','Feldspar Ternary');
plot([0 1],[0 0],'k')
hold on
plot([0 0.5],[0 sqrt(3)/2],'k')
hold on 
plot([0.5 1],[sqrt(3)/2 0],'k')
hold on

%transforms the data to ternary space 
X1=0.5.*(Endmembers(:,1))+(Endmembers(:,3));
Y1=(Endmembers(:,1))*(cos(30*pi()/180));

xlim([0 1])
hold on
ylim([0 1])
hold on

plot(X1(:),Y1(:),'o')
hold off 