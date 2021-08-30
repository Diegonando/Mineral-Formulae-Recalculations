% Sulfide recalculation

%Copyright (c) 2021 Jesse Walters
%The DOI for this release is: 10.5281/zenodo.5110201

%Normalization parameters
cat=2; %number of cations
an=2; %number of anions

D=load('sulfide.txt'); %EPMA data text file goes here

%order of elements (wt %) in input file:
%column 1: S
%column 2: Co
%column 3: Cu
%column 4: As
%column 5: Fe
%column 6: Ni

%% Convert to moles

%average atomic masses
S=32.065;
Co=58.933195;
Cu=63.546;
As=74.9216;
Fe=55.845;
Ni=58.6934;

[m,n]=size(D); %finds the x and y size of the input data matrix
mol=zeros(size(D)); %creates a matrix of zeroes the size of the input data

mol(:,1)=D(:,1)./S; %moles of S
mol(:,2)=D(:,2)./Co; %moles of Co
mol(:,3)=D(:,3)./Cu; %moles of Cu
mol(:,4)=D(:,4)./As; %moles of As
mol(:,5)=D(:,5)./Fe; %moles of Fe
mol(:,6)=D(:,6)./Ni; %moles of Ni

%% normalize to cations

NF=cat./(mol(:,2)+mol(:,3)+mol(:,4)+mol(:,5)+mol(:,6)); %normalization factor

[m,n]=size(D); %finds the x and y size of the input data matrix
CAPFU=zeros(size(D)); %creates a matrix of zeroes the size of the input data
CAPFU(:,1)=mol(:,1).*NF; %S
CAPFU(:,2)=mol(:,2).*NF; %Co
CAPFU(:,3)=mol(:,3).*NF; %Cu
CAPFU(:,4)=mol(:,4).*NF; %As
CAPFU(:,5)=mol(:,5).*NF; %Fe
CAPFU(:,6)=mol(:,6).*NF; %Ni

%% normalize to anions

NF=an./(mol(:,1)); %normalization factor

[m,n]=size(D); %finds the x and y size of the input data matrix
AAPFU=zeros(size(D)); %creates a matrix of zeroes the size of the input data
AAPFU(:,1)=mol(:,1).*NF; %S
AAPFU(:,2)=mol(:,2).*NF; %Co
AAPFU(:,3)=mol(:,3).*NF; %Cu
AAPFU(:,4)=mol(:,4).*NF; %As
AAPFU(:,5)=mol(:,5).*NF; %Fe
AAPFU(:,6)=mol(:,6).*NF; %Ni