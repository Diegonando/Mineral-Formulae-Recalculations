function APFU=Amph_Fe2O3(D,W)

% All ferrous iron, 23 O2 normalization

[m,n]=size(D); %finds the x and y size of the input data matrix

%moles of cations
D_Fe2O3=D;
D_Fe2O3(:,5)=D(:,5)*(159.688./(2*W(:,5))); %convert FeO to Fe2O3

W(:,5)=159.688; %change to Fe2O3 

%moles of cations
MC=D_Fe2O3./(W); 

%Moles of O2
O2=zeros(size(D)); %creates a matrix of zeroes the size of the input data
O2(:,1)=MC(:,1)*2; %SiO2
O2(:,2)=MC(:,2)*2; %TiO2
O2(:,3)=MC(:,3)*3; %Al2O3
O2(:,4)=MC(:,4)*3; %Cr2O3
O2(:,5)=MC(:,5)*3; %Fe2O3
O2(:,6)=MC(:,6); %MnO
O2(:,7)=MC(:,7); %MgO
O2(:,8)=MC(:,8); %CaO
O2(:,9)=MC(:,9); %Na2O
O2(:,10)=MC(:,10); %K2O

O2total=sum(O2,2); %O2 totals
O2_N=23./O2total; %O2 normalization factor

%moles of cations
N_Ox=O2.*O2_N;

%atoms per formula unit
APFU=zeros(size(D));
APFU(:,1)=N_Ox(:,1)./2; %SiO2
APFU(:,2)=N_Ox(:,2)./2; %TiO2
APFU(:,3)=N_Ox(:,3).*(2/3); %Al2O3
APFU(:,4)=N_Ox(:,4).*(2/3); %Cr2O3
APFU(:,5)=N_Ox(:,5)*(2/3); %Fe2O3
APFU(:,6)=0; 
APFU(:,7)=N_Ox(:,6); %MnO
APFU(:,8)=N_Ox(:,7); %MgO
APFU(:,9)=N_Ox(:,8); %CaO
APFU(:,10)=N_Ox(:,9).*2; %Na2O
APFU(:,11)=N_Ox(:,10).*2; %K2O 

for c=1:m
    APFU(c,12)=sum(APFU(c,1:1:11)); %Criteria 2-1: Should be <= 16 APFU
end

for c=1:m
    APFU(c,13)=sum(APFU(c,1:1:9)); %Criteria 3-1: Should be <=15 APFU
end

for c=1:m
    APFU(c,14)=APFU(c,1)+APFU(c,3); % Criteria 1-2: Si + Al = 8
end

for c=1:m
    APFU(c,15)=sum(APFU(c,1:1:10)); %Criteria 2-2: Should be equal to 15 APFU
end

for c=1:m
    APFU(c,16)=sum(APFU(c,1:1:8)); %Criteria 3-2: Should be equal to 13 APFU
end

for c=1:m
    APFU(c,17)=O2_N(c); 
end

end

