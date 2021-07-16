function APFU=Amph_KA(D,W)

% Normalize to [Si + Al + Ti + Fe + Mn + Mg + Ca + Na] = 15 cations

[m,n]=size(D); %finds the x and y size of the input data matrix

%moles of cations
MC=D./(W); 

%moles of cation units
Cat(:,1)=MC(:,1); %Si
Cat(:,2)=MC(:,2); %Ti
Cat(:,3)=MC(:,3).*2; %Al
Cat(:,4)=MC(:,4).*2; %Cr
Cat(:,5)=0; %Fe3+ (blank)
Cat(:,6)=MC(:,5); %Fe2+
Cat(:,7)=MC(:,6); %Mn
Cat(:,8)=MC(:,7); %Mg
Cat(:,9)=MC(:,8); %Ca
Cat(:,10)=MC(:,9).*2; %Na
Cat(:,11)=MC(:,10).*2; %K

for c=1:m
    Norm(c,1)=15./sum(Cat(c,1:1:10)); %Normalize to 15 cations
end

for c=1:m
    N_Cat(c,:)=Cat(c,:).*Norm(c,:); %normalized cation units
end


%Oxygen units
O2(:,1)=N_Cat(:,1).*2; %Si
O2(:,2)=N_Cat(:,2).*2; %Ti
O2(:,3)=N_Cat(:,3).*(3/2); %Al
O2(:,4)=N_Cat(:,4).*(3/2); %Cr
O2(:,5)=N_Cat(:,5).*(3/2); %Fe3+ (blank)
O2(:,6)=N_Cat(:,6); %Fe2+
O2(:,7)=N_Cat(:,7); %Mn
O2(:,8)=N_Cat(:,8); %Mg
O2(:,9)=N_Cat(:,9); %Ca
O2(:,10)=N_Cat(:,10)./2; %Na
O2(:,11)=N_Cat(:,11)./2; %K

O2total=sum(O2,2);

APFU=N_Cat; 

%correction for electroneutrality 
for c=1:m
    if (23-O2total(c,1)) > 0;
        if O2(c,6) > 2.*(23-O2total(c,1));
            APFU(c,5)=2.*(23-O2total(c,1));
            APFU(c,6)=N_Cat(c,6)-APFU(c,5);
        else
            APFU(c,5)=N_Cat(c,6);
            APFU(c,6)=0;
        end
    else
        APFU(c,5)=APFU(c,5);
        APFU(c,6)=APFU(c,6);
    end
end

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

end