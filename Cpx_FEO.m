%%calculates CPX formula, all Fe is FeO

clear all;
close all;

D=load('CPX_EPMA.txt'); %EPMA data is loaded from a text file here

%input wt % oxide in the following order
%column1: SiO2
%column2: TiO2
%column3: Al2O3
%column4: Cr2O3
%column5: FeO
%column6: MnO
%column7: MgO
%column8: CaO
%column9: Na2O
%column10: K2O

%OUTPUT: cations (APFU)
%column1: Si
%column2: Ti
%column3: Al
%column4: Fe3+
%column5: Cr
%column6: Fe2+
%column7: Mn
%column8: Mg
%column9: Ca
%column10: Na
%column11: K
%column12: total
%column13: O2 deficiency 

%OUTPUT: Structural formula (APFU)
%column1: Si (T)
%column2: Al (T)
%column3: Fe3+ (T)
%column4: sum of T site
%column5: Ti (M1)
%column6: Al (M1)
%column7: Fe3+ (M1)
%column8: Cr (M1)
%column9: Mn (M1)
%column10: Mg (M1)
%column11: Fe2+ (M1)
%column12: sum of M1 site
%column13: Mg (M2)
%column14: Fe2+ (M2)
%column15: Ca (M2)
%column16: Na (M2)
%column17: K (M2)
%column18: Sum of M2 site

%Output: endmembers
%column1: XWo
%column2: XFs
%column3: XEn
%column4: XJd
%column5: XAcm
%column6: XKos
%column7: XQuad

cat=4.0; %cations per formula unit
Opfu=6.0; %oxygens per formula unit


%% Molecular weights

SiO2=60.084;
TiO2=79.866;
Al2O3=101.961;
Cr2O3=151.99;
FeO=71.844;
MnO=70.937;
MgO=40.304;
CaO=56.077;
Na2O=61.979;
K2O=94.196;

W=[SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O];

%% Calculate cations units

[m,n]=size(D); %finds the x and y size of the input data matrix
MC=zeros(size(D)); %creates a matrix of zeroes the size of the input data
MC(:,1)=D(:,1)./W(:,1); %for SiO2
MC(:,2)=D(:,2)./W(:,2); %for TiO2
MC(:,3)=(D(:,3)./W(:,3)).*2; %for Al2O3
MC(:,4)=(D(:,4)./W(:,4)).*2; %for Cr2O3
MC(:,5)=D(:,5)./W(:,5); %for FeO
MC(:,6)=D(:,6)./W(:,6); %for MnO
MC(:,7)=D(:,7)./W(:,7); %for MgO
MC(:,8)=D(:,8)./W(:,8); %for CaO
MC(:,9)=(D(:,9)./W(:,9)).*2; %for Na2O
MC(:,10)=(D(:,10)./W(:,10)).*2; %for K2O

MCnormfact=zeros(m,1); %creates a zeromatrix for the totals of the mole cation matrix 
MCnormfact=cat./sum(MC,2); %normalization factor

MCnorm=MCnormfact.*MC; %creates a matrix of normalized cations

%% Calculate Oxygen Units

O2=zeros(size(D));
O2(:,1)=MCnorm(:,1).*2; %for SiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4).*(3/2); %for Cr2O3
O2(:,5)=MCnorm(:,5); %for FeO
O2(:,6)=MCnorm(:,6); %for MnO
O2(:,7)=MCnorm(:,7); %for MgO
O2(:,8)=MCnorm(:,8); %for CaO
O2(:,9)=MCnorm(:,9)./2; %for Na2O
O2(:,10)=MCnorm(:,10)./2; %for K2O

O2total=sum(O2,2); %O2 totals

%% Atoms pfu

APFU=zeros(m,n+3); %matrix of zeroes to be filled, n+2 for Fe3+ and total

APFU(:,1)=MCnorm(:,1); %for Si
APFU(:,2)=MCnorm(:,2); %for Ti
APFU(:,3)=MCnorm(:,3); %for Al
APFU(:,5)=MCnorm(:,4); %for Cr
APFU(:,6)=MCnorm(:,5); %for Fe
APFU(:,7)=MCnorm(:,6); %for Mn
APFU(:,8)=MCnorm(:,7); %for Mg
APFU(:,9)=MCnorm(:,8); %for Ca
APFU(:,10)=MCnorm(:,9); %for Na
APFU(:,11)=MCnorm(:,10); %for K

APFU(:,12)=sum(APFU,2); %calculations the total, which should be 4

%XMg
XMg=APFU(:,8)./(APFU(:,8)+APFU(:,6)); 

%% structural formula calculation

StrctFrm=zeros(m,n+8);%creates a matrix to be filled

%T SITE
%Si 
for c=1:m
    if APFU(c,1)<2.000
    StrctFrm(c,1)=APFU(c,1);
    else
    StrctFrm(c,1)=2;
    end
end

%Al(T)
for c=1:m
    if 2-StrctFrm(c,1)<APFU(c,1)
        StrctFrm(c,2)=2-StrctFrm(c,1);
    else
        StrctFrm(c,2)=APFU(c,3);
    end
end

%Fe3+(T)
for c=1:m
    if 2-StrctFrm(c,1)-StrctFrm(c,2)>0
         StrctFrm(c,3)=2-StrctFrm(c,1)-StrctFrm(c,2);
    else
        StrctFrm(c,3)=0;
    end
end

%Sum of T site
StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3);

%M1 SITE
%Ti
StrctFrm(:,5)=APFU(:,2);

%Al(M1)
for c=1:m
    if StrctFrm(c,2)<APFU(c,3)
        StrctFrm(c,6)=APFU(c,3)-StrctFrm(c,2);
    else
        StrctFrm(c,6)=0;
    end
end

%Fe3+ (M1)
for c=1:m
    if StrctFrm(c,3)<APFU(c,4)
        StrctFrm(c,7)=APFU(c,4)-StrctFrm(c,3);
    else
        StrctFrm(c,7)=0;
    end
end

%Cr3+ (M1)
StrctFrm(:,8)=APFU(:,5);

%Mn (M1)
StrctFrm(:,9)=APFU(:,7);

%Mg (M1)
for c=1:m
    if XMg(c).*(1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,5))<APFU(c,8)
        StrctFrm(c,10)=XMg(c).*(1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,5));
    else
        StrctFrm(c,10)=APFU(c,8);
    end
end

%Fe2+ (M1)
for c=1:m
    if 1-(StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10)+StrctFrm(c,5))<=APFU(c,6)
        StrctFrm(c,11)=1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,10)-StrctFrm(c,5);
    else
        StrctFrm(c,11)=APFU(c,6);
    end
end

%Sum of M1 site
StrctFrm(:,12)=StrctFrm(:,5)+StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)+StrctFrm(:,9)+StrctFrm(:,10)+StrctFrm(:,11);

%M2 SITE
%Mg (M2)
for c=1:m
    if APFU(c,8)-StrctFrm(c,9)>0
        StrctFrm(c,13)=APFU(c,8)-StrctFrm(c,10);
    else
        StrctFrm(c,13)=0;
    end
end

%Fe2+ (M2)
for c=1:m
    if APFU(c,6)-StrctFrm(c,11)>0
        StrctFrm(c,14)=APFU(c,6)-StrctFrm(c,11);
    else
        StrctFrm(c,14)=0;
    end
end

%Ca (M2)
StrctFrm(:,15)=APFU(:,9);

%Na (M2)
StrctFrm(:,16)=APFU(:,10);

%K (M2)
StrctFrm(:,17)=APFU(:,11);

%Sum of M2 site
StrctFrm(:,18)=StrctFrm(:,13)+StrctFrm(:,14)+StrctFrm(:,15)+StrctFrm(:,16)+StrctFrm(:,17);

%% end member calculations

%ortho and calcic pyroxene
%Normalization procedure follows Morimoto et al. (1988)
EnWoFs=zeros(m,3);
EnWoFs(:,1)=(APFU(:,8)./(APFU(:,9)+APFU(:,6)+APFU(:,8)+APFU(:,7)+APFU(:,4))); % XEn
EnWoFs(:,2)=(APFU(:,9)./(APFU(:,9)+APFU(:,6)+APFU(:,8)+APFU(:,7)+APFU(:,4))); % XWo
EnWoFs(:,3)=((APFU(:,6)+APFU(:,7)+APFU(:,4))./(APFU(:,9)+APFU(:,6)+APFU(:,8)+APFU(:,7)+APFU(:,4))); % XFs

%Na-Ca endmembers 
Endmembers=zeros(m,7); %matrix to be filled with endmember data
A=zeros(m,6); %A vector to be filled
A(:,1)=APFU(:,9); %Ca
A(:,2)=APFU(:,3); %Al
A(:,3)=APFU(:,4); %Fe3+
A(:,4)=APFU(:,8); %Mg
A(:,5)=APFU(:,6); %Fe2+
A(:,6)=APFU(:,7); %Cr
AT=transpose(A); %transpose of A

M=[2 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 2 0 0 0; 0 2 0 0 0 0; 0 0 0 0 0 1];

X=zeros(6,m);
for c=1:m
    X(:,c)=inv(M)*AT(:,c); %calculates endmembers
end

Xtot=sum(X);%sum of endmembers
Xnorm=X./sum(X); %normalizes the endmembers to 1
XnormT=transpose(Xnorm); %transposes back 

Endmembers(:,1)=XnormT(:,1); % XWo
Endmembers(:,2)=XnormT(:,2); % XFs
Endmembers(:,3)=XnormT(:,3); % XEn 
Endmembers(:,4)=XnormT(:,4); % XJd
Endmembers(:,5)=XnormT(:,5); % XAcm 
Endmembers(:,6)=XnormT(:,6); % XKos 
Endmembers(:,7)=(XnormT(:,1)+XnormT(:,2)+XnormT(:,3))./(XnormT(:,1)+XnormT(:,2)+XnormT(:,3)+XnormT(:,4)+XnormT(:,5)+XnormT(:,6)); % XQuad

%% Plot ternary diagram

%Q-J diagram
figure('Name','Q-J Diagram')
plot([0 2],[2 0],'k')
hold on
plot([0 1.5],[1.5 0],'k')
hold on
plot([0.3 0.4],[1.2 1.6],'k')
hold on
plot([1.2 1.6],[0.3 0.4],'k')
hold on

Y3=APFU(:,8)+APFU(:,6)+APFU(:,9);
X3=APFU(:,10).*2;

plot(X3(:),Y3(:),'o')
hold off
text(0.25,1.8,'Quad','FontSize',14)
text(1.05,1,'Ca-Na','FontSize',14)
text(1.85,0.2,'Na','FontSize',14)

%plots a ternary for Ca-Mg-Fe pyroxenes
figure('Name','Ca-Mg-Fe Pyroxenes');
plot([0 1],[0 0],'k') 
hold on
plot([0 0.5],[0 sqrt(3)/2],'k')
hold on 
plot([0.5 1],[sqrt(3)/2 0],'k')
hold on
plot([0.25 0.75],[0.4330127 0.4330127],'k')
hold on
plot([0.225 0.775],[0.38971143 0.38971143],'k')
hold on
plot([0.1 0.9],[0.17320508 0.17320508],'k')
hold on
plot([0.025 0.975],[0.04330127 0.04330127],'k')
hold on
plot([0.5 0.5],[0.38971143 0.4330127],'k')
hold on
plot([0.5 0.5],[0 0.04330127],'k')
hold on
text(0,0.1,'En','FontSize',14)
text(0.4,0.85,'Wo','FontSize',14)
text(0.95,0.1,'Fs','FontSize',14)

%transforms the data to ternary space
X2=0.5.*(EnWoFs(:,2))+(EnWoFs(:,3));
Y2=(EnWoFs(:,2))*(cos(30*pi()/180));

plot(X2(:),Y2(:),'o')
hold off

%plots a ternary figure for sodic-calcic pyroxenes 
figure('Name','Sodic-Calcic Pyroxenes');
plot([0 1],[0 0],'k')
hold on
plot([0 0.5],[0 sqrt(3)/2],'k')
hold on 
plot([0.5 1],[sqrt(3)/2 0],'k')
hold on
plot([0.4 0.6],[0.69282032 0.69282032],'k')
hold on
plot([0.1 0.9],[0.17320508 0.17320508],'k')
hold on
plot([0.5 0.5],[0 0.69282032],'k')
hold on
text(0,0.1,'Jd','FontSize',14)
text(0.4,0.85,'Quad','FontSize',14)
text(0.95,0.1,'Acm','FontSize',14)

%transforms the data to ternary space 
X1=0.5.*(Endmembers(:,7))+(Endmembers(:,5));
Y1=(Endmembers(:,7))*(cos(30*pi()/180));

plot(X1(:),Y1(:),'o')
hold off 
