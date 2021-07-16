%%calculates amphibole formula following Hawthorne et al. (2012) and
%Leake et al., (1997)

clear all
close all

%IMPORTANT!!!! 
%1. Assumes OH = 2 and Cl, F, and Li are not included.
%2. DO NOT USE on amphibole with a significant richterite component (See max
%   Fe3+ criteria below)

%Criteria for minimum Fe3+ measurement:
%Criteria 1-1: Amphibole cannot have more than 8 APFU Si
%Criteria 2-1: Amphibole cannot accomodate more than 16 cations
%Criteria 3-1: This criteria assumes that Ca cannot go into the A site. This
%criteria is NOT valid for all amphiboles. For example, amphiboles in
%marbles may contnain Ca in the A site. Apply with caution!!!

%Critera for maximum Fe3+ measurement:
%Criteria 1-2: Not always correct, richterite may contain Ti as a T cation,
%and if Ti is suspected in in the T site then this criteria should not be
%used. 
%Criteria 2-2: Assumes that K does not occur as a B cation, this is not
%always correct. K can occur as a B cation in richterite. 
%Criteria 3-2: This criteria can be wrong if there is Li in the structure if
%significant Fe, Mn, and Mg occur in the B site. 

clear all;
close all;

D=load('amphibole_comp.txt'); %EPMA data is loaded from a text file here

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

%Output (Strct_Frm)
%column1: Si (T)
%column2: Al (T)
%column3: T total
%column4: Al (C)
%column5: Ti (C)
%column6: Cr (C)
%column7: Fe3+ (C)
%column8: Mg (C)
%column9: Fe2+ (C)
%column10: Mn (C)
%column11: C total
%column12: Mg (B)
%column13: Fe2+ (B)
%column14: Mn (B)
%column15: Ca (B)
%column16: Na (B)
%column17: B total
%column18: Ca (A)
%column19: Na (A)
%column20: K (A)
%column21: cation total

% Molecular weights

SiO2=60.084;
TiO2=79.866;
Al2O3=101.961;
Cr2O3=151.99;
FeO=71.844;
MnO=70.9374;
MgO=40.304;
CaO=56.0774;
Na2O=61.979;
K2O=94.196;

W=[SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O];

[m,n]=size(D); %finds the x and y size of the input data matrix
%% Lower Limits on Fe3+ 

%Atoms per formula unit (APFU) assuming that all Fe is FeO
APFU_Fe=Amph_Fe(D,W); 
StrctFrm_Fe=StrctFrm(APFU_Fe);
Fe(1:m,1)={'All Fe is FeO'};

%APFU assuming 8 Si cations
APFU_Si=Amph_Si(D,W);
StrctFrm_Si=StrctFrm(APFU_Si);
Si(1:m,1)={'Criteria 1-1: 8 Si cations on T'};

%APFU assuming 16 cations (no vacancies on A)
APFU_Afull=Amph_Afull(D,W);
StrctFrm_Afull=StrctFrm(APFU_Afull);
Afull(1:m,1)={'Criteria 2-1: 16 total cations (no vac. on A)'};

%APFU assuming Na in A site only 
APFU_NaA=Amph_NaA(D,W);
StrctFrm_NaA=StrctFrm(APFU_NaA);
NaA(1:m,1)={'Criteria 3-1: Na on A site only'};



%Script to decide which lower limit amphibole compositions to choose
%selection is based on the lowest value normalization factor for each of
%the three criteria (1-1 to 3-1, see above) 

for c=1:m
    if (APFU_Fe(c,12) > 1) && (APFU_Fe(c,13) > 1) && (APFU_Fe(c,14) > 1) %if the normalization factors are all greater than 1, then all Fe is assumed to be FeO
        StrctFrm_low(c,:)=StrctFrm_Fe(c,:);
        APFU_low(c,:)=APFU_Fe(c,1:1:11);
        low_check(c,:)=Fe(c,1);
    else
        if (APFU_Fe(c,12) < APFU_Fe(c,13)) && (APFU_Fe(c,12) < APFU_Fe(c,14)) 
            StrctFrm_low(c,:)=StrctFrm_Si(c,:);
            APFU_low(c,:)=APFU_Si(c,1:1:11);
            low_check(c,:)=Si(c,1);
        else
            if (APFU_Fe(c,13) < APFU_Fe(c,12)) && (APFU_Fe(c,13) < APFU_Fe(c,14)) 
                StrctFrm_low(c,:)=StrctFrm_Afull(c,:);
                APFU_low(c,:)=APFU_Afull(c,1:1:11);
                low_check(c,:)=Afull(c,1);
            else
                if (APFU_Fe(c,14) < APFU_Fe(c,12)) && (APFU_Fe(c,14) < APFU_Fe(c,13)) 
                    StrctFrm_low(c,:)=StrctFrm_NaA(c,:);
                    APFU_low(c,:)=APFU_NaA(c,1:1:11);
                    low_check(c,:)=NaA(c,1);
                end
            end
        end
    end
end
 
 
%% Upper Limits on Fe3+ 

%APFU assuming all Fe is Fe2O3
APFU_Fe2O3=Amph_Fe2O3(D,W);
StrctFrm_Fe2O3=StrctFrm(APFU_Fe2O3);
Fe2O3(1:m,1)={'All Fe is Fe2O3'};
O2_Nfact=APFU_Fe2O3(:,17)./APFU_Fe(:,18); %Normalization factor for all Fe as Fe2O3 relative to all Fe as FeO

%APFU assuming Al + Si in T site only 
APFU_SiAlT=Amph_SiAlT(D,W);
StrctFrm_SiAlT=StrctFrm(APFU_SiAlT);
SiAlT(1:m,1)={'Al + Si normalized to 8 in T site'};

%APFU assuming K only on the A site
APFU_KA=Amph_KA(D,W);
StrctFrm_KA=StrctFrm(APFU_KA);
KA(1:m,1)={'K only on the A site, fully occupied T, B, & C sites'};

%APFU assuming fully occupied T & B sites
APFU_CaNaB=Amph_CaNaB(D,W);
StrctFrm_CaNaB=StrctFrm(APFU_CaNaB);
CaNaB(1:m,1)={'Fully occupied T & B sites'};


 %Choses the upper limit for Fe3+
 %The upper limit is chosen using the lowest value of normalization scheme
 %for the three critera (1-2 to 3-2, see above)

 
 for c=1:m
     if (O2_Nfact(c,1) > APFU_Fe(c,15)) && (O2_Nfact(c,1) > APFU_Fe(c,16)) && (O2_Nfact(c,1) > APFU_Fe(c,17))
         StrctFrm_hi(c,:)=StrctFrm_Fe2O3(c,:); 
         APFU_hi(c,:)=APFU_Fe2O3(c,1:1:11);
         hi_check(c,:)=Fe2O3(c,1);
     else 
         if (APFU_Fe(c,15) > O2_Nfact(c,1)) && (APFU_Fe(c,15) > APFU_Fe(c,16)) && (APFU_Fe(c,15) > APFU_Fe(c,17))
             StrctFrm_hi(c,:)=StrctFrm_SiAlT(c,:); 
             APFU_hi(c,:)=APFU_SiAlT(c,1:1:11);
             hi_check(c,:)=SiAlT(c,1);
         else
             if (APFU_Fe(c,16) > O2_Nfact(c,1)) && (APFU_Fe(c,16) > APFU_Fe(c,15)) && (APFU_Fe(c,16) > APFU_Fe(c,17))
                 StrctFrm_hi(c,:)=StrctFrm_KA(c,:); 
                 APFU_hi(c,:)=APFU_KA(c,1:1:11);
                 hi_check(c,:)=KA(c,1);
             else
                 if (APFU_Fe(c,17) > O2_Nfact(c,1)) && (APFU_Fe(c,17) > APFU_Fe(c,15)) && (APFU_Fe(c,17) > APFU_Fe(c,16))
                     StrctFrm_hi(c,:)=StrctFrm_CaNaB(c,:);
                     APFU_hi(c,:)=APFU_CaNaB(c,1:1:11);
                     hi_check(c,:)=CaNaB(c,1);
                 end
             end
         end
     end
 end
    
%% Calculates Mean of Min and Max Fe3+ estimates

Strct_Frm=(StrctFrm_hi+StrctFrm_low)./2;
APFU=(APFU_hi+APFU_low)./2;

%uncomment these if you want to output and plot the maximum values 
%Strct_Frm=(StrctFrm_hi);
%APFU=(APFU_hi);

%% Data Plotting 

Amp_Plot(:,1)=APFU(:,5)./(APFU(:,5)+APFU(:,6)); %Fe3+/Fetotal
Amp_Plot(:,2)=sum(APFU(:,5:1:6),2); %Fetotal
Amp_Plot(:,3)=Strct_Frm(:,15)./(Strct_Frm(:,15)+Strct_Frm(:,16)); %Ca/(Ca+Na) in B
Amp_Plot(:,4)=Strct_Frm(:,19)+Strct_Frm(:,20)+2.*Strct_Frm(:,18); %A(Na + K + 2Ca), Ca=0 
Amp_Plot(:,5)=Strct_Frm(:,4)+Strct_Frm(:,7)+2.*Strct_Frm(:,5); %C(Al + Fe3+ +2Ti)
Amp_Plot(:,6)=Strct_Frm(:,8)./(Strct_Frm(:,8)+Strct_Frm(:,9)+Strct_Frm(:,13)); %XMg
Amp_Plot(:,7)=Strct_Frm(:,9)./(Strct_Frm(:,9)+Strct_Frm(:,8)+Strct_Frm(:,10)); %Fe2+/(Fe2+ + Mg + Mn) in C
Amp_Plot(:,8)=Strct_Frm(:,7)./(Strct_Frm(:,7)+Strct_Frm(:,4)+Strct_Frm(:,5)); %Fe3+/(Fe3+ + Al + Ti) in C
Amp_Plot(:,9)=Strct_Frm(:,1); %Si (T)

figure('Name','Fe3+/Fetotal vs Fetotal');
xlim([0 5])
hold on
ylim([0 1])
hold on
xlabel('Fetotal (apfu)')
ylabel('Fe3+/Fetotal')
hold on
plot(Amp_Plot(:,2),Amp_Plot(:,1),'o')

%Plots for Ca amphiboles
figure('Name','Calcium Amphiboles');
xlim([0 2])
hold on
ylim([0 1])
hold on
plot([0 2],[0.5 0.5],'k')
hold on
plot([0.5 0.5],[0.0 1.0],'k')
hold on
plot([1.5 1.5],[0.0 1.0],'k')
hold on
text(0.25,0.25,'Tremolite','FontSize',14,'HorizontalAlignment','center')
text(0.25,0.75,'Edenite','FontSize',14,'HorizontalAlignment','center')
text(1,0.25,'Magnesio-Hornblende','FontSize',14,'HorizontalAlignment','center')
text(1,0.75,'Pargasite','FontSize',14,'HorizontalAlignment','center')
text(1.75,0.75,'Sadanagaite','FontSize',14,'HorizontalAlignment','center')
text(1.75,0.25,'Tschermakite','FontSize',14,'HorizontalAlignment','center')
xlabel('C(Al + Fe3+ + 2Ti) apfu')
ylabel('A(Na + K + 2Ca) apfu')
hold on

for c=1:m
    if (Strct_Frm(c,16)./Strct_Frm(c,17)) <= 0.25
        plot(Amp_Plot(c,5),Amp_Plot(c,4),'o')
        hold on
    end
end

figure('Name','Calcium Amphiboles 2');
xlim([5.5 8])
hold on
ylim([0 1])
hold on
plot([5.5 8],[0.5 0.5],'k')
hold on
plot([6.5 6.5],[0.0 1.0],'k')
hold on
plot([7.5 7.5],[0.0 1.0],'k')
hold on
plot([7.5 8],[0.9 0.9],'k')
hold on
text(7.0,0.25,'Ferro-hornblende','FontSize',14,'HorizontalAlignment','center')
text(6.0,0.25,'Ferro-tschermakite','FontSize',14,'HorizontalAlignment','center')
text(7.0,0.75,'Magnesio-Hornblende','FontSize',14,'HorizontalAlignment','center')
text(6.0,0.75,'Tschermakite','FontSize',14,'HorizontalAlignment','center')
text(7.75,0.25,'Ferro-actinolite','FontSize',14,'HorizontalAlignment','center')
text(7.75,0.75,'Actinolite','FontSize',14,'HorizontalAlignment','center')
text(7.75,0.95,'Tremolite','FontSize',14,'HorizontalAlignment','center')
xlabel('Si apfu')
ylabel('Mg/(Mg + Fe2+)')
hold on

for c=1:m
    if (Strct_Frm(c,16)./Strct_Frm(c,17)) <= 0.25

        plot(Amp_Plot(c,9),Amp_Plot(c,6),'o')
        hold on
    end
end


%Na-Ca amphiboles
figure('Name','Sodium-Calcium Amphiboles');
xlim([0 2])
hold on
ylim([0 1])
hold on
plot([0 2],[0.5 0.5],'k')
hold on
plot([0.5 0],[0 0.5],'k')
hold on
plot([0.5 0.5],[0.5 1.0],'k')
hold on
plot([1.5 1.5],[0.0 1.0],'k')
hold on
text(0.25,0.75,'Richterite','FontSize',14,'HorizontalAlignment','center')
text(1,0.25,'Winchite','FontSize',14,'HorizontalAlignment','center')
text(1,0.75,'Katophorite','FontSize',14,'HorizontalAlignment','center')
text(1.75,0.75,'Taramite','FontSize',14,'HorizontalAlignment','center')
text(1.75,0.25,'Barroisite','FontSize',14,'HorizontalAlignment','center')
xlabel('C(Al + Fe3+ + 2Ti) apfu')
ylabel('A(Na + K + 2Ca) apfu')
hold on

for c=1:m
    if (Strct_Frm(c,16)./Strct_Frm(c,17)) < 0.75 && (Strct_Frm(c,16)./Strct_Frm(c,17)) > 0.25
       plot(Amp_Plot(c,5),Amp_Plot(c,4),'o')
       hold on
    end
end

%Na amphiboles

figure('Name','Sodium Amphiboles');
xlim([0 2])
hold on
ylim([0 1])
hold on
plot([1 2],[0.5 0.5],'k')
hold on
plot([1.5 1.5],[0.5 1.0],'k')
hold on
plot([0.5 1.5],[1 0],'k')
hold on
text(1.15,0.75,'Eckermannite','FontSize',14,'HorizontalAlignment','center')
text(1.75,0.75,'Nyboite','FontSize',14,'HorizontalAlignment','center')
text(1.75,0.25,'Glaucophane','FontSize',14,'HorizontalAlignment','center')
xlabel('C(Al + Fe3+ + 2Ti) apfu')
ylabel('A(Na + K + 2Ca) apfu')
hold on

for c=1:m
    if (Strct_Frm(c,16)./Strct_Frm(c,17)) >= 0.75
        plot(Amp_Plot(c,5),Amp_Plot(c,4),'o')
        hold on
    end
end

figure('Name','Sodium Amphiboles 2');
xlim([0 1])
hold on
ylim([0 1])
hold on
plot([0 1],[0.5 0.5],'k')
hold on
plot([0.5 0.5],[0 1],'k')
hold on
text(0.25,0.25,'Glaucophane','FontSize',14,'HorizontalAlignment','center')
text(0.25,0.75,'Ferro-Glaucophane','FontSize',14,'HorizontalAlignment','center')
text(0.75,0.25,'Magnesio-Riebeckite','FontSize',14,'HorizontalAlignment','center')
text(0.75,0.75,'Riebeckite','FontSize',14,'HorizontalAlignment','center')
xlabel('Fe3+/(Fe3+ + Al + Ti)')
ylabel('Fe2+/(Fe2+ + Mg + Mn)')
hold on

for c=1:m
    if (Strct_Frm(c,16)./Strct_Frm(c,17)) >= 0.75
        plot(Amp_Plot(c,8),Amp_Plot(c,7),'o')
        hold on
    end
end

