function Frm=StrctFrm(APFU)

[m,n]=size(APFU); %finds the x and y size of the input data matrix

%T SITE

%Si 
for c=1:m
    if APFU(c,1)<8.000
    Frm(c,1)=APFU(c,1);
    else
    Frm(c,1)=8;
    end
end

%Al(T)
for c=1:m
    if (8-Frm(c,1))<APFU(c,1)
        Frm(c,2)=8-Frm(c,1);
    else
        Frm(c,2)=APFU(c,3);
    end
end


%Sum of T site
Frm(:,3)=Frm(:,1)+Frm(:,2);

%C site 

%Al(C)
for c=1:m
    if Frm(c,2)<APFU(c,3)
        Frm(c,4)=APFU(c,3)-Frm(c,2);
    else
        Frm(c,4)=0;
    end
end


Frm(:,5)=APFU(:,2); %Ti (C)
Frm(:,6)=APFU(:,4); %Cr (C)
Frm(:,7)=APFU(:,5); %Fe3+ (C)

%Mg (C)
for c=1:m
    if sum(Frm(c,4:1:7))+APFU(c,8)<5
        Frm(c,8)=APFU(c,8); 
    else
        Frm(c,8)=5-sum(Frm(c,4:1:7));
    end
end
        

%Fe2+ (C)
for c=1:m
    if sum(Frm(c,4:1:8),2)+APFU(c,6)<5
        Frm(c,9)=APFU(c,6);
    else
        Frm(c,9)=5-sum(Frm(c,4:1:8));
    end
end


%Mn (C)
for c=1:m
    if sum(Frm(c,4:1:9))+APFU(c,7)<5
        Frm(c,10)=APFU(c,7);
    else
        Frm(c,10)=5-sum(Frm(c,4:1:9));
    end
end

%Sum of C site
for c=1:m
    Frm(c,11)=sum(Frm(c,4:1:10));
end

%B Site 

Frm(:,12)=APFU(:,8)-Frm(:,8); %Mg (B)
Frm(:,13)=APFU(:,6)-Frm(:,9); %Fe2+ (B)
Frm(:,14)=APFU(:,7)-Frm(:,10); %Mn(B)
Frm(:,15)=APFU(:,9); %Ca (B)

%Ca (B)
for c=1:m
    if sum(Frm(c,12:1:14))+APFU(c,9)>2
        Frm(c,15)=2-sum(Frm(c,12:1:14));
    else
        Frm(c,15)=APFU(c,9);
    end
end


%Na (B)
for c=1:m
    if sum(Frm(c,12:1:15))+APFU(c,10)>2
        Frm(c,16)=2-sum(Frm(c,12:1:15));
    else
        Frm(c,16)=APFU(c,10);
    end
end

Frm(:,17)=Frm(:,12)+Frm(:,13)+Frm(:,14)+Frm(:,15)+Frm(:,16);

%A site 
Frm(:,18)=APFU(:,9)-Frm(:,15); %Ca (A)
Frm(:,19)=APFU(:,10)-Frm(:,16); %Na (A)
Frm(:,20)=APFU(:,11); %K (B)

%Cation total
Frm(:,21)=Frm(:,20)+Frm(:,19)+Frm(:,18)+Frm(:,17)+Frm(:,11)+Frm(:,3);

end