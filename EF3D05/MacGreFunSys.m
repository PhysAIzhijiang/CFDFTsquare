function [A0temp, A1temp]=MacGreFunSys(key,E,H00,invH01,H10,sigmaL,wid,layer,ksi,dense,rd)
%This is the lead sysmtem Green function on the square lattice of quantum well
%    iteration of the center regions
%layer=disorder;%ksi=8;
%key idea try to avoid  matrix inverse!!!!!
 %initial random
    x1=clock;
    %x1=rand+clock;
    x=round(x1(6)*10000);
    for i=1:x
        rand();
    end
if key==1
   Hw=rd*(rand(2*wid*layer,1)-0.5);
   Hw([1:4:4*wid*layer-3 3:4:4*wid*layer-1  2:4:4*wid*layer-2 4:4:4*wid*layer ])=[Hw(1:wid*layer)  Hw(1:wid*layer)  Hw(wid*layer+1:2*wid*layer)  Hw(wid*layer+1:2*wid*layer)];
    %{
    for i=1:wid*layer
        h=rd*rand(1,1)-rd/2;
        Hw(2*i-1:2*i)=[h,h];
    end
    %}
elseif key==2
    %[1,1,1,1] is relative onsite energy of each orbital
    %Hw=rd*kron(DisorderSi( wid,layer,ksi,dense ),[1,1]);
    Hw=rd*kron(DisorderSi2( wid,layer,ksi,dense ),(rand(2,1)-0.5)*2);
elseif key==3
    Hw=zeros(2*wid*layer,1);
end
A0temp=E*eye(4*wid)-H00-diag(Hw(1:4*wid))-sigmaL;%imaginary part neglected since self-energy plugged in
A1temp=(E*eye(4*wid)-H00-diag(Hw(4*wid+1:8*wid)))*invH01*A0temp-H10;
Astemp=1;%initialization of As
S_num=4;
for i=3:layer    
 A3temp=A1temp;      
A1temp =(E*eye(4*wid)-H00-diag( Hw( (i-1)*4*wid+ (1:4*wid)) ) )*invH01*A1temp-H10*invH01*A0temp;
A0temp=A3temp;%effective propogator form the first layer to Mth layer
if mod(i,S_num)==0 & (layer-i)>=2
    A3temp=inv(A1temp);
    Astemp=Astemp*A3temp;
    A1temp=A1temp*A3temp;
    A0temp=A0temp*A3temp;
end %record of As and iteration of A    
end
A3temp=inv(A1temp);
A0temp=invH01*A0temp*A3temp; %avoid defining a new variable, so borrow
A1temp=Astemp*A3temp; %avoid defining a new variable, so borrow
end

