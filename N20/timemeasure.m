clear all;
%maxNumCompThreads(8);
%delete(gcp('nocreate'));
%parpool(8);
%H_k=m1*(cos(k0)-cos(kz))*sigmaz...cosk0 to tune m0, not written here thus
%cosk0<1 since m0<0
%    +m2*(2-cos(kx)-cos(ky))*sigmaz...
%    +tx*sin(kx)*sigmax/2+ty*sin(ky)*sigmay/2  %here consider half of the total hamiltonian
tic;eta=1i*10^(-12);
a=1;
N=20;Ny=N;Nx=N;Nz=N;
c=1;EF_Lead=3+eta;
tx=1;ty=1;m1=1;m2=1;cosk0=0.6;
%TD=[0  0.8  3 ];
TD=[0.5];
EF=[0.01 0.5 0.7 0.9 1.1 1.2 1.3 1.7 2 3];lenE =length(EF);
random=[1:15];average=200;%caution:average=1 for random=0
mz=0;
delta=2*pi/a^2/Ny;
%phi_n=[0,0.00001,0.0001,0.001,0.01,1,2,6]*delta;%[0:delta*0.5:5*delta];
phi_n=[0]*delta;
T10=zeros(length(TD),length(phi_n),length(EF),length(random),average);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i6=1:length(TD)
    td=TD(i6);
    i6
sigma0=eye(2);sigmax=[0,1;1,0];
sigmay=[0,-1i;1i,0];sigmaz=[1,0;0,-1];
Vspx= -1i*tx/(2*a)/2;Vspy= -1i*ty/(2*a)/2;Vspd=-1i*td/(2*a)/2;
%Vspx= -1i*tx/(2*a);Vspy= -1i*ty/(2*a);Vspd=-1i*td/(2*a);
V=m1*cosk0/c^2*sigmaz+2*m2/a^2*sigmaz;
%define the hopping matrix on each lattice
Tx=sigmax*Vspx-sigmaz*m2/(2*a^2);
Ty=sigmay*Vspy-sigmaz*m2/(2*a^2);
Tz=-sigmaz*m1/(2*c^2);
DTx=zeros(4,4);DTy=zeros(4,4);DTz=zeros(4,4);DV=zeros(4,4);
%delta=2*pi/a^2/Ny;
DV(1:2,1:2)=V;DV(3:4,3:4)=V;
DTx(1:2,1:2)=Tx;DTx(3:4,3:4)=-sigmax*Vspx-sigmaz*m2/(2*a^2);
DTy(1:2,1:2)=Ty;DTy(3:4,3:4)=Ty;
DTz(1:2,1:2)=Tz;DTz(3:4,3:4)=Tz;
DTx(1,3)=Vspd;DTx(3,1)=Vspd;DTy(1,3)=-1i*Vspd;DTy(3,1)=1i*Vspd;
for i5=1:length(phi_n)   
    i5
    phi0=phi_n(i5);
%DV(1,1)=DV(1,1)+mz(i5);DV(2,2)=DV(2,2)+mz(i5);DV(3,3)=DV(3,3)-mz(i5);DV(4,4)=DV(4,4)-mz(i5);
%DV(1,1)=DV(1,1);DV(2,2)=DV(2,2)+mz(i5);DV(3,3)=DV(3,3)-mz(i5);DV(4,4)=DV(4,4)-mz(i5);
%h00 hopping matrix of intra-line in y-direction;h01x hopping matrix of inter-line in x-direction;
h00=zeros(4*Ny,4*Ny);h01x=zeros(4*Ny,4*Ny);
h01z=zeros(4*Ny,4*Ny);%inter-line in x-direction;
for i=1:Ny
    h00 ( 4*i-3 : 4*i, 4*i-3 : 4*i )=DV;
    h01x ( 4*i-3 : 4*i, 4*i-3 : 4*i )=DTx*exp(-i*phi0*1i*a^2);%magnetic field
    h01z ( 4*i-3 : 4*i, 4*i-3 : 4*i )=DTz;%exp(-i*phi0*1i);
end
for i=1:Ny-1
    h00 ( 4*i-3+4 : 4*i+4, 4*i-3 : 4*i )=DTy; %lower diagnal
    h00 ( 4*i-3 : 4*i, 4*i-3+4 : 4*i+4 )=DTy';%upper diagnal
end
%impose period boundary condition in y-direction
h00(1:4,4*Ny-3:4*Ny)=DTy;h00(4*Ny-3:4*Ny,1:4)=DTy';
%the "hxy_plane" is total Hamiltonian of xy-plane without considering hopping 
%among different planes in z-diretion
%Hz among xy-plane hopping matrix in z-direction
hxy=zeros(4*Nx*Ny,4*Nx*Ny); Hz=zeros(4*Nx*Ny,4*Nx*Ny);
for i=1:Nx
    indx=4*Ny*(i-1);
    hxy( 1+indx:4*Ny+indx,  1+indx:4*Ny+indx )=h00;
    Hz( 1+indx:4*Ny+indx,  1+indx:4*Ny+indx )=h01z;%*exp(-i*1i*phi0);%the magnetic field 
end
for i=1:Nx-1
    indx=4*Ny*(i-1);
    hxy( 1+indx:4*Ny+indx , 4*Ny+1+indx:8*Ny+indx )=h01x;%upper diagnal
    hxy(4*Ny+1+indx:8*Ny+indx,1+indx:4*Ny+indx  )=h01x'; %lower diagnal
end
%impose period boundary condition in x-direction
hxy( 4*Ny*(Nx-1)+1:4*Ny*Nx , 1:4*Ny )=h01x;
hxy(  1:4*Ny,4*Ny*(Nx-1)+1:4*Ny*Nx )=h01x';
H00=hxy;H01=Hz;
%**************************��ϵ�ĵ��ڲ��� **********************
layer=Nz;%size in z direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i3=1:length(EF)
    E=EF(i3)+eta;%the Fermi energy 'E'd    
    tic;
    sigmaL= H01 * SqGreFunLead(Nx*Ny,EF_Lead,H00,H01 ) * H01'; %change here ,so that no change in freenfunc
    H10=H01';
    sigmaR= H10 * SqGreFunLead(Nx*Ny,EF_Lead,H00,H10 ) * H10'; %change here ,so that no change in freenfunc
    for Nrd=1:length(random)
        rd=random(Nrd);
        tic;
        parfor avr=1:average            
            ksi=1;dense=0.05;
            %[gg,G1M] =SqGreFunSys(1,E,H00,H01,sigmaL,Nx*Ny,layer,ksi,dense,rd);
            [gg,G1M]=MacGreFunSys(1,E,H00,inv(H01),H10,sigmaL,Nx*Ny,layer,ksi,dense,rd);
            gg  = inv(inv(gg)-sigmaR);
            G1M = G1M+G1M*sigmaR*gg;
            T10(i6,i5,i3,Nrd,avr)=-real(trace((sigmaL-sigmaL')*G1M*(sigmaR-sigmaR')*G1M'));
        end
        toc
        save T10_6EF T10 Nx Ny Nz  EF EF_Lead mz random average cosk0 TD phi_n
    end
%end of loop on EF
end
%end of loop on phi_n
end
%end of loop on TD
end
save T10_6EFall T10 Nx Ny Nz  EF EF_Lead mz random average cosk0 TD phi_n
%phi_n=Nsample;
%T10aver=zeros(length(random),length(mz));
%T10std=zeros(length(random),length(mz));
[T10aver,T10std]=drawtdphiW(T10,phi_n,random,average,TD);
%figure(5);title(['N_x' 'x' 'N_y' '='  num2str(Nx) 'x' num2str(Ny) 'EF' '=' num2str(EF)...
 %   'random' '=' num2str(random) 'average' '=' num2str(average) ],...
  %  'Units','normalized','HorizontalAlignment','center','FontSize',12);
%figure(6);title(['N_x' 'x' 'N_y' '='  num2str(Nx) 'x' num2str(Ny) 'EF' '=' num2str(EF)...
  %  'random' '=' num2str(random) 'average' '=' num2str(average) ],...
  %  'Units','normalized','HorizontalAlignment','center','FontSize',12);
%figure(7);title(['N_x' 'x' 'N_y' '='  num2str(Nx) 'x' num2str(Ny) 'EF' '=' num2str(EF)...
 %   'random' '=' num2str(random) 'average' '=' num2str(average) ],...
 %   'Units','normalized','HorizontalAlignment','center','FontSize',12);
%seeL=1;seeR=2;
%figure(8); plot(phi_n(seeL:seeR),T10aver(seeL:seeR));xlabel('phi');ylabel('R');
%figure(8); plot(phi_n(1:see),T10std(1:see));xlabel('phi');ylabel('R');
matlabpool close;