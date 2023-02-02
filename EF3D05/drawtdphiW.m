function [T10aver,T10std]=drawtdphiW(T10,phi_n,random,average,TD)
%T10=zeros(length(TD),length(phi_n),length(random),average);
Lrd=length(random);Nphi=length(phi_n);NTD=length(TD);
T10aver1=zeros(Nphi,Lrd,NTD);T10std1=zeros(Nphi,Lrd,NTD);
T10aver=zeros(Nphi*NTD,Lrd);T10std=zeros(Nphi*NTD,Lrd);
for k1=1:Nphi
    for k2=1:Lrd
        for k3=1:NTD
        T10aver1(k1,k2,k3)=mean( T10(k3,k1,k2,:) );
        T10std1(k1,k2,k3)=std(T10(k3,k1,k2,:));
        T10aver((k1-1)*NTD+k3,k2)=mean( T10(k3,k1,k2,:) );
        T10std((k1-1)*NTD+k3,k2)=std(T10(k3,k1,k2,:));
        %Raver(k1,k2)=mean( 1./T10(k1,k2,:) );
        %Rstd(k1,k2)=std(1./T10(k1,k2,:));
        end
    end
end
%save T10aver T10std
%[Tatemp,Tstemp]= array2cell(T10aver1,T10std1,T10)
  %T10aver=Tatemp{1};
  %T10std=Tstemp{1};
  %phi_n=[0,0.00001,0.0001,0.001,0.01,1,2,6]/10;
  size(T10std)
  FirstPick=1
  seeList=[FirstPick:NTD:(Nphi-1)*NTD+FirstPick]
if 1==2
figure(10);hold off;
plot(random,T10aver(seeList,:),'o-','MarkerFaceColor','g','MarkerSize',10);set(gca,'FontSize',24);hold on;
xlabel('random','fontsize',24);ylabel('Taver','fontsize',24);hold on;axis tight
end
figure(11);hold off;
plot(random,T10std(seeList,:),'o-','MarkerSize',10);set(gca,'FontSize',24);hold on;
xlabel('W','fontsize',24);ylabel('\delta G','fontsize',24);hold on;

axis tight
gcalinewidth=2;
set(gca,'linewidth',gcalinewidth);
text(1.5,0.6,['$ \bf L_x=10\;L_y=10\;L_z=100\;; E_F=3\;away DP\;;D=0 $'],'Interpreter','latex','fontsize',15);
OrangeColor=[0.9290 0.6940 0.1250];MaroonColor=[128 0 0]/256;
plot([random(1)-1 random(end)+1],0.365*[1 1]/sqrt(2),'--','linewidth',2,'color',OrangeColor);
plot([random(1)-1 random(end)+1],0.365*[1 1],'--','linewidth',2,'color',MaroonColor);
plot([random(1)-1 random(end)+1],0.516*[1 1],'--k','linewidth',2);
legend('0(\phi/\pi)','10^{-4}','0.0013','0.1257','0.2513','0.7540');
%legend Box off;