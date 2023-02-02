function [T10aver,T10std]=drawT10all(T10,phi_n,random,average,TD,EF)
%T10=zeros(length(TD),length(phi_n),length(EF),length(random),average);
Lrd=length(random);Nphi=length(phi_n);NTD=length(TD);lenE =length(EF);
T10aver=zeros(Nphi,Lrd,NTD,lenE);T10std=zeros(Nphi,Lrd,NTD,lenE);
%T10aver=zeros(Nphi*NTD,Lrd);T10std=zeros(Nphi*NTD,Lrd);
for k1=1:Nphi
    for k2=1:Lrd
        for k3=1:NTD
            for k4=1:lenE
                %hint: transfer the k labels from the left to the right;
                % new labels are inserted without moving other labels
                T10aver(k1,k2,k3,k4)=mean( T10(k3,k1,k4,k2,:) );
                T10std(k1,k2,k3,k4)=std(T10(k3,k1,k4,k2,:));
            end
        end
    end
end
k1_NphiDefault=1;
k2_LrdDefault=1;
k3_NTDDefault=1;
k4_lenEDefault=1;
  size(T10std)
  figure(202302021);hold off;
  figure(202302022);hold off;
for k4_lenEDefault=1:lenE
    if 1==1
    figure(202302021);hold on;
    plot(random,T10aver(k1_NphiDefault,:,k3_NTDDefault,k4_lenEDefault),'o-','MarkerFaceColor','g','MarkerSize',10);set(gca,'FontSize',24);hold on;
    xlabel('random','fontsize',24);ylabel('Taver','fontsize',24);hold on;axis tight
    end
    figure(202302022);hold on;
    plot(random,T10std(k1_NphiDefault,:,k3_NTDDefault,k4_lenEDefault),'o-','MarkerSize',10);set(gca,'FontSize',24);hold on;
    xlabel('W','fontsize',24);ylabel('\delta G','fontsize',24);hold on;
end

axis tight
gcalinewidth=2;
set(gca,'linewidth',gcalinewidth);
%text(1.5,0.6,['$ \bf L_x=10\;L_y=10\;L_z=100\;; E_F=3\;away DP\;;D=0 $'],'Interpreter','latex','fontsize',15);
OrangeColor=[0.9290 0.6940 0.1250];MaroonColor=[128 0 0]/256;
plot([random(1)-1 random(end)+1],0.365*[1 1]/sqrt(2),'--','linewidth',2,'color',OrangeColor);
plot([random(1)-1 random(end)+1],0.365*[1 1],'--','linewidth',2,'color',MaroonColor);
plot([random(1)-1 random(end)+1],0.516*[1 1],'--k','linewidth',2);
legend('0(\phi/\pi)','10^{-4}','0.0013','0.1257','0.2513','0.7540');
%legend Box off;