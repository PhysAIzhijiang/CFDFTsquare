function [ gr ] = SqGreFunLead(N, E,H00,H01 )
%This is the lead surface Green function on the square lattice of quantum well 
%   Detailed explanation goes here
    invH01=inv( H01);
    TT=[invH01*(E*eye(N*4)-H00),-invH01*H01';eye(N*4),zeros(N*4)];
    [eigenvector,eigenvalue]=eig(TT);   
    e_abs=abs(diag(eigenvalue));  
    [a,b]=sort(e_abs);  
    v2=eigenvector(:,b);
    S1=v2(1:N*4,1:N*4);
    S2=v2(N*4+1:2*N*4,1:N*4); 
    gr=inv(E*eye(N*4)-H00-H01*S1*inv(S2));

end

