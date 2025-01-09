function [A,Phi,m1]=pod_svd_train(U)
% Y0=mean(U,2);
% Y=U-repmat(Y0,1,size(U,2));
% Yde=Y/sqrt(size(U,2)-1);
[~,S,V]=svd(U');
Sv=diag(S);
ss=sum(Sv);
% m1=sum(cumsum(Sv/ss)<0.9999);
 m1=12;
% Sr=Sv(1:m1);
Vr=V(:,1:m1);
Phi=Vr;
A=U'*Phi;
% Lambda=Sr.^2;

% Phi=zeros(Nxo*Nyo,m1);
% for i=1:m1
%     Phi(:,i)=1/Sdr(i)*sum(repmat(Vr(:,i)',Nxo*Nyo,1).*Y,2);
% end 
% 
% %%给一个xi我们就有一个新的yy
% % xi=randn(m1,1);
% % Phi=V;
% % C1=Phi'*Phi;
% yy=@(xi) 1/sqrt(P)*Phi*diag(Sdr)*xi+Y0;