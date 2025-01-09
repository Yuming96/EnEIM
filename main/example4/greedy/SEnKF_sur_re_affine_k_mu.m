function [Q,Usim,err]=SEnKF_sur_re_affine_k_mu(Uin,sig,M,q,p,Uobs,k,co,KQVo)
 m=25;
%   Z=zeros(length(p),M);
  CD=q;
   Q{1}=CD;
%    error=zeros(m,1);
   l=1;
   err=1;
while l<m && err(l)>1e-3
    l
       W1=CD;
       W=W1;
           r1=Uin(W); 
           r=KQVo*r1';
            Z=r(p,:);
Usim{l}=Z;
      mean_q=mean(CD,2);
    mean_z=mean(Z,2);
    Z_obs=repmat(Uobs(:),1,M)+sig.*randn(length(p),M);
    delta_z=(Z-repmat(mean_z,1,M))./sqrt(M-1)./sig;
    delta_q=(CD-repmat(mean_q,1,M))./sqrt(M-1);
    [u,s,v]=svd(delta_z,'econ');
%     [lambda,order]=sort(diag(s),'descend');
%      v=v(:,order);
%      u=u(:,order);
%     ss=lambda/sum(lambda);
%     truc=sum(cumsum(ss)<0.95);
truc=6;

%    s=diag(lambda(1:truc));
   K=delta_q*v(:,1:truc)*s(1:truc,1:truc)*pinv(1+s(1:truc,1:truc).^2)*u(:,1:truc)';
   CD=CD+1/sig.*K*(Z_obs-Z);
   Q{l+1}=CD;
   k_mean=mean(Q{l+1},2);
   k_mean1=mean(Q{l},2);
error=norm(k(co,k_mean')-k(co,k_mean1'))/norm(k(co,k_mean1'));
err=[err error];
l=l+1;
end