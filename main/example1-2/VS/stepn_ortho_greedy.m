function [KQv,VS,UU,w1,uvpon,el,iex]=stepn_ortho_greedy(w1,Atem,Ax,KQv,active,Fx,UU,VS,fs,ks,G,iex,K0)
ELt=zeros(size(w1,2),1);
for ji=1:size(w1,2)
%     avt=KQv*VS(ji,:)';
    avt=UU(:,ji)-KQv*VS(ji,:)';
    ELt(ji)=avt'*Atem*avt;
end
[el,ind]=max(ELt);
% av=UU(:,ind)-KQv*VS(ind,:)';
iex=[iex;ind];
av = UU(:,ind);
av=othorg(av,KQv,Atem);
KQv=[KQv av];
nFG=KQv'*Fx;%%BD
A1g=zeros(size(KQv,2),size(KQv,2),size(ks,2));nFg=zeros(size(KQv,2),size(ks,2));
nFg0=zeros(size(KQv,2),1);A1g0=zeros(size(KQv,2),size(KQv,2),1);
nFg0(:,1)=KQv(active,:)'*K0(active,:)*G;A1g0(:,:,1)=KQv(active,:)'*K0(active,active)*KQv(active,:);
for i=1:size(ks,2)
    nFg(:,i)=KQv(active,:)'*Ax{i}(active,:)*G;
    A1g(:,:,i)=KQv(active,:)'*Ax{i}(active,active)*KQv(active,:);
end
uvpon={nFG,nFg,A1g,nFg0,A1g0};
[VS]=online_RBgreedy(w1',nFG,nFg,A1g,nFg0,A1g0,fs,ks);
%%
w1(:,ind)=[];VS(ind,:)=[];UU(:,ind)=[];

