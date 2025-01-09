function [av,VS,UU,w1,uvpon,el,iex]=stepn_ortho_matrix(w1,Atem,Ax,KQv,active,Nt,Fx,UU,VS,fs,ks,G,iex)
ELt=zeros(size(w1,2),1);
for ji=1:size(w1,2)
%     Ubi=KQv*VS(ji,:)';
    avt=UU(:,ji)-KQv*VS(ji,:)';
    ELt(ji)=avt'*Atem*avt;
end
[el,ind]=max(ELt);
[av] = UU(:,ind)-KQv*VS(ind,:)';
iex=[iex;ind];
% % a1=zeros(1,size(ks,2));%a2c=zeros(size(a1c));
% % A1=zeros(Nt,size(ks,2));A2=0;
% % for i=1:size(ks,2)
% %     a1(i)=av'*Ax{i}*av;%a2c(:,i)=Ax{i}*G(active);
% %     A1(:,i)=KQv'*Ax{i}*av;
% % end
% [av] = residsotion1(active,w1(ind,:),ks,fs,G,Ax,Fx,KQv*VS(ind,:)');
av=othorg(av,KQv,Atem);%av:正交化后的基函数
nF=Fx'*av;
a1=zeros(1,size(ks,2));a2c=zeros(size(a1));
A1=zeros(Nt,size(ks,2));%A2=zeros(Nt+1,size(ks,2));
for i=1:size(ks,2)
    a1(i)=av(active)'*Ax{i}(active,active)*av(active);
    a2c(:,i)=G'*Ax{i}(:,active)*av(active);
    A1(:,i)=KQv(active,:)'*Ax{i}(active,active)*av((active));
%     A2(:,i)=[KQv(T.cnodes,:) av(T.cnodes)]'*Ax{i}(:,active)*av((active));
end
A2=a2c;
uvpon={nF,A1,A2,a1};
[uv]=ups_online(w1',VS,nF,A1,A2,a1,fs,ks);
%%
VS=[VS uv];w1(:,ind)=[];VS(ind,:)=[];UU(:,ind)=[];

