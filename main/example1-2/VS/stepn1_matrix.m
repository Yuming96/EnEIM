function [av,VS,UU,w1,uvpon,el,kk]=stepn1_matrix(w1,Atem,Ax,KQv,active,Nt,Fx,UU,VS,fs,ks,G,K0)
%Ñ¡Ñù±¾
ELt=zeros(size(w1,2),1);
for ji=1:size(w1,2)
%     avt=KQv(:,end)*VS(ji,end)';
%     avt=KQv*VS(ji,:)';
    avt=UU(:,ji)-KQv*VS(ji,:)';
    ELt(ji)=avt'*Atem*avt;
end
[el,ind]=max(ELt);
kk=w1(:,ind);
[av] = residsotion1(active,w1(:,ind)',ks,fs,G,Ax,Fx,KQv*VS(ind,:)',K0);
nF=Fx'*av;
a1=zeros(1,size(ks,2));a1g=zeros(size(a1));
Aq=zeros(Nt,size(ks,2));%A2=zeros(Nt+1,size(ks,2));
a1g0=av(active)'*K0(active,:)*G;a10=av(active)'*K0(active,active)*av(active);
Aq0=KQv(active,:)'*K0(active,active)*av((active));
for i=1:size(ks,2)
    a1(i)=av(active)'*Ax{i}(active,active)*av(active);
    a1g(:,i)=G'*Ax{i}(:,active)*av(active);
    Aq(:,i)=KQv(active,:)'*Ax{i}(active,active)*av((active));
%     A2(:,i)=[KQv(T.cnodes,:) av(T.cnodes)]'*Ax{i}(:,active)*av((active));
end
% A2=a2c;
uvpon={nF,a1,a1g,a1g0,a10,Aq,Aq0};
[uv]=ups_online(w1',VS,nF,a1,a1g,a1g0,a10,Aq,Aq0,fs,ks);
%%
VS=[VS uv];UU(:,ind)=[];VS(ind,:)=[];w1(:,ind)=[];

