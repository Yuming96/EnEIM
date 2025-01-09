function [U,err]=ROM_surro_Noaffine(k,fss,T,W,w1,Fx,Fg,G0,co,I,J,Val_k,K00)
Iter=10;
Fs=fss(w1');
m=length(T.FNodePtrs);
% Ks0=zeros(m,size(T.Nodes,1),size(w1,2));
Ks0=cell(1,size(w1,2));
% tic
for i=1:size(w1,2)
%     aa=k(co,w1(:,i)');
   K0= sparse(I,J,Val_k.*repmat(W(:,i),16,1));
%     Ks0(:,:,i)=K00(T.FNodePtrs,T.FNodePtrs)\K0(T.FNodePtrs,:);
Ks0{i}=K0;
end
% toc
% Ax0=pagemtimes(Kf,Ks0);

U0=Fx(T.FNodePtrs,:)*Fs'-repmat(Fg,1,size(w1,2));
% U0=repmat(Fx-Fg,1,size(w1,2));
% U=cell(1,Iter);
U{1}=K00(T.FNodePtrs,T.FNodePtrs)\U0;
count=1;
% err=zeros(Iter,size(w1,2));
err=max(mean(abs(U{1}./U{1}),2));
while count<Iter&&err(count)>1e-6
         G0(T.FNodePtrs,:)=U{count};
%          bb=repmat(G0,1,1,size(w1,1));
%          GG=dot(pagemtimes(Ax0,repmat(G0,1,1,size(w1,1))),Ks0,3);
%          UU=Urom_online(G0,U0,Ax0,Ks0);
%          UU=U0-dot(GG,Ks0,3);
         GG=zeros(size(U0));
         for i=1:size(w1,2)
%              aa=Ks0{i}(T.FNodePtrs,:);
              GG(:,i)=Ks0{i}(T.FNodePtrs,:)*G0(:,i);
         end
% %          GG=pagemtimes(Ks0,G0);
%          tic
%          GG=E0*La0*E0'*GG;
%          toc
%          tic
%          GG00=pinv(full(K00(T.FNodePtrs,T.FNodePtrs)));
%          GG1=K00(T.FNodePtrs,T.FNodePtrs)\GG;
%          toc
         UU=K00(T.FNodePtrs,T.FNodePtrs)\(U0-reshape(GG,m,size(w1,2)));
         U{count+1}=UU;
        count=count+1
       err=[err;max(mean(abs((U{count}-U{count-1})./U{count-1}),2))];
end

