function [U,err,count]=ROM_surro_affine_f(Ks0,Fs,T,w1,Fx,Fg,G0,U1,M1,ht)
Iter=20;
m=length(T.FNodePtrs);
U0=U1+Fx*Fs'-repmat(Fg,1,size(w1,2));
% U{1}=U0;
count=1;
Urom=zeros(size(G0));
Urom(T.CNodePtrs,:)=G0(T.CNodePtrs,:);
Urom(T.FNodePtrs,:)=M1(T.FNodePtrs,T.FNodePtrs)\U0;
U{count}=Urom;
err=max(mean(abs(U{1}./U{1}),2));
while count<Iter&&err(count)>1e-6
         G0(T.FNodePtrs,:)=U{count}(T.FNodePtrs,:);
         GG=zeros(size(U0));
         for i=1:size(w1,2)
              GG(:,i)=Ks0{i}(T.FNodePtrs,:)*(G0(:,i));
         end
%          GG1=GG.*ht;
         UU=M1(T.FNodePtrs,T.FNodePtrs)\(U0-reshape(GG.*ht,m,size(w1,2)));
%          U{count+1}=UU;
         Urom(T.CNodePtrs,:)=G0(T.CNodePtrs,:);
        Urom(T.FNodePtrs,:)=UU;
        count=count+1;
        U{count}=Urom;
       err=[err;max(mean(abs((U{count}-U{count-1})./U{count-1}),2))];
%        Urom=zeros(size(G0));

end

