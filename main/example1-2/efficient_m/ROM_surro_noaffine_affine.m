function [U,err]=ROM_surro_noaffine_affine(ks,fs,T,w1,Ax,Fx,Fg,G0,K0)
Iter=6;
% Ite=[];
Fs=zeros(size(fs,2),size(w1',1));
for j=1:size(w1',1)
for i=1:size(fs,2)
    Fs(i,j)=fs{i}(w1(:,j)');
end
end
Ks=zeros(size(ks,2),size(w1',1));
% Ks(1,:)=ks{1}(w1');
% Ks0=cell(1,size(ks,2));
% bb=repmatsize(w1',1)(ks{1}(w1'),1,length(T.Nodes))';
% Ks0{1}=Ax{1}*repmat(ks{1}(w1'),1,length(T.Nodes))';
% parfor i=2:size(ks,2)
% %     Ks(i,:)=ks{i}(w1');
%     Ks0{i}=Ax{i}*repmat(ks{i}(w1'),1,length(T.Nodes))';
% end
% Ks=ks{1}(w1');
for j=1:size(w1',1)
for i=1:size(ks,2)
    Ks(i,j)=ks{i}(w1(:,j)');
end
end
U0=Fx(T.FNodePtrs,:)*Fs-repmat(Fg,1,size(w1,2));U=cell(1,Iter);
U{1}=K0(T.FNodePtrs,T.FNodePtrs)\U0;count=1;err=zeros(Iter,size(w1,2));
err(1,:)=norm(U{1},1)./norm(U{1},1);
while count<Iter
         G0(T.FNodePtrs,:)=U{count};
         UU=Urom_online1(G0,U0,Ax,Ks,K0,T);
         U{count+1}=UU;
        count=count+1
       err(count,:)=norm(U{count}-U{count-1},1)./norm(U{count-1},1);
end

