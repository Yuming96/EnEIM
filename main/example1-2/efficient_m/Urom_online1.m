function U=Urom_online1(G0,U0,Ax,Ks,K0,T)
% U=U0;
U1=Ax{1}(T.FNodePtrs,:)*G0.*repmat(Ks(1,:),size(U0,1),1);
for i=2:size(Ks,1)
  U1=U1+Ax{i}(T.FNodePtrs,:)*G0.*repmat(Ks(i,:),size(U0,1),1);
end
U=K0(T.FNodePtrs,T.FNodePtrs)\(U0-U1);