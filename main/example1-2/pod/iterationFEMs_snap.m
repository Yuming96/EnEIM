function [U]=iterationFEMs_snap(kai,f,G,T,W_train,Ntr,If,Jf,Val_f,I,J,Value_K)
U=zeros(size(T.FNodePtrs,1),Ntr);
for j=1:Ntr
   Val_k=repmat(kai(W_train(:,j)),16,1).*Value_K;
   KK=sparse(I,J,Val_k);
   f1=@(x) f(x,W_train(:,j)');
   F1=repmat(f1(T.centriod),4,1).*Val_f;
   F=sparse(If,Jf,F1);
   U0=KK(T.FNodePtrs,T.FNodePtrs)\(F(T.FNodePtrs)-KK(T.FNodePtrs,T.CNodePtrs)*G);
   U(:,j)=U0;
end