function [U,Ite,err1]=iterationFEMs_surro(kai,f,G,T,W_train,Ntr,If,Jf,Val_f,K0,I,J,Value_K)
Iter=20;Tol=1e-4;
U=zeros(size(T.Nodes,1),Ntr);
U(T.CNodePtrs,:)=repmat(G,1,Ntr);
Ite=[];err1=[];
for j=1:Ntr
    j
   Val_k=repmat(kai(W_train(:,j)),16,1).*Value_K;
   KK=sparse(I,J,Val_k);
   f1=@(x) f(x,W_train(:,j)');
   F1=repmat(f1(T.centriod),4,1).*Val_f;
   F=sparse(If,Jf,F1);
   U0=K0(T.FNodePtrs,T.FNodePtrs)\(F(T.FNodePtrs)-K0(T.FNodePtrs,T.CNodePtrs)*G);
   count=1;
   err=max(abs(U0));
    U1=U0;
%     KK0=K0{1}\KK(T.FNodePtrs,T.FNodePtrs);
%     KK0_G=K0{1}\KK(T.FNodePtrs,T.CNodePtrs)*G;
    while count<Iter && err>Tol
%             U2=U0-KK0*U1-KK0_G; 
            U2=U0-K0(T.FNodePtrs,T.FNodePtrs)\(KK(T.FNodePtrs,T.FNodePtrs)*U1+KK(T.FNodePtrs,T.CNodePtrs)*G);
            err=max(abs(U2-U1));
            count=count+1
            U1=U2;
    end
    Ite=[Ite;count];
    err1=[err1;err];
    U(T.FNodePtrs,j)=U2;
end