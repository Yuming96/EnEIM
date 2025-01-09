function [U,Ite,err1]=iterationFEMs_POD(kai,f,G,T,W_train,Ntr,If,Jf,Val_f,K0,I,J,Value_K)
% A0=exp(mu);
% K0=stiff1(T,A0,Nx,Ny);

Iter=20;

Tol=1e-4;
U=zeros(size(T.FNodePtrs,1),Ntr);
% F=zeros(size(T.FNodePtrs,1),Ntr);
Ite=[];err1=[];
for j=1:Ntr
%     a=Kai(:,j); 
%    K=stiff1(T,a,Nx,Ny);
%    KK=0;
%    KK_g=0;
%    for i=1:P
%          KK=KK+W_train(i,j).*Kai{i};
%          KK_g=KK_g+W_train(i,j).*Kai_g{i};
%    end
   Val_k=T.area.*repmat(kai(W_train(:,j)),16,1).*Value_K/4;
   KK=sparse(I,J,Val_k);
   f1=@(x) f(x,W_train(:,j)');
%    F=loaf(T,f1,Nx,Ny);
   F1=repmat(f1(T.centriod),4,1).*Val_f/4;
   F=sparse(If,Jf,F1);
%    F(T.CNodePtrs)=[];
   U0=K0{1}\(F(T.FNodePtrs)-(K0{2})*G);
   count=1;
   err=max(abs(U0));
    U1=U0;
%     KK0=K0{1}\KK(T.FNodePtrs,T.FNodePtrs);
%     KK0_G=K0{1}\KK(T.FNodePtrs,T.CNodePtrs)*G;
    while count<Iter && err>Tol
%             U2=U0-KK0*U1-KK0_G; 
            U2=U0-K0{1}\(KK(T.FNodePtrs,T.FNodePtrs)*U1+KK(T.FNodePtrs,T.CNodePtrs)*G);
            err=max(abs(U2-U1));
            count=count+1;
            U1=U2;
    end
    Ite=[Ite;count];
    err1=[err1;err];
    U(:,j)=U2;
end