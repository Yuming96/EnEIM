function U=iterationFEMs(Te,B,K,t,f1,G,U0,Nxo,Nyo,T)
Nt=length(t)-1;
ht=Te/Nt;
% deta=gamma(2-gama)*ht^gama;
deta=ht;
% U=zeros(size(K,1),Nt+1);
U=repmat(G(:,1),1,Nt+1);
U(:,1)=U0;
for j=2:Nt+1
    f2=@(x) f1(x,t(j));
    F=loaf_m(T,f2,Nxo,Nyo);
%     [c]=coeff_bc(j,gama);
%     U(:,j)=(B{1}+deta*(K{1}+B2{1}))\(B{1}*U(:,1:j-1)*c+deta*F(:,j)-(deta*(K{2}+B2{2}))*G); 
    U(T.FNodePtrs,j)=(B(T.FNodePtrs,T.FNodePtrs)+deta*K(T.FNodePtrs,T.FNodePtrs))\(B(T.FNodePtrs,T.FNodePtrs)*U(T.FNodePtrs,j-1)+...
        deta*F(T.FNodePtrs)-deta*K(T.FNodePtrs,T.CNodePtrs)*G(T.CNodePtrs,1)); 
end
