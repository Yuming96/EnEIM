function Ufem=ufem(W_train,T,Nxo,Nyo,f,g,kai)
Ntrain=size(W_train,2);
Ufem=zeros(size(T.Nodes,1),Ntrain);
for i=1:Ntrain
      k1=@(x) kai(x,W_train(:,i)');
       a=k1(T.centriod);
%       a=kai(W_train(2:end,i));
       K=stiff_m(T,a,Nxo,Nyo);
       f1=@(x) f(x,W_train(:,i)');
        F=loaf_m(T,f1,Nxo,Nyo);
      Ufem(:,i)=iterationFEMs(T,K,F,g);
end