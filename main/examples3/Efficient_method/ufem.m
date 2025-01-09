function Ufem=ufem(W_train,T,Nxo,Nyo,f,G,kai,U0,t,B,Te)
Ntrain=size(W_train,2);
% Ufem=zeros(size(T.Nodes,1),Nt+1,Ntrain);
for i=1:Ntrain
      k1=@(x) kai(x,W_train(:,i)');
       a=k1(T.centriod);
%       a=kai(W_train(:,i));
       K=stiff_m(T,a,Nxo,Nyo);
       f1=@(x,t) f(x,t,W_train(:,i)');
%         F=loaf_m(T,f1,Nxo,Nyo);
%       Ufem(:,i)=iterationFEMs(T,K,F,g);
       Ufem{i}=iterationFEMs(Te,B,K,t,f1,G,U0,Nxo,Nyo,T);
end
