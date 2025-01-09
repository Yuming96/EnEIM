function Ufem=ufem_online_pod(w,T,Nxo,Nyo,f,g,kai,Phi)
N=size(w,2);
Ufem=zeros(N,size(Phi,2));
for i=1:N
      k1=@(x) kai(x,w(:,i)');
       a=k1(T.centriod);
%       a=kai(w(2:end,i));
       K=stiff_m(T,a,Nxo,Nyo);
       f1=@(x) f(x,w(:,i)');
        F=loaf_m(T,f1,Nxo,Nyo);
      Ufem(i,:)=iterationFEMs_online(T,K,F,g,Phi);
end