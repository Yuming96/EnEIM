function U_new=knn_rom(W_train,U_train,n,w_test,kai)
Ntrain=size(W_train,2);
res= sum((kai(W_train)-repmat(kai(w_test),1,Ntrain)).^2,1);
res=res.^0.5;
[d,ind]=sort(res,'ascend');
dr=1./(d(1:n));

dr=dr./sum(dr);
UUt=U_train(:,ind(1:n));
U_new=zeros(size(UUt,1),1);
for i=1:n
    U_new=U_new+dr(i).*UUt(:,i);
end
