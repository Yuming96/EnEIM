%POD在线计算函数 

function su_solu=fast_pod1(theta,U_pod,A_hat,B_hat,ht,tn,u0_hat,f0,X,F_hat)
U_hat=[u0_hat];
delt1=theta(1,2)^2;delt2=theta(1,1);%1是theta,2是mu
u=u0_hat;
u_old=u0_hat;
for j=2:length(tn)
         U_hat=(B_hat+ht*A_hat)\(ht*F_hat*f0(X,tn(j),theta')+B_hat*u_old-ht.*A_hat*G);
       U_hat=[U_hat,u];
        u_old=u;
end
U=U_pod*U_hat;
%%
% b_u=zeros(1,size(U,2));
% su_solu=[b_u;U;b_u];
su_solu=U;