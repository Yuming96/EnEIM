function [kx,ks,k]=Affine1(ki,s0,s01,cof,xo,yo,Nxm,Nym)

Nw=30; 
xm=xo(1):1/Nxm:xo(end);
ym=yo(1):1/Nym:yo(end); 
[X,Y]=meshgrid(xm,ym);
co=[X(:) Y(:)];
elsonk=0.0001;

% % % % %k函数
% % P1=dim;kx1=cell(1,P1);
% % kx0=@(x) 1+0.*x(:,1);
% % % kx0=@(x) exp(x(:,1)+x(:,2));
% % kx1{1}=@(x,s) (1.*cos(pi.*s(:,1).*(x(:,1)+x(:,2)))).*(1.*sin(pi.*s(:,1).*(x(:,2)-3.*x(:,1))));
% % for i=2:P1
% %     kx1{i}=@(x,s) (1.*cos(pi.*s(:,i).*(x(:,1)+x(:,2)))).*(1.*sin(pi.*s(:,i).*(x(:,2)-3.*x(:,1))));
% % end
% % 
% % k1=@(x,s) kx1{1}(x,s);
% % for i=2:P1
% %     k1=@(x,s) k1(x,s)+kx1{i}(x,s);
% % end
% % ki=@(x,s) kx0(x)+k1(x,s);
% Ntr=Nes;
% P=15;
% s0=randn(Ntr,P);
% % for i=1:P
% %     s0(:,i)=2.^(i-1).*(1+0.5.*rand(Ntr,1));
% % end

inda=zeros(1,size(co,1)); 
% KQ =cell(1,Nw);Qq =cell(1,Nw);
count=1;
ke=@(s) ki(co,s);
% ke(s0) 
[a,u] = stepk1(co,ke,s0(1,:),inda);
% KQ{count}=a;KQa=zeros(size(a,1),Nw);KQa(:,count)=a;
% Qq{count}=@(s) u(s);
ks{count}=@(s) u(s);
% [a,u] = stepkxs(hx,Nn,elson,t,ke,w1);
kx{count}=@(x) interp2(X,Y,reshape(a,Nym+1,Nxm+1),x(:,1),x(:,2));
% aa=kx{count}(cof);
k=@(x,s) kx{count}(x).*ks{count}(s);
Ki=@(s) a*u(s);
ke=@(s) ke(s)-Ki(s);
% s01=randn(1,size(s0,2));

el=norm(ke(s01))/norm(ki(co,s01))
while el> elsonk & count <Nw
    count=count+1
    inda=zeros(1,size(co,1)); 
%     w1{1}=randn(2*Nes,dim);
%     w1{2}=randn(2*Nes,dim);
%     [a,u] = stepk(t,ke,w1{mod(count+1,4)+1},inda,KQa(:,1:count-1));
   [a,u] = stepk1(co,ke,s0(count,:),inda);
% %     [a,u] = stepk1_orth(t,ke,w1{mod(count+1,4)+1},inda,KQa,count);
%     KQ{count} =a;KQa(:,count)=a;
%     Qq{count}=@(s) u(s);
    ks{count}=@(s) u(s);
%     [a,u] = stepkxs(hx,Nn,elson,t,ke,w1);
    kx{count}=@(x) interp2(X,Y,reshape(a,Nym+1,Nxm+1),x(:,1),x(:,2));
%     aa=kx{count}(cof);
    k=@(x,s) k(x,s)+kx{count}(x).*ks{count}(s);
    Ki=@(s) Ki(s)+a*u(s);
    ke=@(s) ki(co,s)-Ki(s);
    el=norm(ke(s01))/norm(ki(co,s01))
end
% %%
% wo=2*rand(5,dim)-1;
% coefk=zeros(size(wo,1),1);
% coefs=zeros(size(wo,1),1);
% for j=1:size(wo,1)
%     s=wo(j,:);
%     aa=k(co,s)-ki(co,s);
%     coefk(j)=norm(ke(s))/norm(ki(co,s));
%     coefs(j)=norm(aa)/norm(ki(co,s));
% end
% mean(coefk)
% mean(coefs)