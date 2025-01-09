clear,clc,close all  
tic
[sig,P,Te]=deal(0.01,5,0.1);
Nt=100;
Nxo=100;Nyo=100;
%  Nxo=20;Nyo=20;
bd=[0 1;0 1];
hx=(bd(1,2)-bd(1,1))/Nxo;
hy=(bd(2,2)-bd(2,1))/Nyo;
ht=Te/Nt;
t=0:ht:Te;


xo=bd(1,1):hx:bd(1,2);yo=bd(2,1):hy:bd(2,2);
%  f=@(x,w) 1+0.*x(:,1).*w(:,1);
T=struc(xo,yo);
g=@(x) x(:,2).*x(:,1).*(x(:,1)-1).*(x(:,2)-1)+1;
% g=@(x) 0.*(x(:,1)-1).*(x(:,2)-1)+1;
u0=@(x) sin(2*pi.*x(:,1)).*sin(2*pi.*x(:,2))+0;
% u0=@(x) 1.*x(:,1);
[X,Y]=meshgrid(xo,yo);
No=[X(:) Y(:)];
ke=@(x)  0.*x(:,1).*x(:,2)+1;
U0=u0(T.Nodes);

%%%%use for inverse problem, the observation locations 
[~,Nodes,inter,co,~]=point_numberf(Nxo,Nyo,T);
bdy=T.CNodePtrs;

%%%%KLE expasion for permeability
xmo=1/Nxo/2:1/Nxo:1;
ymo=1/Nyo/2:1/Nyo:1;
[Nxm,Nym]=deal(20,20);

tic
% % %k函数
k0=@(x) exp(x(:,1)+x(:,2));
ks{1}=@(s) s(:,1);
kx{1}=@(x) x(:,1);
ks{2}=@(s) s(:,2);
kx{2}=@(x) x(:,2);
ks{3}=@(s) s(:,3);
kx{3}=@(x) x(:,1).*x(:,2);

kx{1}=@(x) kx{1}(x);
k=@(x,s) k0(x)+ks{1}(s).*kx{1}(x);
k1=@(x,s)ks{1}(s).*kx{1}(x);
kss=@(s) ks{1}(s);
for i=2:size(ks,2)
    kx{i}=@(x) kx{i}(x);
    k=@(x,s) k(x,s)+ks{i}(s)*kx{i}(x);
    k1=@(x,s) k1(x,s)+ks{i}(s)*kx{i}(x);
    kss=@(s) [kss(s) ks{i}(s)];
end

fs{1}=@(s) sin(s(:,1).*s(:,end));
fx{1}=@(x,t) 0.1.*exp(t).*(x(:,2)+3);
f=@(x,t,s) fs{1}(s)*fx{1}(x,t);
fss=@(s) fs{1}(s);
for i=2:size(fs,2)
    f=@(x,t,s) f(x,s)+fs{i}(s)*fx{i}(x,t);
     fss=@(s) [fss(s) fs{i}(s)];
end


com=Nodes(bdy,:);
%% greedy offline
m=length(T.FNodePtrs);
[If,Jf,Val_f]=loaf(T,Nxo,Nyo);
[I,J,Val_k,Val_b]=stiff(T,Nxo,Nyo);
% A0=exp(mu);
 K0=sparse(I,J,Val_k.*repmat(k0(co),16,1));
% % K0=sparse(I,J,Val_k.*repmat(exp(mu),16,1));
%获取列数，含空间的刚度矩阵
Am=sparse(I,J,Val_k);
M=sparse(I,J,Val_b);
%获取列数，含空间的刚度矩阵
G=zeros(size(T.Nodes,1),1);%G边界条件
G(T.CNodePtrs)=g(T.Nodes(T.CNodePtrs,:));
% Kf=K0(T.FNodePtrs,T.FNodePtrs);
M1=M+ht.*K0;
Fg=K0(T.FNodePtrs,:)*G.*ht;

Fx=zeros(m,size(fs,2),Nt+1);Fxq=zeros(size(K0,1),Nt+1,size(fs,2));
%源项的向量
termf=M*fx{1}(T.Nodes,t);
Fx(:,1,:)=ht.*termf(T.FNodePtrs,:);
Fxq(:,:,1)=termf;
f1=@(s) fs{1}(s)*fx{1}(T.Nodes,t);

for i=2:size(fs,2) 
    termf=M*fx{i}(T.Nodes,t);Fxq(:,:,i)=termf;
    Fx(:,i,:)=ht.*termf(T.FNodePtrs,:);
    f1=@(s) f(s)+fs{i}(s)*fx{i}(T.Nodes,t);
end

%% ROM offline
 Nte=1e2;
 load data_w W Ufem
% W_train=2*rand(P,Nte)-1;
W_train=W;
[Fs,Ks0]=Matrix_affine(k1,fss,W_train,co,I,J,Val_k);
Mss=M(T.FNodePtrs,T.FNodePtrs);
G0=zeros(size(T.Nodes,1),Nte);
G0(T.CNodePtrs,:)=repmat(g(T.Nodes(T.CNodePtrs,:)),1,1,Nte);
U00=repmat(U0,1,Nte);
Ueff{1}=U00;
toc
iter=[];
tic
for j=2:Nt+1
    U1=Mss*U00(T.FNodePtrs,:);
% [Urom,ite_tr,err_tra]=iterationFEMs_surro(k1,f,G,T,W_train,Ntrain,M,K0,I,J,Val_k);
Fx1=Fx(:,:,j);
% tic
[Uromo,err_ef,count]=ROM_surro_affine_f(Ks0,Fs,T,W_train,Fx1,Fg,G0,U1,M1,ht);
% toc
Ueff{j}=Uromo;
U00=Uromo{end};
err{j}=err_ef;
iter=[iter count];
end
toc
% % 
% tic
% Ufem=ufem(W_train,T,Nxo,Nyo,f,G0,k,U0,t,M,Te);
% toc

% count=length(n);
Nw=20;Tru=Nw;
coefko=zeros(Nte,Nt+1,Nw);
for ii=2:Nt+1
for i=1:iter(ii-1)
 for j=1:Nte
            Vtemo=Ueff{1,ii}{1,i}(T.FNodePtrs,j)-Ufem{j}(T.FNodePtrs,ii);
        coefko(j,ii,i)=sqrt((Vtemo'*M(T.FNodePtrs,T.FNodePtrs)*Vtemo)/...
            (Ufem{j}(T.FNodePtrs,ii)'*M(T.FNodePtrs,T.FNodePtrs)*Ufem{j}(T.FNodePtrs,ii)));
 end
end
end

eer=mean(coefko,1);

y = reshape(mean(coefko,1),Nt+1,Nw);
e=zeros(Nt+1,Nw);
for j=2:Nt+1
    e_std=reshape(coefko(:,j,:),Nte,Tru);
e(j,:) = std(e_std,1,1);
end
figure(1)
subplot(2,2,1)
errorbar(y(1,:),e(1,:),'xr')
% xlabel('Numbers of the separated terms');  
ylabel('Mean error with error bars');
grid on
subplot(2,2,2)
errorbar(y(Nt/20,:),e(Nt/20,:),'xr')
% xlabel('Numbers of the separated terms');  
% ylabel('Mean error with error bars');
grid on
subplot(2,2,3)
errorbar(y(Nt/2,:),e(Nt/2,:),'xr')
xlabel('Numbers of the separated terms');  
ylabel('Mean error with error bars');
grid on
subplot(2,2,4)
errorbar(y(Nt+1,:),e(Nt+1,:),'xr')
xlabel('Numbers of the separated terms');  
% ylabel('Mean error with error bars');
grid on

figure(2)
plot(1:Tru,y(:,1:end),'r-')
xlabel('Numbers of the iteration steps');  
ylabel('Relative error ');
% legend('N=1','N_p=3','N_p=5') ;
grid on

figure(3)
subplot(2,2,1)
semilogy(1:size(coefko,1),coefko(:,Nt/20,1),'g-')
subplot(2,2,2)
semilogy(1:size(coefko,1),coefko(:,Nt/10,3),'r-')
subplot(2,2,3)
semilogy(1:size(coefko,1),coefko(:,Nt/2,4),'b-')
subplot(2,2,4)
semilogy(1:size(coefko,1),coefko(:,end,5),'b-')
xlabel('Numbers of the separated terms');  
ylabel('Mean error with error bars');
grid on
title ('u')

l=10;
figure(4)
subplot(3,3,1)
mesh(xo,yo,reshape(Ueff{1,Nt/20+1}{1,1}(:,l),Nyo+1,Nxo+1))
title('ROM solution','fontsize',12)
colorbar
subplot(3,3,2)
mesh(xo,yo,reshape(Ueff{1,Nt/20+1}{1,end}(:,l),Nyo+1,Nxo+1))
title('ROM solution','fontsize',12)
colorbar
subplot(3,3,3)
mesh(xo,yo,reshape(Ufem{l}(:,Nt/20+1),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)
colorbar
subplot(3,3,4)
mesh(xo,yo,reshape(Ueff{1,Nt/2+1}{1,1}(:,l),Nyo+1,Nxo+1))
title('ROM solution','fontsize',12)
colorbar
subplot(3,3,5)
mesh(xo,yo,reshape(Ueff{1,Nt/2+1}{1,end}(:,l),Nyo+1,Nxo+1))
title('ROM solution','fontsize',12)
colorbar
subplot(3,3,6)
mesh(xo,yo,reshape(Ufem{l}(:,Nt/2+1),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)
colorbar
subplot(3,3,7)
mesh(xo,yo,reshape(Ueff{1,end}{1,1}(:,l),Nyo+1,Nxo+1))
title('ROM solution','fontsize',12)
colorbar
subplot(3,3,8)
mesh(xo,yo,reshape(Ueff{1,end}{1,end}(:,l),Nyo+1,Nxo+1))
title('ROM solution','fontsize',12)
colorbar
subplot(3,3,9)
mesh(xo,yo,reshape(Ufem{l}(:,end),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)
colorbar
% U_eff{1}=Ueff{1,6};
% U_eff{2}=Ueff{1,51};
% U_eff{3}=Ueff{1,101};

% save data_para_eff_affine U_eff iter err coefko tim_effi W Ufem

