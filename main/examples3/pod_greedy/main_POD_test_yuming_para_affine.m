clear,clc,close all  
tic
[sig,P,Te,Nte]=deal(0.01,4,0.1,1e2);
Nt=100;
Nxo=100;Nyo=100;
%  Nxo=40;Nyo=40;
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

tic
%k函数
k0=@(x) exp(x(:,1)+x(:,2));
ks{1}=@(s) s(:,1);
kx{1}=@(x) x(:,1);
ks{2}=@(s) s(:,2);
kx{2}=@(x) x(:,2);
ks{3}=@(s) s(:,3);
kx{3}=@(x) x(:,1).*x(:,2);

k=@(x,s) k0(x)+ks{1}(s)*kx{1}(x);
for i=2:size(ks,2)
%     kx{i}=@(x) kx{i}(x);
    k=@(x,s) k(x,s)+ks{i}(s)*kx{i}(x);
%     k=@(x,s) k(x,s)+ks{i}(s)*Fai1(:,i);
end

% %f函数
fs{1}=@(s) sin(s(:,1).*s(:,end));
fx{1}=@(x,t) 0.1.*exp(t).*(x(:,2)+3);
f=@(x,t,s) fs{1}(s)*fx{1}(x,t);
for i=2:size(fs,2)
    f=@(x,t,s) f(x,s)+fs{i}(s)*fx{i}(x,t);
end



com=T.Nodes(bdy,:);
% G=g(com);
% F=loaf(T,f,Nxo,Nyo);
m=length(T.Nodes);
[If,Jf,Val_f]=loaf(T,Nxo,Nyo);
[I,J,Val_k,Val_b]=stiff(T,Nxo,Nyo);
K0=sparse(I,J,Val_k.*repmat(k0(co),16,1));
Am=sparse(I,J,Val_k);
M=sparse(I,J,Val_b);
%获取列数，含空间的刚度矩阵
Ax=cell(1,size(ks,2));Axq=cell(1,size(ks,2));
for i=1:size(ks,2)
%     tem=sparse(I,J,Val_k.*repmat(Fai1(:,i),16,1));
   tem=sparse(I,J,Val_k.*repmat(kx{i}(co),16,1));
    Ax{i}=tem(T.FNodePtrs,T.FNodePtrs);Axq{i}=tem;
end
Fx=zeros(m,Nt+1,size(fs,2));
%源项的向量
Fx(:,:,1)=M*fx{1}(T.Nodes,t);
% f=@(s) fs{1}(s)*fx{1}(T.Nodes);
for i=2:size(fs,2)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    Fx(:,:,i)=M*fx{i}(T.Nodes,t);
%     f=@(s) f(s)+fs{i}(s)*fx{i}(T.Nodes);
end

%% greedy offline
% load data_gree wo W_train1
 Ntrain=5e1;
W_train=2*rand(size(ks,2),Ntrain)-1;
G=zeros(size(T.Nodes,1),1);
G(T.CNodePtrs)=g(T.Nodes(T.CNodePtrs,:));

%  tic
% % [Usnap]=snapsolution(W_train,ks,fs,T,g,Axq,Fx,K0);
% Usnap=ufem(W_train,T,Nxo,Nyo,f,g,k);
Usnap=ufem(W_train,T,Nxo,Nyo,f,G,k,U0,t,M,Te);
r=20;
sumU=[];
for i=1:Ntrain
     X2=Usnap{i};
     XX0=X2(T.FNodePtrs,:);
    [U,~,~]=svd(XX0, 'econ');
    Ur = U(:, 1:r);
    sumU=[sumU,Ur];
end
Tru=15;
[U, ~, ~] = svd(sumU, 'econ');
Phi= U(:, 1:Tru);

% B_hat=U_pod'*M*U_pod;
u0_hat=Phi'*U0(T.FNodePtrs);
% F_hat=U_pod'*M;
u0_hat1=Phi*u0_hat;
error=norm(u0_hat1-U0(T.FNodePtrs));

uvp_ono=cell(5,1);
% av=Usnap(:,ind);
av=Phi;


% POD
nFG=av'*Fx(T.FNodePtrs,:,:);%%BD
A1g=zeros(Tru,Tru,size(ks,2));nFg=zeros(Tru,size(ks,2));
% nFg0=zeros(Tru,1);A1g0=zeros(Tru,Tru,1);
nFg0=av'*K0(T.FNodePtrs,:)*G;A1g0=av'*K0(T.FNodePtrs,T.FNodePtrs)*av;
for i=1:size(ks,2)
    nFg(:,i)=av'*Axq{i}(T.FNodePtrs,:)*G;
    A1g(:,:,i)=av'*Axq{i}(T.FNodePtrs,T.FNodePtrs)*av;
end
uvp_ono(:,1)={nFG,nFg,A1g,nFg0,A1g0};
% uvp_ono(:,1)={nFG,nFg,A1g};
KQVo=av; UUo=Usnap;
M_hat=av'*M(T.FNodePtrs,T.FNodePtrs)*av;
toc

% %% ROM online]
Ntest=1e2;
% wo=2*rand(P,Ntest)-1;
 load data_w W Ufem
 wo=W;
 tic
% [nFG,nFg,A1g,nFg0,A1g0]=deal(uvp_ono{:,count});
Uscofo=online_RBpod(wo',nFG,nFg,A1g,nFg0,A1g0,fs,ks,M_hat,ht,Nt,u0_hat);
toc
tic
% % V=snapsolution(wo,ks,fs,T,g,Axq,Fx,K0);
% V=ufem(wo,T,Nxo,Nyo,f,g,k);
V=ufem(wo,T,Nxo,Nyo,f,G,k,U0,t,M,Te);
% V=Ufem;
toc
Nw=Tru;
coefko=zeros(size(wo',1),Nt+1,Nw);
Uvo=[];
uu=repmat(G,1,Nt+1,Nw);
uu(:,1,:)=repmat(U0,1,1,Nw);
for j=1:size(wo',1)
% uu=zeros(m,Nt+1,Nw);
for jj=2:Nt+1
        for i=1:Nw
            uu(T.FNodePtrs,jj,i)=KQVo(:,1:i)*Uscofo{j}(1:i,jj);
        Vtemo=uu(:,jj,i)-V{j}(:,jj);
        coefko(j,jj,i)=sqrt((Vtemo(T.FNodePtrs)'*M(T.FNodePtrs,T.FNodePtrs)*Vtemo(T.FNodePtrs))/...
            (V{j}(T.FNodePtrs,jj)'*M(T.FNodePtrs,T.FNodePtrs)*V{j}(T.FNodePtrs,jj)));
        end
end
Uvo{j}=uu;
end
eer=mean(coefko,1);

y = reshape(mean(coefko,1),Nt+1,Nw);
e=zeros(Nt+1,Nw);
for j=2:Nt+1
    e_std=reshape(coefko(:,j,:),Ntest,Tru);
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
semilogy(1:size(coefko,1),coefko(:,end,1),'g-')
subplot(2,2,2)
semilogy(1:size(coefko,1),coefko(:,end,2),'r-')
subplot(2,2,3)
semilogy(1:size(coefko,1),coefko(:,end,3),'b-')
subplot(2,2,4)
semilogy(1:size(coefko,1),coefko(:,end,Nw),'b-')
xlabel('Numbers of the separated terms');  
ylabel('Mean error with error bars');
grid on
title ('u')

l=10;
figure(4)
subplot(3,3,1)
mesh(xo,yo,reshape(Uvo{l}(:,1,1),Nyo+1,Nxo+1))
title('POD solution','fontsize',12)
colorbar
subplot(3,3,2)
mesh(xo,yo,reshape(Uvo{l}(:,1,Nw),Nyo+1,Nxo+1))
title('POD solution','fontsize',12)
colorbar
subplot(3,3,3)
mesh(xo,yo,reshape(V{l}(:,1),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)
colorbar
subplot(3,3,4)
mesh(xo,yo,reshape(Uvo{l}(:,Nt/5,1),Nyo+1,Nxo+1))
title('POD solution','fontsize',12)
colorbar
subplot(3,3,5)
mesh(xo,yo,reshape(Uvo{l}(:,Nt/5,Nw),Nyo+1,Nxo+1))
title('POD solution','fontsize',12)
colorbar
subplot(3,3,6)
mesh(xo,yo,reshape(V{l}(:,Nt/5),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)
colorbar
subplot(3,3,7)
mesh(xo,yo,reshape(Uvo{l}(:,end,1),Nyo+1,Nxo+1))
title('POD solution','fontsize',12)
colorbar
subplot(3,3,8)
mesh(xo,yo,reshape(Uvo{l}(:,end,Nw),Nyo+1,Nxo+1))
title('POD solution','fontsize',12)
colorbar
subplot(3,3,9)
mesh(xo,yo,reshape(V{l}(:,end),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)
colorbar


