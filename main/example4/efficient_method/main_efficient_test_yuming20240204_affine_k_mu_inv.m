  clear,clc,close all  
% tic
[sig,P,Nw,Nte]=deal(0.01,15,60,1e2);

Nxo=100;Nyo=100;
 Nx=Nxo*2;Ny=Nyo*2;
bd=[0 1;0 1];
hx=(bd(1,2)-bd(1,1))/Nxo;
hy=(bd(2,2)-bd(2,1))/Nyo;
hxc=(bd(1,2)-bd(1,1))/Nx;
hyc=(bd(2,2)-bd(2,1))/Ny;

xo=bd(1,1):hx:bd(1,2);yo=bd(2,1):hy:bd(2,2);
xoc=bd(1,1):hxc:bd(1,2);yoc=bd(2,1):hyc:bd(2,2);
T=struc(xo,yo);
Tc=struc(xoc,yoc);
g=@(x) x(:,2).*x(:,1).*(x(:,1)-1).*(x(:,2)-1);
% g=@(x) 0.*(x(:,1)-1).*(x(:,2)-1)+1;
% u0=@(x) sin(2*pi.*x(:,1)).*sin(2*pi.*x(:,2))+0;
u0=@(x) 1.*x(:,1);
[X,Y]=meshgrid(xo,yo);
No=[X(:) Y(:)];
ke=@(x)  0.*x(:,1).*x(:,2)+1;

%%%%use for inverse problem, the observation locations 
[p,Nodes,inter,co,~]=point_numberf(Nxo,Nyo,T);
% co=T.centriod;
bdy=T.CNodePtrs;
figure(1)
subplot(1,2,1)
plot(Nodes(inter(p),1),Nodes(inter(p),2),'pr')
axis([0 1 0 1])

% %%%%use for inverse problem, the observation locations 
[pc,Nodesc,interc,co1,~]=point_numberf(Nx,Ny,Tc);
tic
%%%%KLE expasion for permeability
xmo=1/Nxo/2:1/Nxo:1;
ymo=1/Nyo/2:1/Nyo:1;
[Nxm,Nym]=deal(20,20);
c1=[0.1 0.2]; sigma1=1;
kx=KL_expansion(xo,yo,P,co,sigma1,c1,Nxm,Nym);

% % % %%%%Choose a proper mu1 to ensure k>0
mu=1.*ones(size(co,1),1);
muc=1.*ones(4*size(co,1),1);
w_true=randn(P,1);
% load data_eff_affine_k_mu_inv Uobs p w_true Q kx
P=length(w_true);

% % % %k函数
k0=@(x) 0.*x(:,1)+exp(mu);
for i=1:P
ks{i}=@(s) s(:,i);
% kx{i}=@(x) 0.*x(:,1)+Fai1(:,i);
end

kx{1}=@(x) kx{1}(x);
k=@(x,s) k0(x)+ks{1}(s).*kx{1}(x);
% k=@(x,s) k0(x)+ks{1}(s)*Fai1(:,1);
kss=@(s) ks{1}(s);
for i=2:size(ks,2)
    kx{i}=@(x) kx{i}(x);
    k=@(x,s) k(x,s)+ks{i}(s)*kx{i}(x);
%     k=@(x,s) k(x,s)+ks{i}(s)*Fai1(:,i);
    kss=@(s) [kss(s) ks{i}(s)];
end
% 
k0c=@(x) 0.*x(:,1)+exp(muc);
kx{1}=@(x) kx{1}(x);
kc=@(x,s) k0c(x)+ks{1}(s).*kx{1}(x);
for i=2:size(ks,2)
    kx{i}=@(x) kx{i}(x);
    kc=@(x,s) kc(x,s)+ks{i}(s)*kx{i}(x);
end

% %f函数
fs{1}=@(s) 10*sin(s(:,1).*s(:,end));
fx{1}=@(x) exp(x(:,1)+x(:,2)+3);
fs{2}=@(s) 10*cos(s(:,1).*s(:,end));
fx{2}=@(x) exp(x(:,1).*x(:,2));
f=@(x,s) fs{1}(s)*fx{1}(x);
fss=@(s) fs{1}(s);
for i=2:size(fs,2)
    f=@(x,s) f(x,s)+fs{i}(s)*fx{i}(x);
     fss=@(s) [fss(s) fs{i}(s)];
end

figure(2)
% subplot(1,2,1)
mesh(ymo,xmo,reshape(k(co,w_true'),Nyo,Nxo))
title('Reference permeability field')
colorbar

% p=[];pc=[];
Ufc=ufem(w_true,Tc,Nx,Ny,f,g,kc,pc);
Uf=ufem(w_true,T,Nxo,Nyo,f,g,k,p);
Uobs=Ufc+sig.*randn(length(pc),1);
tic
com=T.Nodes(bdy,:);
% G=g(com);
% F=loaf(T,f,Nxo,Nyo);
m=length(T.FNodePtrs);
% [If,Jf,Val_f]=loaf(T,Nxo,Nyo);
[I,J,Val_k,Val_b]=stiff(T,Nxo,Nyo);
% A0=exp(mu);
 K0=sparse(I,J,Val_k.*repmat(k0(co),16,1));
% K0=sparse(I,J,Val_k.*repmat(exp(mu),16,1));
Am=sparse(I,J,Val_k);
M=sparse(I,J,Val_b);
%获取列数，含空间的刚度矩阵
G=zeros(size(T.Nodes,1),1);%G边界条件
G(T.CNodePtrs)=g(T.Nodes(T.CNodePtrs,:));
% Kf=K0(T.FNodePtrs,T.FNodePtrs);
% Fg=Kf\K0(T.FNodePtrs,:)*G;
Fg=K0(T.FNodePtrs,:)*G;
Axq=cell(1,size(ks,2));
for i=1:size(ks,2)
%     tem=sparse(I,J,Val_k.*repmat(Fai1(:,i),16,1));
   tem=sparse(I,J,Val_k.*repmat(kx{i}(co),16,1));
    Axq{i}=tem;
end

% Fx=zeros(m,size(fs,2));
Fxq=zeros(size(K0,1),size(fs,2));
%源项的向量
termf=M*fx{1}(T.Nodes);
% Fx(:,1)=Kf\termf(T.FNodePtrs);
Fxq(:,1)=termf;
f1=@(s) fs{1}(s)*fx{1}(T.Nodes);
for i=2:size(fs,2) 
    termf=M*fx{i}(T.Nodes);Fxq(:,i)=termf;
%     Fx(:,i)=Kf\termf(T.FNodePtrs);
    f1=@(s) f(s)+fs{i}(s)*fx{i}(T.Nodes);
end
toc
% %% ROM offline
% %  load data wo
% load data_elliptic_eff_affine wo V
% W=2*rand(size(ks,2),Nte)-1;
% W=wo;
W=1*randn(P,Nte);
G0=zeros(size(T.Nodes,1),Nte);
G0(T.CNodePtrs,:)=repmat(g(T.Nodes(T.CNodePtrs,:)),1,Nte);
% U(T.CNodePtrs,:)=repmat(G(T.CNodePtrs),1,Ntr);

%% 测试集
tic
Ntest=Nte;
W_test=randn(P,Ntest);
Iter=8;
fw=@(w)  ROM_surro_affine(ks,fs,T,w,Axq,Fxq,Fg,G0,K0,Iter);
[Q,Usim,error_w] = SEnKF_sur_re_affine_k_mu(fw,sig,Ntest,W_test,p,Uobs,k,co);

for i=1:length(Q)
w_mean=mean((Q{i}),2);
er_w(i)=norm(w_true-w_mean)/norm(w_true);
er_k(i)=norm(k(co,w_true')-k(co,w_mean'))/norm(k(co,w_true'));
end
toc
% end

mm=length(Q);
w_mean=mean(Q{mm},2);
figure(2)
colormap(jet)
subplot(1,2,1)
mesh(ymo,xmo,reshape(k(co,w_true'),Nyo,Nxo))
colorbar
title('Reference')
hold on
subplot(1,2,2)
mesh(ymo,xmo,reshape(k(co,w_mean'),Nyo,Nxo))
colorbar
title('Estimate')

for ii=1:P
% ii=1;
figure(ii+2)
% ii=1;
% a0=atan(Q{end}(1,:))*4/pi;
 ah0=(Q{mm}(ii,:));
 if ii==1
[fh,xhi] = ksdensity(ah0,'width',0.0051);
 elseif ii==2
     [fh,xhi] = ksdensity(ah0,'width',0.0067);
      elseif ii==3
     [fh,xhi] = ksdensity(ah0,'width',0.0092);
      elseif ii==P
     [fh,xhi] = ksdensity(ah0,'width',0.013);
 end
[nh,xh]=hist(ah0,10);
nh=nh/size(ah0,2);
bar(xh,nh,1)
% hold on
% plot(xhi,fh*(xh(2)-xh(1)),'r--','Linewidth',2)
hold on
plot(w_true(ii),0,'r*','Markersize',8)
hold off
end

figure(P+3)
semilogy(1:length(Q),er_w)
hold on
semilogy(1:length(Q),er_k)

% save data_eff_affine_k_mu_inv tim_eff Q p w_true Nte Uobs Usim Uf kx