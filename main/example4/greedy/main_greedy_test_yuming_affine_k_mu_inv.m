clear,clc,close all  
[sig,P,Nw]=deal(0.01,25,60);
Nxo=100;Nyo=100;
%  Nxo=40;Nyo=40;
bd=[0 1;0 1];
hx=(bd(1,2)-bd(1,1))/Nxo;
hy=(bd(2,2)-bd(2,1))/Nyo;

xo=bd(1,1):hx:bd(1,2);yo=bd(2,1):hy:bd(2,2);
%  f=@(x,w) 1+0.*x(:,1).*w(:,1);
T=struc(xo,yo);
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
figure(1)
bdy=T.CNodePtrs;
plot(Nodes(inter(p),1),Nodes(inter(p),2),'pr')
axis([0 1 0 1])
U0=u0(Nodes(T.FNodePtrs,:));
tic
% %%%%KLE expasion for permeability
xmo=1/Nxo/2:1/Nxo:1;
ymo=1/Nyo/2:1/Nyo:1;
[Nxm,Nym]=deal(20,20);
c1=[0.1 0.2]; sigma1=1;
kx=KL_expansion(xo,yo,P,co,sigma1,c1,Nxm,Nym);

mu=1.*ones(size(co,1),1);
load data_ell_pod_affine_inv_k_mu2 Uobs p w_true Q kx W_train
P=length(w_true);
% tic
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
mesh(ymo,xmo,reshape(k(co,w_true'),Nyo,Nxo))
title('Reference permeability field')
colorbar

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
Fx=zeros(m,size(fs,2));
%源项的向量
Fx(:,1)=M*fx{1}(T.Nodes);
% f=@(s) fs{1}(s)*fx{1}(T.Nodes);
for i=2:size(fs,2)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    Fx(:,i)=M*fx{i}(T.Nodes);
%     f=@(s) f(s)+fs{i}(s)*fx{i}(T.Nodes);
end

%% greedy offline
% load data_ell_pod_affine_inv Uobs p w_true Q W_train
 Ntrain=2e2;count=1;
% W_train=2*rand(size(ks,2),Ntrain)-1;
% % W_train=1*randn(size(ks,2),Ntrain);
W_train0=W_train(:,1:Ntrain);
P=length(w_true);

% tic
% [Usnap]=snapsolution(W_train,ks,fs,T,g,Axq,Fx,K0);
p1=[];
Usnap=ufem(W_train0,T,Nxo,Nyo,f,g,k,p1);
% toc

%所有样本下H的半模
ELt=zeros( Ntrain,1);
for ji=1: Ntrain
    ELt(ji)=Usnap(T.FNodePtrs,ji)'*Am(T.FNodePtrs,T.FNodePtrs)*Usnap(T.FNodePtrs,ji);
end
%残差最大范数下的g1
[el,ind]=max(ELt);EL(count)=el;
uvp_ono=cell(5,Nw);
av=Usnap(:,ind);
G=zeros(size(T.Nodes,1),1);
G(T.CNodePtrs)=g(T.Nodes(T.CNodePtrs,:));
%载荷向量中v取g1
nF=(Fx'*av);
fs1=@(s) fs{1}(s)*nF(1);
for i=2:size(fs,2)
   fs1 =@(s) fs1(s)+fs{i}(s)*nF(i);
end
kss=@(s) 0;
a1=zeros(1,size(ks,2));A2=zeros(1,size(ks,2));
% %VS
for i=1:size(ks,2)
    %刚度阵V取给g1
    a1(i)=av(T.FNodePtrs)'*Ax{i}*av(T.FNodePtrs);A2(i)=av(T.FNodePtrs)'*Axq{i}(T.FNodePtrs,:)*G;
    kss=@(s) kss(s)+a1(i)*ks{i}(s);
    fs1=@(s) fs1(s)-A2(i)*ks{i}(s);
end
% 贪婪
nFG=av'*Fx;%%BD
A1g=zeros(count,count,size(ks,2));nFg=zeros(count,size(ks,2));nFg0=zeros(count,1);A1g0=zeros(count,count,1);
nFg0(count)=av(T.FNodePtrs)'*K0(T.FNodePtrs,:)*G;A1g0(:,:,1)=av(T.FNodePtrs)'*K0(T.FNodePtrs,T.FNodePtrs)*av(T.FNodePtrs);
for i=1:size(ks,2)
    nFg(:,i)=av(T.FNodePtrs)'*Axq{i}(T.FNodePtrs,:)*G;
    A1g(:,:,i)=av(T.FNodePtrs)'*Axq{i}(T.FNodePtrs,T.FNodePtrs)*av(T.FNodePtrs);
end
%含随机变量基函数
kk=sum(W_train0);
u=@(s) fs1(s)./kss(s);uvp_ono(:,1)={nFG,nFg,A1g,nFg0,A1g0};
KQVo(:,count)=av; UUo=Usnap;
W_train0(:,ind)=[]; UUo(:,ind)=[];
%随机基函数
VSo=u(W_train0');
iex=ind;w1o=W_train0;

elsonk=1e-6;
while el> elsonk & count <Nw
    count=count+1;
    %贪婪算法
    [KQVo,VSo,UUo,w1o,uvpono,el,iex] = stepn_ortho_greedy(w1o,M,Axq,KQVo,T.FNodePtrs,Fx,UUo,VSo,fs,ks,G,iex,K0);
    uvp_ono(:,count)=uvpono;
    EL(count)=el;
end
toc

[nFG,nFg,A1g,nFg0,A1g0]=deal(uvp_ono{:,count});
%% 测试集
tic
Ntest=size(Q{1},2);
W_test=Q{1};
fw=@(w)  online_RBgreedy(w',nFG,nFg,A1g,nFg0,A1g0,fs,ks);
[Q,Usim,error_w] = SEnKF_sur_re_affine_k_mu(fw,sig,Ntest,W_test,p,Uobs,k,co,KQVo);

for i=1:length(Q)
w_mean=mean(Q{i},2);
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
 ah0=Q{mm}(ii,:);
 if ii==1
[fh,xhi] = ksdensity(ah0,'width',0.0051);
 elseif ii==2
     [fh,xhi] = ksdensity(ah0,'width',0.0047);
      elseif ii==3
     [fh,xhi] = ksdensity(ah0,'width',0.012);
      elseif ii==P
     [fh,xhi] = ksdensity(ah0,'width',0.007);
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

% save data_ell_greedy_affine_inv_k_mu_0319 tim_greedy KQVo Q p w_true Uobs W_train0 Usim Nw