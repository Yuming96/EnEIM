  clear,clc,close all  
% tic
[sig,P,Nw,Nte]=deal(0.01,25,60,1e3);

Nxo=100;Nyo=100;
bd=[0 1;0 1];
hx=(bd(1,2)-bd(1,1))/Nxo;
hy=(bd(2,2)-bd(2,1))/Nyo;

xo=bd(1,1):hx:bd(1,2);yo=bd(2,1):hy:bd(2,2);
T=struc(xo,yo);
g=@(x) x(:,2).*x(:,1).*(x(:,1)-1).*(x(:,2)-1);
u0=@(x) 1.*x(:,1);
[X,Y]=meshgrid(xo,yo);
No=[X(:) Y(:)];
ke=@(x)  0.*x(:,1).*x(:,2)+1;

%%%%use for inverse problem, the observation locations 
co=T.centriod;
bdy=T.CNodePtrs;

xmo=1/Nxo/2:1/Nxo:1;
ymo=1/Nyo/2:1/Nyo:1;
[Nxm,Nym]=deal(20,20);
c1=[0.1 0.2]; sigma1=1;

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

com=T.Nodes(bdy,:);
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
Fg=K0(T.FNodePtrs,:)*G;
Axq=cell(1,size(ks,2));
for i=1:size(ks,2)
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

%% ROM offline
%  load data wo
load data_elliptic_eff_affine1 wo V
% W=2*rand(size(ks,2),Nte)-1;
W=wo;
% W_train=1*randn(P,Ntrain);
G0=zeros(size(T.Nodes,1),Nte);
G0(T.CNodePtrs,:)=repmat(g(T.Nodes(T.CNodePtrs,:)),1,Nte);
% U(T.CNodePtrs,:)=repmat(G(T.CNodePtrs),1,Ntr);

tic
[Uromo,err_ef]=ROM_surro_affine(ks,fs,T,W,Axq,Fxq,Fg,G0,K0);
toc


Urom=zeros(size(G0));
Urom(T.CNodePtrs,:)=G0(T.CNodePtrs,:);
Urom(T.FNodePtrs,:)=Uromo{end};

% tic
% Ufem=ufem(W,T,Nxo,Nyo,f,g,k);
% % Ufem=snapsolution(W,ks,fs,T,g,Axq,Fxq,K0);
% toc
Ufem=V;

error_te=zeros(Nte,1);
for i=1:Nte
    error_te(i)=norm(Urom(:,i)-Ufem(:,i))/norm(Ufem(:,i));
end

mean(error_te,1)
figure(1)
semilogy(1:Nte,error_te)

count=size(err_ef,1);
coefko=zeros(size(W',1),count);
for i=1:count
 for j=1:size(W',1)
            Vtemo=Uromo{i}(:,j)-Ufem(T.FNodePtrs,j);
        coefko(j,i)=sqrt((Vtemo'*M(T.FNodePtrs,T.FNodePtrs)*Vtemo)/...
            (Ufem(T.FNodePtrs,j)'*M(T.FNodePtrs,T.FNodePtrs)*Ufem(T.FNodePtrs,j)));
 end
end

eer=mean(coefko,1)
y = mean(coefko,1);e = std(coefko,1,1);
close all
figure(2)
errorbar(y,e,'xr')
xlabel('Numbers of the iteration terms');  
ylabel('Mean error with error bars');
grid on


figure(3)
plot(1:count,eer,'r-')
xlabel('Numbers of the iteration terms');  
ylabel('Relative error ');
% legend('N=1','N_p=3','N_p=5') ;
grid on

figure(4)
% subplot(3,1,1)
semilogy(1:size(coefko,1),coefko(:,1),'g-')
% ylabel('Errors (1th iteration)');
hold on
% subplot(3,1,2)
semilogy(1:size(coefko,1),coefko(:,3),'r-')
% ylabel('Errors (3th iteration)');
hold on
% subplot(3,1,3)
semilogy(1:size(coefko,1),coefko(:,count),'b-')
xlabel('Numbers of the samples');  
ylabel('Errors (7th iteration)');
grid on
% title ('u')

l=4;
figure(5)
subplot(1,2,1)
mesh(xo,yo,reshape(Urom(:,l),Nyo+1,Nxo+1))
title('ROM solution','fontsize',12)
hold on
subplot(1,2,2)
mesh(xo,yo,reshape(Ufem(:,l),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)

% save data_elliptic_eff_affine_20240207 tim_effi Uromo V coefko wo count