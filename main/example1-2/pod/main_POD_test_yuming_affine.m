clear,clc,close all  
[sig,P,Nw]=deal(0.01,25,100);
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
% [~,Nodes,inter,co,~]=point_numberf(Nxo,Nyo,T);
co=T.centriod;
bdy=T.CNodePtrs;
% plot(Nodes(inter(p),1),Nodes(inter(p),2),'pr')
% axis([0 1 0 1])
% U0=u0(Nodes(T.FNodePtrs,:));

% %%%%KLE expasion for permeability
% xmo=1/Nxo/2:1/Nxo:1;
% ymo=1/Nyo/2:1/Nyo:1;
% [Nxm,Nym]=deal(20,20);
% c1=[0.1 0.2]; sigma1=1;
% Fai1=KL_expansion(xo,yo,P,co,sigma1,c1,Nxm,Nym);
% % 
% % %%%%Choose a proper mu1 to ensure k>0
% % % kai_mu=@(x) exp(x(:,1)+x(:,2));
% % % kai=@(x,w) x(:,1).*w(:,1)+x(:,2).*w(:,2)+x(:,1).*x(:,2).*w(:,3)+exp(x(:,1)+x(:,2));
% % % kai1=@(x,w) x(:,1).*w(:,1)+x(:,2).*w(:,2)+x(:,1).*x(:,2).*w(:,3);
% mu=1.*ones(size(co,1),1);
% kai=@(w) Fai1*w+exp(mu);
% kai1=@(w) Fai1*w;
% % % % % w=2*rand(P,1)-1;
% % w=randn(P1,1);
% % P=P1+1;
% % Fai=[exp(mu),Fai1];

tic
%k函数
% % k0=@(x) exp(x(:,1)+x(:,2));
% % ks{1}=@(s) sin(s(:,1))+1.5;
% % kx{1}=@(x) x(:,2);
% % ks{2}=@(s) exp(s(:,2)/10);
% % kx{2}=@(x) x(:,1).*x(:,2);
% % ks{3}=@(s) s(:,3)+1.5;
% % kx{3}=@(x) x(:,1);
% % % ks{1}=@(s) 1;
% k0=@(x) 0.*x(:,1)+exp(mu);
% for i=1:P
% ks{i}=@(s) s(:,i);
% kx{1}=@(x) 0.*x(:,1)+Fai1(:,i);
% end
k0=@(x) exp(x(:,1)+x(:,2));
ks{1}=@(s) s(:,1);
kx{1}=@(x) x(:,1);
ks{2}=@(s) s(:,2);
kx{2}=@(x) x(:,2);
ks{3}=@(s) s(:,3);
kx{3}=@(x) x(:,1).*x(:,2);

% % kx{1}=@(x) kx{1}(x);
% k=@(x,s) k0(x)+ks{1}(s)*Fai1(:,1);
k=@(x,s) k0(x)+ks{1}(s)*kx{1}(x);
for i=2:size(ks,2)
%     kx{i}=@(x) kx{i}(x);
    k=@(x,s) k(x,s)+ks{i}(s)*kx{i}(x);
%     k=@(x,s) k(x,s)+ks{i}(s)*Fai1(:,i);
end

% %f函数
fs{1}=@(s) 10*sin(s(:,1).*s(:,end));
fx{1}=@(x) exp(x(:,1)+x(:,2)+3);
fs{2}=@(s) 10*cos(s(:,1).*s(:,end));
fx{2}=@(x) exp(x(:,1).*x(:,2));
f=@(x,s) fs{1}(s)*fx{1}(x);
for i=2:size(fs,2)
    f=@(x,s) f(x,s)+fs{i}(s)*fx{i}(x);
end
% fs{1}=@(s) sin(s(:,1).*s(:,end));
% fx{1}=@(x) x(:,2)+3;
% f=@(x,s) fs{1}(s)*fx{1}(x);
% for i=2:size(fs,2)
%     f=@(x,s) f(x,s)+fs{i}(s)*fx{i}(x);
% end


% % % KA1=cat(2,mu1,Fai1);
% f=@(x,w) (x(:,2)+3).*sin(w(:,1).*w(:,end));
% 
% figure
% mesh(ymo,xmo,reshape(kai(w),Nyo,Nxo))
% title('Reference permeability field')
% colorbar

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
load data_elliptic_eff_affine V wo
% load data_gree wo W_train1
 Ntrain=1e2;
W_train=2*rand(size(ks,2),Ntrain)-1;
% % W_train=1*randn(size(ks,2),Ntrain);
% W_train=W_train1;

%  tic
% [Usnap]=snapsolution(W_train,ks,fs,T,g,Axq,Fx,K0);
Usnap=ufem(W_train,T,Nxo,Nyo,f,g,k);
[A,Phi,Tru]=pod_svd_train(Usnap);
Upod=Phi*A';
error_train=zeros(Ntrain,1);
for i=1:Ntrain
    error_train(i)=norm(Upod(:,i)-Usnap(:,i))/norm(Usnap(:,i));
end
mean(error_train,1)
figure
semilogy(1:Ntrain,error_train)


% %所有样本下H的半模
% ELt=zeros( Ntrain,1);
% for ji=1: Ntrain
%     ELt(ji)=Usnap(T.FNodePtrs,ji)'*Am(T.FNodePtrs,T.FNodePtrs)*Usnap(T.FNodePtrs,ji);
% end
% %残差最大范数下的g1
% [el,ind]=max(ELt);EL(count)=el;
uvp_ono=cell(5,1);
% av=Usnap(:,ind);
av=Phi;
G=zeros(size(T.Nodes,1),1);
G(T.CNodePtrs)=g(T.Nodes(T.CNodePtrs,:));

% POD
nFG=av'*Fx;%%BD
A1g=zeros(Tru,Tru,size(ks,2));nFg=zeros(Tru,size(ks,2));
% nFg0=zeros(Tru,1);A1g0=zeros(Tru,Tru,1);
nFg0=av(T.FNodePtrs,:)'*K0(T.FNodePtrs,:)*G;A1g0=av(T.FNodePtrs,:)'*K0(T.FNodePtrs,T.FNodePtrs)*av(T.FNodePtrs,:);
for i=1:size(ks,2)
    nFg(:,i)=av(T.FNodePtrs,:)'*Axq{i}(T.FNodePtrs,:)*G;
    A1g(:,:,i)=av(T.FNodePtrs,:)'*Axq{i}(T.FNodePtrs,T.FNodePtrs)*av(T.FNodePtrs,:);
end
uvp_ono(:,1)={nFG,nFg,A1g,nFg0,A1g0};
% uvp_ono(:,1)={nFG,nFg,A1g};
KQVo=av; UUo=Usnap;
toc

%% ROM online]
Ntest=1e3;
% wo=2*rand(size(ks,2),Ntest)-1;
 tic
% % [nFG,nFg,A1g,nFg0,A1g0]=deal(uvp_ono{:,count});
Uscofo=online_RBpod(wo',nFG,nFg,A1g,nFg0,A1g0,fs,ks);
% toc
toc
% tic
% % % V=snapsolution(wo,ks,fs,T,g,Axq,Fx,K0);
% V=ufem(wo,T,Nxo,Nyo,f,g,k);
% toc
Nw=Tru;
coefko=zeros(size(wo',1),Nw);
for i=1:Nw
    Uvo=KQVo(:,1:i)*Uscofo(:,1:i)';
    for j=1:size(wo',1)
        Vtemo=Uvo(:,j)-V(:,j);
        coefko(j,i)=sqrt((Vtemo(T.FNodePtrs)'*M(T.FNodePtrs,T.FNodePtrs)*Vtemo(T.FNodePtrs))/...
            (V(T.FNodePtrs,j)'*M(T.FNodePtrs,T.FNodePtrs)*V(T.FNodePtrs,j)));
    end
end
eer=[mean(coefko,1)]

y = mean(coefko,1);e = std(coefko,1,1);
figure(1)
errorbar(y,e,'xr')
xlabel('Numbers of the separated terms');  
ylabel('Mean error with error bars');
grid on

figure(2)
semilogy(1:Tru,eer,'r-')
xlabel('Numbers of the iteration steps');  
ylabel('Relative error ');
% legend('N=1','N_p=3','N_p=5') ;
grid on


figure(3)
semilogy(1:size(coefko,1),coefko(:,1),'g-')
hold on
semilogy(1:size(coefko,1),coefko(:,2),'r-')
hold on
semilogy(1:size(coefko,1),coefko(:,Nw),'b-')
xlabel('Numbers of the separated terms');  
ylabel('Mean error with error bars');
grid on
title ('u')

l=10;
figure(4)
subplot(1,2,1)
mesh(xo,yo,reshape(Uvo(:,l),Nyo+1,Nxo+1))
title('POD solution','fontsize',12)
hold on
subplot(1,2,2)
mesh(xo,yo,reshape(V(:,l),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)

% save data_elliptic_pod_affine tim_fem tim_pod Phi KQVo Uscofo V coefko wo W_train