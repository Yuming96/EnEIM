clear,clc,close all  
[sig,P,Nw]=deal(0.01,25,100);
Nxo=100;Nyo=100;
%  Nxo=40;Nyo=40;
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
co=T.centriod;
bdy=T.CNodePtrs;
[Nxm,Nym]=deal(20,20);

% % %k函数
 tic
P1=2;kx1=cell(1,P1);
% kx0=@(x) 1+0.*x(:,1);
kx0=@(x) exp(x(:,1)+x(:,2));
kx1{1}=@(x,s) (1+1.*cos(pi.*s(:,1).*(x(:,1)+x(:,2)))).*(1+sin(1.*pi.*s(:,1).*(x(:,2)-3.*x(:,1))));
for i=2:P1
    kx1{i}=@(x,s) (1+1.*cos(pi.*s(:,i).*(x(:,1)+x(:,2)))).*(1+1.*sin(pi.*s(:,i).*(x(:,2)-3.*x(:,1))));
end
k1=@(x,s) kx1{1}(x,s);
for i=2:P1
    k1=@(x,s) k1(x,s)+kx1{i}(x,s);
end
ki=@(x,s) kx0(x)+k1(x,s);
Ntr=50;
P=P1;
% s0=2*rand(Ntr,P)-1;
% s01=2*rand(1,P)-1;
load data_elliptic_pod_noaffine W_train wo s0 s01
load data_elliptic_eff_noaffine Ufem
% 
[kx,ks,k]=Affine1(ki,s0,s01,co,xo,yo,Nxm,Nym);
k0=@(x) 0.*x(:,1);

fi=@(x,s) 10*sin(s(:,1)*x(:,2)+3.*s(:,end));
[fx,fs,f]=Affine1(fi,s0,s01,co,xo,yo,Nxm,Nym);
com=T.Nodes(bdy,:);
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
 Ntrain=14;count=1;
W_train=W_train(:,1:14);
Usnap=ufem(W_train,T,Nxo,Nyo,f,g,ki);
Usnap=real(Usnap);
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
kk=sum(W_train);
u=@(s) fs1(s)./kss(s);uvp_ono(:,1)={nFG,nFg,A1g,nFg0,A1g0};
KQVo(:,count)=av; UUo=Usnap;
W_train(:,ind)=[]; UUo(:,ind)=[];
%随机基函数
VSo=zeros(size(W_train,2),1);
for i=1:size(W_train,2)
VSo(i)=u(W_train(:,i)');
end
iex=ind;w1o=W_train;

elsonk=1e-5;
while el> elsonk & count <Nw
    count=count+1;
    %贪婪算法
    [KQVo,VSo,UUo,w1o,uvpono,el,iex] = stepn_ortho_greedy(w1o,M,Axq,KQVo,T.FNodePtrs,Fx,UUo,VSo,fs,ks,G,iex,K0);
    uvp_ono(:,count)=uvpono;
    EL(count)=el;
end
toc

%% ROM online]
Ntest=1e3;
% % load data_effi wo
% wo=1*rand(size(ks,2),Ntest);
tic
[nFG,nFg,A1g,nFg0,A1g0]=deal(uvp_ono{:,count});
Uscofo=online_RBgreedy(wo',nFG,nFg,A1g,nFg0,A1g0,fs,ks);
Uscofo=real(Uscofo);
toc
% tic
% % % V=snapsolution(wo,ks,fs,T,g,Axq,Fx,K0);
% % V=ufem(wo,T,Nxo,Nyo,f,g,ki);
V=Ufem;
% toc
Nw=count;
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
semilogy(1:count,eer,'r-')
xlabel('Numbers of the separated terms');  
ylabel('Relative error ');
% legend('N=1','N_p=3','N_p=5') ;
grid on

figure(3)
semilogy(1:size(coefko,1),coefko(:,1),'g-')
hold on
semilogy(1:size(coefko,1),coefko(:,4),'r-')
hold on
semilogy(1:size(coefko,1),coefko(:,count),'b-')
xlabel('Numbers of the separated terms');  
ylabel('Mean error with error bars');
grid on
title ('u')

l=10;
figure(4)
subplot(1,2,1)
mesh(xo,yo,reshape(Uvo(:,l),Nyo+1,Nxo+1))
title('GA solution','fontsize',12)
hold on
subplot(1,2,2)
mesh(xo,yo,reshape(V(:,l),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)

% save data_elliptic_greedy_noaffine tim_greedy KQVo Uscofo V coefko wo W_train count s0 s01