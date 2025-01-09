  clear,clc,close all  
tic
[sig,P,Nw,Nte]=deal(0.01,25,60,1e3);

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

%%%%KLE expasion for permeability
xmo=1/Nxo/2:1/Nxo:1;
ymo=1/Nyo/2:1/Nyo:1;
[Nxm,Nym]=deal(20,20);
c1=[0.1 0.2]; sigma1=1;
% Fai1=KL_expansion(xo,yo,P,co,sigma1,c1,Nxm,Nym);
% % % 
% % % %%%%Choose a proper mu1 to ensure k>0
% % % % kai_mu=@(x) exp(x(:,1)+x(:,2));
% % % % kai=@(x,w) x(:,1).*w(:,1)+x(:,2).*w(:,2)+x(:,1).*x(:,2).*w(:,3)+exp(x(:,1)+x(:,2));
% % % % kai1=@(x,w) x(:,1).*w(:,1)+x(:,2).*w(:,2)+x(:,1).*x(:,2).*w(:,3);
% mu=1.*ones(size(co,1),1);
% kai=@(w) Fai1*w+exp(mu);
% kai1=@(w) Fai1*w;
% % % % w=2*rand(P,1)-1;
% w=randn(P1,1);
% P=P1+1;
% % Fai=[exp(mu),Fai1];
% Fai=[exp(mu),Fai1];

% % %k函数
tic
P1=2;kx1=cell(1,P1);
% kx0=@(x) 1+0.*x(:,1);
kx0=@(x) exp(x(:,1)+x(:,2));
kx1{1}=@(x,s) (1+1.*cos(pi.*s(:,1).*(x(:,1)+x(:,2)))).*(1+1.*sin(1.*pi.*s(:,1).*(x(:,2)-3.*x(:,1))));
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
load data_elliptic_eff_noaffine_20240222 wo Ufem
k0=kx0(co);
% f=@(x) 1+0.*x(:,1);

% % %f函数
% fs{1}=@(s) 10*sin(s(:,1).*s(:,end));
% fx{1}=@(x) exp(x(:,1)+x(:,2)+3);
% fs{2}=@(s) 10*cos(s(:,1).*s(:,end));
% fx{2}=@(x) exp(x(:,1).*x(:,2));
% f=@(x,s) fs{1}(s)*fx{1}(x);
% fss=@(s) fs{1}(s);
% for i=2:size(fs,2)
%     f=@(x,s) f(x,s)+fs{i}(s)*fx{i}(x);
%      fss=@(s) [fss(s) fs{i}(s)];
% end

% fs{1}=@(s) sin(s(:,1).*s(:,end));
% fx{1}=@(x) x(:,2)+3;
% f=@(x,s) fs{1}(s)*fx{1}(x);
% fss=@(s) fs{1}(s);
% for i=2:size(fs,2)
%     f=@(x,s) f(x,s)+fs{i}(s)*fx{i}(x);
%     fss=@(s) [fss(s) fs{i}(s)];
% end
f=@(x,s) 10*sin(s(:,1)*x(:,2)+3.*s(:,end));


% wo=2*rand(P,Nte)-1;W=wo;W_mu=mean(W,2);
Nte=size(wo,2);
% % W=2*rand(P,Nte)-1;
Nte1=1e3;
Wo1=2*rand(P,Nte1)-1;
Wt=zeros(size(co,1),Nte1);
for i=1:Nte1
Wt(:,i)=ki(co,Wo1(:,i)');
end
W_mu=mean(Wt,2);
W_train=Wt-repmat(W_mu,1,Nte1);
W=zeros(size(co,1),Nte);
for i=1:Nte
W(:,i)=ki(co,wo(:,i)');
end
W_test=W-repmat(W_mu,1,Nte);


com=T.Nodes(bdy,:);
m=length(T.FNodePtrs);
[I,J,Val_k,Val_b]=stiff(T,Nxo,Nyo);

 K0=sparse(I,J,Val_k.*repmat(W_mu,16,1));
% K0=sparse(I,J,Val_k.*repmat(ki(co,W_mu'),16,1));
Am=sparse(I,J,Val_k);
M=sparse(I,J,Val_b);
%获取列数，含空间的刚度矩阵
G=zeros(size(T.Nodes,1),1);%G边界条件
G(T.CNodePtrs)=g(T.Nodes(T.CNodePtrs,:));
Fg=K0(T.FNodePtrs,:)*G;

G0=zeros(size(T.Nodes,1),1,Nte);
G0(T.CNodePtrs,:,:)=repmat(g(T.Nodes(T.CNodePtrs,:)),1,1,Nte);
toc

tic
[Uromo,err_ef]=ROM_surro_Noaffine_20240205(f,T,W_test,wo,Fg,G0,M,I,J,Val_k,K0,ki,co);
toc
Urom=zeros(size(G0));
Urom(T.CNodePtrs,:)=G0(T.CNodePtrs,:);
Urom(T.FNodePtrs,:)=Uromo{end};

% tic
% Ufem=ufem(wo,T,Nxo,Nyo,f,g,ki);
% % % Ufem=snapsolution(W,ks,fs,T,g,Axq,Fxq,K0);
% % Ufem=V;
% toc

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
figure(2)
errorbar(y,e,'xr')
xlabel('Numbers of the iteration terms');  
ylabel('Mean error with error bars');
grid on

% figure(2)
% plot(1:size(err_ef,1),mean(err_ef,2),'r-')
% xlabel('Numbers of the iteration steps');  
% ylabel('L_max error ');
% % legend('N=1','N_p=3','N_p=5') ;
% grid on

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
% ylabel('Errors (7th iteration)');
grid on
% title ('u')

close 
l=10;
figure(5)
subplot(1,2,1)
mesh(xo,yo,reshape(Urom(:,l),Nyo+1,Nxo+1))
title('ROM solution','fontsize',12)
hold on
subplot(1,2,2)
mesh(xo,yo,reshape(Ufem(:,l),Nyo+1,Nxo+1))
title('FEM solution','fontsize',12)

% save data_elliptic_eff_noaffine_20240222 tim_effi Uromo Ufem coefko wo count