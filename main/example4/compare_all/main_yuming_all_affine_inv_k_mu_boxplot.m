  clear,clc,close all  
tic
[sig,P1,Nw,Nte]=deal(0.01,4,60,1e3);

Nxo=100;Nyo=100;
%  Nxo=20;Nyo=20;
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
Mesh = TProd_Mesh(xo,yo);

%%%%use for inverse problem, the observation locations 
[p,Nodes,inter,co,~]=point_numberf1(Nxo,Nyo,T);
% co=T.centriod;
figure(1)
subplot(1,2,1)
bdy=T.CNodePtrs;
plot(Nodes(inter(p),1),Nodes(inter(p),2),'pr')
axis([0 1 0 1])
xlabel('x')
ylabel('y')
title('Observation locations')
set(gca,'FontSize',22)  %是设置刻度字体大小

% %%%%KLE expasion for permeability
xmo=1/Nxo/2:1/Nxo:1;
ymo=1/Nyo/2:1/Nyo:1;
[Nxm,Nym]=deal(20,20);
c1=[0.1 0.2]; sigma1=1;
% kx=KL_expansion(xo,yo,P,co,sigma1,c1,Nxm,Nym);
% % % 
% % % %%%%Choose a proper mu1 to ensure k>0
% % % % kai_mu=@(x) exp(x(:,1)+x(:,2));
% % % % kai=@(x,w) x(:,1).*w(:,1)+x(:,2).*w(:,2)+x(:,1).*x(:,2).*w(:,3)+exp(x(:,1)+x(:,2));
% % % % kai1=@(x,w) x(:,1).*w(:,1)+x(:,2).*w(:,2)+x(:,1).*x(:,2).*w(:,3);
mu=1.*ones(size(co,1),1);
% kai=@(w) Fai1*w+exp(mu);
% kai1=@(w) Fai1*w;
% % % % w=2*rand(P,1)-1;
% w_true=randn(P,1);
% P=P1+1;
% % Fai=[exp(mu),Fai1];
% Fai=[exp(mu),Fai1];

load data_ell_pod_affine_inv_k_mu2 Uobs p w_true Q kx W_train
P=length(w_true);
tic
% % % %k函数
% k0=@(x) exp(x(:,1)+x(:,2));
% ks{1}=@(s) s(:,1);
% kx{1}=@(x) x(:,1);
% ks{2}=@(s) s(:,2);
% kx{2}=@(x) x(:,2);
% ks{3}=@(s) s(:,3);
% kx{3}=@(x) x(:,1).*x(:,2);
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
% fs{1}=@(s) sin(s(:,1).*s(:,end));
% fx{1}=@(x) x(:,2)+3;
% f=@(x,s) fs{1}(s)*fx{1}(x);
% fss=@(s) fs{1}(s);
% for i=2:size(fs,2)
%     f=@(x,s) f(x,s)+fs{i}(s)*fx{i}(x);
%     fss=@(s) [fss(s) fs{i}(s)];
% end

% % A0=k0(co);
% % % A0=exp(mu);
% % K0=stiff_m(T,A0,Nxo,Nyo);
% % aa=k(co,w');
% subplot(1,2,2)
% plot_BFE(k(co,w_true'),Mesh)
% title('Reference permeability field')
% colorbar
% set(gca,'FontSize',22)  %是设置刻度字体大小

load data_ell_pod_affine_inv_k_mu2 Uobs p w_true Q W_train tim_pod Usim
Usi_pod=Usim;
P=length(w_true);
for i=1:length(Q)
w_mean=mean(Q{i},2);
er_wpod(i)=norm(w_true-w_mean)/norm(w_true);
er_kpod(i)=norm(k(co,w_true')-k(co,w_mean'))/norm(k(co,w_true'));
end
Q_pod=Q;
% toc

load data_ell_greedy_affine_inv_k_mu_0319 Q tim_greedy Usim
Usi_gre=Usim;
for i=1:length(Q)
w_mean=mean(Q{i},2);
er_wgre(i)=norm(w_true-w_mean)/norm(w_true);
er_kgre(i)=norm(k(co,w_true')-k(co,w_mean'))/norm(k(co,w_true'));
end
Q_gre=Q;

load data_ell_vs_affine_k_mu_inv Q tim_vs Usim
Usi_vs=Usim;
for i=1:length(Q)
w_mean=mean(Q{i},2);
er_wvs(i)=norm(w_true-w_mean)/norm(w_true);
er_kvs(i)=norm(k(co,w_true')-k(co,w_mean'))/norm(k(co,w_true'));
end
Q_vs=Q;

load data_eff_affine_k_mu_inv Q tim_eff Usim Uf
Usi_eff=Usim;
for i=1:length(Q)
w_mean=mean(Q{i},2);
er_weff(i)=norm(w_true-w_mean)/norm(w_true);
er_keff(i)=norm(k(co,w_true')-k(co,w_mean'))/norm(k(co,w_true'));
end
Q_eff=Q;

subplot(1,2,2)
semilogy(1:length(Q_pod),er_kpod,'k-','Linewidth',2)
hold on
semilogy(1:length(Q_gre),er_kgre,'b-','Linewidth',2)
hold on
semilogy(1:length(Q_vs),er_kvs,'g-','Linewidth',2)
hold on
semilogy(1:length(Q_eff),er_keff,'r-','Linewidth',2)
legend('POD','Greedy','VS','EnEIM','fontsize',16)
grid on
title('Reletive error of a(x,\omega) ')
xlabel('Number of the inversion iterations')
set(gca,'FontSize',22)  %是设置刻度字体大小

% % 计算方差
% % % % P=2*N_pc+3;
% s_ga=zeros(size(k(co,w_true'),1),1);
% s_pod=zeros(size(k(co,w_true'),1),1);
% s_vs=zeros(size(k(co,w_true'),1),1);
% s_rom=zeros(size(k(co,w_true'),1),1);
% % s_fem=zeros(size(k(co,w_true'),1),1);
Q_pod1=Q_pod{end};
Q_gre1=Q_gre{end};
Q_vs1=Q_vs{end};
Q_eff1=Q_eff{end};
for i=1:size(Q{end},2)
    k_pod(:,i)=k(co,Q_pod1(:,i)');
    k_gre(:,i)=k(co,Q_gre1(:,i)');
     k_vs(:,i)=k(co,Q_vs1(:,i)');
      k_eff(:,i)=k(co,Q_eff1(:,i)');
end
for j=1:size(k(co,w_true'),1)
      W= k_gre(j,:);
      s_gre(j)=std(W');
      W= k_pod(j,:);
      s_pod(j)=std(W');
      W=k_vs(j,:);
      s_vs(j)=std(W');
      W=k_eff(j,:);
      s_eff(j)=std(W');
 end
  Varc_gre=s_gre.^2;
Varc_pod=s_pod.^2;
Varc_vs=s_vs.^2;
Varc_eff=s_eff.^2;


figure(2)
tiledlayout(3,2)
nexttile([1 2])
plot_BFE(k(co,w_true'),Mesh)
title('Reference for a(x,\omega*)')
colorbar
set(gca,'FontSize',22)  %是设置刻度字体大小
nexttile(3)
w_mean=mean(Q_pod{end},2);
plot_BFE(k(co,w_mean'),Mesh)
colorbar
title('POD estimate')
set(gca,'FontSize',22)  %是设置刻度字体大小
nexttile(4)
w_mean=mean(Q_gre{end},2);
plot_BFE(k(co,w_mean'),Mesh)
colorbar
title('Greedy  estimate')
set(gca,'FontSize',22)  %是设置刻度字体大小
nexttile(5)
w_mean=mean(Q_vs{end},2);
plot_BFE(k(co,w_mean'),Mesh)
colorbar
title('VS  estimate')
set(gca,'FontSize',22)  %是设置刻度字体大小
nexttile(6)
w_mean=mean(Q_eff{end},2);
plot_BFE(k(co,w_mean'),Mesh)
colorbar
title('EnEIM  estimate')
set(gca,'FontSize',22)  %是设置刻度字体大小

figure(3)
subplot(2,2,1)
plot_BFE(Varc_pod',Mesh)
colorbar
title('POD')
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,2,2)
plot_BFE(Varc_gre',Mesh)
colorbar
title('Greedy')
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,2,3)
plot_BFE(Varc_vs',Mesh)
colorbar
title('VS')
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,2,4)
plot_BFE(Varc_eff',Mesh)
colorbar
title('EnEIM')
set(gca,'FontSize',22)  %是设置刻度字体大小


figure(4)
ind=[1 3 6 9 12 15];
for i=1:length(ind) 
QQ(1,:)=Q_pod1(ind(i),:);
QQ(2,:)=Q_gre1(ind(i),:);
QQ(3,:)=Q_vs1(ind(i),:);
QQ(4,:)=Q_eff1(ind(i),:);
subplot(2,3,i)
boxplot(QQ','OutlierSize',0.8,'Widths',0.5,'Labels',{'POD','Greedy','VS','EnEIM'},'sym','')
hold on
 plot(1:4,w_true(ind(i)):w_true(ind(i)),'r*','markersize',12)
%  legend('POD','Greedy','VS','EnEIM','fontsize',16)
 if i==1
title('\omega_1')
 elseif i==2
title('\omega_3')
 elseif i==3
title('\omega_6')
 elseif i==4
title('\omega_9')
 elseif i==5
title('\omega_{12}')
 else
title('\omega_{15}')
 end
 set(gca,'FontSize',22)  %是设置刻度字体大小
 end

% figure(5)
% subplot(1,2,1)
% semilogy(1:length(Q_pod),er_wpod,'k-','Linewidth',2)
% hold on
% semilogy(1:length(Q_gre),er_wgre,'b-','Linewidth',2)
% hold on
% semilogy(1:length(Q_vs),er_wvs,'g-','Linewidth',2)
% hold on
% semilogy(1:length(Q_eff),er_weff,'r-','Linewidth',2)
% legend('POD','Greedy','VS','EnEIM','fontsize',16)
% ylim([5e-3 5e0])
% grid on
% title('Reletive error of \omega ')
% xlabel('Number of the inversion iterations')
% set(gca,'FontSize',22)  %是设置刻度字体大小
% 
% subplot(1,2,2)
% semilogy(1:length(Q_pod),er_kpod,'k-','Linewidth',2)
% hold on
% semilogy(1:length(Q_gre),er_kgre,'b-','Linewidth',2)
% hold on
% semilogy(1:length(Q_vs),er_kvs,'g-','Linewidth',2)
% hold on
% semilogy(1:length(Q_eff),er_keff,'r-','Linewidth',2)
% legend('POD','Greedy','VS','EnEIM','fontsize',16)
% grid on
% title('Reletive error of k(x,\omega) ')
% xlabel('Number of the inversion iterations')
% set(gca,'FontSize',22)  %是设置刻度字体大小


for i=1:length(Usi_pod)
En_pod(i)=norm(mean(Uobs-Usi_pod{i},2))^2;
end
for i=1:length(Usi_gre)
En_gre(i)=norm(mean(Uobs-Usi_gre{i},2))^2;
end
for i=1:length(Usi_vs)
En_vs(i)=norm(mean(Uobs-Usi_vs{i},2))^2;
end
for i=1:length(Usi_eff)
En_eff(i)=norm(mean(Uobs-Usi_eff{i},2))^2;
end

figure(6)
ss=length(Uf)*sig.^2;
semilogy(1:length(Usi_pod),En_pod,'k-','Linewidth',2)
hold on
semilogy(1:length(Usi_gre),En_gre,'b-','Linewidth',2)
hold on
semilogy(1:length(Usi_vs),En_vs,'g-','Linewidth',2)
hold on
semilogy(1:length(Usi_eff),En_eff,'r-','Linewidth',2)
hold on
semilogy(1:20,ss.*ones(20,1),'r--','Linewidth',2)
legend('POD','Greedy','VS','EnEIM','n_d\sigma^2','fontsize',16)
grid on
title('Discrepancy principle')
xlabel('Number of the inversion iterations')
set(gca,'FontSize',22)  %是设置刻度字体大小

% colorList=[0.9300    0.9500    0.9700
%     0.7900    0.8400    0.9100
%     0.6500    0.7300    0.8500
%     0.5100    0.6200    0.7900
%     0.3700    0.5100    0.7300
%     0.2700    0.4100    0.6300
%     0.2100    0.3200    0.4900
%     0.1500    0.2200    0.3500
%     0.0900    0.1300    0.2100
%     0.0300    0.0400    0.0700];
% 
%  figure(7)
% for i=1:6
%   
%  QQ=Q_pod1;
%  QQ=outlier(QQ);
% if i==1
%   subplot(6,6,1)
% plot(w_true(i),0,'r*','Markersize',8)
% hold on
%  ah0=Q_pod1(ind(i)); 
% [fh,xhi] = ksdensity(ah0,'width',0.0055);
% [nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
% hold on
% plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
%  elseif i==2
%      subplot(6,6,8)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%    ah0=Q_pod1(ind(i)); 
% [fh,xhi] = ksdensity(ah0,'width',0.088);
% [nh,xh]=hist(ah0,10);
% % nh=nh/size(ah0,2);
% % bar(xh,nh,1)
% hold on
% plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
% elseif i==3
%      subplot(6,6,15)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%   ah0=Q_pod1(ind(i)); 
% [fh,xhi] = ksdensity(ah0,'width',0.004);
% [nh,xh]=hist(ah0,10);
% % nh=nh/size(ah0,2);
% % bar(xh,nh,1)
% hold on
% plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
% elseif i==4
%  subplot(6,6,22)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%  ah0=Q_pod1(ind(i)); 
% [fh,xhi] = ksdensity(ah0,'width',0.004);
% [nh,xh]=hist(ah0,10);
% % nh=nh/size(ah0,2);
% % bar(xh,nh,1)
% hold on
% plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
% elseif i==5
%  subplot(6,6,29)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%  ah0=Q_pod1(ind(i)); 
% [fh,xhi] = ksdensity(ah0,'width',0.0058);
% [nh,xh]=hist(ah0,10);
% % nh=nh/size(ah0,2);
% % bar(xh,nh,1)
% hold on
% plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
% elseif i==6
%  subplot(6,6,36)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%  ah0=Q_pod1(ind(i)); 
% [fh,xhi] = ksdensity(ah0,'width',0.0079);
% [nh,xh]=hist(ah0,10);
% % nh=nh/size(ah0,2);
% % bar(xh,nh,1)
% hold on
% plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
%  end
%  if i>1
% for j=1:i-1
%  
% subplot(6,6,j+6*(i-1))
% % h2=(max(QQ(i,:)',1)-min(QQ(i,:)',1))/5;
% % h1=(max(QQ(j,:)',1)-min(QQ(j,:)',1))/5;
% % yli=min(QQ(i,:)',1):h2:max(QQ(i,:)',1);
% % xli=min(QQ(j,:)',1):h1:max(QQ(j,:)',1);
% % xli=-2:0.01:2;
% % yli=-2:0.01:2;
% if j==1&& i==2
% xli=-0.7:0.006:-0.6;
% yli=4.12:0.011:4.34;
% elseif j==1&&i==3
% xli=0.11:0.07:0.25;
% yli=0.07:0.003:0.13;  
% elseif j==2&&i==3
% xli=4.12:0.011:4.34;
% yli=0.07:0.003:0.13;  
% elseif j==1&&i==4
% xli=0.12:0.006:0.24;
% yli=0.79:0.0005:0.89;  
% elseif j==2&&i==4
% xli=4.12:0.011:4.34;
% yli=0.79:0.0005:0.89;  
% elseif j==3&&i==4
% xli=0.065:0.003:0.125;
% yli=0.79:0.0005:0.89;  
% elseif j==1&&i==5
% xli=0.12:0.006:0.24;
% yli=0.775:0.004:0.855;  
% elseif j==2&&i==5
% xli=4.1:0.012:4.34;
% yli=0.775:0.004:0.855;  
% elseif j==3&&i==5
% xli=0.065:0.003:0.125;
% yli=0.78:0.004:0.86;
% elseif j==4&&i==5
% xli=0.79:0.0005:0.89;
% yli=0.78:0.004:0.86;  
% elseif j==1&&i==6
% xli=0.12:0.006:0.24;
% yli=0.12:0.009:0.3;  
% elseif j==2&&i==6
% xli=4.1:0.012:4.34;
% yli=0.12:0.009:0.3;  
% elseif j==3&&i==6
% xli=0.07:0.003:0.13;
% yli=0.12:0.009:0.3;
% elseif j==4&&i==6
% xli=0.79:0.0005:0.89;
% yli=0.12:0.009:0.3; 
% elseif j==5&&i==6
% xli=0.78:0.004:0.86;
% yli=0.12:0.009:0.3; 
% end
% [~,~,XMesh,YMesh,ZMesh,colorList]=density2C(QQ(j,:)',QQ(i,:)',xli,yli,colorList);
% set(gcf,'Color',[1 1 1]);
% contourf(XMesh,YMesh,ZMesh)
% end
%  end
% end
% figure(8)
% for i=1:6
%   
%  QQ=Q_pod1;
% %  QQ=outlier(QQ);
% if i==1
%   subplot(3,2,i)
% plot(w_true(i),0,'r*','Markersize',8)
% hold on
%  ah0=Q_pod1(ind(i),:); 
% [fh,xhi] = ksdensity(ah0,'width',0.000055);
% [nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
% hold on
% plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
%  elseif i==2
%      subplot(3,2,i)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%    ah0=Q_pod1(ind(i),:); 
% [fh,xhi] = ksdensity(ah0,'width',0.00000088);
% [nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
% % hold on
% % plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
% elseif i==3
%      subplot(3,2,i)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%   ah0=Q_pod1(ind(i),:); 
% [fh,xhi] = ksdensity(ah0,'width',0.0000004);
% [nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
% % hold on
% % plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
% elseif i==4
%  subplot(3,2,i)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%  ah0=Q_pod1(ind(i),:); 
% [fh,xhi] = ksdensity(ah0,'width',0.0000004);
% [nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
% % hold on
% % plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
% elseif i==5
%  subplot(3,2,i)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%  ah0=Q_pod1(ind(i),:); 
% [fh,xhi] = ksdensity(ah0,'width',0.0000058);
% [nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
% % hold on
% % plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
% elseif i==6
%  subplot(3,2,i)
% plot(w_true(ind(i)),0,'r*','Markersize',8)
% hold on
%  ah0=Q_pod1(ind(i),:); 
% [fh,xhi] = ksdensity(ah0,'width',0.00000079);
% [nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
% % hold on
% % plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
%  end
%  end


