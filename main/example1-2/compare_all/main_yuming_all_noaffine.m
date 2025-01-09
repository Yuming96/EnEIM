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


co=T.centriod;
bdy=T.CNodePtrs;
com=T.Nodes(bdy,:);
G=g(com);

load data_elliptic_greedy_noaffine coefko KQVo Uscofo tim_greedy count V
Ufem=V;
Nw=count;
% for i=1:Nw
    Uvo=KQVo*Uscofo';
% end
eer=[mean(coefko,1)]
Uvo_ga=Uvo;coefko_ga=coefko;
eer_ga=mean(coefko_ga,1)
% y_ga = mean(coefko_ga,1);e_ga = std(coefko_ga,1,1);
y_ga = mean((coefko_ga),1);e_ga = std((coefko_ga),1,1);
y_ga1=mean(coefko_ga,2);e_ga1 = std(coefko_ga,0,2);

load data_elliptic_pod_noaffine coefko KQVo Uscofo tim_pod
% Nw=count;
% for i=1:size(KQVo,2)
    Uvo=KQVo*Uscofo';
% end
Uvo_pod=Uvo;
coefko_pod=coefko;
eer_pod=mean(coefko_pod,1)
% y_pod = mean(coefko_pod,1);e_pod = std(coefko_pod,1,1);
y_pod = mean((coefko_pod),1);e_pod = std((coefko_pod),1,1);
y_pod1=mean(coefko_pod,2);e_pod1 = std(coefko_pod,0,2);

load data_elliptic_eff_noaffine_20240222 coefko Uromo tim_effi
% load data_elliptic_eff_noaffine coefko Uromo tim_effi
G0=repmat(g(T.Nodes),1,Nte);
Urom=zeros(size(G0));
Urom(T.CNodePtrs,:)=G0(T.CNodePtrs,:);
Urom(T.FNodePtrs,:)=Uromo{end};
coefko_rom=coefko;
eer_rom=mean(coefko_rom,1)
% y_rom = mean(coefko_rom,1);e_rom = std(coefko_rom,1,1);
y_rom = mean((coefko_rom),1);e_rom = std((coefko_rom),1,1);
y_rom1 = mean(coefko_rom,2);e_rom1 = std(coefko_rom,0,2);

load data_elliptic_vs_noaffine coefko Uscof KQV tim_vs
    Uvo=KQV*Uscof';
% end
Uvo_vs=Uvo;
coefko_vs=coefko;
eer_vs=mean(coefko_vs,1)
% y_vs = mean(coefko_vs,1);e_vs = std(coefko_vs,1,1);
y_vs = mean((coefko_vs),1);e_vs = std((coefko_vs),1,1);
y_vs1 = mean(coefko_vs,2);e_vs1 = std(coefko_vs,0,2);

figure(1)
subplot(1,2,1)
errorbar(y_ga,e_ga,'xg','linewidth',1.3)
hold on
errorbar(y_pod,e_pod,'xb','linewidth',1.3)
hold on
errorbar(y_vs,e_vs,'xk','linewidth',1.3); 
ylabel('Mean error with error bars','fontsize',16);
hold on
errorbar(y_rom,e_rom,'xr','linewidth',1.3)
grid on
legend('Greedy','POD','VS','EnEIM','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(1,2,2)
semilogy(1:size(eer_ga,2),eer_ga,'g-','linewidth',2)
hold on
semilogy(1:size(eer_pod,2),eer_pod,'b-','linewidth',2)
hold on
semilogy(1:size(eer_vs,2),eer_vs,'k-','linewidth',2)
hold on
semilogy(1:size(eer_rom,2),eer_rom,'r-','linewidth',2)
ylabel('Relative error ','fontsize',16);
grid on
legend('Greedy','POD','VS','EnEIM','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小

% close
figure(2)
subplot(2,2,1)
semilogy(1:size(coefko_ga,1),coefko_ga(:,1),'b-','linewidth',1)
% ylabel('Relative error','fontsize',16);
hold on
semilogy(1:size(coefko_ga,1),coefko_ga(:,size(eer_ga,2)),'r-','linewidth',1)
legend('1th iteration','11th iteration','fontsize',16)
xlabel('The index of samples','fontsize',12);  
title('Greedy','fontsize',16)
ylim([1e-6 1e1])
set(gca,'FontSize',22)  %是设置刻度字体大小

subplot(2,2,2)
semilogy(1:size(coefko_pod,1),coefko_pod(:,1),'b-','linewidth',1)
hold on
semilogy(1:size(coefko_pod,1),coefko_pod(:,size(eer_pod,2)),'r-','linewidth',1)
legend('r=1','r=26','fontsize',16)
xlabel('The index of samples','fontsize',12);  
title('POD','fontsize',16)
ylim([1e-6 1e1])
set(gca,'FontSize',22)  %是设置刻度字体大小

subplot(2,2,3)
semilogy(1:size(coefko_vs,1),coefko_vs(:,1),'b-','linewidth',1)
hold on
semilogy(1:size(coefko_vs,1),coefko_vs(:,size(eer_vs,2)),'r-','linewidth',1)
legend('1th iteration','14th iteration','fontsize',16)
xlabel('The index of samples','fontsize',12);  
title('VS','fontsize',16)
ylim([1e-6 1e1])
set(gca,'FontSize',22)  %是设置刻度字体大小

subplot(2,2,4)
semilogy(1:size(coefko_rom,1),coefko_rom(:,1),'b-','linewidth',1)
hold on
semilogy(1:size(coefko_rom,1),coefko_rom(:,size(eer_rom,2)),'r-','linewidth',1)
xlabel('The index of samples','fontsize',12);  
% ylabel('Relative error','fontsize',16);
% grid on
legend('1th iteration','10th iteration','fontsize',16)
ylim([1e-6 1e1])
title('EnEIM','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小


% close
figure(3)
% subplot(3,2,1)
Uf=uFDM(Nxo,Nyo,G,mean(Ufem(T.FNodePtrs,:),2),T);
Uga=uFDM(Nxo,Nyo,G,mean(Uvo_ga(T.FNodePtrs,:),2),T);
Upod=uFDM(Nxo,Nyo,G,mean(Uvo_pod(T.FNodePtrs,:),2),T);
Uro=uFDM(Nxo,Nyo,G,mean(Urom(T.FNodePtrs,:),2),T);
Uvs=uFDM(Nxo,Nyo,G,mean(Uvo_vs(T.FNodePtrs,:),2),T);
tiledlayout(3,2)
nexttile([1 2])
plot_BFE(Uf,Mesh)
colorbar
% colormap(jet)
title('Reference','fontsize',16)
% xlabel('The mean of u(\omega)','Fontsize',16)
set(gca,'FontSize',18) %是设置刻度字体大小
nexttile(3)
% mesh(xo,yo,reshape(Uvo_ga(:,l),Nyo+1,Nxo+1))
plot_BFE(Uga,Mesh)
colorbar
% colormap(jet)
title('Greedy','fontsize',16)
set(gca,'FontSize',20) %是设置刻度字体大小
nexttile(4)
% mesh(xo,yo,reshape(Uvo_pod(:,l),Nyo+1,Nxo+1))
plot_BFE(Upod,Mesh)
colorbar
% colormap(jet)
title('POD','fontsize',16)
set(gca,'FontSize',18) %是设置刻度字体大小
nexttile(6)
% mesh(xo,yo,reshape(Urom(:,l),Nyo+1,Nxo+1))
plot_BFE(Uro,Mesh)
colorbar
% colormap(jet)
title('EnEIM','fontsize',16)
% xlabel('The mean of u(\omega)','Fontsize',16)
set(gca,'FontSize',18) %是设置刻度字体大小
nexttile(5)
% mesh(xo,yo,reshape(Urom(:,l),Nyo+1,Nxo+1))
plot_BFE(Uvs,Mesh)
colorbar
% colormap(jet)
title('VS','fontsize',16)
% xlabel('The mean of u(\omega)','Fontsize',16)
set(gca,'FontSize',18) %是设置刻度字体大小


% close
figure(4)
% subplot(1,2,1)
Uf=Ufem(152,:);
Uga=Uvo_ga(152,:);
Upod=Uvo_pod(152,:);
Uvs=Uvo_vs(152,:);
Uro=Urom(152,:);
 ah0=Uf;
[fh,xhi] = ksdensity(ah0,'width',0.21);
[nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
hold on
plot(xhi,fh*(xh(2)-xh(1)),'k-*','Linewidth',2)
 ah0=Uga;
[fh,xhi] = ksdensity(ah0,'width',0.23);
[nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
hold on
plot(xhi,fh*(xh(2)-xh(1)),'g-','Linewidth',2)
hold on
 ah0=Upod;
[fh,xhi] = ksdensity(ah0,'width',0.2);
[nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
hold on
plot(xhi,fh*(xh(2)-xh(1)),'k-','Linewidth',2)
% colorbar
hold on
 ah0=Uvs;
[fh,xhi] = ksdensity(ah0,'width',0.25);
[nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
hold on
plot(xhi,fh*(xh(2)-xh(1)),'b-','Linewidth',2)
hold on
 ah0=Uro;
[fh,xhi] = ksdensity(ah0,'width',0.21);
[nh,xh]=hist(ah0,10);
% nh=nh/size(ah0,2);
% bar(xh,nh,1)
hold on
plot(xhi,fh*(xh(2)-xh(1)),'r-','Linewidth',2)
legend('Reference','Greedy','POD','VS','EnEIM','fontsize',16)
xlabel('Solution u at x_0=(0.01,0.5)','fontsize',16);  
ylabel('The probability of u','fontsize',16);
set(gca,'FontSize',22) %是设置刻度字体大小

