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
Mesh = TProd_Mesh(xo,yo);


co=T.centriod;
bdy=T.CNodePtrs;
com=T.Nodes(bdy,:);
G=g(com);

 load data_para_fem_affine Ufem tim_fem 
 Uvo_fem=zeros(size(Ufem{1},1),Nte,Nt+1);
 for ii=1:Nte
    Uvo_fem(:,ii,:)=Ufem{ii};
 end

load data_para_pod_affine coefko KQVo Uscofo tim_pod Nw
Uvo=zeros(size(T.Nodes,1),Nte,3);
n1=[1 100];nn=[Nt/20+1 Nt/2+1 Nt+1];
for ii=1:Nte
for j=1:length(nn)
    Uvo(T.FNodePtrs,ii,j)=KQVo*Uscofo{ii}(:,nn(j));
    Uvo(T.CNodePtrs,ii,j)=G;
end
end
Uvo_pod=Uvo;
coefko_pod=coefko;
eer_pod= reshape(mean(coefko,1),Nt+1,Nw);
y_pod = reshape(mean(log(coefko),1),Nt+1,Nw);
e_pod=zeros(Nt+1,Nw);
for j=2:Nt+1
    e_podstd=reshape(log(coefko(:,j,:)),Nte,Nw);
e_pod(j,:) = std(e_podstd,1,1);
end

load data_para_eff_affine coefko U_eff tim_effi
Uvo=zeros(size(T.Nodes,1),Nte,3);
n1=[1 100];nn=[Nt/20+1 Nt/2+1 Nt+1];
for i=1:length(nn)
for ii=1:Nte
    Uvo(:,ii,i)=U_eff{1,i}{1,end}(:,ii);
end
end
Uvo_rom=Uvo;
coefko_rom=coefko;

Nw=20;
y_rom = reshape(mean(log(coefko),1),Nt+1,Nw);
eer_rom= reshape(mean(coefko,1),Nt+1,Nw);
e_rom=zeros(Nt+1,Nw);
for j=2:Nt+1
    e_romstd=reshape(log(coefko(:,j,:)),Nte,Nw);
e_rom(j,:) = std(e_romstd,1,1);
end

% % 计算方差
% % % % P=2*N_pc+3;
s_pod1=zeros(size(Uvo_pod,1),1);
s_rom1=zeros(size(Uvo_pod,1),1);
s_fem1=zeros(size(Uvo_pod,1),1);
s_pod2=zeros(size(Uvo_pod,1),1);
s_rom2=zeros(size(Uvo_pod,1),1);
s_fem2=zeros(size(Uvo_pod,1),1);
s_pod3=zeros(size(Uvo_pod,1),1);
s_rom3=zeros(size(Uvo_pod,1),1);
s_fem3=zeros(size(Uvo_pod,1),1);
for j=1:size(Uvo_pod,1)
      W=Uvo_pod(j,:,1);
      s_pod1(j)=std(W');
      W=Uvo_rom(j,:,1);
      s_rom1(j)=std(W');
       W=Uvo_fem(j,:,Nt/20+1);
      s_fem1(j)=std(W');
       W=Uvo_pod(j,:,2);
      s_pod2(j)=std(W');
      W=Uvo_rom(j,:,2);
      s_rom2(j)=std(W');
        W=Uvo_fem(j,:,Nt/2+1);
       s_fem2(j)=std(W');
       W=Uvo_pod(j,:,3);
      s_pod3(j)=std(W');
      W=Uvo_rom(j,:,3);
      s_rom3(j)=std(W');
      W=Uvo_fem(j,:,Nt+1);
      s_fem3(j)=std(W');
end
Varc_pod1=s_pod1.^2;
Varc_rom1=s_rom1.^2;
Varc_fem1=s_fem1.^2;
Varc_pod2=s_pod2.^2;
Varc_rom2=s_rom2.^2;
Varc_fem2=s_fem2.^2;
Varc_pod3=s_pod3.^2;
Varc_rom3=s_rom3.^2;
Varc_fem3=s_fem3.^2;


figure(1)
subplot(2,2,1)
% errorbar(y_ga(1,:),e_ga(1,:),'xr','linewidth',2)
% hold on
errorbar(y_pod(2,:),e_pod(2,:),'xb','linewidth',2)
hold on
errorbar(y_rom(2,1:16),e_rom(2,1:16),'xr','linewidth',2)
% xlabel('Numbers of the iteration terms','fontsize',16);  
ylabel('Mean error with error bars','fontsize',16);
ylim([-40 -0.7])
grid on
set(gca,'yscale','log')  %是设置刻度字体大小
legend('POD-Greedy','EnEIM','fontsize',16)
title('$t=0.01T_e$','Interpreter',"latex",'fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,2,2)
% errorbar(y_ga(Nt/20,:),e_ga(Nt/20,:),'xr','linewidth',2)
hold on
errorbar(y_pod(Nt/20+1,:),e_pod(Nt/20+1,:),'xb','linewidth',2)
hold on
errorbar(y_rom(Nt/20+1,1:12),e_rom(Nt/20+1,1:12),'xr','linewidth',2)
% xlabel('Numbers of the iteration terms','fontsize',16);  
ylabel('Mean error with error bars','fontsize',16);
ylim([-40 -0.1])
grid on
set(gca,'yscale','log')  %是设置刻度字体大小
legend('POD-Greedy','EnEIM','fontsize',16)
title('$t=0.06T_e$','Interpreter',"latex",'fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,2,3)
% errorbar(y_ga(Nt/2,:),e_ga(Nt/2,:),'xr','linewidth',2)
% hold on
errorbar(y_pod(Nt/2+1,2:end),e_pod(Nt/2+1,2:end),'xb','linewidth',2)
hold on
errorbar(y_rom(Nt/2+1,1:6),e_rom(Nt/2+1,1:6),'xr','linewidth',2)
% xlabel('Numbers of the iteration terms','fontsize',16);  
ylabel('Mean error with error bars','fontsize',16);
ylim([-30 -1])
grid on
set(gca,'yscale','log')  %是设置刻度字体大小
legend('POD-Greedy','EnEIM','fontsize',16)
title('$t=0.5T_e$','Interpreter',"latex",'fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,2,4)
% errorbar(y_ga(Nt/2,:),e_ga(Nt/2,:),'xr','linewidth',2)
% hold on
errorbar(y_pod(Nt+1,1:end),e_pod(Nt+1,1:end),'xb','linewidth',2)
hold on
errorbar(y_rom(Nt+1,1:5),e_rom(Nt/2+1,1:5),'xr','linewidth',2)
% xlabel('Numbers of the iteration terms','fontsize',16);  
ylabel('Mean error with error bars','fontsize',16);
ylim([-30 -1])
grid on
set(gca,'yscale','log')  %是设置刻度字体大小
legend('POD-Greedy','EnEIM','fontsize',16)
title('$t=T_e$','Interpreter',"latex",'fontsize',22)
set(gca,'FontSize',22)  %是设置刻度字体大小


% close 
figure(2)
subplot(1,2,1)
errorbar(y_pod(:,end),e_pod(:,end),'xb','linewidth',2)
hold on
errorbar(y_rom(:,5),e_rom(:,5),'xr','linewidth',2)
xlabel('N_t','fontsize',16);  
ylabel('Mean error with error bars','fontsize',16);
grid on
legend('POD-Greedy','EnEIM','fontsize',16)
% xlim([0 0.1])
set(gca,'FontSize',22)  %是设置刻度字体大小

subplot(1,2,2)
% semilogy(t,e_ga(:,end),'r-','linewidth',2)
% hold on
semilogy(t,e_pod(:,end),'b-*','linewidth',2)
hold on
semilogy(t,e_rom(:,5),'r-*','linewidth',2)
xlabel('t','fontsize',16);  
ylabel('Relative error','fontsize',16);
grid on
legend('POD-Greedy','EnEIM','fontsize',16)
ylim([1e-8 1e-1])
set(gca,'FontSize',22)  %是设置刻度字体大小

figure(3)
subplot(1,2,1)
semilogy(1:size(coefko_pod,1),coefko_pod(:,end,1),'k-','linewidth',2)
% ylabel('1th iteration','fontsize',16);
hold on
semilogy(1:size(coefko_pod,1),coefko_pod(:,end,4),'b-','linewidth',2)
hold on
semilogy(1:size(coefko_pod,1),coefko_pod(:,end,size(eer_pod,2)),'r-','linewidth',2)
title('{\bf POD-Greedy} ($t=T_e$)','Interpreter',"latex",'fontsize',16)
legend('r=1','r=4','r=15','fontsize',16)
xlabel('The index of samples','fontsize',16);  
ylabel('Relative error ','fontsize',16);
ylim([1e-4 1e1])
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(1,2,2)
semilogy(1:size(coefko_rom,1),coefko_rom(:,end,1),'k-','linewidth',2)
hold on
semilogy(1:size(coefko_rom,1),coefko_rom(:,end,3),'b-','linewidth',2)
hold on
semilogy(1:size(coefko_rom,1),coefko_rom(:,end,5),'r-','linewidth',2)
title('{\bf EnEIM} ($t=T_e$)','Interpreter',"latex",'fontsize',16)
legend('1th iteration','3th iteration','5th iteration','fontsize',16)
xlabel('The index of samples','fontsize',16);  
ylabel('Relative error ','fontsize',16);
ylim([1e-12 1e-1])
set(gca,'FontSize',22)  %是设置刻度字体大小

figure(4)
Ufem1=uFDM(Nxo,Nyo,G,mean(Uvo_fem(T.FNodePtrs,:,Nt/20+1),2),T);
Ufem3=uFDM(Nxo,Nyo,G,mean(Uvo_fem(T.FNodePtrs,:,Nt+1),2),T);
Upod1=uFDM(Nxo,Nyo,G,mean(Uvo_pod(T.FNodePtrs,:,1),2),T);
Uro1=uFDM(Nxo,Nyo,G,mean(Uvo_rom(T.FNodePtrs,:,1),2),T);
% Upod2=uFDM(Nxo,Nyo,G,mean(Uvo_pod(T.FNodePtrs,:,2),2),T);
% Uro2=uFDM(Nxo,Nyo,G,mean(Uvo_rom(T.FNodePtrs,:,2),2),T);
Upod3=uFDM(Nxo,Nyo,G,mean(Uvo_pod(T.FNodePtrs,:,3),2),T);
Uro3=uFDM(Nxo,Nyo,G,mean(Uvo_rom(T.FNodePtrs,:,3),2),T);
subplot(2,3,1)
plot_BFE(Ufem1,Mesh)
colorbar
title('{\bf Reference} ($t=0.05T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The mean of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,2)
plot_BFE(Upod1,Mesh)
colorbar
title('{\bf POD-Greedy} ($t=0.05T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The mean of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,3)
plot_BFE(Uro1,Mesh)
colorbar
title('{\bf EnEIM} ($t=0.05T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The mean of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,4)
plot_BFE(Ufem3,Mesh)
colorbar
title('{\bf Reference} ($t=T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The mean of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,5)
plot_BFE(Upod3,Mesh)
colorbar
title('{\bf POD-Greedy} ($t=T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The mean of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,6)
plot_BFE(Uro3,Mesh)
colorbar
title('{\bf EnEIM} ($t=T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The mean of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小

% l=2;
% % Uf=uFDM(Nxo,Nyo,G,Varc_pod1(T.FNodePtrs),T);
% % Uga=uFDM(Nxo,Nyo,G,Uvo_ga(T.FNodePtrs,3,l),T);
Upod1=uFDM(Nxo,Nyo,Varc_pod1(T.CNodePtrs),Varc_pod1(T.FNodePtrs),T);
Uro1=uFDM(Nxo,Nyo,Varc_rom1(T.CNodePtrs),Varc_rom1(T.FNodePtrs),T);
Ufem1=uFDM(Nxo,Nyo,Varc_fem1(T.CNodePtrs),Varc_fem1(T.FNodePtrs),T);
Upod2=uFDM(Nxo,Nyo,Varc_pod2(T.CNodePtrs),Varc_pod2(T.FNodePtrs),T);
Uro2=uFDM(Nxo,Nyo,Varc_rom2(T.CNodePtrs),Varc_rom2(T.FNodePtrs),T);
Ufem2=uFDM(Nxo,Nyo,Varc_fem2(T.CNodePtrs),Varc_fem2(T.FNodePtrs),T);
Upod3=uFDM(Nxo,Nyo,Varc_pod3(T.CNodePtrs),Varc_pod3(T.FNodePtrs),T);
Uro3=uFDM(Nxo,Nyo,Varc_rom3(T.CNodePtrs),Varc_rom3(T.FNodePtrs),T);
Ufem3=uFDM(Nxo,Nyo,Varc_fem3(T.CNodePtrs),Varc_fem3(T.FNodePtrs),T);
figure(5)
subplot(2,3,1)
plot_BFE(Ufem2,Mesh)
colorbar
title('{\bf Reference} ($t=0.5T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The variance of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,2)
plot_BFE(Upod2,Mesh)
colorbar
title('{\bf POD-Greedy} ($t=0.5T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The variance of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,3)
plot_BFE(Uro2,Mesh)
colorbar
title('{\bf EnEIM} ($t=0.5T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The variance of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,4)
plot_BFE(Ufem3,Mesh)
colorbar
title('{\bf Reference} ($t=T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The variance of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,5)
plot_BFE(Upod3,Mesh)
colorbar
title('{\bf POD-Greedy} ($t=T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The variance of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
subplot(2,3,6)
plot_BFE(Uro3,Mesh)
colorbar
title('{\bf EnEIM} ($t=T_e$)','Interpreter',"latex",'fontsize',16)
% xlabel('The variance of u(\omega)','fontsize',16)
set(gca,'FontSize',22)  %是设置刻度字体大小
