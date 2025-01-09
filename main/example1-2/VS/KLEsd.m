function [KQ,Qq,P_kle0,W]=KLEsd(LM,Nx,kdata)
% A two-dimensional Log-Gaussian process. 

% Nx = 40; % K-L interplation mesh, use this mesh Nx X Nx for computation
% trunc=6; % truncation of K-L expansion. This means random dim=trunc
% load Thetfile Theta
%%
% Theta=-1+2*rand(trunc,1);
% Build a 2D grid on a coarse mesh 
nx=sqrt(size(kdata,1));
x=linspace(0,100,nx)'; % spatial domain is [0,1]^2
[X,Y]=meshgrid(x,x);
%mesh=[X(:) Y(:)];
P_kle0=mean(kdata,2);
Ptrain=kdata-repmat(P_kle0,1,size(kdata,2));
% [Vp,Dp] = eigs(Ptrain'*Ptrain/size(Ptrain,2),LM);
[Vp,Dp] = eig(Ptrain'*Ptrain/size(Ptrain,2));
dd= diag(Dp);
W=zeros(size(Ptrain,1),LM);
for kc=1:size(Dp,1)
    W(:,kc)=1/sqrt(Dp(kc,kc)*size(Ptrain,2))*(Ptrain*Vp(:,kc));
end

%% Interplate Eigenfunctions W to Nx*Nx mesh.
% x2 = linspace(0,1,Nx)'; % spatial domain is [0,1]^2
% [X2,Y2] = meshgrid(x2,x2);

% KQ =cell(1,LM+1);Qq =cell(1,LM+1);
% KQ{1} =@(x) interp2(X,Y,reshape(P_kle0,nx,nx),x(:,1),x(:,2));
% Qq{1} =@(s) 1;
% ki=@(x,s) KQ{1}(x)*Qq{1}(s);
% for i = 1:LM
%     KQ{i+1} =@(x) interp2(X,Y,reshape(W(:,i),nx,nx),x(:,1),x(:,2));
%     Qq{i+1} =@(s) s(:,i);
%     ki=@(x,s) ki(x,s)+KQ{i+1}(x)*Qq{i+1}(s);
% end

KQ =cell(1,LM);Qq =cell(1,LM);
KQ{1} =@(x) interp2(X,Y,reshape(W(:,1),nx,nx),x(:,1),x(:,2));
Qq{1} =@(s) s(:,1);
ki=@(x,s) KQ{1}(x)*Qq{1}(s);
for i = 2:LM
    KQ{i} =@(x) interp2(X,Y,reshape(W(:,i),nx,nx),x(:,1),x(:,2));
    Qq{i} =@(s) s(:,i);
    ki=@(x,s) ki(x,s)+KQ{i}(x)*Qq{i}(s);
end

% Theta=-1+2*rand(LM,1);%Nx=length(x);
% x2 = linspace(0,100,Nx)'; % spatial domain is [0,1]^2
% [X2,Y2] = meshgrid(x2,x2);
% kle=reshape(ki([X2(:),Y2(:)],Theta'),Nx,Nx);
% figure (2);
% pcolor(kle);shading flat;colorbar;%caxis([-4 11]);
end
