function T=struc(x,y)
Nx=length(x)-1;
Ny=length(y)-1;
T.Degree=2;
T.Index=[1:(Nx+1)*(Ny+1)]';

[X,Y]=meshgrid(x,y);
Nod=[X(:) Y(:)];
T.Nodes=Nod;

T.FNodePtrs=[1:(Ny+1)*(Nx+1)]';
idx = reshape(T.FNodePtrs,Ny+1,Nx+1);
% T.CNodePtrs = [idx(2:end-1,1)',idx(1,:),idx(Ny+1,:),idx(2:end-1,Nx+1)']';   
% T.CNodePtrs= [ idx(1:end,1)',idx(Ny+1,2:end-1), idx(1:end,Nx+1)'];
T.CNodePtrs=[1,(Ny+1)*(Nx+1)];
% T.CNodePtrs=[];
T.FNodePtrs(T.CNodePtrs)=[];


T.Elements=zeros(Nx*Ny,4);
idx1= idx(1:end-1,1:end-1);
T.Elements(:,1)=idx1(:);
idx2 = idx1 + 1;
T.Elements(:,2)=idx2(:);
idx3 = idx1 + (Ny+1);
T.Elements(:,3) = idx3(:);
 idx4=idx1 + 1 + (Ny+1);
T.Elements(:,4)=idx4(:);

xm=x(1)+1/Nx/2:1/Nx:x(end)-1/Nx/2;
ym=y(1)+1/Ny/2:1/Ny:y(end)-1/Ny/2;
[Xm,Ym]=meshgrid(xm,ym);
T.centriod=[Xm(:) Ym(:)];
T.area=1/(Nx*Ny);

T.basis=zeros(4,4,size(T.Elements,1));E=eye(4,4);
for i=1:size(T.Elements,1)
    T.basis(:,:,i)=[ones(4,1) T.Nodes(T.Elements(i,:),1) T.Nodes(T.Elements(i,:),2) T.Nodes(T.Elements(i,:),1).*T.Nodes(T.Elements(i,:),2)]\E;
    %basis functions coefficient a+bx+cy+dxy;
end

T.Gauss_x=zeros(size(T.Elements,1),4);T.Gauss_y=zeros(size(T.Elements,1),4);
% T.Gauss_frax=zeros(size(T.Elements,1),2);T.Gauss_fray=zeros(size(T.Elements,1),2);
for i=1:size(T.Elements,1)
       mid=T.centriod(i,:)-T.Nodes(T.Elements(i,1),:);
       xh=[T.centriod(i,1)-mid(1)/sqrt(3) T.centriod(i,1)+mid(1)/sqrt(3)];
%        xxh=[T.Nodes(T.Elements(i,1),1) T.Nodes(T.Elements(i,1),1)];
       yh=[T.centriod(i,2)-mid(2)/sqrt(3) T.centriod(i,2)+mid(2)/sqrt(3)];
        hx=[xh;xh];hy=[yh' yh'];
       T.Gauss_x(i,:)=hx(:)';
       T.Gauss_y(i,:)=hy(:)';
%        T.Gauss_frax(i,:)=xxh(:)';
%        T.Gauss_fray(i,:)=yh(:)';
end

