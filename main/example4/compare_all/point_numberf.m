function [p,Nodes,inter,co,bdy]=point_numberf(Nxo,Nyo,T)
xo=0:1/Nxo:1;yo=0:1/Nyo:1;
[Xo,Yo]=meshgrid(xo,yo);
Nodes=[Xo(:) Yo(:)];
inter=1:(Nyo+1)*(Nxo+1);
bdy=T.CNodePtrs;
% inter(bdy)=[];
[X,Y] = meshgrid(1/Nxo/2:1/Nxo:1,1/Nyo/2:1/Nyo:1);
co=[X(:) Y(:)];
%  p=reshape(1:1:(Nyo-1)*(Nxo-1),Nyo-1,Nxo-1);
 p=reshape(1:(Nyo+1)*(Nxo+1),Nyo+1,Nxo+1);
% %  p=reshape(inter,Nyo+1,Nxo+1);
% % %  m=Nxo/10;
m=Nxo/5;
 p=p(1:m:end,:);
p=p(:,m/2+1:m:end);

p=p(:);
%  p=[p(3:4:end-2,1)',p(1,1:2:end),p(Nyo+1,1:2:end),p(3:4:end-2,Nxo+1)']';   
% p= [p(1,m/4:m/4:end)';p(Nyo+1,m/4:m/4:end)'];