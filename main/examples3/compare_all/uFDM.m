function [uu]=uFDM(Nx,Ny,G,U,T)
% hx=1/Nx;hy=1/Ny;
U1=zeros((Ny+1)*(Nx+1),1);
U1(T.FNodePtrs,:)=U;
U1(T.CNodePtrs,:)=G;
% u1=reshape(U1,Ny+1,Nx+1);
uu=[U1(T.Elements(:,1)) U1(T.Elements(:,2)) U1(T.Elements(:,3)) U1(T.Elements(:,4))];
% % u1=1/2.*[uu(:,1)+uu(:,2) uu(:,1)+uu(:,3) uu(:,2)+uu(:,4) uu(:,3)+uu(:,4)];
% % ux=(u1(:,4)-u1(:,1))/hx;
% % uy=(u1(:,3)-u1(:,2))/hy;
% % uu=1/4.*(u1(:,1)+u1(:,4)+u1(:,2)+u1(:,3));
% ux=((u1(1:end-1,2:end)-u1(1:end-1,1:end-1))/hx+(u1(2:end,2:end)-u1(2:end,1:end-1))/hx)/2;
% uy=((u1(2:end,1:end-1)-u1(1:end-1,1:end-1))/hy+(u1(2:end,2:end)-u1(1:end-1,2:end))/hy)/2;
uu=sum(uu,2)/4;