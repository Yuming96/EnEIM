function [K_s]=stiff_m(T,kai,Nx,Ny)
el=T.Elements;

I=[kron(ones(4,1),el(:,1));kron(ones(4,1),el(:,2));kron(ones(4,1),el(:,3));kron(ones(4,1),el(:,4))];
J=kron(ones(4,1),[el(:,1);el(:,2);el(:,3);el(:,4)]);

A1=zeros(Nx*Ny,4);A2=zeros(Nx*Ny,4);
A3=zeros(Nx*Ny,4);A4=zeros(Nx*Ny,4);
for i=1:Nx*Ny
    A1(i,:)=T.basis(:,1,i);
    A2(i,:)=T.basis(:,2,i);
    A3(i,:)=T.basis(:,3,i);
    A4(i,:)=T.basis(:,4,i);
end

Ar=[kron(ones(4,1),A1);kron(ones(4,1),A2);kron(ones(4,1),A3);kron(ones(4,1),A4)];
Br=kron(ones(4,1),[A1;A2;A3;A4]);
Gaussx=repmat(T.Gauss_x,16,1);
Gaussy=repmat(T.Gauss_y,16,1);

Value_K=((Ar(:,2)+Ar(:,4).*Gaussy(:,1)).*(Br(:,2)+Br(:,4).*Gaussy(:,1))+(Ar(:,3)+Ar(:,4).*Gaussx(:,1)).*(Br(:,3)+Br(:,4).*Gaussx(:,1))+...
            (Ar(:,2)+Ar(:,4).*Gaussy(:,2)).*(Br(:,2)+Br(:,4).*Gaussy(:,2))+(Ar(:,3)+Ar(:,4).*Gaussx(:,2)).*(Br(:,3)+Br(:,4).*Gaussx(:,2))+...
            (Ar(:,2)+Ar(:,4).*Gaussy(:,3)).*(Br(:,2)+Br(:,4).*Gaussy(:,3))+(Ar(:,3)+Ar(:,4).*Gaussx(:,3)).*(Br(:,3)+Br(:,4).*Gaussx(:,3))+...
            (Ar(:,2)+Ar(:,4).*Gaussy(:,4)).*(Br(:,2)+Br(:,4).*Gaussy(:,4))+(Ar(:,3)+Ar(:,4).*Gaussx(:,4)).*(Br(:,3)+Br(:,4).*Gaussx(:,4)));
Val_k=T.area.*repmat(kai,16,1).*Value_K/4;
K_s=sparse(I,J,Val_k);

% Value_B=(Ar(:,1)+Ar(:,2).*Gaussx(:,1)+Ar(:,3).*Gaussy(:,1)+Ar(:,4).*Gaussx(:,1).*Gaussy(:,1)).*(Br(:,1)+Br(:,2).*Gaussx(:,1)+Br(:,3).*Gaussy(:,1)+Br(:,4).*Gaussx(:,1).*Gaussy(:,1))+...
%            (Ar(:,1)+Ar(:,2).*Gaussx(:,2)+Ar(:,3).*Gaussy(:,2)+Ar(:,4).*Gaussx(:,2).*Gaussy(:,2)).*(Br(:,1)+Br(:,2).*Gaussx(:,2)+Br(:,3).*Gaussy(:,2)+Br(:,4).*Gaussx(:,2).*Gaussy(:,2))+...
%            (Ar(:,1)+Ar(:,2).*Gaussx(:,3)+Ar(:,3).*Gaussy(:,3)+Ar(:,4).*Gaussx(:,3).*Gaussy(:,3)).*(Br(:,1)+Br(:,2).*Gaussx(:,3)+Br(:,3).*Gaussy(:,3)+Br(:,4).*Gaussx(:,3).*Gaussy(:,3))+...
%             (Ar(:,1)+Ar(:,2).*Gaussx(:,4)+Ar(:,3).*Gaussy(:,4)+Ar(:,4).*Gaussx(:,4).*Gaussy(:,4)).*(Br(:,1)+Br(:,2).*Gaussx(:,4)+Br(:,3).*Gaussy(:,4)+Br(:,4).*Gaussx(:,4).*Gaussy(:,4));
% % Val_b=T.area.*ones(16*Nx*Ny,1).*Value_B/4;
% Val_b=T.area.*repmat(a0,16,1).*Value_B/4;
% B_m=sparse(I,J,Val_b);
% % Val_p=T.area.*repmat(p,16,1).*Value_B/4;
% % B_p=sparse(I,J,Val_p);
% 
% % KK=K_s;
% K_s(T.CNodePtrs,:)=[];
% K{2}=K_s(:,T.CNodePtrs);
% K_s(:,T.CNodePtrs)=[];
% K{1}=K_s;
% % B_m(T.CNodePtrs,:)=[];
% % B{2}=B_m(:,T.CNodePtrs);
% % B_m(:,T.CNodePtrs)=[];
% % B{1}=B_m;
% % % B_p(T.CNodePtrs,:)=[];
% % % B2{2}=B_p(:,T.CNodePtrs);
% % % B_p(:,T.CNodePtrs)=[];
% % % B2{1}=B_p;


