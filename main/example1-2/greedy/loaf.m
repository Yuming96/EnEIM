function [If,Jf,Value_f]=loaf(T,Nx,Ny)
el=T.Elements;
If=[el(:,1);el(:,2);el(:,3);el(:,4)];
Jf=ones(4*Nx*Ny,1);

A1=zeros(Nx*Ny,4);A2=zeros(Nx*Ny,4);
A3=zeros(Nx*Ny,4);A4=zeros(Nx*Ny,4);
for i=1:Nx*Ny
    A1(i,:)=T.basis(:,1,i);
    A2(i,:)=T.basis(:,2,i);
    A3(i,:)=T.basis(:,3,i);
    A4(i,:)=T.basis(:,4,i);
end

Br=[A1;A2;A3;A4];

Gaussx=repmat(T.Gauss_x,4,1);
Gaussy=repmat(T.Gauss_y,4,1);
% a=f([Gaussx(:,1),Gaussy(:,1)]);
Value_f=(Br(:,1)+Br(:,2).*Gaussx(:,1)+Br(:,3).*Gaussy(:,1)+Br(:,4).*Gaussx(:,1).*Gaussy(:,1))+...
                 (Br(:,1)+Br(:,2).*Gaussx(:,2)+Br(:,3).*Gaussy(:,2)+Br(:,4).*Gaussx(:,2).*Gaussy(:,2))+...
           (Br(:,1)+Br(:,2).*Gaussx(:,3)+Br(:,3).*Gaussy(:,3)+Br(:,4).*Gaussx(:,3).*Gaussy(:,3))+...
            (Br(:,1)+Br(:,2).*Gaussx(:,4)+Br(:,3).*Gaussy(:,4)+Br(:,4).*Gaussx(:,4).*Gaussy(:,4));
Value_f=T.area.*Value_f/4;
        
% Val_b=T.area.*Value_B/4;
% F=sparse(I,J,Val_b);
% F(T.CNodePtrs)=[];

