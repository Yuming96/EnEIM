function [VS]=online_RBpod(w1,nF,nFg,A1,nFg0,A0,fs,ks,M,ht,Nt,U0)
VS=cell(1,size(w1,1));
    VS1=zeros(size(A1,1),Nt+1);
    VS1(:,1)=U0;
for j=1:size(w1,1)

    A=A0+ks{1}(w1(j,:))*A1(:,:,1);
        Fg=nFg0+ks{1}(w1(j,:))*nFg(:,1);
        for i=2:size(ks,2)
            A=A+ks{i}(w1(j,:))*A1(:,:,i);
            Fg=Fg+ks{i}(w1(j,:))*nFg(:,i);
        end
%         u=U0;
%         u_old=U0;
    for jj=2:Nt+1
        F=fs{1}(w1(j,:))*nF(:,jj,1)-Fg;
        for i=2:size(fs,2)
           F=F+fs{i}(w1(j,:))*nF(:,jj,i);
        end
%         err=norm(u,2);
%         while err>1e-6
%             u1=u;
%              F1=(M+ht.*A)*u-ht*F-M*u_old;
%              Df=M+ht.*A;
%              u=u-Df\F1;
%              err=norm(u-u1,2)/norm(u,2)
%         end
%         VS1(:,jj)=u;
%         u_old=u;
        VS1(:,jj)=(M+ht.*A)\(ht*F+M*VS1(:,jj-1));
    end
    VS{j}=VS1;
end
% VS=zeros(size(w1,1),size(A1,1));
% for j=1:size(w1,1)
%     for j=2:Nt+1
%         A=A0+ks{1}(w1(j,:))*A1(:,:,1);
%         Fg=nFg0+ks{1}(w1(j,:))*nFg(:,1);
%         for i=2:size(ks,2)
%             A=A+ks{i}(w1(j,:))*A1(:,:,i);
%             Fg=Fg+ks{i}(w1(j,:))*nFg(:,i);
%         end
%         F=fs{1}(w(j,:))*nF(:,1)-Fg;
%         for i=2:size(fs,2)
%            F=F+fs{i}(w(j,:))*nF(:,i);
%         end
%         VS(j,:)=(A\F)';
%     end
% end