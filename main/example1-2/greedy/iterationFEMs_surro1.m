function [U,Ite]=iterationFEMs_surro1(ks,fs,T,w1,Ax,Fx,Fg,G0)
Iter=20;
Ite=[];
Fs=zeros(size(fs,2),size(w1',1));Fs(1,:)=fs{1}(w1');
for i=2:size(fs,2)
    Fs(i,:)=fs{i}(w1');
end
Ks=zeros(size(ks,2),size(w1',1));Ks(1,:)=ks{1}(w1');
for i=2:size(ks,2)
    Ks(i,:)=ks{i}(w1');
end
% for j=1:Ntr
%     j
%      A=ks{1}(w1(:,j)')*Ax{2};
% %      A=ks{1}(w1(:,j))
%     Ag=ks{1}(w1(:,j)')*Axg{1};
%     %带随机的刚度矩阵带边界
%     for i=2:size(ks,2)
%         A=A+ks{i}(w1(:,j)')*Ax{i};
%         Ag=Ag+ks{i}(w1(:,j)')*Axg{i};
%     end
%     %带随机的质量矩阵
%     F=fs{1}(w1(:,j)')*Fx(:,1);
%     for i=2:size(fs,2)
%         F=F+fs{i}(w1(:,j)')*Fx(:,i);
%     end
%     U0=K0(T.FNodePtrs,T.FNodePtrs)\(F(T.FNodePtrs)-Kg(T.FNodePtrs));
    U0=Fx*Fs-repmat(Fg,1,size(w1,2));
%    err=max(abs(U0),1);
    U=U0;count=1;
    while count<Iter
             G0(T.FNodePtrs,:)=U;
             U=Urom_online(G0,U0,Ax,Ks);
            count=count+1
%             U1=U2;
    end
    Ite=[Ite;count];
% G0(T.FNodePtrs,:)=U;