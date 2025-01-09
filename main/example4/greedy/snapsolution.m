%ÿ�������½�g1
function [U]=snapsolution(w1,ks,fs,T,g,Ax,Fx,K0)
U=zeros(size(T.Nodes,1),size(w1,2));
for j=1:size(w1,2)
    A=K0+ks{1}(w1(:,j)')*Ax{1};
    %������ĸնȾ�����߽�
    for i=2:size(ks,2)
        A=A+ks{i}(w1(:,j)')*Ax{i};
    end
    %���������������
    F=fs{1}(w1(:,j)')*Fx(:,1);
    for i=2:size(fs,2)
        F=F+fs{i}(w1(:,j)')*Fx(:,i);
    end
    G=zeros(size(T.Nodes,1),1);
    G(T.CNodePtrs)=g(T.Nodes(T.CNodePtrs,:));
    Fg=A*G;
    G(T.FNodePtrs)=A(T.FNodePtrs,T.FNodePtrs)\(F(T.FNodePtrs)-Fg(T.FNodePtrs));
    U(:,j)=G;
end


