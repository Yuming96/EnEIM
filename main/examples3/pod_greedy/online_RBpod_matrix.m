function [A,F]=online_RBpod_matrix(w1,nF,nFg,A1,nFg0,A0,fs,ks)
% VS=zeros(size(w1,1),size(A1,1));
% for j=1:size(w1,1)
    A=A0+ks{1}(w1)*A1(:,:,1);
    Fg=nFg0+ks{1}(w1)*nFg(:,1);
    for i=2:size(ks,2)
        A=A+ks{i}(w1)*A1(:,:,i);
        Fg=Fg+ks{i}(w1)*nFg(:,i);
    end
    F=fs{1}(w1)*nF(:,1)-Fg;
    for i=2:size(fs,2)
       F=F+fs{i}(w1)*nF(:,i);
    end
%     VS(j,:)=(A\F)';
% end