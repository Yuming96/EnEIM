function [VS]=online_RBgreedy(w1,nF,nFg,A1,nF0,A10,fs,ks)
VS=zeros(size(w1,1),size(A1,1));
for j=1:size(w1,1)
    A=A10+ks{1}(w1(j,:))*A1(:,:,1);
    Fg=nF0+ks{1}(w1(j,:))*nFg(:,1);
    for i=2:size(ks,2)
        A=A+ks{i}(w1(j,:))*A1(:,:,i);
        Fg=Fg+ks{i}(w1(j,:))*nFg(:,i);
    end
    F=fs{1}(w1(j,:))*nF(:,1)-Fg;
    for i=2:size(fs,2)
       F=F+fs{i}(w1(j,:))*nF(:,i);
    end
    VS(j,:)=(A\F)';
end