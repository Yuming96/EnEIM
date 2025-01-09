function [uv]=ups_online(w1,VS,nF,A1,A2,a1,fs,ks)
% fd1=fs{1}(w1)*nF(:,1);
% for i=2:size(fs,2)
%    fd1=fd1+fs{i}(w1)*nF(:,i);
% end
% K1=zeros(size(w1,1),size(ks,2));K1(:,1)=ks{1}(w1);
% for i=2:size(ks,2)
%     K1(:,i)=ks{i}(w1);
% end
% kss=K1*a1';fd1=fd1-sum(VS.*(A1*K1')',2);uv=fd1./kss;

fd1=fs{1}(w1)*nF(1);
for i=2:size(fs,2)
   fd1=fd1+fs{i}(w1)*nF(i);
end
K1=zeros(size(w1,1),size(ks,2));K1(:,1)=ks{1}(w1);
for i=2:size(ks,2)
    K1(:,i)=ks{i}(w1);
end
kss=K1*reshape(a1,1,size(ks,2))';
fd1=fd1-sum(VS.*(A1*K1')',2)-K1*A2';uv=fd1./kss;



