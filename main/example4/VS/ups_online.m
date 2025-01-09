function [uv]=ups_online(w1,VS,nF,a1,a1g,a1g0,a10,Aq,Aq0,fs,ks)
% fd1=fs{1}(w1)*nF(:,1);
% for i=2:size(fs,2)
%    fd1=fd1+fs{i}(w1)*nF(:,i);
% end
% K1=zeros(size(w1,1),size(ks,2));K1(:,1)=ks{1}(w1);
% for i=2:size(ks,2)
%     K1(:,i)=ks{i}(w1);
% end
% kss=K1*a1';fd1=fd1-sum(VS.*(A1*K1')',2);uv=fd1./kss;
fd1=zeros(size(w1,1),1);K1=zeros(size(w1,1),size(ks,2));
for j=1:size(w1,1)
    fd=fs{1}(w1(j,:))*nF(1);
for i=2:size(fs,2)
   fd=fd+fs{i}(w1(j,:))*nF(i);
end
fd1(j)=fd;
K1(j,1)=ks{1}(w1(j,:));
for i=2:size(ks,2)
    K1(j,i)=ks{i}(w1(j,:));
end
end
kss=K1*a1'+a10;
fd1=fd1-sum(VS.*(Aq*K1'+Aq0)',2)-K1*a1g'-a1g0;
uv=fd1./kss;



