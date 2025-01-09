function [av]=othorg(av0,KQv,z,Atem)
B=zeros(size(KQv,1),size(KQv,2));
for j=1:size(KQv,2) 
% B(:,j)=(av0'*Atem*KQv(:,j))/(KQv(:,j)'*Atem*KQv(:,j))*KQv(:,j);
% B(:,j)=(av0'*Atem*KQv(:,j))/(KQv(:,j)'*Atem*KQv(:,j))*(:,j);
B(:,j)=(av0'*KQv(:,j))*z(:,j);
end
av=av0-sum(B,2);
end