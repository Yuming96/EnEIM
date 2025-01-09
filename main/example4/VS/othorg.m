function [av]=othorg(av0,KQv,Atem)
B=zeros(size(KQv,1),size(KQv,2));
for j=1:size(KQv,2) 
B(:,j)=(av0'*Atem*KQv(:,j))/(KQv(:,j)'*Atem*KQv(:,j))*KQv(:,j);
end
av=av0-sum(B,2);
end