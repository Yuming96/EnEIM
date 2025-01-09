function [a,u ] = stepk1(t,ki,w1,inda)
Ktrain=zeros(size(t,1),size(w1,1));
for i=1:size(w1,1)
    s0=w1(i,:);
    Ktrain(:,i)=ki(s0);
end
% u=@(s) sum(hx.^2*ki(s));
% a=mean(Ktrain,2)/mean(sum(hx.^2*Ktrain,1));
% u=@(s) sum(hx.^2*ki(s));
a=mean(Ktrain,2);
index=find(abs(a)==max(abs(a)));
inda(index(1))=1;
u=@(s) inda*ki(s)./a(index(1));
end

