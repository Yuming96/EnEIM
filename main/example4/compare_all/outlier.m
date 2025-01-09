function q=outlier(theta)
[P,M_EKF]=size(theta);
t_low=quantile(theta',0.25)-1.5*iqr(theta');
t_up=quantile(theta',0.75)+1.5*iqr(theta');
I=[];
for j=1:P
    q=theta(j,:);
    Id=zeros(M_EKF,1);
    z1=t_low(j);
    z2=t_up(j);
    parfor i=1:M_EKF
       if q(i)<z1||q(i)>z2
           Id(i)=i;
       end
    end
    Id(Id==0)=[];
    I=[I ;Id];
end
I=unique(I);
q=theta;
q(:,I)=[];