function T=fem_structnum(x,y,h)
T.area=(x(2)-x(1))*(y(2)-y(1));
[X,Y]=meshgrid(x,y);
T.nodes=[X(:),Y(:)];
[m,n]=size(X);
in=[1:m*n]';
%将原来的矩阵变成m行n列，index是矩阵
index=reshape(1:m*n,m,n);
inde=index;
T.fbnodes=cell(1,length(h));
for i=1:length(h) 
       if i==1
          if h(i)==0
             index(end,:)=[];
          else
             T.fbnodes{1,i}=inde(end,:);
          end
       elseif i==2
           if h(i)==0
              index(:,1)=[];
           else
             T.fbnodes{1,i}=inde(:,1);
           end
       elseif i==3
           if h(i)==0  
             index(:,end)=[];
           else
             T.fbnodes{1,i}=inde(:,end);
           end 
       else
           if h(i)==0  
            index(1,:)=[];
           else
            T.fbnodes{1,i}=inde(1,:);
           end  
       end
end
T.fnodes=index(:);
in(T.fnodes)=[];
T.cnodes=in(:);
T.cnflag=zeros(m*n,1); T.cnflag(T.cnodes)=1;
T.edges=zeros(2*m*n-m-n,2); k=0;
T.nelem=zeros(m*n-m-n+1,4); kn=0;
for j=1:n
    for i=1:m
        if i>1 && j==1
            k=k+1; T.edges(k,:)=[i-1,i];
        elseif i==1 && j>1
            k=k+1; T.edges(k,:)=[j*m-2*m+1, j*m-m+1];
        elseif i>1 && j>1
            k=k+2; T.edges(k-1,:)=[j*m-2*m+i, j*m-m+i]; T.edges(k,:)=[j*m-m+i-1, j*m-m+i];
            kn=kn+1; T.nelem(kn,:)=[j*m-2*m+i-1, j*m-2*m+i, j*m-m+i-1, j*m-m+i];
        end
    end
end
T.midp=zeros(m*n-m-n+1,2);
for i=1:m*n-m-n+1
    T.midp(i,:)=mean(T.nodes(T.nelem(i,:),:));
end
% toc