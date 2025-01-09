function U=Urom_online(G0,U0,Ax,Ks)
U=U0;
for i=1:size(Ks,1)
  U=U-Ax{i}*G0.*repmat(Ks(i,:),size(U0,1),1);
end
