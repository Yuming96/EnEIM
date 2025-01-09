function [G]=residsotion1(active,w1,ks,fs,G,Ax,Fx,Un,K0)
ks{1}(w1)
A=K0+ks{1}(w1)*Ax{1};
for i=2:size(ks,2)
    A=A+ks{i}(w1)*Ax{i};
end
F=fs{1}(w1)*Fx(:,1);
for i=2:size(fs,2)
    F=F+fs{i}(w1)*Fx(:,i);
end
Fg=A*G;
G(active)=A(active,active)\(F(active)-Fg(active)-A(active,active)*Un(active));


