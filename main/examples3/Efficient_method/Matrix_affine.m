function [Fs,Ks0]=Matrix_affine(k,fss,w1,co,I,J,Val_k)

Fs=fss(w1');
% m=length(T.FNodePtrs);
% % Ks0=zeros(m,size(T.Nodes,1),size(w1,2));
Ks0=cell(1,size(w1,2));
% tic
for i=1:size(w1,2)
%     aa=k(co,w1(:,i)');
   K0= sparse(I,J,Val_k.*repmat(k(co,w1(:,i)'),16,1));
%     Ks0(:,:,i)=K00(T.FNodePtrs,T.FNodePtrs)\K0(T.FNodePtrs,:);
Ks0{i}=K0;
end



