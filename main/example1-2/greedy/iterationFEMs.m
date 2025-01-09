function G=iterationFEMs(T,K,F,g)
G=zeros(size(T.Nodes,1),1);
G(T.CNodePtrs)=g(T.Nodes(T.CNodePtrs,:));
U=K(T.FNodePtrs,T.FNodePtrs)\(F(T.FNodePtrs)-K(T.FNodePtrs,T.CNodePtrs)*G(T.CNodePtrs)); 
G(T.FNodePtrs)=U;

% if isempty(p)==0
% U=U(p);
% % U=U(:);
% end