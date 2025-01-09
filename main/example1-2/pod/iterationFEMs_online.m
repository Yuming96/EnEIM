function Ur=iterationFEMs_online(T,K,F,g,Phi)
% G=zeros(size(T.Nodes,1),1);
% G(T.CNodePtrs)=g(T.Nodes(T.CNodePtrs,:));
Gr=Phi'*K(T.FNodePtrs,:)*g(T.Nodes);
Kr=Phi'*K(T.FNodePtrs,T.FNodePtrs)*Phi;
Fr=Phi'*F(T.FNodePtrs);
Ur=Kr\(Fr-Gr); 
% G(T.FNodePtrs)=Phi*Ur;

% if isempty(p)==0
% U=U(p);
% % U=U(:);
% end