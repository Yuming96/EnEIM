function Mesh = TProd_Mesh(x,y)
%TPROD_MESH Create Tensor-Product mesh
%
%   MESH = TPROD_MESH(X,Y) initializes a tensor-product mesh with vertices
%   at (X(i),Y(j)).  The vectors X and Y must be sorted in ascending order.
%
%   MESH = TPROD_MESH(X) is short for MESH = TPROD_MESH(X,X);
%
%   Example :
%
%     Mesh = TProd_Mesh(0:1/n:1);

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Handle input arguments
  if(nargin < 2)
    y = x;
  end
  
  % Create coordinates
  [X,Y] = meshgrid(x,y);
  Mesh.Coordinates = [X(:),Y(:)];
  
  % Create elements
  nX = length(x);
  nY = length(y);
  nElem = (nX-1)*(nY-1);
  Mesh.Elements = zeros(nElem,4);
  for i=1:nX-1
    Mesh.Elements((i-1)*(nY-1)+(1:nY-1),1) = (i-1)*nY + (1:nY-1);
    Mesh.Elements((i-1)*(nY-1)+(1:nY-1),4) = (i-1)*nY + (2:nY);
    Mesh.Elements((i-1)*(nY-1)+(1:nY-1),2) = i*nY + (1:nY-1);
    Mesh.Elements((i-1)*(nY-1)+(1:nY-1),3) = i*nY + (2:nY);
  end
  for i=1:size( Mesh.Elements,1)
    Mesh.midnodes(i,:)=mean(Mesh.Coordinates(Mesh.Elements(i,:),:));
  end
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);

% Numbering of loca vertices and edges
%    4 ___ 3          ___
%     |   |          | 4 |
%     |___|        1 |___| 2
%    1     2           3

el_node = sparse([1:nElements 1:nElements 1:nElements 1:nElements],Mesh.Elements(:),[ones(1,nElements) 2*ones(1,nElements) 3*ones(1,nElements) 4*ones(1,nElements)]);
node_node = el_node'*el_node;
node_node = node_node-diag(diag(node_node));
node_node(node_node==3)=0;
node_node(node_node==8)=0;
node_node = tril(node_node);
% full(node_node)
[In,Jn,~]=find(node_node);
Mesh.nEdges = numel(In);
Mesh.Vert2Edge = sparse(In,Jn,1:numel(In),nCoordinates,nCoordinates);
Mesh.Vert2Edge = Mesh.Vert2Edge + Mesh.Vert2Edge'; 
Mesh.EdgeLoc = [Mesh.Vert2Edge(Mesh.Elements(:,4)+(Mesh.Elements(:,1)-1)*nCoordinates) Mesh.Vert2Edge(Mesh.Elements(:,2)+(Mesh.Elements(:,3)-1)*nCoordinates) ...
Mesh.Vert2Edge(Mesh.Elements(:,1)+(Mesh.Elements(:,2)-1)*nCoordinates) Mesh.Vert2Edge(Mesh.Elements(:,4)+(Mesh.Elements(:,3)-1)*nCoordinates)]';
Mesh.nElements = size(Mesh.Elements,1);
return