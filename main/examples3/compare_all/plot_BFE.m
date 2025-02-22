function varargout = plot_BFE(U_aux,Mesh)
% PLOT_BFE Plot finite element solution.
%
%   PLOT_BFE(U,MESH) generates a plot fo the finite element solution U on the
%   mesh MESH.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where M is
%                equal to the number of vertices contained in the mesh.
%    ELEMENTS    M-by-4 matrix specifying the elements of the mesh, where M is
%                equal to the number of elements contained in the mesh.
%
%   H = PLOT(U,MESH) also returns the handle to the figure.
%
%   Example:
%
%   plot_LFE(U,MESH);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;
  nElements = size(Mesh.Elements,1);
  
  % Compute axes limits
  
  XMin = min(Mesh.Coordinates(:,1));
  XMax = max(Mesh.Coordinates(:,1));
  YMin = min(Mesh.Coordinates(:,2));
  YMax = max(Mesh.Coordinates(:,2));
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1]; 
  YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
  
  % Compute the nodal value of the solution
  
  Coordinates = zeros(nElements*4,2);
  Elements = zeros(nElements,4);
  for i =1:nElements
     vidx = Mesh.Elements(i,:);     
     idx = 4*(i-1)+[1 2 3 4]; 
     Elements(i,:) = idx;
     Coordinates(idx,:) = Mesh.Coordinates(vidx,:);
  end
  U = zeros(4*size(U_aux,1),1);
  for i = 1:size(U_aux,1)
      U(4*i-3:4*i,1)=U_aux(i);
  end
  % Generate figure
    
  if(isreal(U))
  
    % Compute color axes limits 
      
    CMin = min(U);
    CMax = max(U);
    if(CMin < CMax)
      CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    else
      CLim = [1-OFFSET 1+OFFSET]*CMin;   
    end
    
    % Plot real finite element solution  
      
%     fig = figure('Name','Bilinear finite elements');
    patch('faces', Elements, ...
          'vertices', [Coordinates(:,1) Coordinates(:,2) U], ...
          'CData', U, ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    
    if(nargout > 0)
      varargout{1} = fig;
    end
 
  else  
      
    % Compute color axes limits 
      
    CMin = min([real(U); imag(U)]);
    CMax = max([real(U); imag(U)]);
    CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    
    % Plot imaginary finite element solution  
      
    fig_1 = figure('Name','Bi-Linear finite elements');
    patch('faces', Mesh.Elements, ...
          'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) real(U)], ...
          'CData', real(U), ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    fig_2 = figure('Name','Bi-Linear finite elements');
    patch('faces', Mesh.Elements, ...
          'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) imag(U)], ...
          'CData', imag(U), ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');  
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    if(nargout > 0)
      varargout{1} = fig_1;
      varargout{2} = fig_2;
    end
      
  end
  
return