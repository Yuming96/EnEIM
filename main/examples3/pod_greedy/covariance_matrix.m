function C = covariance_matrix(cv,mesh1,mesh2,spthresh)
% COVARIANCE_MATRICES builds the covariance matrices
%
% C = covariance_matrix(cv,mesh1)
% C = covariance_matrix(cv,mesh1,mesh2)
% C = covariance_matrix(cv,mesh1,mesh2,spthresh)
%
% Outputs:
%   C:      The covariance matrix.
%
% Required inputs:
%   cv:     A function handle to a two point covariance function. The
%           function handle must be of the form cv=@(x1,x2) my_cov(x1,x2),
%           where x1 and x2 are row vectors giving the coordinates of two
%           points in the given mesh.
%
%   mesh1:  A matrix of size nx1 by d, where nx1 is the number of points in
%           the mesh and d is the dimension.
%
% Optional inputs:
%   mesh2:      A matrix of nx2 by d, where nx2 is the number of points in
%               the mesh and d is the dimension.
%
%   spthresh:   Threshold for setting an element of the covariance matrix
%               to zero in the sparse matrix. (Default 0 means include all
%               elements).
%
% This function is useful for large meshes where one wishes to build and
% store the covariance matrix for a given covariance function and mesh.
%
% COPYRIGHT 2010 Qiqi Wang (qiqi@mit.edu) and Paul G. Constantine 
% (pconsta@sandia.gov).

symm=0;
if ~exist('mesh2','var') || isempty(mesh2) 
    mesh2=mesh1; 
    symm=1;
end
if ~exist('spthresh','var') || isempty(spthresh), spthresh=0; end

nx1=size(mesh1,1); nx2=size(mesh2,1); 
if spthresh
    nx=max(nx1,nx2);
    ibuffer=ones(nx,1);
    jbuffer=ones(nx,1);
    Cbuffer=zeros(nx,1);
    count=1;
    if symm==1
        for i=1:nx1
            for j=i:nx2
                xi=mesh1(i,:); xj=mesh2(j,:);
                Cij=cv(xi,xj);
                if abs(Cij)>spthresh
                    if i~=j
                        ibuffer(count)=i; ibuffer(count+1)=j;
                        jbuffer(count)=j; jbuffer(count+1)=i;
                        Cbuffer(count)=Cij; Cbuffer(count+1)=Cij; 
                        count=count+2;
                    else
                        ibuffer(count)=i;
                        jbuffer(count)=j;
                        Cbuffer(count)=Cij;
                        count=count+1;
                    end
                    
                    % increase the buffer size
                    if count>length(ibuffer)
                        ibuffer=[ibuffer; ones(nx,1)];
                        jbuffer=[jbuffer; ones(nx,1)];
                        Cbuffer=[Cbuffer; zeros(nx,1)];
                    end
                end
            end
        end
    else
        for i=1:nx1
            for j=1:nx2
                xi=mesh1(i,:); xj=mesh2(j,:);
                Cij=cv(xi,xj);
                if abs(Cij)>spthresh
                    ibuffer(count)=i;
                    jbuffer(count)=j;
                    Cbuffer(count)=Cij;
                    count=count+1;
                    
                    % increase the buffer size
                    if count>length(ibuffer)
                        ibuffer=[ibuffer; ones(nx,1)];
                        jbuffer=[jbuffer; ones(nx,1)];
                        Cbuffer=[Cbuffer; zeros(nx,1)];
                    end
                end
            end
        end
    end
    C=sparse(ibuffer,jbuffer,Cbuffer);
else
    C=zeros(nx1,nx2);
    if symm
        for i=1:nx1
            for j=i:nx2
                xi=mesh1(i,:); xj=mesh2(j,:);
                Cij=cv(xi,xj);
                C(i,j)=Cij; C(j,i)=Cij;
            end
        end
    else
        for i=1:nx1
            for j=1:nx2
                xi=mesh1(i,:); xj=mesh2(j,:);
                C(i,j)=cv(xi,xj);
            end
        end
    end
end



