function [Am,I,J,Mm,Ib,Jb] = Mifemofof(T,hx,hy)
I=[T.nelem(:); T.nelem(:,1);T.nelem(:,1);T.nelem(:,2);T.nelem(:,3);T.nelem(:,2);T.nelem(:,3);T.nelem(:,4);T.nelem(:,4);T.nelem(:,1);T.nelem(:,2);T.nelem(:,4);T.nelem(:,3)];
J=[T.nelem(:); T.nelem(:,2);T.nelem(:,3);T.nelem(:,4);T.nelem(:,4);T.nelem(:,1);T.nelem(:,1);T.nelem(:,2);T.nelem(:,3);T.nelem(:,4);T.nelem(:,3);T.nelem(:,1);T.nelem(:,2)];
Am=[2/3*ones(size(T.nelem(:))); -1/6*ones(size([ T.nelem(:,1);T.nelem(:,1);T.nelem(:,2);T.nelem(:,3);T.nelem(:,2);T.nelem(:,3);T.nelem(:,4);T.nelem(:,4)]));...
    -1/3*ones(size([T.nelem(:,1);T.nelem(:,2);T.nelem(:,4);T.nelem(:,3)]))];
Ib=[T.nelem(:,1); T.nelem(:,1);T.nelem(:,1);T.nelem(:,1);T.nelem(:,2);T.nelem(:,2);T.nelem(:,2);T.nelem(:,2);T.nelem(:,3);T.nelem(:,3);T.nelem(:,3);T.nelem(:,3);...
    T.nelem(:,4);T.nelem(:,4);T.nelem(:,4);T.nelem(:,4)];
Jb=[T.nelem(:,1); T.nelem(:,2);T.nelem(:,3);T.nelem(:,4);T.nelem(:,1);T.nelem(:,2);T.nelem(:,3);T.nelem(:,4);T.nelem(:,1);T.nelem(:,2);T.nelem(:,3);T.nelem(:,4);...
    T.nelem(:,1);T.nelem(:,2);T.nelem(:,3);T.nelem(:,4)];
Mm=[hx*hy/9*ones(size(T.nelem(:,1))); hx*hy/18*ones(size(T.nelem(:,1))); hx*hy/18*ones(size(T.nelem(:,1))); hx*hy/36*ones(size(T.nelem(:,1))); hx*hy/18*ones(size(T.nelem(:,1))); ...
    hx*hy/9*ones(size(T.nelem(:,1))); hx*hy/36*ones(size(T.nelem(:,1))); hx*hy/18*ones(size(T.nelem(:,1))); hx*hy/18*ones(size(T.nelem(:,1))); hx*hy/36*ones(size(T.nelem(:,1)));...
    hx*hy/9*ones(size(T.nelem(:,1))); hx*hy/18*ones(size(T.nelem(:,1))); hx*hy/36*ones(size(T.nelem(:,1))); hx*hy/18*ones(size(T.nelem(:,1))); hx*hy/18*ones(size(T.nelem(:,1))); ...
    hx*hy/9*ones(size(T.nelem(:,1)))];
end

