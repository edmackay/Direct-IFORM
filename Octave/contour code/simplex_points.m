function w = simplex_points(ndim,npoints)

% creates array of points on unit simplex that is surface of L1 sphere in ndim-dimensions
% inputs:
%   ndim = number of dimensions
%   npoints = number of points along edge of simplex
% outputs:
%   w = [M x ndim] array of coordinates of regularly spaced points on L1 sphere

% create grid of points with spacing dx
dx=1/(npoints-1); % spacing along axis
xin=0:dx:1;
sizearray(1:ndim) = npoints;
s=sizearray;
s(1)=1;
grids = cell(1,ndim);
gridsum=repmat(zeros(npoints,1),s);
for i=1:ndim
    x = xin;
    s = ones(1,ndim);
    s(i) = numel(x);
    x = reshape(x,s);
    s = sizearray;
    s(i) = 1;
    grids{i} = repmat(x,s);
    gridsum=gridsum+grids{i};
end

% keep only points which lie on L1 unit sphere
good = gridsum>1-dx/1e4 & gridsum<1+dx/1e4;
for i=1:ndim
    x=grids{i};
    grids{i} = x(good);
end

% output points
M=nchoosek(npoints+ndim-2,ndim-1); % total number of points
w=zeros(M,ndim);
for i=1:ndim
    w(:,i)=grids{i};
end
