function [U, W] = direction_vectors(ndim,npoints,posneg)

% posneg = 1 x ndim array, with entry:
% posneg(i)=1 if X(:,i) is positive
% posneg(i)=-1 if X(:,i) is negative
% posneg(i)=0 if X(:,i) is both positive and negative

% set L1 directions
W = simplex_points(ndim,npoints);

% add negative directions
for i=1:ndim
    if posneg(i)==1
        % add extra single negative value
        line=zeros(1,ndim);
        line(i)=-1;
        W=[W;line];
    elseif posneg(i)==-1
        % flip and add extra single positive value
        W(:,i)=-W(:,i);
        line=zeros(1,ndim);
        line(i)=1;
        W=[W;line];
    elseif posneg(i)==0
        % add negative directions
        Wneg=W;
        remove=Wneg(:,i)==0;
        Wneg(remove,:)=[];
        Wneg(:,i)=-Wneg(:,i);
        W=[W;Wneg];
    end
end

% normalise
R=sqrt(sum(W.^2,2));
U=W./repmat(R,1,ndim);