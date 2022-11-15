function [xmax,imax,winmax] = decluster(x,r)

% [xmax,imax] = decluster(x,r)
%
% finds maxima of series x defined as maxima within sliding window of +/- r samples
% assumes time series x is continuous

% find maximum in window of +/- r samples
winmax = movmax(x,[r r]);

% define peaks
imax = find(x==winmax);
xmax = x(imax);

% go through peaks to check if there are any which are too close
% (this occurs if there are equal maxima within a window)
while 1
    bad=find(diff(imax)<=r,1,'first');
    if isempty(bad)
        break
    end
    % remove second peak
    xmax(bad+1)=[];
    imax(bad+1)=[];
end
