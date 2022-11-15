function OUTPUTS = fit_contour(INPUTS)

% INPUTS data structure with fields:
%   INPUTS.Data               - Ndat x Ndim array of observations of Ndim-dimensional random vector X
%   INPUTS.TimeStep           - time step of original data in hours (assumed that there are no missing time steps)
%   INPUTS.PeakSepTime        - separation of storm peaks in hours
%   INPUTS.ThreshExProb       - GP threshold exceedance probability
%   INPUTS.DirectionVectors   - Ndir x Ndim array of unit vectors, specifying directions for analysis
%   INPUTS.ReturnPeriods      - NRP x 1 array of return periods (years) to calculate contours for
%   
% OUTPUTS data structure with fields:
%   OUTPUTS.NORMALISATION.STD         - 1 x Ndim array of STD of each dimension of data
%   OUTPUTS.NORMALISATION.Median      - 1 x Ndim array of STD of each dimension of data
%   OUTPUTS.DIRECTIONAL.PEAKS         - Ndir x 1 cell array of peaks of projected data
%   OUTPUTS.DIRECTIONAL.GP.Threshold  - Ndir x 1 array of GP threshold parameter
%   OUTPUTS.DIRECTIONAL.GP.Scale      - Ndir x 1 array of GP scale parameter
%   OUTPUTS.DIRECTIONAL.GP.Shape      - Ndir x 1 array of GP shape parameter
%   OUTPUTS.DIRECTIONAL.ReturnValues  - Ndir x NRP array of angular return values
%   OUTPUTS.CONTOUR                   - NRP x 1 cell array of coordinates of contours

% Parse inputs
X = INPUTS.Data;
[Ndat, Ndim] = size(X);
zeta = INPUTS.ThreshExProb;
dt_data = INPUTS.TimeStep;
dt_peaks = INPUTS.PeakSepTime;
RP = INPUTS.ReturnPeriods;
NRP = length(RP);
U = INPUTS.DirectionVectors;
Ndir = size(U,1);
Nyears = length(X)/((24/dt_data)*365.25); % dataset length in years

% Normalise data
Xstd = std(X,'omitnan');
Xmed = median(X,'omitnan');
Y = (X-repmat(Xmed,Ndat,1))./repmat(Xstd,Ndat,1);

% Define output structures
OUTPUTS.NORMALISATION.STD = Xstd;
OUTPUTS.NORMALISATION.Median = Xmed;
OUTPUTS.DIRECTIONAL.PEAKS = cell(Ndir,1);
OUTPUTS.DIRECTIONAL.GP.Threshold = zeros(Ndir,1);
OUTPUTS.DIRECTIONAL.GP.Scale = zeros(Ndir,1);
OUTPUTS.DIRECTIONAL.GP.Shape = zeros(Ndir,1);
OUTPUTS.DIRECTIONAL.ReturnValues = zeros(Ndir,NRP);
OUTPUTS.CONTOUR = cell(NRP,1);

% Peaks-over-threshold analysis for each direction
disp('Conducting POT analysis for:')
opts=optimset('display','off');
for i=1:Ndir
    disp(['   Direction ' num2str(i) '/' num2str(Ndir)])
    
    % define projected data
    R = Y * (U(i,:))';
    
    % find peaks in projected data time series
    Rpeaks = decluster(R,dt_peaks/dt_data);
    
    % set threshold 
    n = length(Rpeaks);
    P = (1:n)/(n+1);
    mu = interp1(P,sort(Rpeaks),1-zeta);
    
    % define threshold exceedances
    Z = Rpeaks - mu;
    Z = Z(Z>0);
    
    % fit generalised Pareto distribution to threshold exceedances
    [sigma, xi]=gp_fit(Z,-1,0,opts);
    
    % calculate return values
    lambda = length(Z)/Nyears;
    if xi~=0
        RV = mu + sigma/xi*((RP*lambda).^xi - 1);
    else
        RV = mu + sigma*log(RP*lambda);
    end
    
    % set output variables
    OUTPUTS.DIRECTIONAL.PEAKS{i} = Rpeaks;
    OUTPUTS.DIRECTIONAL.GP.Threshold(i) = mu;
    OUTPUTS.DIRECTIONAL.GP.Scale(i) = sigma;
    OUTPUTS.DIRECTIONAL.GP.Shape(i) = xi;
    OUTPUTS.DIRECTIONAL.ReturnValues(i,:) = RV;
end

% find intersection of half-space regions
for i=1:NRP
    % create reflection points
    S=2*repmat(OUTPUTS.DIRECTIONAL.ReturnValues(:,i),1,Ndim).*U;
    
    % calculate contour
    [verts,cells]=voronoin([zeros(1,Ndim); S]);
    verts=verts(cells{1},:); % select cell containing origin
    
    % transform back to original margins
    OUTPUTS.CONTOUR{i}=verts.*Xstd + Xmed;
end

