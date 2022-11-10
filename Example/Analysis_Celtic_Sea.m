clear

% addpath for contour functions
addpath('Contour code')

% load data
data=load('Data_Celtic_Sea','Hs','Tm','dir_rel','U10');

% set variables
U10=data.U10;
Tm=data.Tm;
HL=sqrt(data.Hs).*cos(data.dir_rel*pi/180);
HT=sqrt(data.Hs).*sin(data.dir_rel*pi/180);
clear data

% set input data structure
DATA = [U10, HL, HT, Tm];

% set directions
ndim=4;
npoints=11;
posneg=[0 0 0 0];
U = direction_vectors(ndim,npoints,posneg);

% define INPUT data structure
INPUTS.Data = DATA;
INPUTS.TimeStep = 1;       
INPUTS.PeakSepTime = 48;   
INPUTS.ThreshExProb = 0.1; 
INPUTS.DirectionVectors = U; 
INPUTS.ReturnPeriods = [1 5 50]; % return periods to calculate contours

% fit contours
OUTPUTS = fit_contour(INPUTS);
save('Celtic_Sea_contour_data','INPUTS','OUTPUTS')

