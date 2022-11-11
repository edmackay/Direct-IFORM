clear

% load data
load('Celtic_Sea_contour_data.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check projected return values
dims=[2 3];
plot_directional_return_values(INPUTS,OUTPUTS,dims)
ylim([-3 6])

dims=[1 4];
plot_directional_return_values(INPUTS,OUTPUTS,dims)
ylim([-3 6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set direction vectors to check GP fit
U=[1 0 0 0;
   0 1 0 0;
   0 0 1 0;
   0 0 0 1;
   1 1 0 0;
   1 0 1 0;
   1 0 0 1;
   0 1 1 0;
   0 1 0 1;
   0 0 1 1;
   0 -1 0 0;
   1 -1 0 0];
r=sqrt(sum(U.^2,2));
U=U./repmat(r,1,4);

% plot fit
plot_directional_GP_fit(INPUTS,OUTPUTS,U)
