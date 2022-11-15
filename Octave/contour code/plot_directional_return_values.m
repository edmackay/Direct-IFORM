function plot_directional_return_values(INPUTS,OUTPUTS,dims)

% set return periods
T = INPUTS.ReturnPeriods;
NRP = length(T);

% direction vectors
U = INPUTS.DirectionVectors;

% find direction vectors that are in desired plane
Ndim = size(U,2);
k = 1:Ndim;
for i = 1:length(dims)
    k(k==dims(i))=[];
end
inds = find(sum(abs(U(:,k)),2)==0);
Ud = U(inds,dims);

% corresponding angles
theta = mod(atan2(Ud(:,2),Ud(:,1))*180/pi,360);
[theta,I]= sort(theta);
inds=inds(I);

% plot angular peaks and return values
figure
hold on; box on
for i=1:length(inds)
    y=OUTPUTS.DIRECTIONAL.PEAKS{inds(i)};
    x=0*y+theta(i);
    h=plot(x,y,'.','color',[1 1 1]*0.7);
    %h.DisplayName='Declustered peaks';
    %if i>1
    %    h.HandleVisibility='off';
    %end

    % sose: octave's version
    set(h, 'DisplayName', 'Declustered peaks');
    if i>1
        set(h, 'HandleVisibility', 'off');
    end
end
h=plot(theta,OUTPUTS.DIRECTIONAL.GP.Threshold(inds),'k','LineWidth',2);
%h.DisplayName='Threshold';
set(h, 'DisplayName', 'Threshold');
for i=1:NRP
%     h=plotcol(theta,OUTPUTS.DIRECTIONAL.ReturnValues(inds,i),i,1,NRP+1,'-',2);
    h=plot(theta,OUTPUTS.DIRECTIONAL.ReturnValues(inds,i),'linewidth',2);
    %h.DisplayName=[num2str(T(i)) '-year RV'];
    set(h, 'DisplayName', [num2str(T(i)) '-year RV']);
end
xlim([0 360])
xticks(0:45:360)
xlabel('Projection angle [deg]')
ylabel('Projected value')
L=legend;
%L.Box='off'; L.Location='eastoutside';
set(L, 'Box', 'off'); set(L, 'Location', 'eastoutside');

