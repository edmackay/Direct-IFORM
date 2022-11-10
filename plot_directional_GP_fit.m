function plot_directional_GP_fit(INPUTS,OUTPUTS,U)

% direction vectors
Ucalc = INPUTS.DirectionVectors;
inds=find(ismember(Ucalc,U,'rows'));

% plot directional POT
dfws=get(0,'defaultfigurewindowstyle');
set(0,'defaultfigurewindowstyle','docked');
for i=1:length(inds)
    Zp = OUTPUTS.DIRECTIONAL.PEAKS{inds(i)};
    u = OUTPUTS.DIRECTIONAL.GP.Threshold(inds(i));
    sig = OUTPUTS.DIRECTIONAL.GP.Scale(inds(i));
    xi = OUTPUTS.DIRECTIONAL.GP.Shape(inds(i));
    
    figure
    gp_exceedance_plot(Zp,u,sig,xi,1)
    xlabel('Projected value')
    U1=num2str(Ucalc(inds(i),1),3);
    U2=num2str(Ucalc(inds(i),2),3);
    U3=num2str(Ucalc(inds(i),3),3);
    U4=num2str(Ucalc(inds(i),4),3);
%     title(['$\mathbf{u}=(' U1 ',' U2 ',' U3 ',' U4 ')$'])
    title(['u=(' U1 ',' U2 ',' U3 ',' U4 ')'])
    set(gca,'yminortick','off')
end
set(0,'defaultfigurewindowstyle',dfws);
