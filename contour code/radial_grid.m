function radial_grid(radii,angles)

ax=axis;
hold on
t=linspace(0,2*pi,1000);
for r=radii
    h=plot(r*cos(t),r*sin(t),'color',[1 1 1]*0.8,'LineWidth',0.5);
    h.HandleVisibility='off';
    uistack(h,'bottom')
end
r=[0 max(radii)];
for t=angles
    h=plot(r*cos(t),r*sin(t),'color',[1 1 1]*0.8,'LineWidth',0.5);
    h.HandleVisibility='off';
    uistack(h,'bottom')
end
axis([-1 1 -1 1]*max(abs(ax)));
axis square
