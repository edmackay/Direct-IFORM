clear

load('Celtic_Sea_contour_data.mat')
load('Data_Celtic_Sea.mat')

Hs_long=Hs.*cos(dir_rel*pi/180);
Hs_trans=Hs.*sin(dir_rel*pi/180);

% wind speed for slice
U10_slice=11;

% take slice through 3D contour
dim=1;
v3D = slice_contour(OUTPUTS.CONTOUR,dim,U10_slice);

% calculate bin for wind speeds
binU=U10>=U10_slice-1 & U10<=U10_slice+1;

% calculate 2D slices
sp=0;
figure
for t=4:12
    bin = binU & t-0.5<=Tm & Tm<t+0.5;

    sp=sp+1;
    ax=subplot(3,3,sp);
    hold on; box on
    title(['T_m = ' num2str(t) ' s'])
    radial_grid(2:2:10,0:pi/8:pi)
    axis([-6 8 0 8])
    set(ax, 'PlotBoxAspectRatio', [14 8 1]);
    xlabel('H_{s,long}')
    ylabel('H_{s,trans}')

    % plot observations
    plot(Hs_long(bin),Hs_trans(bin),'.','color',[1 1 1]*0.6)

    % calculate 2D slices
    v2D = slice_contour(v3D,3,t);
    for i=1:3
        v=v2D{i};
        if ~isempty(v)
            % convert to Hs and rel dir
            rel_dir=atan2(v(:,2),v(:,1));
            vHs=v(:,1).^2+v(:,2).^2;
            v1=vHs.*cos(rel_dir);
            v2=vHs.*sin(rel_dir);
            TRI=convhull(v(:,1),v(:,2));
            plot(v1(TRI),v2(TRI),'LineWidth',2)
        end
    end
end

