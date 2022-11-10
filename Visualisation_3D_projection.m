clear

load('Celtic_Sea_contour_data.mat')
load('Data_Celtic_Sea.mat')

HL=sqrt(Hs).*cos(dir_rel*pi/180);
HT=sqrt(Hs).*sin(dir_rel*pi/180);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D surface: U10-Hlong-Htrans
c=OUTPUTS.CONTOUR{3}; % select 50-year contour
TRI=convhulln(c(:,1:3)); % calculate triangulation for dimensions 1-3

figure
ax=axes;
hold on; grid on
scatter3(U10,HL,HT,'k.')
p=trisurf(TRI,c(:,1),c(:,2),c(:,3));
p.EdgeColor="none";
p.FaceColor="interp";
p.FaceAlpha=0.7;
rotate3d on
L=light;
lighting gouraud
view([-40 35])
axis([0 25 -2 4 0 3])
xticks(0:5:25)
yticks(-3:4)
zticks(0:3)
ax.PlotBoxAspectRatio=[2 2 1];
axis vis3d
xlabel('U_{10} [m/s]')
ylabel('H_{L} [m^{1/2}]')
zlabel('H_{T} [m^{1/2}]')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D surface: U10-Hs-Tm
c=OUTPUTS.CONTOUR{3}; % select 50-year contour
cU=c(:,1);
cHs=c(:,2).^2+c(:,3).^2; % calculate Hs
cTm=c(:,4);
TRI=convhulln([cU cTm sqrt(cHs)]); % calculate triangulation

figure
ax=axes;
hold on; grid on
scatter3(U10,Tm,Hs,'k.')
p=trisurf(TRI,cU,cTm,cHs);
p.EdgeColor="none";
p.FaceColor="interp";
p.FaceAlpha=0.7;
rotate3d on
L=light;
L.Position=[9.4 -15.8   18.7];
lighting gouraud
view([-40 35])
axis([0 25 0 15 0 10])
ax.PlotBoxAspectRatio=[2 1.5 1];
xlabel('U_{10} [m/s]')
ylabel('T_m [s]')
zlabel('H_s [m]')

