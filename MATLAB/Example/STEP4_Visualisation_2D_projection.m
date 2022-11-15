clear

load('Celtic_Sea_contour_data.mat')
load('Data_Celtic_Sea.mat')

Hs_long=Hs.*cos(dir_rel*pi/180);
Hs_trans=Hs.*sin(dir_rel*pi/180);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2D projections: Hs-Tm
figure
hold on; grid on; box on
plot(Tm,Hs,'.','color',[1 1 1]*0.6)
for irp=1:length(INPUTS.ReturnPeriods)
    c=OUTPUTS.CONTOUR{irp};
    Hx=c(:,2);
    Hy=c(:,3);
    c_H=sqrt(Hx.^2 + Hy.^2);
    c_Hs=c_H.^2;
    c_Tm=c(:,4);
    ind=convhull(c_Tm,c_H);
    plot(c_Tm(ind),c_Hs(ind),'LineWidth',2)
end
xlabel('$T_m$ [s]','Interpreter','latex')
ylabel('$H_s$ [m]','Interpreter','latex')
axis([0 16 0 10])
xticks(0:2:16)
yticks(0:2:10)
L=legend('Data','1-year','5-year','50-year');
L.Box='off';
L.Location="northwest";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D projections: Hs_long-Hs_trans
figure
ax=axes;
hold on; box on
plot(Hs_long,Hs_trans,'.','color',[1 1 1]*0.6)
for irp=1:length(INPUTS.ReturnPeriods)
    c=OUTPUTS.CONTOUR{irp};
    % convert back to Hs
    c_Hx=c(:,2);
    c_Hy=c(:,3);
    c_Hs=c_Hx.^2+c_Hy.^2;
    c_dir=atan2(c_Hy,c_Hx);
    c_Hsx=c_Hs.*cos(c_dir);
    c_Hsy=c_Hs.*sin(c_dir);
    
    ind=convhull(c_Hsx,c_Hsy);
    plot(c_Hsx(ind),c_Hsy(ind),'LineWidth',2)
end
radial_grid(2:2:20,pi/8:pi/8:pi-pi/8)
axis([-6 10 0 10])
xticks(-10:2:10)
ax.PlotBoxAspectRatio=[1.6 1 1];
xlabel('$H_{s,long}$ [m]','Interpreter','latex')
ylabel('$H_{s,trans}$ [m]','Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D projections: U10-Hs
figure
axis([0 25 0 10])
xticks(0:5:30)
yticks(0:2:30)
hold on; box on; grid on
plot(U10,Hs,'.','color',[1 1 1]*0.6)
for irp=1:length(INPUTS.ReturnPeriods)
    c=OUTPUTS.CONTOUR{irp};
    c_U=c(:,1);
    c_Hx=c(:,2);
    c_Hy=c(:,3);
    c_H=sqrt(c_Hx.^2 + c_Hy.^2);
    c_Hs=c_H.^2;
    ind=convhull(c_U,c_H);
    plot(c_U(ind),c_Hs(ind),'LineWidth',2)
end
xlabel('$U_{10}$ [m/s]','Interpreter','latex')
ylabel('$H_s$ [m]','Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D projections: U10-Tm
figure
axis([0 25 0 16])
xticks(0:5:30)
yticks(0:2:30)
hold on; box on; grid on
plot(U10,Tm,'.','color',[1 1 1]*0.6)
for irp=1:length(INPUTS.ReturnPeriods)
    c=OUTPUTS.CONTOUR{irp};
    c_U=c(:,1);
    c_T=c(:,4);
    ind=convhull(c_U,c_T);
    plot(c_U(ind),c_T(ind),'LineWidth',2)
end
xlabel('$U_{10}$ [m/s]','Interpreter','latex')
ylabel('$T_m$ [s]','Interpreter','latex')
