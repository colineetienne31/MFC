function plot_airfoil(p)
%PLOT_AIRFOIL Draw the diamond airfoil only as a boundary marker.
x1=p.xAirfoil; x2=x1+p.c; xm=x1+p.c/2;
yt=p.tHalf;
plot([x1 xm x2],[0 yt 0],'k-','LineWidth',1.5);
end
