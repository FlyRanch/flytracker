% Make plots to illustrate the roll correction

p = sol;
pnew = fixRollTPlane(p,PAR);

%p = pnew;

[x,y,z] = flymod(p,PAR.params,PAR);
dorsalpts = [x{1}(:,10),y{1}(:,10),z{1}(:,10)];


[Gbody,GL,GR] = GTransform(p,PAR);

Vr = -GR(1:3,2);
Vl = GL(1:3,2);
%Coordinate frame attached to body
Zb = Gbody(1:3,3);
Yb = Gbody(1:3,2);
Xb = Gbody(1:3,1);
Tb = Gbody(1:3,4);

TL = GL(1:3,4);
TR = GR(1:3,4);
%% Plotting
figure; 
for i = 1:length(x)
    surf(x{i},y{i},z{i},'edgecolor','k','facecolor','r','facelighting','phong');
    hold on;
end

%% Vectors plot
hold on; quiver3(Tb(1),Tb(2),Tb(3),Zb(1),Zb(2),Zb(3),'k','linewidth',3)
hold on; quiver3(Tb(1),Tb(2),Tb(3),Yb(1),Yb(2),Yb(3),'k','linewidth',3)

hold on; quiver3(TL(1),TL(2),TL(3),Vl(1),Vl(2),Vl(3),'b','linewidth',3)
hold on; quiver3(TR(1),TR(2),TR(3),Vr(1),Vr(2),Vr(3),'b','linewidth',3)


%% Transverse plane plot
w = 1;
u = [-1 1];
v = [-1 1];
%Make transverse plane surface 
for i=1:length(u)
    for j = 1:length(v)
        pp = Tb + u(i)*w.*Yb + v(j)*w.*Zb;
        xp(i,j) = pp(1);
        yp(i,j) = pp(2);
        zp(i,j) = pp(3);
    end
end
hold on; 
surf(xp,yp,zp,'facecolor','k','facealpha',.3,'edgecolor','none')
