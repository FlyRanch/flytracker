function IMout = renderfly(fignum,imgres)

%This function will create a surface plot of the points in xx,yy,zz and
%will plot them on the screen white with a black background and capture the
%image at the specified resolution.
figure(fignum);
ax1 = gca;
surfh = get(ax1,'children');
for i = 1:length(surfh)
    set(surfh(i),'facecolor','w','edgecolor','none');
end
set(fignum,'color','k');

%set(ax1,'color','none','visible','off');
set(ax1,'units','pixels','position',[1 1 imgres(2) imgres(1)],'color','none',...
    'visible','off');

figure(fignum);
set(fignum,'units','pixels','position',[50 200 imgres(2) imgres(1)]);
 
IM = getframe(fignum,[1 1 imgres(2) imgres(1)]);
%close(1);

IMout = double(rgb2gray(frame2im(IM)));