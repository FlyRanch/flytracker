function [] = zoom2(fig_H,subfig_H,type)

%This is a special zoom function to handle the figures produced by
%'paste_image'.  Because they are two axes layered on top of each
%other (image, then model), we had to zoom in like this in a
%sequential manner.
%
%[] = zoom2(subax,cam,type)
% subax - the vector of axis handles defined by 'paste_image.m'
% cam - the number of the camera view to zoom in on
% type  - either the string 'in' or 'out'
% e.g - zoom2(1,'in') zooms in on figure 1 after using paste_image_auto 
% the 'out' command needs to be fixed, I think.

% ax_H = get(fig_H,'children');
% ax_H = flipud(ax_H);

ax_H = fig_H;

ax1 = ax_H(2*subfig_H-1); 
ax2 = ax_H(2*subfig_H);

switch lower(type)
 case 'in'
  axes(ax1)
  org_x = get(ax1,'xlim');
  org_y = get(ax1,'ylim');
  zoom on
  pause
  zoom off
  xrange = get(ax1,'xlim');
  yrange = get(ax1,'ylim');
  
  
  %yrange = fliplr(org_y(2) - yrange); %coordinate change to keep origin in lower left.
  
  set(ax2,'color','none','xlim',xrange,'ylim',yrange);
  axes(ax2);
 case 'out'
  org_x = [0 PAR.imgres(2)] +.5;
  org_y = [0 PAR.imgres(1)] +.5;
  set(ax1,'color','none','xlim',org_x,'ylim',org_y); 
  %set(ax2,'color','none','xlim',org_x,'ylim',fliplr(org_y(2) - org_y));
  set(ax2,'color','none','xlim',org_x,'ylim',org_y);
end
