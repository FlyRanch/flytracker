%Make Radius updated plot

PAR.videopath = 'video/';
PAR.filetag = 'exp098';


PAR.solutionpath = [PAR.videopath '/solutions/'];
PAR.stub = [PAR.filetag];

%Load initial condition
load([PAR.solutionpath 'fly_' PAR.stub '/' 'ManualFit_' PAR.stub]);


params = ManualFit.params;
%% Plot the Canonical model
%body
len = params.bodylen;
radius = params.bodyrad_old;
sbody = linspace(-len/2,len/2,100);
Rbody =  radiusspline1(sbody',radius,length(sbody),0,len);

%head
len = params.headlen;
radius = params.headrad_old;
shead = linspace(-len/2,len/2,20);
Rhead =  radiusspline1(shead',radius,length(shead),0,len);

%plotting
%shift head profile function parameter to start at the end of the body one
shead = shead + len/2 + sbody(end);

figure; 
hh(1) = plot([sbody shead],[Rbody ; Rhead],'-');
axis equal

%% Plot the Updated model
%body
len = params.bodylen;
radius = params.bodyrad;
sbody = linspace(-len/2,len/2,100);
Rbody =  radiusspline1(sbody',radius,length(sbody),0,len);

%head
len = params.headlen;
radius = params.headrad;
shead = linspace(-len/2,len/2,20);
Rhead =  radiusspline1(shead',radius,length(shead),0,len);

%%plotting
%shift head profile function parameter to start at the end of the body one
shead = shead + len/2 + sbody(end);

hold on;
%hh(2) = plot([sbody shead].*params.bodyscale,[Rbody ; Rhead].*params.bodyscale,'r.--');
hh(2) = plot([sbody shead],[Rbody ; Rhead],'r--');
axis equal
legend(hh,'Original','Updated')
