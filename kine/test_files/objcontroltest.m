function objcontroltest(num_params)
% Create the controller figure window
controller = figure(...
    'Tag','controller',...
    'Units','normalized',...
    'Position',[0 0.05 0.3 .91],...
    'Name','KINE CONTROLLER',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'CloseRequestFcn','controller_closereq'...
    );

%---------------------------------------------------------------------
% Create color palate (could make this part of setup later)

data.colors.red     = [ 1.000 0.000 0.000 ];
data.colors.orange  = [ 1.000 0.502 0.000 ];
data.colors.yellow  = [ 1.000 1.000 0.000 ];
data.colors.green   = [ 0.000 0.502 0.000 ];
data.colors.cyan    = [ 0.000 1.000 1.000 ];
data.colors.blue    = [ 0.000 0.000 1.000 ];
data.colors.purple  = [ 0.502 0.000 0.502 ];
data.colors.magenta = [ 1.000 0.000 1.000 ];
data.colors.lemon   = [ 0.910 0.950 0.650 ];     
data.colors.lime    = [ 0.000 1.000 0.000 ];
data.colors.aqua    = [ 0.000 0.502 0.502 ];
data.colors.sky     = [ 0.000 0.502 1.000 ];
data.colors.pink    = [ 1.000 0.502 1.000 ];
data.colors.black   = [ 0.000 0.000 0.000 ]; 
data.colors.white   = [ 1.000 1.000 1.000 ];
data.colors.gray    = [ 0.850 0.850 0.850 ];       

background = get(controller,'Color');   % This should get the default background color
data.colors.background = background;    % which may be different on different systems
          
button_color = data.colors.lemon;
obj_color = data.colors.lemon;
data.colors.button_color = button_color;
data.colors.obj_color = obj_color;


obj_frame = uicontrol(...
    'Style','frame',...
    'Tag','obj_frame',...
    'Units','normalized',...
    'Position',[0.05 0.02 0.9 0.18],...
    'BackgroundColor',background,...
    'Visible','on'...
    );



obj_header = uicontrol(...
    'Style','text',...
    'Tag','obj_header',...
    'Units','normalized',...
    'Position',[.13 0.19 0.43 0.023],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',12,...
    'FontWeight','bold',...
    'String','Object Controls',...
    'Visible','on'...
    );


obj_frame_pos = get(obj_frame,'Position');
xmin = obj_frame_pos(1);
xmax = obj_frame_pos(1) + obj_frame_pos(3);
ymin = obj_frame_pos(2);
ymax = obj_frame_pos(2) + obj_frame_pos(4);

horiz_marg = xmin + 0.03;
vert_marg = ymin + 0.0075;
horiz_btwn = 0.05;
vert_btwn = 0.0075;
l_edit = 0.09;
h_edit = 0.025;
l_slider = 2*l_edit;
h_slider = h_edit;
text_gap = 0.02;
l_text = 0.10;
h_text = h_edit;

for i = 1:num_params
    
    if i <= 5
        x = horiz_marg;
        y = vert_marg + (5 - i)*(h_edit + vert_btwn);
        
    elseif i > 5 & i <= 10
        x = horiz_marg + 1*(l_slider + l_edit + text_gap + l_text + horiz_btwn);
        y = vert_marg + (10 - i)*(h_edit + vert_btwn);
        
    elseif i > 10 & i <= 15
        x = horiz_marg + 2*(l_slider + l_edit + text_gap + l_text + horiz_btwn);
        y = vert_marg + (15 - i)*(h_edit + vert_btwn);
        
    end
    
    figure(controller)
    
    s_max = 360;
    s_min = 0;
    s_step(1) = 5/(s_max - s_min );
    s_step(2) = 1/(s_max - s_min );
    
    uicontrol(...
        'Style','slider',...
        'Tag','obj_control',...
        'Units','normalized',...
        'Position',[x y l_slider h_slider],...
        'BackgroundColor',data.colors.background,...
        'Max',s_max,...
        'Min',s_min,....
        'SliderStep',s_step,...
        'Value',0,...
        'UserData',i,...
        'Callback','param_slider_callback'...
        );

    uicontrol(...
        'Style','edit',...
        'Tag','obj_control',...
        'Units','normalized',...
        'Position',[x+l_slider y l_edit h_edit],...
        'BackgroundColor','white',...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',12,...
        'FontWeight','normal',...
        'String',0,...
        'UserData',i,...
        'Callback','param_edit_callback'...
        );

    uicontrol(...
        'Style','text',...
        'Tag','obj_control',...
        'Units','normalized',...
        'Position',[x+l_slider+l_edit+text_gap y-.005 l_text h_text],...
        'BackgroundColor',data.colors.background,...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',11,...
        'FontWeight','normal',...
        'HorizontalAlignment','left',...
        'UserData',i,...
        'String','alpha'....
        );    
    
end
