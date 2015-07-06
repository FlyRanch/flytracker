% Function GET_ARC_3D
% 
% CALLING FUNCTION: get_arc, mark_points_second_view
% ACTIONS: checks to see if all arc info entered and then resolves
% i        into 3d coordinates using dlt coefficients and polyfit
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 5, 2004 by gwyneth (modified from resolve_3d by Will
%                Dickson)

function get_arc_3d

controller = findobj('Tag','controller');
data = guidata(controller);

cur_obj = get_string(data.handles.cur_object_menu);
frame = str2num(get(data.handles.frame_box,'String'));
%poly_order = data.kine.(cur_obj).info.poly_order;  %%%%CHANGE ONCE EXISTS
poly_order = 3;
x = 1;
y = 2;
z = 3;

if data.kine.(cur_obj).data.coords(1,1,frame) ~= 0 & data.kine.(cur_obj).data.coords(end,1,frame) ~= 0 ...
        & isempty(data.kine.(cur_obj).info.temp) == 0
    
    % Fit a polynomial to points in each view
    for view = 1:2
        
        % Get 2d projection of start and end points
        cam = data.kine.(cur_obj).info.temp.(['cam',num2str(view)])(end);
        DLT = data.cal.coeff.(['DLT_',cam]);_
        
        [s_u s_v] = dlt_3D_to_2D(DLT, data.kine.(cur_obj).data.coords(1,x,frame),...
                                      data.kine.(cur_obj).data.coords(1,y,frame),...
                                      data.kine.(cur_obj).data.coords(1,z,frame));
                                  
        [e_u e_v] = dlt_3D_to_2D(DLT, data.kine.(cur_obj).data.coords(end,x,frame),...
                                      data.kine.(cur_obj).data.coords(end,y,frame),...
                                      data.kine.(cur_obj).data.coords(end,z,frame));
                                  
        % Retrieve u,v ginput arc points for this view
        u = data.kine.(cur_obj).info.temp.(['u',num2str(view)]);
        v = data.kine.(cur_obj).info.temp.(['v',num2str(view)]);
        
        % Arrange in point array
        pt_array = [s_u s_v;...
                      u   v;...
                    e_u e_v];
        
        % Fit polynomial to this view
        [ poly(view).u, poly(view).v ] = fit_poly( pt_array, poly_order );
        
    end % should now have structure 'poly' with polynomial values for arcs for both views stored
    
    
    