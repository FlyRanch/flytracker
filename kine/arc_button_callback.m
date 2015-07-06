% Function ARC_BUTTON_CALLBACK
% 
% CALLING FUNCTION: Callback for arc_button
% ACTIONS: a program to get points for arc and resolve in 3d
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 5, 2004 by gwyneth (modified from arc programs
%                by Will Dickson)

function arc_button_callback

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

% Set pointer function to zoom
toggle_pointer_fcn('Zoom')

cur_obj = get_string(data.handles.cur_object_menu);
frame = str2num(get(data.handles.frame_box,'String'));

poly_order = 3;        % poly_order = data.kine.(cur_obj).info.poly_order;  %%%%CHANGE ONCE EXISTS
num_pts = 10;          % number of points to click on along arc
                       % should get from data.kine.(cur_obj) structure
num_poly_pts = 10;     % number of points plotted for line and arc to be compared, 
                       % want to be high because will find closest match numerically by comparison
                        
pts_3d = [];

% First mark start and end points using mark_point scripts
for pt = 1:2 % this for loop is for 1 = starting point and 2 = ending point

    redo = 1;
    while redo == 1
        % Set pointer function to zoom%%%DEBUG find out why need to inside loop
        toggle_pointer_fcn('Zoom')

        dig_data = guidata(dig_fig);
            
        if pt == 1
            msgbox('Press RETURN, then click on the arc STARTING POINT in two camera views')
            mark = 'r^';
        elseif pt == 2
            msgbox('Press RETURN, then click on the arc ENDING POINT in the same two camera views')
            mark = 'ro';
        end
        pause
        
        % First point
        [u v] = ginput(1);
        cam = get(gca,'Tag');
        
        if pt == 2                                              % make sure first view for start point and end points are the same cam
            while strcmp(endpts(1).cam,cam) == 0
                warndlg(['That is not the image of the first camera view!  You want ',endpts(1).cam],'Incorrect Camera View')
                [u v] = ginput(1);
                cam = get(gca,'Tag');
            end
        end
        
        endpts(1).cam = cam;
        cam = str2num(cam(end));
        others = find([ 1 2 3] ~= cam);
        
        endpts(1).u(pt) = u;
        endpts(1).v(pt) = v;
        
        plot( u, v, mark, 'MarkerFaceColor','r', 'Tag', 'pt')                          % Plot marked point
        plot_lines( u, v, cam, others )   % Plot corresponding lines in the other views
        
        % Second point 
        [u v] = ginput(1);
        cam = get(gca,'Tag');
        
        if pt == 2                                              % make sure second view for start point and end points are the same cam
            while strcmp(endpts(2).cam,cam) == 0
                warndlg(['That is not the image of the second camera view!  You want ',endpts(2).cam],'Incorrect Camera View')
                [u v] = ginput(1);
                cam = get(gca,'Tag');
            end
        end
        
        endpts(2).cam = cam;
        
        dig_data = guidata(dig_fig);                            % Delete lines
        points_to_delete = fieldnames(dig_data.temp);           
        for i = 1:length(points_to_delete)
            delete(dig_data.temp.(points_to_delete{i}))
            dig_data.temp = rmfield(dig_data.temp,points_to_delete{i});
        end
        guidata(dig_fig,dig_data)
        
        endpts(2).u(pt) = u;
        endpts(2).v(pt) = v;

        plot( u, v, mark, 'MarkerFaceColor','r', 'Tag', 'pt');   % Plot marked point
        
        % Redo if these points are not acceptable
        proceed = questdlg('Accept this point?','Accept Point','Yes','Redo','Cancel','Yes');
        switch proceed
            case 'Yes'
                redo = 0;
            case 'Redo'
                to_clear = findobj('Tag','pt');
                delete(to_clear)
                redo == 1;
            case 'Cancel'
                to_clear = findobj('Tag','pt');
                delete(to_clear)
                redo == 1;
                return
        end                                                     % end proceed switch
    end                                                         % redo loop
end                                                             % for start/end pts


% Then get points along arcs in two views
for view = 1:2   % this for loop is for 1 = cam view 1, 2 = cam view 2
    
    if view == 1, txt = 'FIRST'; elseif view == 2, txt = 'SECOND'; end
    
    redo = 1;
    while redo == 1
        
        question = ['Click on ',num2str(num_pts),' points along the arc (excluding the arc start and end points) in order from the start (triangle) to the end (circle) points, in the ',txt,' camera view'];
        title = ['Get Arc Points:  ',txt,' view'];
        proceed = questdlg(question,title,'OK','Cancel','OK');
        
        switch proceed
            
            case 'OK'
                [u,v] = ginput(num_pts);                                        % Ginput points for arc view 1
                cam = get(gca,'Tag');
                
                while strcmp(endpts(view).cam,cam) == 0                         % Do again if the wrong camera view was used
                    warndlg(['That is not the image of the ',txt,' camera view!  You want ',endpts(view).cam],'Incorrect Camera View')
                    [u v] = ginput(num_pts);
                    cam = get(gca,'Tag');
                end
                
                plot(u,v,'bo','Tag',['arc_points_',num2str(view)])             % plot the inputed points to verify they're okay
                
                question = 'Are these points okay?';
                title = ['Get Arc Points:  ',txt,' view'];
                proceed = questdlg(question,title,'Yes','Redo','Cancel','Yes');
                
                switch proceed
                    case 'Yes'
                        redo = 0;
                        pts(view).u = u;
                        pts(view).v = v;
                        
                    case 'Redo'
                        redo == 1;
                        to_clear = findobj('Tag',['arc_points_',num2str(view)]);
                        delete(to_clear)
                        
                    case 'Cancel'
                        return
                end
                
            case 'Cancel'
                return
                
        end                                                                     % end switch to proceed with getting points
    end                                                                         % end while loop for redo
end                                                                             % end cycling through both views


% Now fit polynomial to each 2D arc
for view = 1:2
    
    % Make point array
    pt_array = [endpts(view).u(1) endpts(view).v(1);...
                   pts(view).u       pts(view).v;...
                endpts(view).u(2) endpts(view).v(2)];
        

        
    % Get polynomial parameters
    [ poly(view).u, poly(view).v ] = fit_poly( pt_array, poly_order );
    
    % Plot the fit
    if view == 1
        num_plot_pts = num_pts;     % This is the number of points we will step through to compare arcs
    elseif view == 2
        num_plot_pts = num_poly_pts;% This is large to get high density of points for identifying intersection
    end
    
    t = linspace( 0, 1, num_plot_pts )';
    
    fit_pts(view).u = polyval( poly(view).u, t );
    fit_pts(view).v = polyval( poly(view).v, t );
    
    target_axes = findobj('Tag',endpts(view).cam);
    axes(target_axes)
    plot( fit_pts(view).u, fit_pts(view).v, 'r' ,'Tag','arc_2d');

end

% If these fits are okay, get 3D points for the arc
proceed = questdlg('Proceed with 3D resolution of arc?','Resolve 3D','Yes','Cancel','Yes');
switch proceed
    
    case 'Yes'
        % project each point in fit_pts(1) into second view and find
        % intersection with second arc
        DLT_1 = data.cal.coeff.(['DLT_',endpts(1).cam(end)]);
        DLT_2 = data.cal.coeff.(['DLT_',endpts(2).cam(end)]);
        
        xres = get(dig_data.handles.(['im',endpts(2).cam(end)]),'XData') + 500; % get image size info for second image where line will be drawn
        yres = get(dig_data.handles.(['im',endpts(2).cam(end)]),'YData') + 500;

        for i = 1:length(fit_pts(1).u)
            
            [ line_u, line_v ] = im_pt_2_im_line( fit_pts(1).u(i), fit_pts(1).v(i), DLT_1, DLT_2, [ -500, xres, -500, yres ], num_poly_pts );
            
            axes(dig_data.handles.(endpts(1).cam))
            plot(fit_pts(1).u(i),fit_pts(1).v(i),'r*','Tag','int_pt')
            
            axes(dig_data.handles.(endpts(2).cam))
            to_clear = findobj('Tag','int_line');
            delete(to_clear)
            
            plot( line_u, line_v, 'b' ,'Tag','int_line');
            
            % Find line, arc intersection
            [int_x int_y] = polyxpoly(line_u,line_v, fit_pts(2).u,fit_pts(2).v);
            
            if isempty(int_x) == 0 | isempty(int_y) == 0
                if size(int_x,1) > 1
                    if size(int_x, 1) == 2 & int_x(1) == int_x(2)
                        int_x = int_x(1);
                        int_y = int_y(1);
                    else
                        % get user to pick which intersection to use
                        plot( int_x, int_y,'y*','Tag','mult_int_pt')
                        
                        questdlg('Multiple intersections found, click on the correct one','Choose Correct Intersection','OK','OK');
                        [pick_x pick_y] = ginput(1);
                        
                        for j = 1:length(int_x)
                            d(j) = pt_dist(pick_x,pick_y,int_x(j),int_y(j));
                        end
                        
                        [d_min ind] = min(d);
                        
                        int_x = int_x(ind);
                        int_y = int_y(ind);
                        
                        to_clear = findobj('Tag','mult_int_pt');
                        delete(to_clear)
                    end
                end
                plot( int_x, int_y,'r*','Tag','int_pt')
                
                % Get 3d coords
                x_cam_1 = fit_pts(1).u(i);
                y_cam_1 = fit_pts(1).v(i);
                x_cam_2 = int_x;
                y_cam_2 = int_y;
                A = [DLT_1, DLT_2];
                L = [x_cam_1, y_cam_1, x_cam_2, y_cam_2];
                
                [H] = reconfu(A,L);
                
                pts_3d = [pts_3d; H(1:3)]; % usually use store(H), first three are the x,y,z
                
                pause
            end
        end
        
        % Test 3d points by plotting on 3rd view
        % Clear all old points
        clear_tags = {'arc_points_1','arc_points_2','pt','int_pt','int_line','arc_2d'};
        for i = 1:length(clear_tags)
            to_clear = findobj('Tag',clear_tags{i});
            delete(to_clear)
        end
                
        % get 2d coords for each cam
        clear u_cam* v_cam*
        for i = 1:size(pts_3d,1)
            x = pts_3d(i,1);
            y = pts_3d(i,2);
            z = pts_3d(i,3);
            
            [ u_cam_1(i), v_cam_1(i) ] = dlt_3D_to_2D( data.cal.coeff.DLT_1, x, y, z );
            [ u_cam_2(i), v_cam_2(i) ] = dlt_3D_to_2D( data.cal.coeff.DLT_2, x, y, z );
            [ u_cam_3(i), v_cam_3(i) ] = dlt_3D_to_2D( data.cal.coeff.DLT_3, x, y, z );
        end
        
        % plot!
        axes(dig_data.handles.cam1)
        plot(u_cam_1,v_cam_1,'c','Tag','arc')
        
        axes(dig_data.handles.cam2)
        plot(u_cam_2,v_cam_2,'c','Tag','arc')

        axes(dig_data.handles.cam3)
        plot(u_cam_3,v_cam_3,'c','Tag','arc')
        
    case 'Cancel'
        to_clear = findobj('Tag','pt');
        to_clear = [to_clear, findobj('Tag',['arc_points_',num2str(view)])];
        delete(to_clear)
        return
        
end