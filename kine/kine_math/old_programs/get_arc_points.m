% Function GET_ARC_POINTS
% 
% CALLING FUNCTION: button down function for arc objects, cur point =
%                   'points'
% ACTIONS: gets user to click on a certain number of points along arc in
%          two images
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 5, 2004 by gwyneth (modified from get_arc by Will
%                Dickson)

function get_arc_points

controller = findobj('Tag','controller');
data = guidata(controller);

cur_obj = get_string(data.handles.cur_object_menu);

% num_pts = data.kine.(cur_obj).info.num_pts
num_pts = 3;  % CHANGE this when add info field

for view = 1:2   % do this twice for two camera views
    
    if view == 1, txt = 'first'; elseif view == 2, txt = 'second'; end
    
    redo = 0;
    while redo == 1
        
        question = ['Click on ',num_pts,' points along the arc (excluding the arc start and end points) in ',txt,'camera view'];
        title = ['Get Arc Points:  ',txt,' view'];
        proceed = questdlg(question,title,'OK','Cancel','OK');
        
        switch proceed
            
            case 'OK'
                [u,v] = ginput(num_pts);                                        % Ginput points for arc view 1
                plot(u,v,'bo','Tag',['arc_points_',num2str(view)],)             % plot the inputed points to verify they're okay
                
                % Fit poly to these points
                
                question = 'Are these points okay?';
                title = 'Get Arc Points:  First View';
                proceed = questdlg(question,title,'Yes','Redo','Cancel','Yes');
                
                switch proceed
                    case 'Yes'
                        redo = 0;
                        data.kine.(cur_obj).info.temp.(['u',num2str(view)]) = u;% save input points into temp structure if okay
                        data.kine.(cur_obj).info.temp.(['v',num2str(view)]) = v;
                        data.kine.(cur_obj).info.temp.(['cam',num2str(view)]) = get(gca,'Tag');
                        guidata(controller,data)
                        
                    case 'Redo'
                        redo == 1;
                        to_clear = findobj('Tag',['arc_points_',num2str(view)]);
                        delete(to_clear)
                        
                    case 'Cancel'
                        data.kine.(cur_obj).info.temp = [];
                        guidata(controller,data)
                        return
                end
                
            case 'Cancel'
                data.kine.(cur_obj).info.temp = [];
                guidata(controller,data)
                return
                
        end % end switch to proceed with getting points
        
    end % end while loop for redo
    
end % end cycling through both views

% check if arc info complete and resolve arc in 3d if so
get_arc_3d