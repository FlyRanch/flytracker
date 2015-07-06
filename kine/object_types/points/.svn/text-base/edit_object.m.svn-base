% Function EDIT_OBJECT
% 
% CALLING FUNCTION: Callback for Edit button in obj_fig for 'POINTS'
%                   objects
% ACTIONS: Opens edit object dialog then sets values of uicontrols
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 10, 2004 by gwyneth

function edit_object

controller = findobj('Tag','controller');
data = guidata(controller);                                 % Retrieve main guidata

addpath(data.setup.mfile_path)

% Set value of uicontrols based on currently selected object
obj_fig = findobj('Tag','obj_fig');
edit_fig = findobj('Tag','edit_fig');

obj_data = guidata(obj_fig);
edit_data = guidata(edit_fig);

obj_name = get(edit_data.handles.name_edit,'String');

% Set number of points
num_pts = obj_data.kine.(obj_name).config.num_pts;
% set(edit_data.handles.num_edit,'String',num2str(num_pts))

% Run num_edit_callback to display correct pt boxes
num_edit_callback   

% Set point names
if num_pts ~= 0
    
    h = findobj('Tag','pt_edit');
    for i = 1:num_pts                   % Get edit handles in correct order for points
        pos = get(h(i),'Position');
        point_order(i) = pos(2);
    end

    [sorted, ind] = sort(point_order);
    h(ind) = h;
    h = flipud(h);

    for i = 1:num_pts
        set(h(i),'String',obj_data.kine.(obj_name).config.points{i})
    end
end
        
guidata(controller, data)


    
    
