% Function ORDER_OK_CALLBACK
% 
% CALLING FUNCTION: Callback for OK button in order_fig
% ACTIONS: Reads edit boxes for new order of objects, calls change_order
%          function to reorder the fields in obj_data.kine
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 21, 2004 by gwyneth

function order_ok_callback

order_data = guidata(gcf);

old_order = order_data.objects;
num_obj = length(order_data.objects);

new_index = 1;
obj_name = 2;

% Get user input new index number
for i = 1:num_obj
    
    order_vector(i) = str2num(get(order_data.handles.(['obj',num2str(i),'_edit']),'String'));
    
end

% Error checking
for i = 1:num_obj
    
    % Check for two of the same index
    check = find(order_vector == i);
    if length(check) > 1
        warndlg('Two objects cannot have the same order number','That''s not gonna work!')
        return
    end
    
    % Check for an index larger than the number of objects
    if order_vector(i) > num_obj
        warndlg('You have entered an order number larger than the number of objects','That''s not gonna work!')
        return
    end
    
end

% Reorder the fieldnames matrix
for i = 1:num_obj
    new_order(i) = old_order(order_vector(i));
end

% Reorder obj_data.kine fields
change_order('obj_fig','kine',new_order)

% Make the current configuration 'custom'
add_custom

% Close the Change Order dialog
delete(gcf)