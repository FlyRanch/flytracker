% Function CHANGE_ORDER(fig_tag, fieldname, new_order)
% 
% CALLING FUNCTION: order_ok_callback, edit_ok_callback
% ACTIONS: Input ("new_order") is a matrix of fieldnames in the new order 
%          desired.  The ordering is done on the fields under 'fieldname' 
%          in the data matrix associated with the figure with tag 'fig_tag'.
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: February 21, 2004 by gwyneth

function change_order(fig_tag, fieldname, new_order)

fig_handle = findobj('Tag',fig_tag);
fig_data = guidata(fig_handle);

num_obj = length(new_order);

% Make new structure with fields in the new order
for i = 1:num_obj
    new_struct.(new_order{i}) = fig_data.(fieldname).(new_order{i});
end

% Make old fieldname equal to correctly ordered new structure
fig_data.(fieldname) = new_struct;

guidata(fig_handle,fig_data)