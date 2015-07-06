% Function send_im_back
% 
% CALLING FUNCTION: draw_it 
% ACTIONS: Sends image to back of objects by setting order of children
%          handles for all three camera axes
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: August 16, 2004 by gwyneth

function send_im_back

for h = 1:3 %%DEBUG for diff # cameras
    ax = findobj('Tag',['cam',num2str(h)]);
    im = findobj('Tag',['im',num2str(h)]);
    
    v = get(ax,'Children');
    
    ind_image_handle = find( v == im );
    ind_other_handles = find( v ~= im);
    
    if isempty(ind_image_handle) == 0
        new_v = v(ind_other_handles);
        new_v = [new_v; im];
    else
        disp('Phooey, no image handle found!')
    end
    
    clear v new_v ind*
end
