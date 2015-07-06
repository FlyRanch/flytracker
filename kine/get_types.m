% Function types = GET_TYPES(mfile_path)
% 
% CALLING FUNCTION: Kine_v2_0
% ACTIONS: Looks into the obj_types folder for which object types exist and
%          returns a cell array of the names of the relevant subfolders
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: May 6, 2004 by gwyneth

function types = get_types(mfile_path)


% Get list of names
if exist([mfile_path,filesep,'object_types']) ~= 0
    
    d = dir([mfile_path,filesep,'object_types']);
    
    for i = 1:length(d)
        types{i} = d(i).name;               % Put type names in cell array
    end
    
    types = types(3:end);                   % Get rid of '.' and '..'
    
else
    
    types = {};
    return
    
end


% Check each type folder to make sure has all necessary files
reqd_files = {'type_menu_callback',...
              'edit_object',...
              'edit_ok_callback',...
              'dig_setup',...
              'obj_program'};
    
for i = 1:length(types)
    for j = 1:length(reqd_files)
        
        ok(i,j) = exist([mfile_path,filesep,'object_types',filesep,types{i},filesep,reqd_files{j},'.m']);
        % Note, exist should return a value of 2 if the file exists as an mfile
        
    end
end

ok = sum(ok,2);

% Only include those types that got a score of 2 (exists) for each required file
include = find(ok==2*length(reqd_files));
types = types(include);
    