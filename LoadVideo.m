function PAR = LoadVideo

%This function allows the user to select any image in the video sequence
%so that the tracker saves the video path and filename information.

% display('Locate the folder where you would like to save the tracker solution');
% PAR.solutionpath = uigetdir('','Locate the folder where you would like to save the tracker solution');
% %Add a forward slash
% PAR.solutionpath = [PAR.solutionpath '/'];

 PAR.solutionpath = ['solutions/'];

% [FileName,PathName] = uigetfile({'*.bmp';'*.tif'},'Click on any image in video sequence',PAR.solutionpath);
FileName = 'exp095000001.bmp';
PathName = '/home/florian/DATA/flighttracker/Takeoff/exp095/cam001/';

temp = imread([PathName FileName],FileName(end-2:end));

PAR.imgres = size(temp);
PAR.imgres_crop = PAR.imgres;

PAR.image_filter = ['*' FileName(end-3:end)];

dir_str = findstr(filesep,PathName);     
PAR.imagepath = PathName(1:dir_str(end-1));  % get path one folder up
% I assume that there are subfolders inside PAR.imagepath that begin with
% "cam*"

PAR.imagefilename = FileName;


filenameCells = makeFileList(PathName, PAR.image_filter );
PAR.stub = findRoot( filenameCells );
PAR.numframes = findFrames( filenameCells, PAR.stub );

%Clear the comand window
clc;

%--------------------------------------------------------------------------
function filenameCells = makeFileList( path, image_filter )

% Get list of all .bmp files in directory
namesStruct = dir( [path , image_filter] );    % struct array :)
namesStruct_length = length( namesStruct );                 % number of files in each folder

% Sort the list
filenameCells = {namesStruct(:).name};
filenameCells = sort( filenameCells );

%--------------------------------------------------------------------------
function file_root = findRoot( filenameCells )

% Find the root of the file name (without frame values appended) - the only
% assumption made is that the name is of the form [root, numbers].bmp
% Names should be in order, so compare first and last to get value range
% All letters 1:n-1 are the root, n:end are the frame numbers

root = 1;
n = 0;
while root == 1

    n = n + 1;
    root = strncmp(filenameCells{1}, filenameCells{end}, n);

end

file_root = filenameCells{1}(1:n-1);

%--------------------------------------------------------------------------
function frames = findFrames( filenameCells, file_root );

% HARDWIRED: length of appended file type is 4
ft_length = 4;
n = length(filenameCells{1}) - length(file_root) - ft_length;

% Assume last four string values are the suffix (usually .bmp)
first_frame = str2num(filenameCells{1}(end-(ft_length+n-1):end-ft_length));
last_frame = str2num(filenameCells{end}(end-(ft_length+n-1):end-ft_length));

frames = last_frame - first_frame + 1;