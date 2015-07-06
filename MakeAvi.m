function elapsed_time = MakeAvi(exp_name, scale, sFrame, nFrame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Program Name]
%  MakeAvi.m
%[Description]
%  Generate avi file from image sequences
%[Usage]
%  MakeAvi('exp035', 1/2, 17, 1775)
%[Programmer]
%  Atsushi Yamashita (Shizuoka University)
%[Development Environment]
%  Matlab R2007a
%[Version]
%  ver.1.0  06/08/2007 Initial Version (exp002)
%  ver.1.1	06/12/2007 exp002, exp035, exp83, exp098, exp101
%[Input]
%  exp_name     experiment number
%  scale        scale of avi file
%  sFrame       first frame number
%  nFrame       final frame number
%[Return]
%  elapsed_time computation time
%[Output]
%  Avi file
%[Comments]
%  This is my first matlab program ... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[variable Name]
%  start_time       Start time
%  elapsed_time     Elasped time
%  exp_name         Experiment name
%  dir_name         Input and output directory name
%  fileNameIn       Input file name
%  fileNameOut      Output file name
%  aviObj           Avi object
%  iFrame           Frame number
%  Image(i,j)       Input image
%  TempImage(i,j)   Temporal image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measure computation time (start)
start_time = cputime;

% Input and output directory name
dir_name = 'result_1_1';

% Display frame number
sprintf('Frame [%04d]-[%04d]', sFrame, nFrame)

% Output file name
if scale ~= 1
    fileNameOut = sprintf('./%s/%s/%s-extract-s.avi', dir_name, exp_name, exp_name);
else
    fileNameOut = sprintf('./%s/%s/%s-extract.avi', dir_name, exp_name, exp_name);
end

clear mex;

% Frame rate setting
frameRate = 30;

% Avi setting (Compression: Indeo5)
aviObj = avifile(fileNameOut, 'FPS', frameRate, 'COMPRESSION', 'Indeo5');

% Iteration from sFrame to nFrame
for iFrame = sFrame:nFrame
    % Display step number
    sprintf('STEP = %04d', iFrame-sFrame+1)    
    
    % Load present frame
    fileNameIn = sprintf('./%s/%s/cam-all/extract%04d.jpg', dir_name, exp_name, iFrame);
    Image = imread(fileNameIn);
    
    % Resize image
    if scale ~= 1
        TempImage = imresize(Image, scale);
        Image = TempImage;
    end

    % Add frame
    aviObj = addframe(aviObj, im2frame(Image));
end

% Close
aviObj = close(aviObj);

% Measure computation time (end)
elapsed_time = cputime - start_time;
