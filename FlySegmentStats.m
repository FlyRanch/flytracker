function [IMBW,flypix] = FlySegmentStats(type, frame1, frame2,PAR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Program Name]
%  FlySegmentStats.m
%[Description]
%  Extract Fly Body and Wings from Image Sequences, the look at the
%  distribution of pixel intensities that belong to the fly
%[Usage]
%  exp002: FlyExtract(0, 1, 1775)
%  exp035: FlyExtract(0, 1, 223)
%  exp083: FlyExtract(0, 1, 640)
%  exp098: FlyExtract(0, 1, 569)
%  exp101: FlyExtract(0, 1, 494)
%[Programmer]
%  Atsushi Yamashita (Shizuoka University)
%  Ebraheem Fontaine
%[Development Environment]
%  Matlab R2007a
%[Version]
%  ver.1.0  06/08/2007 Initial Version (exp002)
%  ver.1.1	06/12/2007 exp002, exp035, exp83, exp098, exp101
%  ver.1.2  11/14/2007 save binary images to disk and load them later
%  during tracker. (Ebraheem)
%[Input]
%  type     =0: All cameras (Camera 1-3)
%           =1: Camera 1
%           =2: Camera 2
%           =3: Camera 3
%  frame1   first file number
%  frame2   final file number
%  PAR      structure containing various parameters of the tracking
%           sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[variable Name]
%  cam                      Camera number
%  frame                    Frame number
%  p_frame                  Frame number for iteration
%  exp_name                 Experiment name
%  input_filename           Input file name
%  dir_name                 Input and output directory name
%  output_dir_name          Output directory name
%  output_filename          Output file name
%  level(cam)               Threshold values for binarize
%  Input(i,j,cam)           Original (input) image
%  Image(i,j,cam)           Image after contrast enhancement
%  Background(i,j,cam)      Background image
%  Subtract(i,j,cam)        Subtract image
%  TempImage(i,j,cam)       Temporal image
%  FlyAll(i,j,cam)          Fly body and wings
%  FlyBody(i,j,cam)         Fly body
%  FlyWing(i,j,cam)         Fly wings
%  Final1(i,j,color,cam)    Image for generating final result 1
%  Final2(i,j,color,cam)    Image for generating final result 2
%  Composite(i,j,color,cam) Composite image for final result
%  PrevImage(i,j,cam,num)   Previous image
%  flypix                   [#frames #cam] cell containing a vector of  fly
%                           intensity values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Setting
% Measure computation time (start)
tic;

% Input and output directory name (Version name)
dir_name = PAR.videopath;
exp_name = PAR.filetag;

% Experiment number
% exp_name = 'exp002';
% exp_name = 'exp035';
% exp_name = 'exp083';
% exp_name = 'exp098';
% exp_name = 'exp101';

% Threshold value for binarize
LEVEL = 0.15;

% Logical multiply number
if exp_name == 'exp002'
    and_number = 16;
elseif exp_name == 'exp035'
    and_number = 0;
elseif exp_name == 'exp083'
    and_number = 0;
elseif exp_name == 'exp098'
    and_number = 0;
elseif exp_name == 'exp101'
    and_number = 0;
else
    error('Incorrect input ... exp_name %s', exp_name);
end

% Initialize
PrevImage = zeros(512,512,and_number,3);
IMBW = zeros(512,512,length(1:frame2),3,'uint8');
TempImage = zeros(512,512,3,'uint8');
Input = zeros(512,512,3,'uint8');
Image = zeros(512,512,3,'uint8');
Subtract = zeros(512,512,3,'uint8');
FlyAll = zeros(512,512,3,'uint8');
FlyBody = zeros(512,512,3,'uint8');
FlyWing = zeros(512,512,3,'uint8');
flypix = cell(length(1:frame2),3);

% All cameras
if type == 0
    cam1 = 1;
    cam2 = 3;
    sprintf('Camera 1-3 ... Generating background images');
    
% Camera 1
elseif type == 1
    cam1 = 1;
    cam2 = 1;
    sprintf('Camera 1 ... Generating background image');
    mkdir(output_dir_name, 'cam1');

% Camera 2
elseif type == 2
    cam1 = 2;
    cam2 = 2;
    sprintf('Camera 2 ... Generating background image');
    mkdir(output_dir_name, 'cam2');

% Camera 3
elseif type == 3
    cam1 = 3;
    cam2 = 3;
    sprintf('Camera 3 ... Generating background image');
    mkdir(output_dir_name, 'cam3');

% Incorrect Input
else
    error('Incorrect Input');
end

% Background image generation
Background = GetBackground(exp_name);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration from frame1 to frame2
for frame=frame1:frame2
    % Load images
    for cam=cam1:cam2
        input_filename = sprintf('%s/%s/cam%03d/%s%06d.bmp', ...
            PAR.videopath,exp_name, cam, exp_name, frame);
        Input(:,:,cam) = imread(input_filename);
    end

    % Contrast enhancement
    for cam=cam1:cam2
        Image(:,:,cam) = imadjust(Input(:,:,cam));
    end
    
    % Image substraction (background removal)
    for cam=cam1:cam2
        Subtract(:,:,cam) = imsubtract(Background(:,:,cam), Image(:,:,cam));
    end

    % Contrast enhancement
    for cam=cam1:cam2
        TempImage(:,:,cam) = imadjust(Subtract(:,:,cam));
        Subtract(:,:,cam) = TempImage(:,:,cam);
    end

    % Convert images to binary images
    for cam=cam1:cam2
        level(cam) = graythresh(Subtract(:,:,cam));
        FlyAll(:,:,cam) = im2bw(Subtract(:,:,cam), .250);
        %FlyAll(:,:,cam) = im2bw(Subtract(:,:,cam), level(cam));
        % FlyAll(:,:,cam) = im2bw(Subtract(:,:,cam), LEVEL);
    end

    % Noise removal (Median filtering)
    for cam=cam1:cam2
        TempImage(:,:,cam) = medfilt2(FlyAll(:,:,cam), [3 3]);
        FlyAll(:,:,cam) = TempImage(:,:,cam);
    end
    
%     % Save each image to disk
%     for cam=cam1:cam2
%         %make new directory if it doesn't exist
%         dirname = sprintf('%s/%s/cam%03d_BG/',dir_name,exp_name,cam);
%         if ~exist(dirname,'dir')
%             mkdir(dirname);
%         end
%         
%         outputimage = sprintf('%s/%s%06d.bmp',dirname, exp_name, frame);
%         imwrite(FlyAll(:,:,cam),outputimage,'bmp');
%     end
    % Noise removal
    if strcmp(exp_name,'exp002')
        % Bottom platform in camera 2 images
        FlyAll(450:512, 1:512, 2) = 0;
    elseif strcmp(exp_name,'exp035')
        % Glass in camera 2 image
        FlyAll(1:512, 365:512, 2) = 0;
        % Noise in camera 3 image
        if frame > 120
            FlyAll(1:512, 300:512, 3) = 0;
        end
    elseif strcmp(exp_name,'exp101')
        % Bottom platform in camera 2 images
        % FlyAll(440:512, 1:512, 2) = 0;
    end
        
    % Fly body and wing separation
    for cam=cam1:cam2
        FlyBody(:,:,cam) = FlyAll(:,:,cam);
        for p_frame=1:and_number
            TempImage(:,:,cam) = FlyBody(:,:,cam) & PrevImage(:,:,p_frame,cam);
            FlyBody(:,:,cam) = TempImage(:,:,cam);
        end
        FlyWing(:,:,cam) = imsubtract(FlyBody(:,:,cam), FlyAll(:,:,cam));
    end

    % Noise removal (Median filtering)
    for cam=cam1:cam2
        TempImage(:,:,cam) = medfilt2(~FlyWing(:,:,cam), [3 3]);
        FlyWing(:,:,cam) = ~TempImage(:,:,cam);
    end
    
    
    % Generation of final results
    for cam=cam1:cam2
        tmp =uint8(FlyBody(:,:,cam) | logical(FlyWing(:,:,cam)));
        IMBW(:,:,frame,cam) = tmp;
        imgray = Input(:,:,cam);
        flypix{frame,cam} = imgray(tmp == 1);
    end
    
    % Reserve results for next frame processing
    for cam=cam1:cam2
        for p_frame =and_number:-1:2
        PrevImage(:,:,p_frame,cam) = PrevImage(:,:,p_frame-1,cam);
        end
        PrevImage(:,:,1,cam) = FlyAll(:,:,cam);
    end
end

% Measure computation time (end)
% toc


