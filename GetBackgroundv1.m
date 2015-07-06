function Background = GetBackgroundv1(frame,PAR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Program Name]
%  GetBackgroundv1.m
%[Description]
%  Generate background images
%[Usage]
%  GetBackground( exp002,frame)
%[Programmer]
%  Ebraheem Fontaine
%[Development Environment]
%  Matlab R2007a
%[Version]
%  ver.1.0  06/08/2007 Initial Version (exp002)
%  ver.1.1	06/12/2007 exp002, exp035, exp83, exp098, exp101
%  ver.1.0  01/14/2008 Allow user to select the region to smooth over and
%  make the background
%[Input]
%  dir_name     Input and output directory name
%  exp_name     Experiment name
%[Return]
%  Background   background image
%[Output]
%  Background images (jpeg files)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[variable Name]
%  cam                      Camera number
%  frame                    Frame number
%  p_frame                  Frame number for iteration
%  input_filename           Input file name
%  output_filename          Output file name
%  Input(i,j,cam)           Original (input) image
%  Image(i,j,cam)           Image after contrast enhancement
%  Background(i,j,cam)      Background image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cam = 1:3

    % Load images for background image generation
    % Camera 1 (Frame 1)
    
    %Number of digit places for the total number of frames
    digits = length(num2str(PAR.numframes));
    
    input_filename = sprintf(['%s/cam%03d/%s%0' num2str(digits) 'd%s'], ...
        PAR.imagepath, cam, PAR.stub, frame,PAR.image_filter(2:end));
    
    Input = imread(input_filename);

    % Contrast enhancement
    Image = imadjust(Input);
    
    imagesc(Image); colormap gray
    
    Background(:,:,cam) = roifill;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
