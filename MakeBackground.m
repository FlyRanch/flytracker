function Background = MakeBackground(dir_name, exp_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Program Name]
%  MakeBackground.m
%[Description]
%  Generate background images
%[Usage]
%  MakeBackground(result_1_1, exp002)
%[Programmer]
%  Atsushi Yamashita (Shizuoka University)
%[Development Environment]
%  Matlab R2007a
%[Version]
%  ver.1.0  06/08/2007 Initial Version (exp002)
%  ver.1.1	06/12/2007 exp002, exp035, exp83, exp098, exp101
%[Input]
%  dir_name     Input and output directory name
%  exp_name     Experiment name
%[Return]
%  Background   background image
%[Output]
%  Background images (jpeg files)
%[Comments]
%  This is my first matlab program ... 
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
% exp002
if exp_name == 'exp002'
    % Load images for background image generation
    % Camera 1 (Frame 1)
    frame = 1;
    input_filename = sprintf('./../0Data/FlyTracker/%s/cam001/%s%06d.bmp', ...
        exp_name, exp_name, frame);
    Input(:,:,1) = imread(input_filename);

    % Camera 2 (Frame 1)
    frame = 1;
    input_filename = sprintf('./../0Data/FlyTracker/%s/cam002/%s%06d.bmp', ...
        exp_name, exp_name, frame);
    Input(:,:,2) = imread(input_filename);

    % Camera 3 (Frame 5)
    frame = 5;
    input_filename = sprintf('./../0Data/FlyTracker/%s/cam003/%s%06d.bmp', ...
        exp_name, exp_name, frame);
    Input(:,:,3) = imread(input_filename);

    % Contrast enhancement
    for cam=1:3
        Image(:,:,cam) = imadjust(Input(:,:,cam));
    end

    % Vertices of polygon for interpolation (Camera 1, Frame 1)
    c(:,1) = [235 235 235 240 215 129 120 160 190 185 193 169 174 250 289 255 247];
    r(:,1) = [  1   2   3 190 212 197 210 241 257 289 306 318 328 315 229 189   1];

    % Vertices of polygon for interpolation (Camera 2, Frame 1)
    c(:,2) = [264 260 244 209 234 246 270 356 360 343 401 408 396 301 281 276 278];
    r(:,2) = [  1 230 235 254 295 295 330 370 330 290 260 229 220 245 244 233   1];

    % Vertices of polygon for interpolation (Camera 3, Frame 5)
    c(:,3) = [250 275 243 151 149 165 237 243 221 222 302 331 382 375 303 288 263];
    r(:,3) = [512 327 295 253 234 228 265 246 241 226 244 302 315 335 320 332 512];

    % Background image generation
    for cam=1:3
        Background(:,:,cam) = roifill(Image(:,:,cam), c(:,cam), r(:,cam));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exp035
elseif exp_name == 'exp035'
    % Load images for background image generation
    frame = 223;
    for cam=1:3
        input_filename = sprintf('./../0Data/FlyTracker/%s/cam%03d/%s%06d.bmp', ...
            exp_name, cam, exp_name, frame);
        Input(:,:,cam) = imread(input_filename);
    end

    % Contrast enhancement
    for cam=1:3
        Image(:,:,cam) = imadjust(Input(:,:,cam));
    end
    
    % Background image generation
    for cam=1:3
        Background(:,:,cam) = Image(:,:,cam);        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exp083
elseif exp_name == 'exp083'
    % Load images for background image generation
    % Camera 1 (Frame 1)
    frame = 1;
    input_filename = sprintf('./../0Data/FlyTracker/%s/cam001/%s%06d.bmp', ...
        exp_name, exp_name, frame);
    Input(:,:,1) = imread(input_filename);

    % Camera 2 (Frame 1)
    frame = 1;
    input_filename = sprintf('./../0Data/FlyTracker/%s/cam002/%s%06d.bmp', ...
        exp_name, exp_name, frame);
    Input(:,:,2) = imread(input_filename);

    % Camera 3 (Frame 1)
    frame = 1;
    input_filename = sprintf('./../0Data/FlyTracker/%s/cam003/%s%06d.bmp', ...
        exp_name, exp_name, frame);
    Input(:,:,3) = imread(input_filename);

    % Contrast enhancement
    for cam=1:3
        Image(:,:,cam) = imadjust(Input(:,:,cam));
    end

    % Vertices of polygon for interpolation (Camera 1, Frame 1)
    c(:,1) = [371 168 216 271 318 293 512 512 371];
    r(:,1) = [512 371 292 285 309 368 501 512 512];

    % Vertices of polygon for interpolation (Camera 2, Frame 1)
    c(:,2) = [  1   1   1   1 264 250 343 382 140];
    r(:,2) = [512 511 510 500 341 291 277 361 512];

    % Vertices of polygon for interpolation (Camera 3, Frame 1)
    c(:,3) = [  1   1   1   1   1 256 309 257  74];
    r(:,3) = [512 511 510 509 469 169 207 308 512];

    % Background image generation
    for cam=1:3
        Background(:,:,cam) = roifill(Image(:,:,cam), c(:,cam), r(:,cam));
    end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exp098
elseif exp_name == 'exp098'
    % Load images for background image generation
    frame = 569;
    for cam=1:3
        input_filename = sprintf('%s/%s/cam%03d/%s%06d.bmp', ...
            dir_name,exp_name, cam, exp_name, frame);
        Input(:,:,cam) = imread(input_filename);
        Background(:,:,cam) = imadjust(Input(:,:,cam));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif exp_name == 'exp101'
    % Load images for background image generation
    % Camera 1 (Frame 1)
    frame = 1;
    input_filename = sprintf('./../0Data/FlyTracker/%s/cam001/%s%06d.bmp', ...
        exp_name, exp_name, frame);
    Input(:,:,1) = imread(input_filename);

    % Camera 2 (Frame 1)
    frame = 1;
    input_filename = sprintf('./../0Data/FlyTracker/%s/cam002/%s%06d.bmp', ...
        exp_name, exp_name, frame);
    Input(:,:,2) = imread(input_filename);

    % Camera 3 (Frame 1)
    frame = 1;
    input_filename = sprintf('./../0Data/FlyTracker/%s/cam003/%s%06d.bmp', ...
        exp_name, exp_name, frame);
    Input(:,:,3) = imread(input_filename);
    
    % Contrast enhancement
    for cam=1:3
        Image(:,:,cam) = imadjust(Input(:,:,cam));
    end

    % Vertices of polygon for interpolation (Camera 1, Frame 1)
    c(:,1) = [300 168 250 311 291 450];
    r(:,1) = [512 420 332 374 398 512];

    % Vertices of polygon for interpolation (Camera 2, Frame 1)
    c(:,2) = [ 32 208 194 302 357 190];
    r(:,2) = [512 406 371 325 411 512];

    % Vertices of polygon for interpolation (Camera 3, Frame 1)
    c(:,3) = [  1   1 252 305 280 110];
    r(:,3) = [512 490 212 249 318 512];

    %{
    % Vertices of polygon for interpolation (Camera 1, Frame 213)
    c(:,1) = [160 294 296 258 242 201];
    r(:,1) = [266 244 281 302 376 365];

    % Vertices of polygon for interpolation (Camera 1, Frame 300)
    c(:,1) = [ 84 189 166 180 148  68];
    r(:,1) = [152 146 181 244 277 225];
    
    % Vertices of polygon for interpolation (Camera 2, Frame 300)
    c(:,2) = [187 261 368 350 338 280];
    r(:,2) = [218 147 146 182 221 275];

    % Vertices of polygon for interpolation (Camera 3, Frame 300)
    c(:,3) = [308 341 407 443 402 323];
    r(:,3) = [207 194 248 346 362 302];
    %}

    % Background image generation
    for cam=1:3
        Background(:,:,cam) = roifill(Image(:,:,cam), c(:,cam), r(:,cam));
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Incorrect input
else
    error('Incorrect Input');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save backgrounds as jpeg files 
for cam=1:3
    output_filename = sprintf('./%s/%s/BG/background%01d.jpg', ...
        dir_name, exp_name, cam);
    imwrite(Background(:,:,cam), output_filename, 'jpg');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
