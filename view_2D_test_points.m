% Function to get the 2D test points from the calibration images

% File name strings chars A-S
file_names = char(65:83);

%number of cameras
Numcam = 3;
path = 'video/calibration/point_set/';

%Initialize point array
%L = zeros([length(file_names),2,Numcam]);

for i = 1:Numcam
    cam(i).pics = zeros(512,512,length(file_names));
    for j = 1:length(file_names)
        im = imread([path file_names(j) '/cam00' num2str(i) '/cal' lower(file_names(j)) '000001.bmp'],'bmp');
        cam(i).pics(:,:,j) = im;
    end
end
        