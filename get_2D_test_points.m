% Function to get the 2D test points from the calibration images

% File name strings chars A-S
file_names = char(65:83);

%number of cameras
Numcam = 3;
path = 'video/calibration/point_set/';

%Initialize point array
load 'all_calib_pts'
%L = zeros([length(file_names),2,Numcam]);

for i = 1:length(file_names)
    for j = 1:Numcam
        im = imread([path file_names(i) '/cam00' num2str(j) '/cal' lower(file_names(i)) '000001.bmp'],'bmp');
        
%         if j == 3             % Flip plotting axes to align with real world axes
%             im = flipud(im);
%         end

        imagesc(im);
        colormap gray; 
        axis image;
        
        title(['CAMERA ',num2str(j),', Point ',num2str(i),' of ',num2str(length(file_names))])
        xlabel('Click on point to enlarge, click again to record point location')
        
        ptapprox = ginput(1);
        xmin = round(ptapprox(1) - 50);
        xmax = round(ptapprox(1) + 50);
        ymin = round(ptapprox(2) - 50);
        ymax = round(ptapprox(2) + 50);
        axis([xmin xmax ymin ymax])
        
        L(i,:,j) = ginput(1);
        
        imagesc(im);
        colormap gray; 
        axis image;
        hold on;
        plot(L(i,1,j),L(i,2,j),'r*')
        pause
    end
end

save 'all_calib_pts' L F
        
L_1 = L(:,:,1);
L_2 = L(:,:,2);
L_3 = L(:,:,3);

save('kine/Calibrations for Kine/cal_for_20040823_done_20070718.mat','L_1','L_2','L_3','F');