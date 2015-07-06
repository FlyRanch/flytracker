% function [L_1, L_2, L_3, F] = get_calibration_points_g
%
% This is a routine for getting the calibration points under the current
% calibration system:
%     - Calibrate by moving pipette tip to set points defined by x,y,z in
%       get_test_points_g.m
%
%     - points are saved with the following path structure:
%       C:\Documents and Settings\gwyneth\My Documents\Research\20031202\calibration\point a\cam001\
%
%     - L_1, L_2, and L_3 are the ginputed points from cam001, cam002, and
%       cam003 respectively; the rows of the L matrices are points a to o

function [L_1, L_2, L_3, F] = get_calibration_points_g

alphabet = ['abcdefghijklmnopqrstuvwxyz'];
count = 0;

cal_fig = figure;

last = inputdlg('How many calibration points are there?', '', 1, {'15'});
last = str2num(last{1});

% Get pixel coordinates of calibration points from all three camera views
[ fileName, pName ] = uigetfile('*.bmp','Please click on the first calibration image file');

for cam = 1:3
    for point = 1:last
        
        count = count + 1;
        
        pName(end - 1) = num2str(cam);
        pName(end - 8) = alphabet(point);
        imag = imread(([pName,fileName]),'bmp');
        
        if cam == 1 | 2
            imag = imag;
        elseif cam == 3             % Flip plotting axes to align with real world axes
            imag = flipud(imag);
        end
        
        imagesc(imag,[min(min(imag)) max(max(imag))]);
        colormap(gray)
        title(['CAMERA ',num2str(cam),', Point ',num2str(count),' of ',num2str(last*3)])
        xlabel('Click on point to enlarge, click again to record point location')
        
        ptapprox = ginput(1);
            xmin = round(ptapprox(1) - 50);
            xmax = round(ptapprox(1) + 50);
            ymin = round(ptapprox(2) - 50);
            ymax = round(ptapprox(2) + 50);
        axis([xmin xmax ymin ymax])
        
        L{cam}(point,:) = ginput(1);
        
    end
end

close(cal_fig)
L_1 = L{1};
L_2 = L{2};
L_3 = L{3};

clear x*

% Get corresponding real world coordinates from user input
form = questdlg('Are the calibration points ordered in the standard format?', '', 'YES!', 'Alas, no','YES!');

switch form
    case 'YES!'
        
        prompts = {'X mid', 'X min', 'X max', 'Y mid', 'Y min', 'Y max', 'Z mid', 'Z min', 'Z max'};
        graph_title = 'Real World Calibration point values';
        
        F_inputs = inputdlg(prompts, graph_title, 1);
        
        xmid = str2num(F_inputs{1});
        xmin = str2num(F_inputs{2});
        xmax = str2num(F_inputs{3});
        ymid = str2num(F_inputs{4});
        ymin = str2num(F_inputs{5});
        ymax = str2num(F_inputs{6});
        zmid = str2num(F_inputs{7});
        zmin = str2num(F_inputs{8});
        zmax = str2num(F_inputs{9});
        
        F = [ xmid ymid zmid; ...
                xmax ymid zmid; ...
                xmin ymid zmid; ...
                xmid ymax zmid; ...
                xmid ymin zmid; ...
                xmid ymid zmax; ...
                xmid ymid zmin; ...
                xmin ymin zmin; ...
                xmax ymin zmin; ...
                xmax ymax zmin; ...
                xmin ymax zmin; ...
                xmin ymin zmax; ...
                xmax ymin zmax; ...
                xmax ymax zmax; ...
                xmin ymax zmax ];
                            
    case 'Alas, no'
          disp('\n Enter x, y, then z coordinate of each calibration point');
          
          for i = 1:last
              for j = 1:3
                  
                  prompts = {['Point ',alphabet(i),', ',alphabet(end - 3 + j),' =  ']};
                  graph_title = 'Real World Calibration point values';
                  
                  F_inputs = inputdlg(prompts, graph_title, 1);
                  
                  F(i,j) = str2num(F_inputs{1});
              end
          end
          
  end
  
  % SAVE
  answer = questdlg('Do you want to save the calibration points?','To save or not to save...');
  
  switch answer
      
      case 'Yes'
          clear f* F_inputs L a* c* g* i* l* p* x* y* z*
          uisave 
          
      case 'No'
          return
          
      case 'Cancel'
          return
          
  end
  

                


        

