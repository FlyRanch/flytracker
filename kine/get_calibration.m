% Function success = GET_CALIBRATION
% 
% CALLING FUNCTION: Callback function for the calibrate button
% ACTIONS: either loads real-world and camera coordinates then performs
%          DLT calibration, OR runs get_test_points_g to get the real-world and
%          camera coordinates then performs DLT calibration
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function success = get_calibration

controller = findobj('Tag','controller');
data = guidata(controller);

% USER INPUT: load calibration points or get them
    question = 'Do you want to load previously saved calibration points or get them from camera data and real-world measurements now?';
    title = 'Calibration points input method';
    button1 = 'Load';
    button2 = 'Get now';
input_method = questdlg(question,title,button1,button2,button1);

switch input_method
    
    case 'Load'
        cd(data.setup.data_path)
        [cal_filename, cal_pathname ] = uigetfile;
        load([cal_pathname, cal_filename])        
        cd(data.setup.mfile_path)
        
    case 'Get now'
        [L_1, L_2, L_3, F] = eval(data.setup.cal_mfile);    % Use mfile defined in setup to get calibration points
        
end

if exist('F') == 1
    success = 1;
    % Write calibration values into guidata structure
    for i = 1:data.setup.cam_num
        
        Lmat = ['L_',num2str(i)];
        data.cal.points.(Lmat) = eval(Lmat);
        
    end
%     data.cal.points.L_1 = L_1;
%     data.cal.points.L_2 = L_2;
%     data.cal.points.L_3 = L_3;
    data.cal.points.F = F;
    
    data.cal.filename = cal_filename;
else
    success = 0;
    warndlg('Unable to find coordinates needed for calibration.  No DLT coefficients will be calculated.','No F matrix')
end

guidata(controller,data)   % Store guidata     