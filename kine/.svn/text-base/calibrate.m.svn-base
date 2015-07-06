% Function CALIBRATE
% 
% CALLING FUNCTION: Callback function for the calibrate button
% ACTIONS: Generates calibration coefficients from values in data.cal.points 
%          and the type of calibration specified in data.cal.type 
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function calibrate

controller = findobj('Tag','controller');
data = guidata(controller);

% % Get calibration points from data structure
% L_1 = data.cal.points.L_1;
% L_2 = data.cal.points.L_2;
% L_3 = data.cal.points.L_3;
% F = data.cal.points.F;
% 
% % Right now only one kind of calibration, so:
% data.setup.cal_type = 'DLT';
switch data.setup.cal_type

    case 'DLT'
        
        % FIND DLT COEFFIECIENTS
        for i = 1:data.setup.cam_num
            
            Lmat = ['L_',num2str(i)];
            DLTmat = ['DLT_',num2str(i)];
            [ data.cal.coeff.(DLTmat), avgres ] = dltfu( data.cal.points.F, data.cal.points.(Lmat), [] );
            
        end
            
%         [ data.cal.coeff.DLT_1, avgres ] = dltfu( F, L_1, [] );
%         [ data.cal.coeff.DLT_2, avgres ] = dltfu( F, L_2, [] );
%         [ data.cal.coeff.DLT_3, avgres ] = dltfu( F, L_3, [] );
        
end

guidata(controller,data)