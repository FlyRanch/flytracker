% Function LOAD_CALIBRATION
% 
% CALLING FUNCTION: Callback function for the calibrate button
% ACTIONS: runs scripts that load or perform the calibration 
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function load_calibration

success = get_calibration;
if success == 1
    calibrate
    update_controller
end


