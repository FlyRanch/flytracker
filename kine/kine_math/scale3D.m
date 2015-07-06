% Function SCALE3D(model_in, s)
% 
% CALLING FUNCTION: obj_function (frame)
% ACTIONS: Scales the model by a factor of s
% PARENT PROGRAM: Kine_v3_0
% LAST MODIFIED: September 27, 2007 by gwyneth
%
% Assumes that model_in is still centered about the origin

function model_out = scale3D(model_in, s)

model_out = s*model_in;