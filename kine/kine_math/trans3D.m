% Function TRANS3D(model_in, vb_anch, vf_anch)
% 
% CALLING FUNCTION: obj_function (frame)
% ACTIONS: Translates model to correct place based on anchor point location
% PARENT PROGRAM: Kine_v3_0
% LAST MODIFIED: September 27, 2007 by gwyneth
%
% model_in is a 3xn vector of model points
% vb_anch and vf_anch should both be 3x1 vectors

function [model_out, v_move] = trans3D(model_in, vb_anch, vf_anch)

v_move = vf_anch - vb_anch;
mat_move = repmat(v_move,1,size(model_in,2));

model_out = model_in + mat_move;