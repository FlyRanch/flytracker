function [ poly_1, poly_2 ] = fit_poly( pt_array, poly_order )
% FIT_POLY: fit polynomial of given order to pt_array

t = linspace( 0, 1, length(pt_array) )'; 

poly_1 = polyfit( t, pt_array(:,1), poly_order );

poly_2 = polyfit( t, pt_array(:,2), poly_order );

