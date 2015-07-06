% DLT test script - modified to work with gwyneth's data 20040129

% clear all
% 
% num_pts = 20;
% 
% line_param_1 = linspace( -40, 5, 100 );
% 
% line_param_2 = linspace( -40, 5, 100 );
% 
% line_param_3 = linspace( -5, 50, 100 );
% 
% [ L_1, L_2, L_3, F ] = get_test_pts;  ----> run get_test_points_g before this

% Changes to run this test with get_test_po* GW Pulford, ''Taxonomy of multiple target tracking methods'', IEE Radar Sonar Navig, 2005
** ([[:Image:Pulford2005_IEERadarSonarNavig.pdf]])ints_g
num_pts = 40;

line_param_1 = linspace( -50, 600, 100 ); %red

line_param_2 = linspace( -60, 500, 100 ); %green

%line_param_3 = linspace( 0, 100, 100 ); %yellow



% Get cube bndry
coords_min = min( F );

coords_max = max( F );

% Get some random points 
rnd_pts = rand( num_pts, 3 );

% Scale random pts into calibration data

a = repmat( ( coords_max - coords_min ), num_pts, 1 );

b = repmat( coords_min, num_pts, 1 );

rnd_pts = a.*rnd_pts + b;

figure%(1)

hold on

plot3( rnd_pts(:,1), rnd_pts(:,2), rnd_pts(:,3), '.b' ) 

% Camera 1 plot ------------------------------------------------


% Get DLT coeff 1
[ DLT_1, avgres ] = dltfu( F, L_1, [] );

% Get image points
[ u, v ] = dlt_3D_to_2D( DLT_1, rnd_pts(:,1), rnd_pts(:,2), rnd_pts(:,3) );

% Get the line params
for i = 1:num_pts

    [ params_x( i, 1:2 ), params_y( i, 1:2 ), params_z( i, 1:2 ) ] = dlt_2D_to_3D( DLT_1, u(i), v(i) );
    
    line(i).x = params_x(i,1) + params_x(i,2)*line_param_1;
    
    line(i).y = params_y(i,1) + params_y(i,2)*line_param_1;
    
    line(i).z = params_z(i,1) + params_z(i,2)*line_param_1;
    
    plot3( line(i).x, line(i).y, line(i).z, 'r' ); 
   
    %pause
    
end

%axis equal

% Camera 2 plot ------------------------------------------------

% Get DLT coeff 1
[ DLT_2, avgres ] = dltfu( F, L_2, [] );

% Get image points
[ u, v ] = dlt_3D_to_2D( DLT_2, rnd_pts(:,1), rnd_pts(:,2), rnd_pts(:,3) );

% Get the line params
for i = 1:num_pts

    [ params_x( i, 1:2 ), params_y( i, 1:2 ), params_z( i, 1:2 ) ] = dlt_2D_to_3D( DLT_2, u(i), v(i) );
    
    line(i).x = params_x(i,1) + params_x(i,2)*line_param_2;
    
    line(i).y = params_y(i,1) + params_y(i,2)*line_param_2;
    
    line(i).z = params_z(i,1) + params_z(i,2)*line_param_2;
    
    plot3( line(i).x, line(i).y, line(i).z, 'g' ); 
   
    %pause
    
end

axis equal

% % Camera 3 plot ------------------------------------------------
% 
% % Get DLT coeff 1
% [ DLT_3, avgres ] = dltfu( F, L_3, [] );
% 
% % Get image points
% [ u, v ] = dlt_3D_to_2D( DLT_3, rnd_pts(:,1), rnd_pts(:,2), rnd_pts(:,3) );
% 
% % Get the line params
% for i = 1:num_pts
% 
%     [ params_x( i, 1:2 ), params_y( i, 1:2 ), params_z( i, 1:2 ) ] = dlt_2D_to_3D( DLT_3, u(i), v(i) );
%     
%     line(i).x = params_x(i,1) + params_x(i,2)*line_param_3;
%     
%     line(i).y = params_y(i,1) + params_y(i,2)*line_param_3;
%     
%     line(i).z = params_z(i,1) + params_z(i,2)*line_param_3;
%     
%     plot3( line(i).x, line(i).y, line(i).z, 'y' ); 
%    
%     %pause
%     
% end

axis tight

view( 30, 45 )

