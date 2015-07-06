function plotquat(pf, s, cur_obj, cx1, cy1, cz1)

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
dig_data = guidata(dig_fig);

% Scale wing to measured length (may be different lengths each frame)
pf = s * pf;

% Translate wing model to right place
m_fit(1,1) = 0;
m_fit(2,1) = 0;
m_fit(3,1) = 0;

m_fit(1,2) = pf(2);
m_fit(2,2) = pf(3);
m_fit(3,2) = pf(4);

m_fit(1,:) = m_fit(1,:) + cx1;
m_fit(2,:) = m_fit(2,:) + cy1;
m_fit(3,:) = m_fit(3,:) + cz1;

% use calibration data to calculate what the wing model should look like in each 2D window
%%%DEBUG when add other calibration methods
x = m_fit(1,:);
y = m_fit(2,:);
z = m_fit(3,:);

for i = 1:data.setup.cam_num
    
    [u(i).pts v(i).pts] = dlt_3D_to_2D( data.cal.coeff.(['DLT_',num2str(i)]), x, y, z );
    axes(dig_data.handles.(['cam',num2str(i)]))
    cam_pt(i) = plot(u(i).pts, v(i).pts, 'Tag', cur_obj);
end

col_rgb = 'm';
set(cam_pt,'Color',col_rgb,'Visible',data.kine.(cur_obj).config.visible)