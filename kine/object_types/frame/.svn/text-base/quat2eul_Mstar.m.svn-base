function [bank, elev, head] = quat2eul_Mstar(q)

e0 = q(1);
ex = q(2);
ey = q(3);
ez = q(4);

% Quaternion to Euler
if round( 10^9*(e0*ey-ex*ez))/10^9 == 0.5
    elev = pi/2;
    head = 0;
    bank = 2*asin(ex/cos(pi/4)) + head;
    
elseif round( 10^9*(e0*ey-ex*ez))/10^9 == -0.5
    elev = -pi/2;
    head = 0;
    bank = 2*asin(ex/cos(pi/4)) - head;
    
else
    % pi/2 rotation
    elev = asin( -e0^2+ex^2+ey^2-ez^2 );
    head = atan2( (ey*ez - ex*e0), (ex*ez + ey*e0) );
    bank = atan2( (ey*ez + ex*e0), (ey*e0 - ex*ez) );
% % - pi/2 rotation
%     elev = asin( e0^2-ex^2-ey^2+ez^2 );
%     head = atan2( (ey*ez - ex*e0), (ex*ez + ey*e0) );
%     bank = atan2( (ey*ez + ex*e0), (ex*ez - ey*e0) );
    
end

% bank = bank*180/pi;
% elev = elev*180/pi;
% head = head*180/pi;
% 
% eul = [bank; elev; head];