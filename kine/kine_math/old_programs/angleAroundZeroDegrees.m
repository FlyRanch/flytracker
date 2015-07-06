% Take an angle in degrees and put it between -180 and 180 degrees

function newAngle=angleAroundZeroDegrees(angle)

while angle > 180
    angle = angle - 360;
end

while angle <= -180
    angle = angle + 360;
end

newAngle = angle;