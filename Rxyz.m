function R = Rxyz(theta)

%theta is a vector with xyz Euler angles ordered theta = [psi theta phi]

psi = theta(3);
the = theta(2);
phi = theta(1);

Cpsi = cos(psi);
Cthe = cos(the);
Cphi = cos(phi);

Spsi = sin(psi);
Sthe = sin(the);
Sphi = sin(phi);

R = [Cthe*Cpsi Cthe*Spsi -Sthe
    Sphi*Sthe*Cpsi-Cphi*Spsi Sphi*Sthe*Spsi+Cphi*Cpsi Cthe*Sphi
    Cphi*Sthe*Cpsi+Sphi*Spsi Cphi*Sthe*Spsi-Sphi*Cpsi Cthe*Cphi];


    