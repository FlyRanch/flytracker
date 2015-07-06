function Qnew = slerpMAT(Q,N)

%Q = slerp(Q,N) performs spherical interpolation between the quaternions
%defined in the matrix Q.  Each time interval between the quaternions is
%subdivided by N.
%
% Q = [4 M] matrix
% N is integer

M = size(Q,2);


B = acos(sum(Q(:,1:end-1).*Q(:,2:end),1));

%t = 0:1/N:1;
Tfull = 0:1/N:(M-1);

%Qnew = zeros(4,N*(M-1)+1);
Qnew = zeros(4,length(Tfull));

%We sample from this initial point
Qnew(:,1) = Q(:,1);

beg = 2;
for n = 1:M-1
    q1 = Q(:,n);
    q2 = Q(:,n+1);
    %Grab the scaling parameters that are in the current interval [n n+1]
    %of the matrix
    t = Tfull((n-1) < Tfull & Tfull <= n);
    
    %Rescale parameter to interval [0 1);
    t = t - (n-1);
%     Qnew(:,beg+(0:length(t)-1)) = ...
%
%   This is the first term in the spherical interpolation A.*AA
%         repmat(q1,1,length(t)) .* repmat([sin(B(n).*(1-t)) ./ sin(B(n))],4,1) +...
    A = [q1(1)*ones(1,length(t))
    q1(2)*ones(1,length(t))
    q1(3)*ones(1,length(t))
    q1(4)*ones(1,length(t))];
    
    AA = zeros(4,length(t));
    CC = zeros(4,length(t));
    
    a = [sin(B(n).*(1-t)) ./ sin(B(n))];
    for j=1:length(t)
        AA(:,j) = a(j)*ones(4,1);
    end
    
%   This is the second term C.*CC
%   repmat(q2,1,length(t)) .* repmat([sin(B(n).*t) ./ sin(B(n))],4,1);
    C = [q2(1)*ones(1,length(t))
    q2(2)*ones(1,length(t))
    q2(3)*ones(1,length(t))
    q2(4)*ones(1,length(t))];
    
    c = sin(B(n).*t) ./ sin(B(n));
    for j=1:length(t)
        CC(:,j) = c(j)*ones(4,1);
    end 
    
    %Computing the interpolated values this way is faster that using the
    %'repmat' function.
    Qnew(:,beg+(0:length(t)-1)) = A.*AA + C.*CC;
        
    beg = beg + length(t);
    
%     Qnew(:,(1:length(t))+(n-1)*(length(t)-1)) = ...
%     repmat(q1,1,length(t)) .* repmat([sin(B(n).*(1-t)) ./ sin(B(n))],4,1) +... 
% 	   repmat(q2,1,length(t)) .* repmat([sin(B(n).*t) ./ sin(B(n))],4,1);


end

