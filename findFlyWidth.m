function X = findFlyWidth(BB,p0,len,Y,PAR)


npts = length(p0);
c = 4;
knotv = oknot(npts,c,len,0);
ss = linspace(-len/2,len/2,200);
AA = dbasis(c,ss,npts,knotv);
B = repmat(AA,2,1);
A = -AA;
b = zeros(size(A,1),1);

options = optimset('Jacobian','off','Gradobj','off','TolX',1e-6,'TolFun',1e-8, ...
    'display','iter','largescale','on', ...
    'LevenbergMarquardt','on','outputFcn',@plotShape);
%keyboard


% [X,resnorm,residual,exitflag,output] = lsqnonlin(@modcurvefunc,p0,[],[],options);
[X] = fmincon(@modcurvefunc,p0,A(2:end-1,:),b(2:end-1,:),A([1 end],:),b([1 end],:),[],[],[],options);
    function [F,G] = modcurvefunc(p)
        F = [];
        [x,y,z,Frenet] = modcurveradiusP(BB,0,p,len,PAR);

        for k = 1:PAR.numworms

            mdlpts = [x(:,1) y(:,1)
                x(:,3) y(:,3)];
            
            Nrml = [-1.*Frenet.N
                Frenet.N];

            Nrml = Nrml(:,1:2); %2D vector

            
            %This is [k 2] matrix.  First transpose,then reshape(..,2*k,1)
            %to column vector with form
            %[X1model - X1data
            % Y1model - Y1data
            %    ..
            %    ..
            % Xnmodel - Xndata
            % Ynmodel - Yndata];

            Data = kdtree(Y,mdlpts);
            
            % point distance projected onto the Normal vector
            Error_Edge =  sum(Nrml.*(mdlpts - Data),2);

            F = sum(Error_Edge.^2);
        end
        if nargout > 1
            G = zeros(1,length(p));
            for k = 1:length(p)
                G(k) = sum(sum(Nrml.*(repmat(B(:,k),1,2).*Nrml),2));
            end
        end
    end %(end of modcurveRENfunc)
    function stop = plotShape(p,optimValues,state)
        if strcmp(state,'done')
            [x,y,z,Frenet] = modcurveradiusP(BB,0,p,len,PAR);
            
            ctrpts = [x(:,2) y(:,2)
                x(:,2) y(:,2)];
            
            mdlpts = [x(:,1) y(:,1)
                x(:,3) y(:,3)];
            Nrml = [-1.*Frenet.N
                Frenet.N];

            Nrml = Nrml(:,1:2); %2D vector

            Data = kdtree(Y,mdlpts);
            datarad = sum(Nrml.*(ctrpts - Data),2);
            datarad1 = datarad(1:length(ss),:)./PAR.pixpermm;
            datarad2 = -datarad(length(ss)+1:end,:)./PAR.pixpermm;
            figure; plot(ss,datarad1,'r--',ss,datarad2,'b-.',[ss fliplr(ss)],[AA*p;flipud(-AA*p)],'k:')
            
            knotv = oknot(npts,c,len,0);
            hold on; plot(knotv(c-2:end-(c-1)),p,'k*')
            stop = false;
        else
            stop = false;
        end
    end
        
end %(end of nested func)