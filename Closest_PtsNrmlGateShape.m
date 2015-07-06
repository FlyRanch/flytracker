function [Features] = Closest_PtsNrmlGateShape(Features,p,i,PAR)

%Features = Closest_PtsNrmlGate(Features,p,i,PAR) finds the closest
%points in Features to the model evaluated at p. It renders the 3D model
%in each camera frame and then searches for the corresponding boundary
%points in the image data (obtained from Features).
%
%The search for corresponding image features in performed by examining the
%locations of high intensity along the normal vectors of the projected 
%model view in the current camera view
% See Blake and Isard "Active Contours", Chap. 5

skew = @(p) [0 -p(3) p(2);p(3) 0 -p(1);-p(2) p(1) 0];
numfeatpts = 30;

%==============================================
%Iterate over each fly
for k = 1:PAR.numfly
    state = p((k-1)*PAR.statedim+1:k*PAR.statedim);

    %Evaluate model at current state
    [x,y,z] = flymod_BodyShapeOnly(PAR.p0,state,PAR.params,PAR);

    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
        % reshape surface matrix structure to N x 3 matrix of all points
        % for ith part
        pts{j} = [reshape(x{j},[],1) reshape(y{j},[],1) reshape(z{j},[],1)];
    end

    for m = 1:PAR.numcam %Iterate over each camera
        %Return camera parameters
        R = PAR.cam(m).R;
        K = PAR.cam(m).K;
        camC = PAR.cam(m).C;
        % First, render the 3D model points to the current camera view
        % Return the corresponding outward normal vectors to perform search
        % along
        [pts2D,idx,IMout,Nrml] = renderflyMOD(pts,PAR,m);

        % Save the length of each cell
        for j = 1:length(pts2D)
            numpts(j) = size(pts2D{j},1);
        end
        %keyboard
        AllModelBndyPts = cell2mat(pts2D);
        AllNrmlVec = cell2mat(Nrml);
        
        %initialize
        Features(i).DataptsFull{m,k} = zeros(size(AllModelBndyPts));
        Features(i).DataRaysFull{m,k} = zeros(size(AllModelBndyPts,1),6);
        
        %now iterate over each model boundary point to search for 1-D image
        %features
        occluded_idx = [];
        for n = 1:size(AllModelBndyPts,1)
            xx = AllModelBndyPts(n,:);
            NN = AllNrmlVec(n,:);
            thresh = 80;
            C = [-0.25 -0.5 0 0.5 0.25];
            %C = fliplr([-0.25 -0.5 0 0.5 0.25]);

            %RR = 2*sqrt(NRy(n );
            RR = 15;
            datapts = find_features(xx,NN,linspace(-RR,RR,numfeatpts),C,...
                thresh,Features(i).IMbodyfull{m});
            Edgepts{n,m} = datapts;
            
            if isempty(datapts)
                occluded_idx(end+1) = n;
            else
                % determine closest point
                cp = kdtree(datapts,xx);
                
                %Store this corresponding point
                Features(i).DataptsFull{m,k}(n,:) = cp;
                
                %Now calculate the projection ray associated with this
                %corresponding data point
                % Add homogenous coordinate
                Allpts = [cp 1];

                % We need to represent these image coordinates in the fixed world frame

                % undo intrinsic camera parameters
                Allpts = Allpts*inv(K)';

                %this is now a vector pointing from camera origin to pixel location
                %Unrotate it to place this vector in the world frame
                % R'*x ==> x'*R
                Allpts = Allpts*R;

                % normalize this vector
                MAG = sqrt(sum(Allpts.^2,2));
                ProjVec = Allpts./repmat(MAG,1,3);

                % We will represent our projection rays in Pluker coordinates
                % For motivation see :
                % Rosenhahn, B.; Brox, T. & Weickert, J. Three-Dimensional Shape
                % Knowledge for Joint Image Segmentation and Pose Tracking IJCV, 2007,
                % 73, 243-262

                % The moment for each of these lines is the cross product of the vector
                % with the camera center;
                Moment = ( skew(camC)*ProjVec' )';
                
                % Store the Plüker coordinates of the line
                Features(i).DataRaysFull{m,k}(n,:) = [ProjVec Moment];
            end
        end
        
        %Store all the values
        Features(i).occluded_idx{m,k} = occluded_idx;
        
        %Get rid of [0 0] place holders for the occluded points
        obs_dumFull = Features(i).DataptsFull{m,k};
        obs_dumFull(occluded_idx,:) = [];
        Features(i).Datapts{m,k} = obs_dumFull;
        
        %Get rid of [0 0] place holders for the occluded rays
        obs_dumFull = Features(i).DataRaysFull{m,k};
        obs_dumFull(occluded_idx,:) = [];
        Features(i).DataRays{m,k} = obs_dumFull;
        
        Features(i).Edgepts{m,k} = Edgepts;

        AllIdxTo3DPts = idx;

        Features(i).Modelpts{m,k} = AllModelBndyPts;
        Features(i).IdxTo3DPts{m,k} = AllIdxTo3DPts;
    end %iterate over 'm', camera
end

