function [Features] = Closest_Pts(Features,p,i,PAR)

%Features = Closest_Pts(Features,p,i,PAR) finds the closest
%points in Features to the model evaluated at p. It renders the 3D model
%in each camera frame and then searches for the corresponding boundary
%points in the image data (obtained from Features).



%==============================================
%Iterate over each fly
for k = 1:PAR.numfly
    state = p((k-1)*PAR.statedim+1:k*PAR.statedim);

    %Evaluate model at current state
    [x,y,z] = flymod(state,PAR.params,PAR);

    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
        % reshape surface matrix structure to N x 3 matrix of all points
        % for ith part
        pts{j} = [reshape(x{j},[],1) reshape(y{j},[],1) reshape(z{j},[],1)];
    end

    for m = 1:PAR.numcam %Iterate over each camera
        %First, render the 3D model points to the current camera view
        [pts2D,idx] = renderflyMOD(pts,PAR,m);

        % Save the length of each cell
        for j = 1:length(pts2D)
            numpts(j) = size(pts2D{j},1);
        end

        % For each camera view, search for the closest points in image
        % boundary to each point in the model
        AllModelBndyPts = cell2mat(pts2D);
        %AllIdxTo3DPts = cell2mat(idx);
        AllIdxTo3DPts = idx;
        
        datapts = Features(i).All_Bndy_Points{m};
        datarays = Features(i).All_Proj_Rays{m};
        
        %these are the indices to points in datapts
        indx = kdtreeidx(datapts,AllModelBndyPts);
        
        Features(i).Datapts{m,k} = datapts(indx,:);
        Features(i).DataRays{m,k} = datarays(indx,:);
        
        Features(i).Modelpts{m,k} = AllModelBndyPts;
        Features(i).IdxTo3DPts{m,k} = AllIdxTo3DPts;
    end
end

