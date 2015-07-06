function [Y,idx,IMout,Nrml] = renderflyMOD(X,PAR,M,X1)

% [Y,IDX,IMout] = renderflyMOD(X,PAR)
% This function takes the 3D surface points of fly model, projects them
% according to camera M, then finds which points in 3D correspond to the 2D
% boundary points.
% X - cell of N x 3 matrices of fly surface points
% PAR.DLT - camera projection
% M - M'th camera
%
% Y - cell of M x 2 matrices of fly boundary points
% idx - cell of M x 1 vectors of indices for points in X; Y{m}(i,:) ==> X{m}(idx{m}(i),:)
% IMout - binary image of Rendered Fly

%amount to subsample the boundary to create less points than that induced
%by the camera resolution
subs = 4;

% Initialize projected planar points
u = cell(length(X),1);
v = cell(length(X),1);
uT = cell(length(X),1);
vT = cell(length(X),1);
Y = cell(length(X),1);
Nrml = cell(length(X),1);
idx = cell(length(X),1);
X2 = cell(length(X),1);
if nargin > 3
    X3 = cell(length(X1),1);
    u1 = cell(length(X1),1);
    v1 = cell(length(X1),1);
    uT1 = cell(length(X1),1);
    vT1 = cell(length(X1),1);
end
IMout = zeros(PAR.imgres(1),PAR.imgres(2));
%%  Create the binary image
for j = 1:length(X)
    
    
    % This is ALL the model points, not just the boundary points
    %X2{j} = dlt_3D_to_2D(PAR.DLT(:,M),X{j});
    X2{j} = PAR.cam(M).P*[X{j} ones(size(X{j},1),1)]';
    
    X2{j} = X2{j}';
    %homogenize
    X2{j} = X2{j}./repmat(X2{j}(:,3),1,3);
    X2{j} = X2{j}(:,1:2);
    
    numpts(j) = size(X2{j},1);
    
    % These arrays are used to create to rendered image using 'roipoly'
    u{j} = reshape(X2{j}(:,1),PAR.modsample(j),[]);
    v{j} = reshape(X2{j}(:,2),PAR.modsample(j),[]);
   
    %Transpose
    uT{j} = reshape(X2{j}(:,1),PAR.modsample(j),[])';
    vT{j} = reshape(X2{j}(:,2),PAR.modsample(j),[])';
        
    %fit mesh points to image grid
    tmp = roipoly(IMout,u{j},v{j});
    tmp1 = roipoly(IMout,uT{j},vT{j});
    %Take the intersection between these two logical image masks.
    tmp(tmp1 == 1) = 1;
    %keyboard
    % For the wings, also render the model points that correspond to their
    % thickness.
    if ((j > 1) && (nargin > 3))
        X3{j} = PAR.cam(M).P*[X1{j} ones(size(X1{j},1),1)]';

        X3{j} = X3{j}';
        %homogenize
        X3{j} = X3{j}./repmat(X3{j}(:,3),1,3);
        X3{j} = X3{j}(:,1:2);

        % These arrays are used to create to rendered image using 'roipoly'
        u1{j} = reshape(X3{j}(:,1),PAR.modsample1(j),[]);
        v1{j} = reshape(X3{j}(:,2),PAR.modsample1(j),[]);
        
        tmp1 = roipoly(IMout,u1{j},v1{j});
        %Take the intersection between these two logical image masks.
        tmp(tmp1 == 1) = 1;
        
        %Transpose
        uT1{j} = reshape(X3{j}(:,1),PAR.modsample1(j),[])';
        vT1{j} = reshape(X3{j}(:,2),PAR.modsample1(j),[])';
        
        tmp1 = roipoly(IMout,uT1{j},vT1{j});
        %Take the intersection between these two logical image masks.
        tmp(tmp1 == 1) = 1;
    end
    
    % Fix for the wings if they don't render correctly because leading edge
    % directly faces the camera
    %
    % the perimeter of the leading edge of a wing at 512x512 resolution is
    % about 115 pixels.  Chord length is ~50 pixels
    %
    % If the euler number is greater than 1 for a wing, then it didn't
    % render properly
    % Let's connect the stray pixels by rendering at the next finer of
    % resolution
    minperm = 40;
    perm = nnz(bwperim(tmp));
    supsample = 0;
    if (j > 1 && perm < minperm && ~isempty(X{j}))
        %keyboard
        utmp = u1{j};
        vtmp = v1{j};
        uTtmp = uT1{j};
        vTtmp = vT1{j};
        while ( (perm(end) < (2^supsample)*minperm) && (supsample < 4))
            supsample = supsample + 1;
            supscale = 2^supsample;
            IM1 = zeros(supscale*PAR.imgres(1),supscale*PAR.imgres(2),'uint8');
            utmp = 2.*utmp - 0.5;
            vtmp = 2.*vtmp - 0.5;
            uTtmp = 2.*uTtmp - 0.5;
            vTtmp = 2.*vTtmp - 0.5;
            tmp = roipoly(IM1,utmp,vtmp);
            tmp1 = roipoly(IM1,uTtmp,vTtmp);
            tmp(tmp1 == 1) = 1; 
            %IM1 = bwmorph(IM1,'bridge','inf');
            perm(end+1) = nnz(bwperim(tmp));
            
        end
        % subsample the high resolution image to get back to the original
        % resolution
        tmp = subsample2(tmp,supsample);
        tmp(tmp > 0) = 1;
        tmp = logical(tmp);
    end
        % 
        
    IMout(tmp) = 1;
    tmp = roipoly(IMout,uT{j},vT{j});
    IMout(tmp) = 1;
end
%keyboard
% Fill in the holes
IMout = bwmorph(IMout,'bridge','inf');
IMout = imfill(IMout,'holes');
%IMout = bwmorph(IMout,'close',Inf);
%% Extract image boundaries
B = bwboundaries(IMout);

% determine which boundary is which
% Pick the boundaries that are not the perimeter of the image.
borderidx = [];
for j = 1:length(B)
    if size(B{j},1) == 2*sum(PAR.imgres - 1)+1
        borderidx = j;
    end
end

% Ignore them for now
% %Find enclosed boundaries
% enclosed_boundaries = [];
% for j = 1:length(B)
%   enclosed_boundaries = [enclosed_boundaries ; find(A(:,j))];
% end

bndy2rid = [];

for k = 1:length(B)
    %I will throw away contours that are smaller than a particular length

    if (size(B{k},1) < round(50) )
        %This boundary is too small, probably noise
        bndy2rid(end+1) = k;
    elseif (all(diff(B{k}(:,1)) == 0) || ...
            all(diff(B{k}(:,2)) == 0))
        %This boundary is the image border
        bndy2rid(end+1) = k;
    elseif (size(B{k},1) == PAR.imgres(1) || size(B{k},1) >= PAR.imgres(2))
        %This boundary is also probably an image border
        bndy2rid(end+1) = k;
    end
    %subsample the boundary
    B{k} = B{k}(1:subs:end,:);
end

idx2keep = setdiff(1:length(B),[borderidx bndy2rid]);

B = B(idx2keep);
nrmltemp = cell(size(B));
%Smooth the boundaries
for m = 1:length(B)
    %Switch to because 'rows' correspond to 'y' pts
    %Origin is in the lower left of the image plane
    B{m} = [B{m}(:,2) B{m}(:,1)];
    
    boundary = B{m};
    extendedboundary = [];
    smoothboundary = [];

    %order of butterworth filter
    filter_order = 4;
    %filter parameters definieeren
    %keyboard
    [b a] = butter( filter_order, .95);

    %Mirror the boundary for filtering
    extendedboundary(1:length(boundary),:)=boundary(1:length(boundary),:);
    extendedboundary(length(boundary)+1:2*length(boundary),:)=boundary(1:length(boundary),:);
    extendedboundary(2*length(boundary)+1:3*length(boundary),:)=boundary(1:length(boundary),:);

    LONGsmoothboundary=filtfilt( b, a, extendedboundary);
    smoothboundary(1:length(boundary),:)=LONGsmoothboundary(length(boundary)+1:2*length(boundary),:);
    B{m} = smoothboundary;
    
    %% Calculate the outward normal vectors
    if nargout > 3
        % Fit a Closed spline to the boundary
        Ytemp = B{m};
        dd = reshape(Ytemp,[],1);

        %-- the number of control points be about one sixth of the total
        %data points
        c = 4; %cubic spline
        numrep = c-1; %this is the number of control points that are repeated
        %at the beginning and end of the vector;
        npts = max(c+numrep,round(size(Ytemp,1)/2));
        T = linspace(-0.5,0.5,size(Ytemp,1));
        knotv = oknotP(npts,c,1,0);
        [N,D1] = dbasisP(c,T,npts,knotv);

        Aeq = zeros(2*numrep,2*npts);
        for i = 1:2
            Aeq((i-1)*numrep+(1:numrep),(i-1)*npts+(1:numrep)) = eye(numrep);
            Aeq((i-1)*numrep+(1:numrep),(i-1)*npts+(npts-(numrep-1):npts)) = -1*eye(numrep);
        end
        options = optimset('largescale','off','display','off');
        bb = lsqlin(blkdiag(N,N),dd,[],[],Aeq,zeros(2*numrep,1),[],[],[],options);

        BB = reshape(bb,[],2);

        B{m} = N*BB;

        skew_z = [0 -1 0
            1 0 0
            0 0 0];
        % --- Perform cross product with skew matrix and make unit length
        Tang = D1*BB;

        Tang = [Tang zeros(size(Tang,1),1)];
        MAG = sqrt(sum(Tang.^2,2));
        Tang = Tang./repmat(MAG,1,3);

        nrmltemp{m,:} = (skew_z*Tang')';
        %Change sign so that it is outward
        nrmltemp{m,:} = -nrmltemp{m,:}(:,1:2);
        
        % Ignore this for now
        %Flip normal direction if fully enclosed
%         if any(m == enclosed_boundaries)
%             nrmltemp{m,:} = -1.*nrmltemp{m,:};
%         end
    end
end
Bpts = cell2mat(B);
Nvec = cell2mat(nrmltemp);



%% Now, determine the correspondence

% Make sure it's a column vector
X2 = X2(:);

%Concatenate the matrices of 2D points corresponding to each body part on
%top of each other.
ALL2DPts = cell2mat(X2);

%Search for closest points in 'ALL2DPts' to each point in Bpts
ALLidx = kdtreeidx(ALL2DPts,Bpts);

%Now transform this index back to one that makes sense for each cell
for j = 1:length(X2)
    % Get the indices to ALL2DPts that correspond to the beginning and end 
    % of part 'j'
    
    if j == 1
        beg = 1;
    else
        beg = sum(numpts(1:j-1)) + 1;
    end
    
    endd = beg + numpts(j) - 1;
    
    Y{j} = Bpts((ALLidx >= beg & ALLidx <= endd),:);
    
    if nargout > 3
        Nrml{j} = Nvec((ALLidx >= beg & ALLidx <= endd),:);
    end
    
    idx{j} = ALLidx((ALLidx >= beg & ALLidx <= endd),:);
    
    %Transform this index back to one that correseponds to the index within
    %the jth cell, (for this specific model part);
    idx{j} = idx{j} - beg + 1;
end

