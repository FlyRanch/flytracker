function IMout = renderflySYN(X,PAR,M,X1)

% [Y,IDX,IMout] = renderflySYN(X,PAR,M,X1)
% This function takes the 3D surface points of fly model, projects them
% according to camera M, and returns a synthetic image of the fly at the
% current camera view
%
% X - cell of N x 3 matrices of fly surface points
% X1 - cell of N x 3 matrices of fly wing thickness points
% PAR.DLT - camera projection
% M - M'th camera
%
% idx - cell of m x 1 vectors of indices for points in X; Y{m}(i,:) ==> X{m}(idx{m}(i),:)
% IMout - image of Rendered Fly

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
IMout = 255.*ones(PAR.imgres(1),PAR.imgres(2));
%%  Create the binary image
for j = length(X):-1:1
    
    
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
    
    tmp = zeros(size(IMout));
    %fit mesh points to image grid
    for k = 1:size(u{j},2)
        ttmp = roipoly(IMout,u{j}(:,k),v{j}(:,k));
        tmp(ttmp == 1) = 1;
    end
    tmp = logical(tmp);
    %tmp1 = roipoly(IMout,uT{j},vT{j});
    %Take the intersection between these two logical image masks.
    %tmp(tmp1 == 1) = 1;
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
        
        %Transpose
        uT1{j} = reshape(X3{j}(:,1),PAR.modsample1(j),[])';
        vT1{j} = reshape(X3{j}(:,2),PAR.modsample1(j),[])';

        %tmp1 = roipoly(IMout,u1{j},v1{j});
        %Take the intersection between these two logical image masks.
        for k = 1:size(u1{j},2)
            ttmp = roipoly(IMout,u1{j}(:,k),v1{j}(:,k));
            tmp(ttmp == 1) = 1;
        end
        %tmp(tmp1 == 1) = 1;
        
        %Transpose
%         uT1{j} = reshape(X3{j}(:,1),PAR.modsample1(j),[])';
%         vT1{j} = reshape(X3{j}(:,2),PAR.modsample1(j),[])';
%         
%         tmp1 = roipoly(IMout,uT1{j},vT1{j});
%         %Take the intersection between these two logical image masks.
%         tmp(tmp1 == 1) = 1;
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
%     if (j > 1 && perm < minperm)
%         %keyboard
%         utmp = u1{j};
%         vtmp = v1{j};
%         uTtmp = uT1{j};
%         vTtmp = vT1{j};
%         while ( (perm(end) < (2^supsample)*minperm) && (supsample < 4))
%             supsample = supsample + 1;
%             supscale = 2^supsample;
%             IM1 = zeros(supscale*PAR.imgres(1),supscale*PAR.imgres(2),'uint8');
%             utmp = 2.*utmp - 0.5;
%             vtmp = 2.*vtmp - 0.5;
%             uTtmp = 2.*uTtmp - 0.5;
%             vTtmp = 2.*vTtmp - 0.5;
%             tmp = roipoly(IM1,utmp,vtmp);
%             tmp1 = roipoly(IM1,uTtmp,vTtmp);
%             tmp(tmp1 == 1) = 1; 
%             %IM1 = bwmorph(IM1,'bridge','inf');
%             perm(end+1) = nnz(bwperim(tmp));
%             
%         end
%         % subsample the high resolution image to get back to the original
%         % resolution
%         tmp = subsample2(tmp,supsample);
%         tmp(tmp > 0) = 1;
%         tmp = logical(tmp);
%     end
        % 
    
    %% Add Gaussian Noise to make it work with EM algorithm in tracker 
    if j > 1 % the wings
        WingIntensity = 33 + 3.*randn(size(IMout)); %Model for wing intensity 
        IMout(tmp) = WingIntensity(tmp);
    else
        BodyIntensity = 15 + 1.*randn(size(IMout)); %Model for wing intensity 
        IMout(tmp) = BodyIntensity(tmp);
    end
end
%keyboard
% Fill in the holes
% IMout = bwmorph(IMout,'bridge','inf');
% IMout = imfill(IMout,'holes');
%IMout = bwmorph(IMout,'close',Inf);
