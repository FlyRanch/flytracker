function [Features,IMBW] = feat_detect(n,Features,PAR)
% n is the frame number
%  Calculate segmented images
IMBW = FlySegment(0,n,n,PAR.BG,PAR);
%Only take the images from the current frame
IMBW = squeeze(IMBW(:,:,n,:));

%% Iterate over each camera
for mm = 1:PAR.numcam
    data = IMBW(:,:,mm);    
   
    %load the raw grayscale image
    %Number of digit places for the total number of frames
    digits = length(num2str(PAR.numframes));
        
    input_filename = sprintf(['%s/cam%03d/%s%0' num2str(digits) 'd%s'], ...
        PAR.imagepath, mm, PAR.stub, n,PAR.image_filter(2:end));
    
    IMgray = imread(input_filename);
    flypix = IMgray(data == 1);
   
    %define initial conditions for EM algorithm
    init.W = [.5 .5];
    init.M = [mode(double(flypix)) 55];
    init.V(:,:,1) = 3;
    init.V(:,:,2) = 100;
    [W,M,V,L,Q] = EM_GM_fast(double(flypix),2,[],[],1,init);
    %Calculate the local minimum between the two gaussian peaks.
    xidx = find(min(M) < Q(:,1) & Q(:,1) < max(M));
    [vel,pidx] = min(Q(xidx,2));
    tmpx = Q(xidx,1);
    thresh = tmpx(pidx);

    IMbody = zeros(size(data));
    IMbody(min(flypix) <= IMgray & IMgray <= thresh) = 1;
    Features(n).IMbodyfull{mm} = IMbody.*255;
    
    %Convert to 256 grayscale image
    %data(logical(IMbody)) = 0;
    
    %Fill in stray pixels
    data = bwmorph(data,'bridge',inf);
    data = bwmorph(data,'fill',inf);
    %data = imfill(data,'holes');
    
    % Add artificial occlusion
    if ~isempty(PAR.OccludeShape{mm})
        keyboard
        fisp('feat detect')
        tmp = data;
        BW = roipoly(tmp,PAR.OccludeShape{mm}(:,1),PAR.OccludeShape{mm}(:,2));
        tmp(BW) = 0;
        data = tmp;
    end

    Features(n).IMfull{mm} = data.*255;
    
%     %% save images FTMmod 20120607
%     imwrite(Features(n).IMfull{mm},    [PAR.solutionpath 'Images_' PAR.stub '/flyimage' num2str(n) '_cam' num2str(mm) '_full.bmp']);
%     imwrite(Features(n).IMbodyfull{mm},[PAR.solutionpath 'Images_' PAR.stub '/flyimage' num2str(n) '_cam' num2str(mm) '_body.bmp']);

end % iterate over 'm', repeat for each camera