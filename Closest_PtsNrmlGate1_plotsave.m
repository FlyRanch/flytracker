    function Closest_PtsNrmlGate1_plotsave(Features,p,i,PAR)

skew = @(p) [0 -p(3) p(2);p(3) 0 -p(1);-p(2) p(1) 0];

%==============================================
%Iterate over each fly
for k = 1:PAR.numfly
    state = p((k-1)*PAR.statedim+1:k*PAR.statedim);
    
    %Evaluate model at current state
    [x,y,z] = flymodQ(state,PAR.params,PAR);

    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
        % reshape surface matrix structure to N x 3 matrix of all points
        % for ith part
        pts{j} = [reshape(x{j},[],1) reshape(y{j},[],1) reshape(z{j},[],1)];
    end
    
    for m = 1:PAR.numcam %Iterate over each camera
        % First, render the 3D model points to the current camera view
        % Return the corresponding outward normal vectors to perform search
        % along
        [pts2D,idx,IMout,Nrml] = renderflyMOD2(pts,PAR,m);    
        
        IMtest = (Features(i).IMfull{m} + Features(i).IMbodyfull{m});
        IMtest=round(IMtest./2);
        hold off; imagesc(IMtest); colormap gray; hold on
        
        if size(pts2D,1)==1
            % body markers
            plot(pts2D{1}(:,1),pts2D{1}(:,2),'b-')
        elseif size(pts2D,1)==2
            % body markers
            plot(pts2D{1}(:,1),pts2D{1}(:,2),'b-')
            % left wing markers
            plot(pts2D{2}(:,1),pts2D{2}(:,2),'r-')
        else
            % body markers
            plot(pts2D{1}(:,1),pts2D{1}(:,2),'b-')
            % left wing markers
            plot(pts2D{2}(:,1),pts2D{2}(:,2),'r-')
            % right wing markers
            plot(pts2D{3}(:,1),pts2D{3}(:,2),'g-')
        end
        
        % save image
%         saveas(gcf,[PAR.solutionpath 'Images_' PAR.stub '/flyimage' num2str(i) '_cam' num2str(m) '.jpg'])
        saveas(gcf,[PAR.solutionpath 'Images_' PAR.stub '/cam' num2str(m) '/flyimage_' num2str(i) '.jpg'])

    end %iterate over 'm', camera
end

