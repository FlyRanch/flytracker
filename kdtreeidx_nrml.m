function idx = kdtreeidx_nrml(R,M)
%idx = kdtreeidx_nrml(R,M)
%
%finds the indices in R that are closest to each point in M
%
%the points in R and M are normalized between [0 1] and then the kdtree is
%applied
%
minx = min(R,[],1);
maxx = max(R,[],1);
newR = (R - repmat(minx,size(R,1),1)) ./ repmat([maxx - minx],size(R,1),1);


minx = min(M,[],1);
maxx = max(M,[],1);
newM = (M - repmat(minx,size(M,1),1)) ./ repmat([maxx - minx],size(M,1),1);

idx = kdtreeidx(newR,newM);

