load WING_Struct.mat
clear TRI
allpts = zeros(3*length(WING),3);
for k = 1:length(WING)
  allpts(3*k-2:3*k,:) = WING(k).triangle;
end

%get rid of repetitions

allpts = unique(allpts,'rows');

for k = 1:length(WING)
  for j = 1:3
    TRI(k,j)  = find( (WING(k).triangle(j,1) == allpts(:,1) &...
                       WING(k).triangle(j,2) == allpts(:,2)&...
                       WING(k).triangle(j,3) == allpts(:,3) ));
  end
end

flywing.TRI = TRI;
flywing.pts = allpts;
save 'flywing.mat' flywing

figure; 
trisurf(TRI,allpts(:,1),allpts(:,2),allpts(:,3),'facecolor','g');
axis equal
