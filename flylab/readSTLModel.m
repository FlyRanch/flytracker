% Matlab File to visualize the output file of trackCorners
%close(fignum)
clear all


datapath =  ['flylab/'];
filename = 'wing.stl';

fid = fopen([datapath filename],'r' );

%readsolid ascii
dum =  fscanf(fid,'%s %s',[1 2]);

k = 0;

%just read until the end of file when you get an error.
while 1 < 2
  k = k + 1;
  %read facet normal
  dum =  fscanf(fid,'%s %s',[1 2]);
  
  WING(k).normal = fscanf(fid,'%lg %lg',[1 3]);
  
  %read outer loop
  dum =  fscanf(fid,'%s %s',[1 2]);
  
  for j = 1:3
    %read vertex
    dum =  fscanf(fid,'%s %s',[1 1]);
    WING(k).triangle(j,:) = fscanf(fid,'%lg %lg',[1 3]);
  end
  
  %read endloop
  dum =  fscanf(fid,'%s %s',[1 1]);
  
  %read endfacet
  dum =  fscanf(fid,'%s %s',[1 1]);
end


    

