function dups = checkdups(x)

% CHECKDUPS  Checks for the presence of duplicate entries in an array
%
%     dups = CHECKDUPS(x)
%
%     dups = number_of_duplicates if there are duplicate entries in x
%     otherwise dups = 0

%   Copyright  (c) Rudolph van der Merwe (2002)
%
%   This file is part of the ReBEL Toolkit. The ReBEL Toolkit is available free for
%   academic use only (see included license file) and can be obtained by contacting
%   rvdmerwe@ece.ogi.edu.  Businesses wishing to obtain a copy of the software should
%   contact ericwan@ece.ogi.edu for commercial licensing information.
%
%   See LICENSE (which should be part of the main toolkit distribution) for more
%   detail.

%=============================================================================================

x=x(:);

if (length(x)>1)

  x = sort(x);
  y = [x(end); x(1:end-1)];
  d = x-y;

  dups = sum(d==0);

else

  dups = 0;

end

