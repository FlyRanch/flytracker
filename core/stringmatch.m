function match = stringmatch(string1,string2)

% STRINGMATCH  Returns match > 0 if string1 and string2 match (string2 can be a cell array of
%              strings). match = 0 is returned if not match is found.
%
%   match = stringmatch(string1,string2)
%

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

  match = sum(strcmp(string1,string2));

