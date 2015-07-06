function errstring = checkstructfields(ds,varargin)

% CHECKSTRUCTFIELDS Checks if a data structure has the required fields.
%
%   errstring = checkstructfields(ds, field1, field2, ...)
%
%   This function checks if data structure 'ds' has the fields specified by the
%   character arrays (strings) 'field1', 'field2', ... e.g.
%
%        missing_fields = checkstructfields(ds, 'field1', 'field2', 'field3');
%
%   The names of the missing fields are returned in 'errstring'.
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

%===============================================================================================

errstring = '';

%-- Check existence of required data structure fields
for j=1:length(varargin),
  if ~ischar(varargin{j})
    error(' [ checkstructfields ]   The field name arguments should be character arrays (strings).');
  end
  if ~isfield(ds,varargin{j}),
    errstring = [varargin{j} ' '];
  end
end

