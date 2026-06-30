function opts = parse_opts(opts, varargin)
if mod(numel(varargin),2)~=0, error('Name/value pairs expected.'); end
for k=1:2:numel(varargin)
 name = varargin{k}; val = varargin{k+1};
 if ~isfield(opts, name), error('Unknown option: %s', name); end
 opts.(name) = val;
end
end