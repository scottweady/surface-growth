function bar = progressBar(frac, width)

  if nargin < 2
    width = 20;
  end

  nfill = round(frac * width);
  bar = sprintf('[%s%s] %3.0f%%', repmat('#', 1, nfill), ...
      repmat('-', 1, width - nfill), 100 * frac);

end
