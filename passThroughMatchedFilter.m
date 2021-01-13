function sig = passThroughMatchedFilter(orgSig, Fs)
  ht = zeros(1, length(orgSig));
  ht(1: Fs) = repelem([1], Fs);
  
  expandedSig = conv2(orgSig, ht);
  
  tempSig = [];
  for i = 1: Fs : length(orgSig)
    if (expandedSig(i) > 0)
      tempSig = [tempSig, 1];
    else
      tempSig = [tempSig, 0];
    end
  end
  
  sig = tempSig;
end
