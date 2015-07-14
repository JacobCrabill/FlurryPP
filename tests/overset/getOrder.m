% function ind = getOrder(list)
%
% returns the indices of list such that list(ind) is in ascending order
function ind = getOrder(list)
  tmp = 1:length(list);
  ind = zeros(length(list),1);
  for i=1:length(list)
    iMin = find(list(tmp)==min(list(tmp)),1,'first');
    ind(i) = tmp(iMin);
    tmp(iMin) = [];
  end
end