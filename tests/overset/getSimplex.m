% function X = getSimplex(nDims,x0)
%
% Generates an n-dimensional simplex with side length L centered at x0

function X = getSimplex(nDims,x0,L)
    dx = ones(nDims+1,1);
    X = zeros(nDims,nDims+1);
    for i=1:nDims
        dx(i+1) = sqrt(1-(i*dx(i)/(i+1))^2);
        ei = zeros(nDims,1); ei(i) = 1;
        X(:,i+1) = 1/(i)*sum(X,2) + dx(i)*ei;
    end
    
    X = X*L;
    for i=1:nDims
      X(i,:) = X(i,:) + x0(i);
    end
end