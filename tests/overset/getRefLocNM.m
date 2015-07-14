%% Reference Location Lookup  -  July 13 2015
%  Use the Nelder-Meade Algorithm to find the reference location within an
%  element of a given physical position

% Nodes of test element (hex)
zmax = .001;

xyz = [0,0,0; 1,.1,0; .9,1.1,0; .1,1,0; ...
       0,0,zmax; .9,.1,zmax; 1,1.1,zmax; .1,1,zmax];

% Physical location to match
pos = [.25,.25,zmax/2];

% Objective Function
f = @(loc) norm(pos - calcPos(loc,xyz));

maxIters = 200;
maxEvals = 1000;

% == Initialize ===========================================================
nDims = 3;
nPts = nDims + 1;
x0 = 0;
x0 = x0*ones(nDims,1);
L = .2;
X = getSimplex(nDims,x0,L);
F = zeros(nPts,1);
numEvals = 0;
for i=1:nPts
  F(i) = f(X(:,i));
  numEvals = numEvals + 1;
end
history = min(F);
Evals = numEvals;
     
% == Iterate ==============================================================
tol = 1e-10;

iter = 1;
while iter<=maxIters && numEvals <= maxEvals - 1 && min(F)>tol
  ind = getOrder(F);
  xn = X(:,ind(nPts));                  % Point with the highets value of F
  x0 = sum(X(:,ind(1:nDims)),2)/nDims;  % Centriod of all other points
  xr = x0 + (x0-xn);                    % Reflected point
  Fr = f(xr);
  numEvals = numEvals + 1;
  
  % Determine what to do with the new point
  if Fr < F(ind(nPts-1)) 
    % We will be keeping this point
    if Fr < F(ind(1))
      % This one's good; keep going!
      xe = xr + x0-xn;
      Fe = f(xe);
      numEvals = numEvals + 1;
      if Fe < Fr
        % This one's even better; use it instead
        X(:,ind(nPts)) = xe;
        F(ind(nPts)) = Fe;
      else
        % No better; stick with xr
        X(:,ind(nPts)) = xr;
        F(ind(nPts)) = Fr;
      end
    else
      % This one's somewhere in the middle; replace xn with xr
      X(:,ind(nPts)) = xr;
      F(ind(nPts)) = Fr;
    end
  else
    % Try reducing the size of the simplex
    xc = x0 - 1/2*(x0-xn);
    Fc = f(xc);
    numEvals = numEvals + 1;
    if Fc < F(ind(nPts))
      % Bringing this point in is better; use it
      X(:,ind(nPts)) = xc;
      F(ind(nPts)) = Fc;
    else
      % Bringing this point in didn't work; shrink the simplex onto the
      % smallest-valued vertex
      x1 = X(:,ind(1));
      for i=2:nPts
        X(:,ind(i)) = x1 + 1/2*(X(:,ind(i))-x1);
        F(ind(i)) = f(X(:,ind(i)));
        numEvals = numEvals + 1;
      end
    end
  end
  history = [history; min(F)];
  Evals = [Evals,numEvals];
  iter = iter + 1;
end

% == Result ===============================================================
X_Final = sum(X,2)/nPts;  % Take centroid of final simplex
F_Final = f(X_Final);     % Evaluate function at X_Final

if F_Final > min(F)
  ind = getOrder(F);
  F_Final = F(ind(1));
  X_Final = X(:,ind(1));
end  

disp('Final Result:');
iter
calcPos(X_Final,xyz)
X_Final
F_Final