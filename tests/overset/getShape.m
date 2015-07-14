function shape = getShape(loc)

xi = loc(1);
eta = loc(2);
mu = loc(3);

shape = zeros(1,8);
shape(1) = 0.125*(1-xi)*(1-eta)*(1-mu);
shape(2) = 0.125*(1+xi)*(1-eta)*(1-mu);
shape(3) = 0.125*(1+xi)*(1+eta)*(1-mu);
shape(4) = 0.125*(1-xi)*(1+eta)*(1-mu);

shape(5) = 0.125*(1-xi)*(1-eta)*(1+mu);
shape(6) = 0.125*(1+xi)*(1-eta)*(1+mu);
shape(7) = 0.125*(1+xi)*(1+eta)*(1+mu);
shape(8) = 0.125*(1-xi)*(1+eta)*(1+mu);

end