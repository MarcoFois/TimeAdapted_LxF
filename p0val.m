function xb = p0val(x)

% This function computes the mean of the object x inside each element

xb = .5 * (x(2:end,:) + x(1:end-1,:));
xb = .5 * (xb(:,2:end)+xb(:,1:end-1));

end
