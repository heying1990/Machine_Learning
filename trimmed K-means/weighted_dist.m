function res = weighted_dist(w, x, c)

%     res = 0;
%     for i=1:size(w,2)
%         res = res + w(i)*(p1(i) - p2(i))^2;
%     end
%     res = sqrt(res);


[ndata, dimx] = size(x);
[ncentres, dimc] = size(c);
if dimx ~= dimc
	error('wdist2.m: Data dimension does not match dimension of centres')
end


W = diag(w);
Wx = repmat(w,[ndata 1]);
Wc = repmat(w,[ncentres 1]);

res = sqrt((ones(ncentres, 1) * sum((x.*Wx.*x)', 1))' + ...
  		    ones(ndata, 1) * sum((c.*Wc.*c)',1) - ...
  		    2.*(x*W*(c')));
    
end