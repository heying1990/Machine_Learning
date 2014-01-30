function res = unweighted_dist(x, c)

%     res = 0;
%     for i=1:size(p1,2)
%         res = res + abs(p1(i) - p2(i));
%     end


[ndata, dimx] = size(x);
[ncentres, dimc] = size(c);
if dimx ~= dimc
	error('wdist2.m: Data dimension does not match dimension of centres')
end


res = sqrt((ones(ncentres, 1) * sum((x.*x)', 1))' + ...
  		    ones(ndata, 1) * sum((c.*c)',1) - ...
  		    2.*(x*(c')));



end