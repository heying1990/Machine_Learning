function res = my_plus(x)

for i=1:size(x,2)
    if(x(1,i) < 0)
        res(1,i) = 0;
    else
        res(1,i) = x(1,i);
    end
end

end