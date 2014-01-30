function res = w_of_delta(del, d, BCSS)




res = sum(abs(my_plus(BCSS-del*ones(1,d))/sqrt(my_plus(BCSS-del*ones(1,d))*my_plus(BCSS-del*ones(1,d))')));




end