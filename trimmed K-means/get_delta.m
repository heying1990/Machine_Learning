function del = get_delta(BCSS, d, lass)
   
    lam1 = 0;
    lam2 = max(abs(BCSS));
    while((lam2-lam1)>(1e-4))
        if(w_of_delta((lam1+lam2)/2, d, BCSS) < lass)
            lam2 = (lam1+lam2)/2;
        else
            lam1 = (lam1+lam2)/2;
        end
    end
    del = (lam1+lam2)/2;
    
end