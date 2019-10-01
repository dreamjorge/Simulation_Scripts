for jj=1:numberpoints
    %vectors of zeros for each x,y-component ray for H1 and H2
    xp = pt(jj,1)
    yp = pt(jj,2)
    zp = 0;
    ray(jj).rxHH11(1) = xp;
    ray(jj).rxHH12(1) = xp;
    ray(jj).rxHH21(1) = xp;
    ray(jj).rxHH22(1) = xp;
    ray(jj).ryHH11(1) = yp;
    ray(jj).ryHH12(1) = yp;
    ray(jj).ryHH21(1) = yp;
    ray(jj).ryHH22(1) = yp;
    
end






for jj = numberpoints
  
    x=ray(jj).rxHH11(1)
    y=ray(jj).ryHH11(1)
    plot(x,y,'*')
end