function HH=hankelH2(a,b,nu,mu,wo,zl,x,y,z)

[HGx,NHGx]=hermitegz(nu,wo,zl,x,z);
[HGy,NHGy]=hermitegz(mu,wo,zl,y,z);

if (a==1)
    
    if (b==1)
    
    HH=((HGy+1i*NHGy)).*(HGx+1i*NHGx);
    
    else
        
    HH=((HGy-1i*NHGy)).*(HGx+1i*NHGx);
    
    end
 
else
    
    if (b==1)
    
    HH=((HGy+1i*NHGy)).*(HGx-1i*NHGx);
    
    else
        
    HH=((HGy-1i*NHGy)).*(HGx-1i*NHGx);
 
    end
 
end

end
  
 