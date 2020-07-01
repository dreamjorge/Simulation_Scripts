function HH=hankelHHermite(a,b,nu,mu,wo,zo,x,y,z)

    HGx  =  hermitePhysicalGaussBeam(nu,wo,zo,x,z);
    NHGx = XhermitePhysicalGaussBeam(nu,wo,zo,x,z);

    HGy  =  hermitePhysicalGaussBeam(mu,wo,zo,y,z);
    NHGy = XhermitePhysicalGaussBeam(mu,wo,zo,y,z);


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
  
 