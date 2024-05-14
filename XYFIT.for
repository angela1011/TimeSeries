
!
!       fitting line using X,Y
!
        subroutine xyfit(x,y,ndata,cept,slop,rms,r_cor,sdv)
        real*4 x(ndata),y(ndata)

        xsum  =0.0
        ysum  =0.0
        xxsum =0.0
        yysum =0.0
        xysum =0.0
        sum   =0.0
        errsum=0.0

        do i=1,ndata
          xp  =x(i)
          yp  =y(i)
          sum =sum+1.0
          xsum=xsum+xp
          ysum=ysum+yp
          xxsum=xxsum+xp*xp
          yysum=yysum+yp*yp
          xysum=xysum+xp*yp
        enddo


        r_cor=(sum*xysum-xsum*ysum)/sqrt( (sum*xxsum-xsum*xsum)*
     +        (sum*yysum-ysum*ysum) )


        xavg=xsum/sum
        yavg=ysum/sum
        su2 =xxsum-xsum*xavg
        sv2 =yysum-ysum*yavg
        suv =xysum-xsum*yavg
        sdx =sqrt( su2/(sum-1.0) )
        sdy =sqrt( sv2/(sum-1.0) )
        r2  =suv/sqrt(su2*sv2)
        suv2=sv2-su2
        b2  =( suv2+sqrt( suv2*suv2+4.0*suv*suv ) )/(2.0*suv)
        b1  =yavg-b2*xavg
        sdb2=b2*sqrt( (1.-r2*r2)/sum )/r2
        part1=(sdy-sdx*b2)**2.0
        part2=2.*sdx*sdy+(b2*(1.+r2)*xavg**2.0)/r2**2.
        sdb1=sqrt( ( part1+(1.-r2)*b2*part2 )/sum )

        cept=b1 !-  a     Y=aX+b
        slop=b2 !-  b

        do i=1,ndata
          err =b2*x(i)+b1-y(i)
          err2=err*err
          errsum=errsum+err2
        enddo

        sdv=sqrt(errsum/float(ndata-1))

        rms=sqrt(abs(errsum))/float(ndata)

        return
        end

