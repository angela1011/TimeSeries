
!
!       fitting line using X,Y
!
        subroutine xylsfit(x,y,ndata,cept,slop,rms,r_cor,sdv)
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

        b=(ndata*xysum-xsum*ysum)/(ndata*xxsum-xsum*xsum)
        a=(ysum-b*xsum)/real(ndata)

        do i=1,ndata
          err =a+b*x(i)-y(i)
          err2=err*err
          errsum=errsum+err2
        enddo
        slop=b
        cept=a
        sdv=sqrt(errsum/float(ndata-1))

        rms=sqrt(abs(errsum))/float(ndata)

        return
        end

