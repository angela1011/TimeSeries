
!-------------------------------------------------------------------
!-- Find matrix inverse from numerical recipes
subroutine MATRIXINV(c,n)

      parameter (nb=3)
      dimension c(nb,nb),indx(nb),y(nb,nb)

      DO I=1,n
       DO J=1,n
         y(i,j)=0.
       ENDDO
       y(i,i)=1.
      ENDDO

      call LUDCMP(C,N,INDX,D)

      do j=1,n
        CALL LUBKSB(C,n,indx,y(1,j))
      enddo

      DO I=1,n
       DO J=1,n
         C(I,J)=y(i,j)
       ENDDO
      ENDDO

      RETURN
END subroutine MATRIXINV


SUBROUTINE LUDCMP(C,N,INDX,D)
      parameter (nb=3,tiny=1.e-20)
      DIMENSION C(nb,nb),INDX(nb),VV(nb)

      D=1.0
      DO I=1,N
        AAMAX=0.0
        DO J=1,N
          IF (ABS(C(I,J)).GT.AAMAX) AAMAX=ABS(C(I,J))
        enddo
        IF (AAMAX.EQ.0.0) PAUSE 'Singular matrix.'
        VV(I)=1.0/AAMAX
      ENDDO

      DO J=1,N
        IF (J.GT.1) THEN
          DO I=1,J-1
            SUM=C(I,J)
            IF (I.GT.1)THEN
              DO K=1,I-1
                SUM=SUM-C(I,K)*C(K,J)
              enddo
              C(I,J)=SUM
            ENDIF
          ENDDO
        ENDIF

        AAMAX=0.0
        DO I=J,N
          SUM=C(I,J)
          IF (J.GT.1)THEN
            DO K=1,J-1
              SUM=SUM-C(I,K)*C(K,J)
            ENDDO
            C(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
        ENDDO

        IF (J.NE.IMAX)THEN
          DO K=1,N
            DUM=C(IMAX,K)
            C(IMAX,K)=C(J,K)
            C(J,K)=DUM
          ENDDO
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX

        IF(J.NE.N)THEN
          IF(C(J,J).EQ.0.0)C(J,J)=TINY
          DUM=1.d0/C(J,J)
          DO I=J+1,N
            C(I,J)=C(I,J)*DUM
          ENDDO
        ENDIF
      enddo

      IF(C(N,N).EQ.0.0)C(N,N)=TINY

      RETURN
END SUBROUTINE LUDCMP


SUBROUTINE LUBKSB(a,n,indx,b)
      parameter (nb=3)
      dimension indx(nb),b(nb),a(nb,nb)
      real sum

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if(ii.ne.0) then
         do j=ii,i-1
          sum=sum-a(i,j)*b(j)
         enddo
        else if(sum.ne.0) then
           ii=i
        endif
        b(i)=sum
      enddo

      do i= n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo

      return
end SUBROUTINE LUBKSB
