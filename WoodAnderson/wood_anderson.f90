program ex1127
 Type :: Header
	character*4  Code
	real*8 orgintime
	integer*2 ncom
	integer*4 ndata
	real :: dt
	character*4 blk
  End type Header 
Type (header) :: wh
real :: f(2048),df
real :: x(4096),y(4096),xr,xi,damping,pi,w0,w,h,x1,y1,freq
real, allocatable :: wd(:,:)
real, allocatable :: Sec(:)
! integer :: istat


open(8,file='seismicdata.bin', status='old', form='unformatted', access='stream')
Read(8)wh
Allocate (wd(wh%ncom,wh%ndata))
Allocate (Sec(wh%ndata))

Read(8, end=99)wd
99 close(8) 

do i=1,wh%ndata
 sec(i)=(i-1)*(wh%dt)
enddo

 x = wd(1,:)

 istat=PGOPEN('wood_anderson.ps/vcps')  !PostScript
 if(istat<=0)stop 'ERR opening for PS file!'

 do i=1,4096
	x(i) = wd(1, i) - (sum(wd(1, :)) / 4096)
	y(i) = 0.0
enddo


 call pgsubp(1,2)
 call pgsci(1)
 call pgslw(4)
 call pgsch(2.0)
 call pgscf(3)
 call pgenv(0.0,46.0,minval(x),maxval(x),0,1) !just,axis
 call pglab('Time (sec)','Amplitude','original data')
 call pgslw(2)
 call pgsci(2)
 call pgline(4096,Sec,x)


call FFT(x,y,4096,12)
df=1./wh%dt/REAL(4096)
damping = 0.0
do i=1,2048
 freq=df*(i-1)
 if(i.eq.1)then
  x(1)=0.0
  y(1)=0.0
  x(4096)=0.0
  y(4096)=0.0
 else
  call woodand_acc(freq,xr,xi,damping)
  x1=x(i)/REAL(4096)
  y1=y(i)/REAL(4096)
  x(i)=x1*xr-y1*xi
  y(i)=-(x1*xi+y1*xr)
  xi=-xi

  x1=x(4097-i)/REAL(4096)
  y1=y(4097-i)/REAL(4096)
  x(4097-i)=x1*xr-y1*xi
  y(4097-i)=-(x1*xi+y1*xr)
 endif
enddo

 call FFT(x,y,4096,12)



 call pgsci(1)
 call pgslw(4)
 call pgsch(2.0)
 call pgscf(3)
 call pgenv(0.0,50.0,minval(x),maxval(x),0,1) !just,axis
 call pglab('Time (sec)','Amplitude','woodanderson data')
 call pgslw(2)
 call pgsci(2)
 call pgline(4096,sec,x)

call pgend
end

!gmt psconvert Recurslve.ps -A -Tg 

SUBROUTINE FFT(XREAL,XIMAG,N,NU)
	DIMENSION XREAL(N),XIMAG(N)
	N2=N/2
	NU1=NU-1
	K=0
	DO 100 L=1,NU
102	DO 101 I=1,N2
	P=IBITR(K/2**NU1,NU)
	ARG=6.283185*P/FLOAT(N)
	C=COS(ARG)
	S=SIN(ARG)
	K1=K+1
	K1N2=K1+N2
	TREAL=XREAL(K1N2)*C+XIMAG(K1N2)*S
	TIMAG=XIMAG(K1N2)*C-XREAL(K1N2)*S
	XREAL(K1N2)=XREAL(K1)-TREAL
	XIMAG(K1N2)=XIMAG(K1)-TIMAG
	XREAL(K1)=XREAL(K1)+TREAL
	XIMAG(K1)=XIMAG(K1)+TIMAG
101	K=K+1
	K=K+N2
	IF(K.LT.N) GOTO 102
	K=0
	NU1=NU1-1
100	N2=N2/2
	DO 103 K=1,N
	I=IBITR(K-1,NU)+1
	IF(I.LE.K) GOTO 103
	TREAL=XREAL(K)
	TIMAG=XIMAG(K)
	XREAL(K)=XREAL(I)
	XIMAG(K)=XIMAG(I)
	XREAL(I)=TREAL
	XIMAG(I)=TIMAG
103	CONTINUE
	RETURN
	END
	FUNCTION IBITR(J,NU)
	J1=J
	IBITR=0
	DO 200 I=1,NU
	J2=J1/2
	IBITR=IBITR*2+(J1-2*J2)
200	J1=J2
	RETURN
	END

SUBROUTINE woodand_acc(freq,xr,xi,damping)
	complex af
	pi=3.141592654
	w0=1.25*2.*pi
	w=2.*pi*freq
	if(damping.eq.0.0)then
		h=0.8
	else
		h=damping
	endif
	xr=-w*w+w0*w0
	xi=2.*w0*w*h
	af=cmplx(xr,xi)
	af=-1.0/af
	xr=real(af)
	xi=aimag(af)
	RETURN
end