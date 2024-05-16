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
real :: fxi(4096),fyi(4096),fzi(4096),fx(4096),fy(4096),fz(4096),f(2048),df
real :: zza(4096),zzb(4096),zzc(4096),aa(4096),bb(4096),cc(4096)
real :: Z(4096),N(4096),E(4096)
real, allocatable :: wd(:,:)
real, allocatable :: Sec(:), a(:)

open(8,file='seismicdata.bin', status='old', form='unformatted', access='stream')
Read(8)wh
Allocate (wd(wh%ncom,wh%ndata))
Allocate (Sec(wh%ndata))
Allocate (a(wh%ndata))

Read(8, end=99)wd
99 close(8) 

do i=1,wh%ndata
 sec(i)=(i-1)*(wh%dt)
enddo

 aa = wd(1,:)

 istat=PGOPEN('butterworth.ps/vcps')  !PostScript
 if(istat<=0)stop 'ERR opening for PS file!'

 call pgsubp(1,2)

do i=1,4096
	aa(i) = wd(1, i) - (sum(wd(1, :)) / 4096)
	zza(i) = 0.0
enddo

call FFT(aa,zza,4096,12)
fx=aa/4096.0
fxi=zza/4096.0
df=1/(4096*0.01)

do i=1,2048
Z(i)=(((fx(i))**2 + (fxi(i))**2))**1/2
f(i)=i*df
end do

 do i=1,2048
  if((1 > f(i)) .or. (f(i) > 3.0))then
   fx(i) = 0.0
   fxi(i) = 0.0
   fx(4097-i)= 0.0
   fxi(4097-i)= 0.0
  endif
 end do

 call FFT(fx,-fxi,4096,12)
 call pgsci(1)
 call pgslw(4)
 call pgsch(2.0)
 call pgscf(3)
 call pgenv(0.0,50.0,-5.0,5.0,0,1) !just,axis
 call pglab('Time (sec)','Amplitude','Plot Z')
 call pgslw(2)
 call pgsci(2)
 call pgline(4096,Sec,fx)
 
 do i=1,wh%ndata
  a(i)=wd(1,i)
 enddo

 call IIRFILT(a,wh%ndata,'BUTTER',2,'BP',1.0,3.0,0.01,1)
 call pgsci(1)
 call pgslw(4)
 call pgsch(2.0)
 call pgscf(3)
 call pgenv(0.0,50.0,-5.0,5.0,0,1) !just,axis
 call pglab('Time (sec)','Amplitude','Plot Z')
 call pgslw(2)
 call pgsci(2)
 call pgline(wh%ndata,Sec,a)


call pgend
end

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