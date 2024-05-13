program ex0429
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
REAL*4 OUTPUT,DATA1(4096)

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
wd(1,:) = DATA1

 istat=PGOPEN('Recurslve.ps/vcps')  !PostScript
 if(istat<=0)stop 'ERR opening for PS file!'

 call pgsubp(1,2)
 call pgsci(1)
 call pgslw(4)
 call pgsch(2.0)
 call pgscf(3)
 call pgenv(0.0,50.0,-5.0,5.0,0,1) !just,axis
 call pglab('Time (sec)','Amplitude','Plot Z')
 call pgslw(2)
 call pgsci(2)
 call pgline(wh%ndata,Sec,aa)

do i=1,4096
	aa(i) = wd(1, i) - (sum(wd(1, :)) / 4096)
	zza(i) = 0.0
enddo

 do i=1,wh%ndata
  a(i)=wd(1,i)
 enddo

X1 = 0.0
X2 = 0.0
Y1 = 0.0
Y2 = 0.0
B0 = 5.91907054E-02
B1 = 0.00000000
B2 = -5.91907054E-02
A1 = -1.87047243
A2 = 0.8816185591 

DO  I=1,4096
	OUTPUT = B0*aa(i) + B1*X1 + B2*X2 -  A1*Y1 - A2*Y2 
	Y2 = Y1
	Y1 = OUTPUT
	X2 = X1
	X1 = aa(i)
	aa(i) = OUTPUT
ENDDO

print *, aa

call pgsci(1)
 call pgslw(4)
 call pgsch(2.0)
 call pgscf(3)
 call pgenv(0.0,50.0,-5.0,5.0,0,1) !just,axis
 call pglab('Time (sec)','Amplitude','Plot Z')
 call pgslw(2)
 call pgsci(2)
 call pgline(4096,Sec,OUTPUT)

call pgend

end