program ex20240520
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
	real :: x(4096),y(4096),xr,xi,damping,pi,w0,w,h,x1,y1,freq
   
	open(8,file='../seismicdata.bin', status='old', form='unformatted', access='stream')
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

! ! call FFT(x,y,4096,12)
! df=1./wh%dt/REAL(4096)
! damping = 0.0
! do i=1,2048
!  freq=df*(i-1)
!  if(i.eq.1)then
!   x(1)=0.0
!   y(1)=0.0
!   x(4096)=0.0
!   y(4096)=0.0
!  else
!   call woodand_acc(freq,xr,xi,damping)
!   x1=x(i)/REAL(4096)
!   y1=y(i)/REAL(4096)
!   x(i)=x1*xr-y1*xi
!   y(i)=-(x1*xi+y1*xr)
!   xi=-xi

!   x1=x(4097-i)/REAL(4096)
!   y1=y(4097-i)/REAL(4096)
!   x(4097-i)=x1*xr-y1*xi
!   y(4097-i)=-(x1*xi+y1*xr)
!  endif
! enddo

!  call pgsci(1)
!  call pgslw(4)
!  call pgsch(2.0)
!  call pgscf(3)
!  call pgenv(0.0,50.0,minval(x),maxval(x),0,1) !just,axis
!  call pglab('Time (sec)','Amplitude','woodanderson data')
!  call pgslw(2)
!  call pgsci(2)
!  call pgline(4096,sec,x)

call pgend
end

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