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
real :: f(2048),df
real, allocatable :: wd(:,:)
real, allocatable :: Sec(:)
real :: x(4096),y(4096),xr,xi,damping,pi,w0,w,h,x1,y1,freq,c1,c2,WA(4096),dt
   
	open(8,file='seismicdata.bin', status='old', form='unformatted', access='stream')
	Read(8)wh
	Allocate (wd(wh%ncom,wh%ndata))
	Allocate (Sec(wh%ndata))


	Read(8, end=99)wd
	99 close(8) 
	
	do i=1,wh%ndata
	 sec(i)=(i-1)*(wh%dt)
	enddo
	
	dt=wh%dt
	x = wd(1,:)
   
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

	pi=3.141592654	
	W0=2*pi/0.8
	h=0.8
	c1=1+h*W0*dt
	c2=1+2*h*W0*dt+(W0*dt)**2
	WA=0.0
	do i=3,4096
		WA(i)=1/c2*(2800*x(i)*(dt**2)+2*c1*WA(i-1)-WA(i-2))
	enddo

	call pgsci(1)
	call pgslw(4)
	call pgsch(2.0)
	call pgscf(3)
	call pgenv(0.0,46.0,minval(WA),maxval(WA),0,1) !just,axis
	call pglab('Time (sec)','Amplitude','Plot Z')
	call pgslw(2)
	call pgsci(2)
	call pgline(4096,Sec,WA)

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