	subroutine nevoidN(x,y,z,nevN,FvN,hitvoid,wvoid)
c returns electron density nevN and fluctuation parameter FvN 
c at position designated by l,b,d,x,y,z c for a set of  
c voids with parameters read in from file  nevoidN.dat

c input:
c	x,y,z	coordinates	(kpc)  (as in TC93)
c
c output:
c	nevN	electron density in void at (x,y,z)
c	FvN	fluctuation parameter
c 	hitvoid =   0:   no void hit
c		  j>0:   j-th void hit
c	wvoid = 0,1:	 void weight	

	implicit none
	real x,y,z,nevN,FvN
	integer hitvoid, wvoid

	integer nvoidsmax
	parameter (nvoidsmax=2000)

c	character*12 losname(nvoidsmax)
	real lv(nvoidsmax), bv(nvoidsmax), dv(nvoidsmax)
        real nev(nvoidsmax), Fv(nvoidsmax)
        real aav(nvoidsmax), bbv(nvoidsmax), ccv(nvoidsmax)
        real thvy(nvoidsmax), thvz(nvoidsmax)

        real xv(nvoidsmax), yv(nvoidsmax), zv(nvoidsmax)
	real c1(nvoidsmax), s1(nvoidsmax), c2(nvoidsmax), s2(nvoidsmax),
     .       cc12(nvoidsmax), ss12(nvoidsmax), 
     .       cs21(nvoidsmax), cs12(nvoidsmax)
        integer edge(nvoidsmax)
	integer nvoids
        integer hitvoidflag
	common /voids/ nvoids, hitvoidflag(nvoidsmax)

c parameters:
c	lv	= galactic longitude of void center
c	bv	= galactic latitude of void center
c	dv	= distance from Sun of void center
c	(xv,yv,zv) = void center location (calculated)
c       nev	= internal peak electron density
c       Fv      = void fluctuation parameter
c	aav	= void major axis at 1/e
c	bbv	= void minor axis at 1/e
c	ccv	= void minor axis at 1/e
c	thvy	= rotation axis of void about y axis
c	thvz	= rotation axis of void about z axis
c	edge    = 0 => use exponential rolloff out to 5rc
c                 1 => uniform and truncated at 1/e


	real radian
	parameter(radian = 57.29577951)

	real rsun
	parameter (rsun=8.5)

	logical first
c	data first/.true./

	real slc, clc, sbc, cbc
	real rgalc
	real q
	real dx, dy, dz
	real th1, th2

c	integer luvoid
c	data luvoid/11/
	integer j
	integer voidflag
        integer just_for_fun
	save

	data first/.true./
c	data luvoid/11/



c first time through, calculate xc, yc, zc

	if(first) then 		!read void parameters
	  j=1
c	  write(6,*) 'reading nevoidN.dat.clean'
c	  open(luvoid, file='nevoidN.NE2001.dat', status='old')
c	  read(luvoid,*)				! label line
c    5     read(luvoid,*,end=99) voidflag, 
c     .      lv(j),bv(j),dv(j),				! deg, deg, kpc
c     .      nev(j),Fv(j),				! cm^{-3}, dimensionless
c     .      aav(j),bbv(j),ccv(j),thvy(j), thvz(j),  	! kpc,kpc,kpc,deg,deg
c     .      edge(j)					! 0 or 1

5         just_for_fun=0         

         if (j .le. 18) then
          voidflag=0
         else
          voidflag=1
          go to 99
         endif

           lv(1) = -81.5
           bv(1) = -0.6
           dv(1) = 0.5
           nev(1) = 0.5
           Fv(1) = 1
           aav(1) = 0.02
           bbv(1) = 0.02
           ccv(1) = 0.04
           thvy(1) = 0
           thvz(1) = 9
           edge(1) = 1
           lv(2) = -60.016
           bv(2) = -1.415
           dv(2) = 1.9
           nev(2) = 0.002
           Fv(2) = 0.1
           aav(2) = 0.5
           bbv(2) = 0.05
           ccv(2) = 0.1
           thvy(2) = 0
           thvz(2) = 30
           edge(2) = 1
           lv(3) = -44.267
           bv(3) = -4.427
           dv(3) = 1.5
           nev(3) = 0.005
           Fv(3) = 0.1
           aav(3) = 0.25
           bbv(3) = 0.1
           ccv(3) = 0.1
           thvy(3) = 0
           thvz(3) = 46
           edge(3) = 1
           lv(4) = -29.31
           bv(4) = 1.631
           dv(4) = 1.5
           nev(4) = 0.001
           Fv(4) = 0.1
           aav(4) = 1
           bbv(4) = 0.1
           ccv(4) = 0.1
           thvy(4) = 0
           thvz(4) = 61
           edge(4) = 1
           lv(5) = -29.31
           bv(5) = 1.631
           dv(5) = 3.6
           nev(5) = 0.001
           Fv(5) = 0.1
           aav(5) = 0.5
           bbv(5) = 0.2
           ccv(5) = 0.2
           thvy(5) = 0
           thvz(5) = 61
           edge(5) = 1
           lv(6) = -25.46
           bv(6) = 6.367
           dv(6) = 1.3
           nev(6) = 0.02
           Fv(6) = 0
           aav(6) = 0.9
           bbv(6) = 0.07
           ccv(6) = 0.07
           thvy(6) = -6.37
           thvz(6) = 64.54
           edge(6) = 1
           lv(7) = -16.901
           bv(7) = -2.683
           dv(7) = 1.5
           nev(7) = 0.01
           Fv(7) = 0.1
           aav(7) = 0.4
           bbv(7) = 0.1
           ccv(7) = 0.1
           thvy(7) = 0
           thvz(7) = 73
           edge(7) = 1
           lv(8) = 5.31
           bv(8) = 0.017
           dv(8) = 3
           nev(8) = 0.035
           Fv(8) = 10
           aav(8) = 1.2
           bbv(8) = 0.03
           ccv(8) = 0.1
           thvy(8) = 0
           thvz(8) = 95.31
           edge(8) = 1
           lv(9) = 7.47
           bv(9) = 0.809
           dv(9) = 2
           nev(9) = 0.055
           Fv(9) = 1.3
           aav(9) = 1.2
           bbv(9) = 0.02
           ccv(9) = 0.1
           thvy(9) = 0
           thvz(9) = 97.47
           edge(9) = 1
           lv(10) = 7.797
           bv(10) = -5.578
           dv(10) = 1.5
           nev(10) = 0.005
           Fv(10) = 0.1
           aav(10) = 0.4
           bbv(10) = 0.2
           ccv(10) = 0.2
           thvy(10) = 0
           thvz(10) = 0
           edge(10) = 1
           lv(11) = 19.85
           bv(11) = 48.34
           dv(11) = 0.35
           nev(11) = 0.004
           Fv(11) = 0
           aav(11) = 0.2
           bbv(11) = 0.2
           ccv(11) = 0.3
           thvy(11) = -41.66
           thvz(11) = -70.15
           edge(11) = 1
           lv(12) = 31
           bv(12) = 0
           dv(12) = 4.2
           nev(12) = 0.01
           Fv(12) = 0.1
           aav(12) = 0.8
           bbv(12) = 0.45
           ccv(12) = 0.2
           thvy(12) = 0
           thvz(12) = -25
           edge(12) = 1
           lv(13) = 37.21
           bv(13) = -0.637
           dv(13) = 6
           nev(13) = 0.1
           Fv(13) = 1
           aav(13) = 0.9
           bbv(13) = 0.02
           ccv(13) = 0.04
           thvy(13) = 0
           thvz(13) = 127.21
           edge(13) = 1
           lv(14) = 45
           bv(14) = 0
           dv(14) = 2
           nev(14) = 0.024
           Fv(14) = 0.1
           aav(14) = 0.8
           bbv(14) = 0.3
           ccv(14) = 0.2
           thvy(14) = 0
           thvz(14) = -25
           edge(14) = 1
           lv(15) = 48.26
           bv(15) = 0.624
           dv(15) = 3
           nev(15) = 0.024
           Fv(15) = 20
           aav(15) = 0.2
           bbv(15) = 0.2
           ccv(15) = 0.2
           thvy(15) = 0
           thvz(15) = -41.74
           edge(15) = 1
           lv(16) = 114.284
           bv(16) = 0.234
           dv(16) = 2.5
           nev(16) = 0.005
           Fv(16) = 0.1
           aav(16) = 0.5
           bbv(16) = 0.1
           ccv(16) = 0.1
           thvy(16) = 0
           thvz(16) = 24
           edge(16) = 1
           lv(17) = 129.15
           bv(17) = -2.105
           dv(17) = 1.5
           nev(17) = 0.017
           Fv(17) = 0.1
           aav(17) = 0.7
           bbv(17) = 0.1
           ccv(17) = 0.1
           thvy(17) = 0
           thvz(17) = 39
           edge(17) = 1
           

          if(voidflag .eq. 0) then
	    slc = sin(lv(j)/radian)
	    clc = cos(lv(j)/radian)
	    sbc = sin(bv(j)/radian)
	    cbc = cos(bv(j)/radian)
	    rgalc = dv(j)*cbc
	    xv(j) = rgalc*slc
	    yv(j) = rsun-rgalc*clc
	    zv(j) = dv(j)*sbc
            th1=thvy(j)
            th2=thvz(j)
            s1(j) = sin(th1/radian) 
            c1(j) = cos(th1/radian) 
            s2(j) = sin(th2/radian) 
            c2(j) = cos(th2/radian) 
            cc12(j) = c1(j)*c2(j)
            ss12(j) = s1(j)*s2(j)
            cs21(j) = c2(j)*s1(j)
            cs12(j) = c1(j)*s2(j)

c	  write(6,"(a12,1x,13(f7.3,1x))") 
c    .           losname(j),lv(j),bv(j),dv(j),
c    .           nev(j),Fv(j),xv(j),yv(j),zv(j),
c    .           aav(j), bbv(j), ccv(j),
c    .           th1, th2
	    j=j+1
	  endif
	  go to 5
   99     continue
	  first = .false.
	  nvoids = j-1
c	  close(luvoid)
	endif


	nevN = 0.
	FvN = 0.
	hitvoid = 0
	wvoid = 0
c note rotation matrix in the 'q = ' statement below
c corresponds to \Lambda_z\Lambda_y
c where \Lambda_y = rotation around y axis
c       \Lambda_z = rotation around z axis
c defined as
c \Lambda_y =  c1  0  s1 
c               0  1   0
c             -s1  0  c1

c \Lambda_z =  c2 s2   0 
c             -s2 c2   0
c               0  0   1
c =>
c \Lambda_z\Lambda_y =  c1*c2   s2   s1*c2
c                      -s2*c1   c2  -s1*s2
c                         -s1    0      c1  

c so the rotation is around the y axis first, then the z axis
	do j=1,nvoids
	  dx = x-xv(j)
	  dy = y-yv(j)
	  dz = z-zv(j)
	  q = ( cc12(j)*dx + s2(j)*dy + cs21(j)*dz)**2 / aav(j)**2 
     .      + (-cs12(j)*dx + c2(j)*dy - ss12(j)*dz)**2 / bbv(j)**2
     .      + (  -s1(j)*dx         +      c1(j)*dz)**2 / ccv(j)**2 
	  if(edge(j) .eq. 0 .and. q .lt. 3.) then
	    nevN = nev(j) * exp(-q)
            FvN = Fv(j)
	    hitvoid = j
            hitvoidflag(j)=1
	  endif
	  if(edge(j) .eq. 1 .and. q .le. 1.) then
	    nevN = nev(j)
            FvN = Fv(j)
	    hitvoid = j
            hitvoidflag(j)=1
	  endif
	enddo
	
	if(hitvoid .ne. 0) wvoid = 1
	

	return
	end




