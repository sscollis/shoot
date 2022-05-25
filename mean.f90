
!.... Use a profile list to generate the mean flow

!=============================================================================!
  module meanflow
!=============================================================================!
  	integer :: nym, ndofm=5, kord=3
  	real, allocatable :: ym(:), vmt(:,:), g2vmt(:,:), g22vmt(:,:)
  	real, allocatable :: vms(:,:), g2vms(:,:), g22vms(:,:), knot(:)

  end module meanflow

!=============================================================================!
	subroutine initmean(imean)
!=============================================================================!
	return
	end subroutine initmean

!=============================================================================!
	subroutine mean(imean, i)
!
! Read a meanflow profile and spline it
!
!=============================================================================!
	use meanflow
	implicit none

	integer :: imean, i

!.... local

	integer :: j, k
	integer :: ier
	real    :: tmp

	character*80 :: base, fname
!=============================================================================!

!.... read the mean field and spline to the new grid

	base = 'profile'
	call makename(base,i,fname)
#if VERBOSE
  write(*,"('Reading profile ',a)") fname
#endif
	open (imean, file=fname, form='formatted', status='old', err=1000)

!.... determine the number of points in the profile

	nym = 0
  20 continue
 	read(imean,*,end=30) tmp
	nym = nym + 1
	goto 20
  30 continue
 	rewind(imean)
#if VERBOSE
  write(*,"('Nym = ',i5)") nym
#endif

!.... allocate room for the profile data

	if ( allocated(ym) ) then
	  deallocate( ym, vmt, vms, g2vmt, g2vms, g22vmt, g22vms )
	end if
	allocate( ym(nym), vmt(nym,ndofm), vms(nym,ndofm), &
            g2vmt(nym,ndofm), g2vms(nym,ndofm), &
            g22vmt(nym,ndofm), g22vms(nym,ndofm), STAT=ier )
	if (ier .ne. 0) then
	  write(*,"('ERROR:  allocating mean field')")
	  call exit(1)
	end if

	do j = 1, nym
	  read(imean,*) ym(j), ( vmt(j,k), k = 1, ndofm )
	end do
	close(imean)

#ifdef USE_BSLIB
	allocate( knot(nym+kord) )
	call BSNAK( nym, ym, kord, knot)
	call BSINT( nym, ym, vmt(1,1), kord, knot, vms(1,1) )
	call BSINT( nym, ym, vmt(1,2), kord, knot, vms(1,2) )
	call BSINT( nym, ym, vmt(1,3), kord, knot, vms(1,3) )
	call BSINT( nym, ym, vmt(1,4), kord, knot, vms(1,4) )
	call BSINT( nym, ym, vmt(1,5), kord, knot, vms(1,5) )
#else
	call SPLINE(nym, ym, vmt(1,1), vms(1,1))
	call SPLINE(nym, ym, vmt(1,2), vms(1,2))
	call SPLINE(nym, ym, vmt(1,3), vms(1,3))
	call SPLINE(nym, ym, vmt(1,4), vms(1,4))
	call SPLINE(nym, ym, vmt(1,5), vms(1,5))
#endif

!.... read the first derivative and spline to the new grid

	base = 'first'
	call makename(base,i,fname)
#if VERBOSE
  write(*,"('Reading profile ',a)") fname
#endif
	open(imean, file=fname, form='formatted', status='old', err=1010)
	do j = 1, nym
	  read(imean,*) ym(j), ( g2vmt(j,k), k = 1, ndofm )
	end do
	close(imean)

#ifdef USE_BSLIB
	call BSINT( nym, ym, g2vmt(1,1), kord, knot, g2vms(1,1) )
	call BSINT( nym, ym, g2vmt(1,2), kord, knot, g2vms(1,2) )
	call BSINT( nym, ym, g2vmt(1,3), kord, knot, g2vms(1,3) )
	call BSINT( nym, ym, g2vmt(1,4), kord, knot, g2vms(1,4) )
	call BSINT( nym, ym, g2vmt(1,5), kord, knot, g2vms(1,5) )
#else
	call SPLINE(nym, ym, g2vmt(1,1), g2vms(1,1))
	call SPLINE(nym, ym, g2vmt(1,2), g2vms(1,2))
	call SPLINE(nym, ym, g2vmt(1,3), g2vms(1,3))
	call SPLINE(nym, ym, g2vmt(1,4), g2vms(1,4))
	call SPLINE(nym, ym, g2vmt(1,5), g2vms(1,5))
#endif

!.... read the second derivative and spline to the new grid

	base = 'second'
	call makename(base,i,fname)
#if VERBOSE
  write(*,"('Reading profile ',a)") fname
#endif
	open(imean, file=fname, form='formatted', status='old', err=1020)
	do j = 1, nym
	  read(imean,*) ym(j), ( g22vmt(j,k), k = 1, ndofm )
	end do
	close(imean)

#if USE_BSLIB
	call BSINT( nym, ym, g22vmt(1,1), kord, knot, g22vms(1,1) )
	call BSINT( nym, ym, g22vmt(1,2), kord, knot, g22vms(1,2) )
	call BSINT( nym, ym, g22vmt(1,3), kord, knot, g22vms(1,3) )
	call BSINT( nym, ym, g22vmt(1,4), kord, knot, g22vms(1,4) )
	call BSINT( nym, ym, g22vmt(1,5), kord, knot, g22vms(1,5) )
#else
	call SPLINE(nym, ym, g22vmt(1,1), g22vms(1,1))
	call SPLINE(nym, ym, g22vmt(1,2), g22vms(1,2))
	call SPLINE(nym, ym, g22vmt(1,3), g22vms(1,3))
	call SPLINE(nym, ym, g22vmt(1,4), g22vms(1,4))
	call SPLINE(nym, ym, g22vmt(1,5), g22vms(1,5))
#endif

	return

  10 format(8(1pe13.6,1x))

  1000 write(*,"('ERROR:  reading profile')")
	call exit(1)
  1010 write(*,"('ERROR:  reading first')")
	call exit(1)
  1020 write(*,"('ERROR:  reading second')")
	call exit(1)

  end subroutine mean

!=============================================================================!
	subroutine getmean( i, y, v, g1v, g2v, g11v, g12v, g22v )
!=============================================================================!
	use meanflow
	implicit none

	integer :: i
	real    :: y, v(ndofm), g1v(ndofm), g2v(ndofm)
  real    :: g11v(ndofm), g12v(ndofm) , g22v(ndofm)

	integer :: idof

  logical, save :: first_time = .true.
!=============================================================================!

#if 1

  if (first_time) then
    first_time = .false.
    write(*,*)
    write(*,"('WARNING:  getmean is not supported for 2d profile list')")
    write(*,"('          Using parallel flow approximation, so that')")
    write(*,"('          nonparallel terms should be identically zero')")
  endif

  call getmeanp( i, y, v, g2v, g22v )

  v(2)    = 0.0
  g1v     = 0.0
  g2v(3)  = 0.0
  g11v    = 0.0
  g12v    = 0.0
  g22v(3) = 0.0

#else

  if ( y .gt. ym(nym) ) then
    do idof = 1, ndofm
      v(idof)    =    vm(nym,idof)
      g1v(idof)  =  g1vm(nym,idof)
      g2v(idof)  =  g2vm(nym,idof)
      g11v(idof) = g11vm(nym,idof)
      g12v(idof) = g12vm(nym,idof)
      g22v(idof) = g22vm(nym,idof)
    end do
    return
  end if
  do idof = 1, ndofm
    call SPEVAL(nym, ym,    vm(1,idof),    vms(1,idof),  y,    v(idof))
    call SPEVAL(nym, ym,  g1vm(1,idof),  g1vms(1,idof),  y,  g1v(idof))
    call SPEVAL(nym, ym,  g2vm(1,idof),  g2vms(1,idof),  y,  g2v(idof))
    call SPEVAL(nym, ym, g11vm(1,idof), g11vms(1,idof),  y, g11v(idof))
    call SPEVAL(nym, ym, g12vm(1,idof), g12vms(1,idof),  y, g12v(idof))
    call SPEVAL(nym, ym, g22vm(1,idof), g22vms(1,idof),  y, g22v(idof))
  enddo

#endif

	return
	end subroutine getmean

!=============================================================================!
	subroutine getmeanp( i, y, vm, g2vm, g22vm )
!=============================================================================!
	use meanflow
	implicit none

	integer :: i
	real    :: y, vm(ndofm), g2vm(ndofm), g22vm(ndofm)

!.... local

	integer :: idof
	real, external :: BSVAL
!=============================================================================!

	if ( y .le. ym(nym) ) then
	  do idof = 1, ndofm
#ifdef USE_BSLIB
  	  vm(idof)    = BSVAL( y, kord, knot, nym,    vms(1,idof) )
  	  g2vm(idof)  = BSVAL( y, kord, knot, nym,  g2vms(1,idof) )
  	  g22vm(idof) = BSVAL( y, kord, knot, nym, g22vms(1,idof) )
#else
  	  call SPEVAL(nym, ym,    vmt(1,idof),    vms(1,idof), y,    vm(idof))
	    call SPEVAL(nym, ym,  g2vmt(1,idof),  g2vms(1,idof), y,  g2vm(idof))
	    call SPEVAL(nym, ym, g22vmt(1,idof), g22vms(1,idof), y, g22vm(idof))
#endif
    end do
	else
	  do idof = 1, ndofm
	    vm(idof)    = vmt(nym,idof)
	    g2vm(idof)  = g2vmt(nym,idof)
	    g22vm(idof) = g22vmt(nym,idof)
	  end do
	end if

	return
  end subroutine getmeanp
