
! \todo Use a Streett bl_sta.out file for the meanflow

! SSC:  this is not implemented yet...

!=============================================================================!
  module meanflow
!=============================================================================!
    integer :: nym, ndofm=5, kord
    real, allocatable :: ym(:), vmt(:,:), g2vmt(:,:), g22vmt(:,:)
    real, allocatable :: vms(:,:), g2vms(:,:), g22vms(:,:), knot(:)
  end module meanflow

!=============================================================================!
  subroutine initmean(imean)
!=============================================================================!
  return
  end

!=============================================================================!
  subroutine mean(imean, i)
!
!  Read a meanflow profile and spline it (often computing using compbl)
!
!=============================================================================!
  use meanflow
  implicit none

  integer :: imean, i

!.... local

  integer :: j, k
  integer :: ier
  real    :: tmp

  character(80) :: base, fname
!=============================================================================!

!.... read the mean field and spline to the new grid

  base = 'profile'
  call makename(base,i,fname)
  open (imean, file=fname, form='formatted', status='old', err=1000)

!.... determine the number of points in the profile

  nym = 0
  20 continue
   read(imean,*,end=30) tmp
  nym = nym + 1
  goto 20
  30 continue
   rewind(imean)

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

  allocate( knot(nym+kord) )
  call BSNAK( nym, ym, kord, knot)
  call BSINT( nym, ym, vmt(1,1), kord, knot, vms(1,1) )
  call BSINT( nym, ym, vmt(1,2), kord, knot, vms(1,2) )
  call BSINT( nym, ym, vmt(1,3), kord, knot, vms(1,3) )
  call BSINT( nym, ym, vmt(1,4), kord, knot, vms(1,4) )
  call BSINT( nym, ym, vmt(1,5), kord, knot, vms(1,5) )

!  call SPLINE(nym, ym, vmt(1,1), vms(1,1))
!  call SPLINE(nym, ym, vmt(1,2), vms(1,2))
!  call SPLINE(nym, ym, vmt(1,3), vms(1,3))
!  call SPLINE(nym, ym, vmt(1,4), vms(1,4))
!  call SPLINE(nym, ym, vmt(1,5), vms(1,5))

!.... read the first derivative and spline to the new grid

  base = 'first'
  call makename(base,i,fname)
  open(imean, file=fname, form='formatted', status='old',err=1010)
  do j = 1, nym
    read(imean,*) ym(j), ( g2vmt(j,k), k = 1, ndofm )
  end do
  close(imean)

  call BSINT( nym, ym, g2vmt(1,1), kord, knot, g2vms(1,1) )
  call BSINT( nym, ym, g2vmt(1,2), kord, knot, g2vms(1,2) )
  call BSINT( nym, ym, g2vmt(1,3), kord, knot, g2vms(1,3) )
  call BSINT( nym, ym, g2vmt(1,4), kord, knot, g2vms(1,4) )
  call BSINT( nym, ym, g2vmt(1,5), kord, knot, g2vms(1,5) )

!  call SPLINE(nym, ym, g2vmt(1,1), g2vms(1,1))
!  call SPLINE(nym, ym, g2vmt(1,2), g2vms(1,2))
!  call SPLINE(nym, ym, g2vmt(1,3), g2vms(1,3))
!  call SPLINE(nym, ym, g2vmt(1,4), g2vms(1,4))
!  call SPLINE(nym, ym, g2vmt(1,5), g2vms(1,5))

!.... read the second derivative and spline to the new grid

  base = 'second'
  call makename(base,i,fname)
  open(imean, file=fname, form='formatted', status='old',err=1020)
  do j = 1, nym
    read(imean,*) ym(j), ( g22vmt(j,k), k = 1, ndofm )
  end do
  close(imean)

  call BSINT( nym, ym, g22vmt(1,1), kord, knot, g22vms(1,1) )
  call BSINT( nym, ym, g22vmt(1,2), kord, knot, g22vms(1,2) )
  call BSINT( nym, ym, g22vmt(1,3), kord, knot, g22vms(1,3) )
  call BSINT( nym, ym, g22vmt(1,4), kord, knot, g22vms(1,4) )
  call BSINT( nym, ym, g22vmt(1,5), kord, knot, g22vms(1,5) )

!  call SPLINE(nym, ym, g22vmt(1,1), g22vms(1,1))
!  call SPLINE(nym, ym, g22vmt(1,2), g22vms(1,2))
!  call SPLINE(nym, ym, g22vmt(1,3), g22vms(1,3))
!  call SPLINE(nym, ym, g22vmt(1,4), g22vms(1,4))
!  call SPLINE(nym, ym, g22vmt(1,5), g22vms(1,5))

  return

  10 format(8(1pe13.6,1x))

  1000 write(*,"('ERROR:  reading profile')")
  call exit(1)
  1010 write(*,"('ERROR:  reading first')")
  call exit(1)
  1020 write(*,"('ERROR:  reading second')")
  call exit(1)
  end

!=============================================================================!
  subroutine getmean( i, y, v, g1v, g2v, g11v, g12v, g22v )
!=============================================================================!
  use meanflow
  implicit none

  integer :: i
  real    :: y, v(ndofm), g1v(ndofm), g2v(ndofm)
  real    :: g11v(ndofm), g12v(ndofm) , g22v(ndofm)

  integer :: idof
!=============================================================================!
  write(*,"('ERROR:  Getmean is not supported in this version')")
  call exit(1)
  return
  end

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
      vm(idof)    = BSVAL( y, kord, knot, nym,    vms(1,idof) )
      g2vm(idof)  = BSVAL( y, kord, knot, nym,  g2vms(1,idof) )
      g22vm(idof) = BSVAL( y, kord, knot, nym, g22vms(1,idof) )
!      call SPEVAL(nym, ym,    vmt(1,idof),    vms(1,idof), y,    vm(idof))
!      call SPEVAL(nym, ym,  g2vmt(1,idof),  g2vms(1,idof), y,  g2vm(idof))
!      call SPEVAL(nym, ym, g22vmt(1,idof), g22vms(1,idof), y, g22vm(idof))
    end do
  else
    do idof = 1, ndofm
      vm(idof)    = vmt(nym,idof)
      g2vm(idof)  = g2vmt(nym,idof)
      g22vm(idof) = g22vmt(nym,idof)
    end do
  end if

  return
  end
