
!.... Use a LNS3D field file for nonparallel stability

!============================================================================!
       module meanflow
!============================================================================!

       integer :: nxm, nym, nzm, ndofm

!.... flow variables and derivatives

       real, allocatable :: vm(:,:,:)
       real, allocatable :: g1vm(:,:,:), g2vm(:,:,:)
       real, allocatable :: g11vm(:,:,:), g12vm(:,:,:), g22vm(:,:,:)
       real, allocatable :: g1vl(:), g2vl(:), g11vl(:), g12vl(:), g22vl(:)

!.... mesh and metrics

       real :: dxi, deta
       real, allocatable :: xm(:,:), ym(:,:)
       real, allocatable :: m1(:,:),  m2(:,:),  n1(:,:),  n2(:,:),   &
                            m11(:,:), m12(:,:), m22(:,:),            &
                            n11(:,:), n12(:,:), n22(:,:)

!.... mean splines

       real, allocatable :: vms(:,:), g1vms(:,:), g2vms(:,:),     &
                            g11vms(:,:), g12vms(:,:), g22vms(:,:)

       end module meanflow

!============================================================================!
       subroutine initmean(imean)
!============================================================================!
       use meanflow
       implicit none

       integer :: imean

!.... derivative parameters

       integer :: optx=-1, opty=-1
       logical :: xper=.false., yper=.false.
       logical :: lsym=.false. , rsym=.false., bsym=.false., tsym=.false.

!.... local

       integer :: lstep
       real    :: time, Ma, Re, Pr, gamma, cv
       integer :: i, j, k, idof
       real    :: tmp
!============================================================================!

!.... read the grid file

  open(imean,file='lstx.dat',form='unformatted',status='old',err=100)
  read(imean) nxm, nym, nzm
  allocate( xm(nym,nxm), ym(nym,nxm) )
  read(imean) (((xm(j,i), i = 1, nxm), j = 1, nym), k = 1, nzm), &
              (((ym(j,i), i = 1, nxm), j = 1, nym), k = 1, nzm), &
              (((    tmp, i = 1, nxm), j = 1, nym), k = 1, nzm)
  close(imean)

  dxi  = 1.0 / float(nxm-1)
  deta = 1.0 / float(nym-1)

!.... read the metrics

  allocate (m1(nym,nxm),  m2(nym,nxm),  n1(nym,nxm),  n2(nym,nxm), &
            m11(nym,nxm), m12(nym,nxm), m22(nym,nxm),            &
            n11(nym,nxm), n12(nym,nxm), n22(nym,nxm) )
  open (imean,file='lstm.dat',form='unformatted', status='old',err=200)
  read(imean) m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
  close(imean)

!.... read the datafile

  open(imean, file='lstq.dat', form='unformatted', status='old',err=300)
  read(imean) lstep, time, nxm, nym, nzm, ndofm, Re, Ma, Pr, gamma, cv
  allocate( vm(nym,nxm,ndofm) )
  read(imean) vm
  close(imean)

!.... Compute first derivatives of field in the mapped space

  allocate( g1vm(nym,nxm,ndofm), g2vm(nym,nxm,ndofm) )
  call grad(ndofm, nxm, nym, vm, g1vm, g2vm, dxi, deta, optx, opty, &
            xper, yper, lsym, rsym, bsym, tsym)

!.... Compute second derivatives of field

  allocate( g11vm(nym,nxm,ndofm), g12vm(nym,nxm,ndofm), &
            g22vm(nym,nxm,ndofm) )
  call grad2(ndofm, nxm, nym, vm, g1vm, g11vm, g12vm, g22vm, dxi, deta, &
             optx, opty, xper, yper, lsym, rsym, bsym, tsym)

!.... transform the gradients to physical space

  allocate(g1vl(nym), g2vl(nym), g11vl(nym), g12vl(nym), g22vl(nym))
  do idof = 1, ndofm
    do i = 1, nxm
      g1vl  = g1vm(:,i,idof)*m1(:,i) + g2vm(:,i,idof)*n1(:,i)
      g2vl  = g1vm(:,i,idof)*m2(:,i) + g2vm(:,i,idof)*n2(:,i)

      g11vl = g11vm(:,i,idof) * m1(:,i)*m1(:,i)  + &
        2.0 * g12vm(:,i,idof) * m1(:,i)*n1(:,i)  + &
        g22vm(:,i,idof)       * n1(:,i)*n1(:,i)  + &
        g1vm(:,i,idof)        * m11(:,i)    + &
        g2vm(:,i,idof)        * n11(:,i)

      g12vl = g11vm(:,i,idof) * m1(:,i)*m2(:,i)  + &
        g12vm(:,i,idof)       * m1(:,i)*n2(:,i)  + &
        g12vm(:,i,idof)       * m2(:,i)*n1(:,i)  + &
        g22vm(:,i,idof)       * n1(:,i)*n2(:,i)  + &
        g1vm(:,i,idof)        * m12(:,i)    + &
        g2vm(:,i,idof)        * n12(:,i)

      g22vl = g11vm(:,i,idof) * m2(:,i)*m2(:,i)  + &
        2.0 * g12vm(:,i,idof) * m2(:,i)*n2(:,i)  + &
        g22vm(:,i,idof)       * n2(:,i)*n2(:,i)  + &
        g1vm(:,i,idof)        * m22(:,i)    + &
        g2vm(:,i,idof)        * n22(:,i)

      g1vm(:,i,idof)  = g1vl
      g2vm(:,i,idof)  = g2vl
      g11vm(:,i,idof) = g11vl
      g12vm(:,i,idof) = g12vl
      g22vm(:,i,idof) = g22vl
    end do
  end do
  deallocate( g1vl, g2vl, g11vl, g12vl, g22vl )

!.... don't need the metrics anymore

  deallocate( m1, m2, n1, n2, m11, m12, m22, n11, n12, n22 )

!.... allocate space for the meanflow splines

  allocate( vms(nym,ndofm), g1vms(nym,ndofm), g2vms(nym,ndofm), &
            g11vms(nym,ndofm), g12vms(nym,ndofm), g22vms(nym,ndofm) )

  return

  100  write(*,"('ERROR:  reading lstx.dat')")
  call exit(1)
  200  write(*,"('ERROR:  reading lstm.dat')")
  call exit(1)
  300  write(*,"('ERROR:  reading lstq.dat')")
  call exit(1)

  end

!============================================================================!
  subroutine mean(imean, i)
!
! Spline the meanflow profile at station i
!
!============================================================================!
  use meanflow
  implicit none

  integer :: imean, i

!.... local variables

  integer :: idof
!============================================================================!
  do idof = 1, ndofm
    call SPLINE(nym, ym(1,i),    vm(1,i,idof),    vms(1,idof))
    call SPLINE(nym, ym(1,i),  g1vm(1,i,idof),  g1vms(1,idof))
    call SPLINE(nym, ym(1,i),  g2vm(1,i,idof),  g2vms(1,idof))
    call SPLINE(nym, ym(1,i), g11vm(1,i,idof), g11vms(1,idof))
    call SPLINE(nym, ym(1,i), g12vm(1,i,idof), g12vms(1,idof))
    call SPLINE(nym, ym(1,i), g22vm(1,i,idof), g22vms(1,idof))
  end do
  return
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
  if ( y .gt. ym(nym,i) ) then
    do idof = 1, ndofm
      v(idof)    = vm(nym,i,idof)
      g1v(idof)  = g1vm(nym,i,idof)
      g2v(idof)  = g2vm(nym,i,idof)
      g11v(idof) = g11vm(nym,i,idof)
      g12v(idof) = g12vm(nym,i,idof)
      g22v(idof) = g22vm(nym,i,idof)
    end do
    return
  end if
  call MSPEVAL(nym, ym(1,i), ndofm, vm(:,i,:),    &
               vms(:,:), y, v(:))
  call MSPEVAL(nym, ym(1,i), ndofm, g1vm(:,i,:),  &
               g1vms(:,:), y, g1v(:))
  call MSPEVAL(nym, ym(1,i), ndofm, g2vm(:,i,:),  &
               g2vms(:,:), y, g2v(:))
  call MSPEVAL(nym, ym(1,i), ndofm, g11vm(:,i,:), &
               g11vms(:,:), y, g11v(:))
  call MSPEVAL(nym, ym(1,i), ndofm, g12vm(:,i,:), &
               g12vms(:,:), y, g12v(:))
  call MSPEVAL(nym, ym(1,i), ndofm, g22vm(:,i,:), &
               g22vms(:,:), y, g22v(:))
  return
  end

!=============================================================================!
  subroutine getmeanp( i, y, v, g2v, g22v )
!=============================================================================!
  use meanflow
  implicit none

  integer :: i
  real    :: y, v(ndofm), g2v(ndofm), g22v(ndofm)

!.... local

  integer :: idof
!=============================================================================!
  if ( y .gt. ym(nym,i) ) then
    do idof = 1, ndofm
      v(idof)    = vm(nym,i,idof)
      g2v(idof)  = g2vm(nym,i,idof)
      g22v(idof) = g22vm(nym,i,idof)
    end do
    return
  end if
  call MSPEVAL(nym, ym(1,i), ndofm, vm(:,i,:),    &
               vms(:,:), y, v(:))
  call MSPEVAL(nym, ym(1,i), ndofm, g2vm(:,i,:),  &
               g2vms(:,:), y, g2v(:))
  call MSPEVAL(nym, ym(1,i), ndofm, g22vm(:,i,:), &
               g22vms(:,:), y, g22v(:))
  return
  end

!=============================================================================!
  subroutine mspeval(n, x, ndof, y, ys, xl, yl)
!=============================================================================!
  integer :: n, ndof
  real :: x(n), y(n,ndof), ys(n,ndof), xl, yl(ndof)
  real, parameter :: onesixth = 0.16666666666666666
!=============================================================================!
  KLO=1
  KHI=N
1 if (KHI-KLO.GT.1) then
    K=(KHI+KLO) / 2
    if(X(K).GT.XL)then
      KHI=K
    else
      KLO=K
    endif
    goto 1
  endif
  DXM    = XL - X(KLO)
  DXP    = X(KHI) - XL
  DEL    = X(KHI) - X(KLO)
  DELINV = 1.0 / DEL
  yl(:) = ( YS(KLO,:)*DXP*(DXP*DXP*DELINV - DEL) + &
            YS(KHI,:)*DXM*(DXM*DXM*DELINV - DEL) ) * onesixth + &
          ( Y(KLO,:)*DXP + Y(KHI,:)*DXM ) * DELINV
  return
  end
