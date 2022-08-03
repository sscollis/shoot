!=============================================================================!
        module bspline
          integer           :: nbs, kyord=5
          real, allocatable :: yknot(:)
          real, allocatable :: bs(:)
          real, external    :: bsder
        end module bspline
!=============================================================================!
        subroutine nonpar
!
! Compute the nonparallel corrections to the LST growth-rate and
! wavenumber.  Updated to second-order boundary derivatives and
! now it computes the growth-rate and wavenumber for the wall-
! normal maximum of (u,v,w).
!
! S. Scott Collis
!
! Revised:  10-15-96
!=============================================================================!
        use bspline
        implicit none

        integer, parameter :: ndof=8
        real, parameter    :: zero=0.0, pt5=0.5, one=1.0, two=2.0
        complex, parameter :: im=(0.0,1.0)

        integer :: nx, ny, i, j, idof
        real    :: dx, dy, ymax, Ma, Re, Pr

        real    :: grdke, gr(ndof), wn(ndof)
        complex :: h1, h2
        complex, external :: inprod
        real, external :: rtsafe, fumax
        real :: yumax, umaxr, umaxi

        real, allocatable    :: x(:), y(:), ke(:), dkedx(:)
        integer, allocatable :: ind(:)
        complex, allocatable :: alpha(:), beta(:), omega(:), emax(:,:)
        complex, allocatable :: q(:,:,:), dqdy(:,:,:), a(:,:,:), demax(:,:)

        complex, allocatable :: c1(:,:), c2(:,:), c3(:,:)
        complex, allocatable :: z1(:,:,:), z2(:,:,:), z3(:,:,:)

        complex, allocatable :: dqdx(:,:,:), dqdxy(:,:,:), dalphadx(:)

        character(80) :: base, fname

        integer :: iin=12, iout=13

        real, allocatable :: ty(:)
        complex, allocatable :: tq(:,:,:)

!=============================================================================!

!.... open the output.dat file and read the parameters

        write(*,'(/,"Shoot::Nonpar:  Reading output.dat...",/)')
        open(iin,file='output.dat',form='unformatted',err=100)
        read(iin) nx, ny, ymax, Ma, Re, Pr

        if (nx.lt.1) then
          write(*,*) 'Nx must be greater than 0'
          call exit(1)
        end if

        allocate( x(nx), y(ny), alpha(nx), beta(nx), omega(nx), dalphadx(nx) )
        allocate( q(ndof,ny,nx), dqdy(ndof,ny,nx), dqdx(ndof,ny,nx), &
                  dqdxy(ndof,ny,nx), a(ndof,ny,nx), ke(nx), ind(nx), &
                  dkedx(nx) )
        allocate( c1(ny,nx), c2(ny,nx), c3(ny,nx), z1(ndof,ny,nx), &
                  z2(ndof,ny,nx), z3(ndof,ny,nx), emax(ndof,nx), &
                  demax(ndof,nx) )

!.... Note that y is decreasing here

        dy = ymax / real(ny-1)
        do j = 1, ny
          y(j) = ymax - dy * real(j-1)
        end do
        do i = 1, nx
          read(iin) ind(i), x(i), alpha(i), beta(i), omega(i), &
                    q(:,:,i), dqdy(:,:,i), a(:,:,i), &
                    c1(:,i), c2(:,i), c3(:,i), &
                    z1(:,:,i), z2(:,:,i), z3(:,:,i)
        end do
        close(iin)

!.... loop over the stations

        write(*,'("Outputing polished eigenvalues and eigenfunctions",/)')

        open(22,file='NPparm.dat')
        write(22,'("# ",i5)') nx
        write(22,'("# ind(i), s(i), real(alpha(i)), aimag(alpha(i))")') 

        parm_new:  do i = 1, nx

!==============================================================================
!.... parm.new (These are the polished parallel eigensolutions on unit=22)
!==============================================================================
!     1.  profile index
!     2.  station location [usually this is arc length, s(i)]
!     3.  alpha_r
!     4.  alpha_i

          write(22,20) ind(i), x(i), real(alpha(i)), aimag(alpha(i))

!.... write-out the regular eigenfunctions

          base = 'efun'
          call makename(base,ind(i),fname)
          open (unit=iout, file=fname)
          do j = 1, ny
            write(iout,"(17(1pe13.6,1x))") y(j), &
              (real(q(idof,j,i)), aimag(q(idof,j,i)), idof=1,ndof)
          end do
          close(iout)  ! efun

!.... write-out the adjoint eigenfunctions

          base = 'adj'
          call makename(base,ind(i),fname)
          open (unit=iout, file=fname)
          do j = 1, ny
            write(iout,"(17(1pe13.6,1x))") y(j), &
              (real(a(idof,j,i)), aimag(a(idof,j,i)), idof=1,ndof)
          end do
          close(iout)

        end do parm_new  ! adj

        close(22)  ! parm.new

!============================================================================!

        if (nx.le.1) then
          write(*,*) 'Nx must be greater than 1 for nonparallel output'
          write(*,*) '  since streamwise derivatives are needed'
          call exit(1)
        end if

        write(*,"('Computing nonparallel correction...',/)")

!.... form the streamwise derivatives of the eigenfunctions

        if (.false.) then
!        dqdx(:,:,1)  = ( q(:,:,2) - q(:,:,1) ) / ( x(2) - x(1) )
        dqdx(:,:,1)  = ( -0.5 * q(:,:,3) + 2.0 * q(:,:,2) - &
                          1.5 * q(:,:,1) ) / ( x(2) - x(1) )
!        dqdxy(:,:,1) = ( dqdy(:,:,2) - dqdy(:,:,1) ) / ( x(2) - x(1) )
        dqdxy(:,:,1) = ( -0.5 * dqdy(:,:,3) + 2.0 * dqdy(:,:,2) - &
                          1.5 * dqdy(:,:,1) ) / ( x(2) - x(1) )
        do i = 2, nx-1
          dqdx(:,:,i) = ( q(:,:,i+1) - q(:,:,i-1) ) / ( x(i+1) - x(i-1) )
          dqdxy(:,:,i) = ( dqdy(:,:,i+1) - dqdy(:,:,i-1) ) / ( x(i+1)-x(i-1) )
        end do
!        dqdx(:,:,nx) = ( q(:,:,nx) - q(:,:,nx-1) ) / ( x(nx) - x(nx-1) )
        dqdx(:,:,nx)  = ( 1.5 * q(:,:,nx) - 2.0 * q(:,:,nx-1) + &
                          0.5 * q(:,:,nx-2) ) / ( x(nx) - x(nx-1) )
!        dqdxy(:,:,nx) = ( dqdy(:,:,nx) - dqdy(:,:,nx-1) ) / ( x(nx) - x(nx-1) )
        dqdxy(:,:,nx) = ( 1.5 * dqdy(:,:,nx) - 2.0 * dqdy(:,:,nx-1) + &
                          0.5 * dqdy(:,:,nx-2) ) / ( x(nx) - x(nx-1) )
        end if

        dx = (x(nx) - x(1)) / real(nx-1)
        call cg1( ndof*ny, q, dqdx, nx, dx, 0, .false. )
        call cg1( ndof*ny, dqdy, dqdxy, nx, dx, 0, .false. )

!.... form the derivative of quasi-parallel growth-rate, alpha

        if (.false.) then
!        dalphadx(1) = ( alpha(2) - alpha(1) ) / ( x(2) - x(1) )
        dalphadx(1) = ( -0.5 * alpha(3) + 2.0 * alpha(2) - &
                         1.5 * alpha(1) ) / ( x(2) - x(1) )
        do i = 2, nx-1
          dalphadx(i) = ( alpha(i+1) - alpha(i-1) ) / ( x(i+1) - x(i-1) )
        end do
!        dalphadx(nx) = ( alpha(nx) - alpha(nx-1) ) / ( x(nx) - x(nx-1) )
        dalphadx(nx) = ( 1.5 * alpha(nx) - 2.0 * alpha(nx-1) + &
                         0.5 * alpha(nx-2) ) / ( x(nx) - x(nx-1) )
        end if

        call cg1( 1, alpha, dalphadx, nx, dx, 0, .false. )

!.... compute the disturbance kinetic energy integral (trapezoid)

        do i = 1, nx
          j = 1
          ke(i) = pt5 * dy * ( abs(q(2,j,i))**2 + abs(q(3,j,i))**2 + &
                               abs(q(4,j,i))**2 )
          do j = 2, ny-1
            ke(i) = ke(i) + dy * ( abs(q(2,j,i))**2 + abs(q(3,j,i))**2 + &
                                   abs(q(4,j,i))**2 )
          end do
          j = ny
          ke(i) = ke(i) + pt5 * dy * ( abs(q(2,j,i))**2 + abs(q(3,j,i))**2 + &
                                       abs(q(4,j,i))**2 )
        end do

!.... find the maximum magnitude of the u-velocity eigenfunction

#ifdef IMSL
!.... need to switch around y and q for SLATEC versions of splines
        allocate(ty(ny), tq(ndof,ny,nx))
        do j = 1, ny
          ty(j) = y(ny-j+1)
          tq(:,j,:) = q(:,ny-j+1,:)
        enddo
        emax = 0
        nbs  = ny
        allocate( yknot(nbs+kyord), bs(nbs) )
        call BSNAK( nbs, ty, kyord, yknot )
        do i = 1, nx
          call BSINT( nbs, ty, abs(tq(2,:,i)), kyord, yknot, bs )
          do j = 1, ny-1 
            if ( abs(tq(2,j+1,i)) .lt. abs(tq(2,j,i)) .and. j.gt.1 ) goto 30
          end do
30        continue
#ifdef DEBUG 
          write(*,*) "i=",i," j=",j,"abs(q(2,1,i))=",abs(tq(2,1,i))
          write(*,*) "abs(q(2,j+1,i))=",abs(tq(2,j+1,i))
          write(*,*) "abs(q(2,j,i))=",abs(tq(2,j,i))
#endif
          if (j.gt.1) then
            yumax = rtsafe( fumax, ty(j-1), ty(j+1), 1.0e-14 )
          else
            yumax = rtsafe( fumax, ty(j), ty(j+1), 1.0e-14 )
          endif
#ifdef DEBUG
          write(*,"(2(1pe13.6,1x))") x(i), yumax
#endif
          do idof = 2, 4
            call BSINT( nbs, ty, real(tq(idof,:,i)), kyord, yknot, bs )
            umaxr = BSDER( 0, yumax, kyord, yknot, nbs, bs )
            call BSINT( nbs, ty, aimag(tq(idof,:,i)), kyord, yknot, bs )
            umaxi = BSDER( 0, yumax, kyord, yknot, nbs, bs )
            emax(idof,i) = cmplx(umaxr,umaxi)
          end do
        end do
        deallocate( yknot, bs )
        deallocate( ty, tq )
#endif

!.... compute the streamwise derivative of the disturbance kinetic energy

        if (.false.) then
!        dkedx(1) = ( ke(2) - ke(1) ) / ( x(2) - x(1) )
        dkedx(i) = ( -0.5 * ke(3) + 2.0 * ke(2) - 1.5 * ke(1) ) / &
                   ( x(2) - x(1) )
        do i = 2, nx
          dkedx(i) = ( ke(i+1) - ke(i-1) ) / ( x(i+1) - x(i-1) )
        end do
!        dkedx(nx) = ( ke(nx) - ke(nx-1) ) / ( x(nx) - x(nx-1) )
        dkedx(nx) = ( 1.5 * ke(nx) - 2.0 * ke(nx-1) + 0.5 * ke(nx-2) ) / &
                    ( x(nx) - x(nx-1) )
        end if

        call g1( 1, ke, dkedx, nx, dx, 0, .false. )

!.... compute the streamwise derivative of the maximum u velocity

        if (.false.) then
        do idof = 2, 4
!          demax(idof,1) = ( emax(idof,2) - emax(idof,1) ) / ( x(2) - x(1) )
          demax(idof,1) = ( -0.5 * emax(idof,3) + 2.0 * emax(idof,2) - &
                             1.5 * emax(idof,1) ) / ( x(2) - x(1) )
          do i = 2, nx-1
            demax(idof,i) = ( emax(idof,i+1) - emax(idof,i-1) ) / &
                            ( x(i+1) - x(i-1) )
          end do
!          demax(idof,nx) = ( emax(idof,nx) - emax(idof,nx-1) ) / &
!                           ( x(nx) - x(nx-1) )
          demax(idof,nx) = ( 1.5 * emax(idof,nx) - 2.0 * emax(idof,nx-1) + &
                             0.5 * emax(idof,nx-2) ) / ( x(nx) - x(nx-1) )
        end do
        end if

        call cg1( ndof, emax, demax, nx, dx, 0, .false. )

!.... open output files for nonparallel corrections

        open(20,file='NPsigma.dat')
        write(20,'("# s(i), -aimag(alpha(i)), -real(h2/h1),  &
                  &real(demax(2,i) / emax(2,i)), gr(2),  &
                  &real(demax(3,i) / emax(3,i)), gr(3),  &
                  &real(demax(4,i) / emax(4,i)), gr(4),  &
                  &pt5*dkedx(i)/ke(i), grdke")')
        open(21,file='NPalpha.dat')
        write(21,'("# s(i), real(alpha(i)), -aimag(h2/h1),  &
                   &aimag(demax(2,i)/emax(2,i)), wn(2),  &
                   &aimag(demax(2,i)/emax(3,i)), wn(3),  &
                   &aimag(demax(4,i)/emax(4,i)), wn(4)")')

!.... \TODO this would be better to use consistent RK4 integration

!.... form h1 and h2 using trapezoid integration at each station

        loop_h12:  do i = 1, nx

#ifdef VERBOSE
          write(*,*) "i = ", i
#endif
          j = 1
          h1 = pt5 * dy * ( c1(j,i) + inprod(ndof, z1(:,j,i), dqdx(:,j,i)) )
          h2 = pt5 * dy * ( dalphadx(i) * c2(j,i) + c3(j,i) +           &
                            inprod(ndof, z2(:,j,i), dqdx(:,j,i) ) +     &
                            inprod(ndof, z3(:,j,i), dqdxy(:,j,i) ) )
          do j = 2, ny-1
            h1 = h1 + dy * ( c1(j,i) + inprod(ndof, z1(:,j,i), dqdx(:,j,i)) )
            h2 = h2 + dy * ( dalphadx(i) * c2(j,i) + c3(j,i) +          &
                             inprod(ndof, z2(:,j,i), dqdx(:,j,i)) +     &
                             inprod(ndof, z3(:,j,i), dqdxy(:,j,i)) )
          end do
          j = ny
          h1 = h1 + pt5 * dy * ( c1(j,i) + inprod(ndof, z1(:,j,i),      &
                                 dqdx(:,j,i)) )
          h2 = h2 + pt5 * dy * ( dalphadx(i) * c2(j,i) + c3(j,i) +      &
                                 inprod(ndof, z2(:,j,i), dqdx(:,j,i)) + &
                                 inprod(ndof, z3(:,j,i), dqdxy(:,j,i)) )

!.... compute the growth-rate and wavenumber with nonparallel corrections

          grdke = -aimag(alpha(i)) - real(h2/h1) + pt5 * dkedx(i) / ke(i)

          do idof = 2, 4
            gr(idof) = -aimag(alpha(i)) - real(h2/h1) +                &
                        real(demax(idof,i) / emax(idof,i))
            wn(idof) =  real(alpha(i))  - aimag(h2/h1) +               &
                        aimag(demax(idof,i) / emax(idof,i))
          end do

!==============================================================================
!....               O u t p u t   t h e   r e s u l t s
!==============================================================================
!.... NPsigma.dat
!==============================================================================
!     1.   x(i) [s(i) for body fitted mesh]
!     2.   parallel growth-rate
!     3.   nonparallel correction
!     4-9. component data
!     10.  eigenfunction growth term
!     11.  Complete nonparallel growth-rate

          write(20,10) x(i), -aimag(alpha(i)), -real(h2/h1), &
                       real(demax(2,i) / emax(2,i)), gr(2),  &
                       real(demax(3,i) / emax(3,i)), gr(3),  &
                       real(demax(4,i) / emax(4,i)), gr(4),  &
                       pt5*dkedx(i)/ke(i), grdke
!==============================================================================
!.... NPalpha.dat
!==============================================================================
!     1.   x(i) [s(i) for body fitted mesh]
!     2.   parallel wavenumber
!     3.   nonparallel correction
!     4-9. component data (usually interested in the streamwise component

          write(21,10) x(i), real(alpha(i)), -aimag(h2/h1), &
                       aimag(demax(2,i)/emax(2,i)), wn(2),  &
                       aimag(demax(2,i)/emax(3,i)), wn(3),  &
                       aimag(demax(4,i)/emax(4,i)), wn(4)

!.... some useful output statements when debugging

#define DEBUG
#ifdef DEBUG
          write(23,10) x(i), real(h1), aimag(h1), abs(h1),  &
                             real(h2), aimag(h2), abs(h2)
          write(24,10) x(i), real(dalphadx(i)), aimag(dalphadx(i))
          j = ny-20
          write(25,10) x(i), real(c1(j,i)), aimag(c1(j,i)), abs(c1(j,i)),    &
                             real(c2(j,i)), aimag(c2(j,i)), abs(c2(j,i)),    &
                             real(c3(j,i)), aimag(c3(j,i)), abs(c3(j,i))
          write(26,10) x(i), real(inprod(ndof, z2(:,j,i), dqdx(:,j,i) )),    &
                             aimag(inprod(ndof, z2(:,j,i), dqdx(:,j,i) )),   &
                             abs(inprod(ndof, z2(:,j,i), dqdx(:,j,i) ))
          write(27,10) x(i), real(q(2,j,i)), aimag(q(2,j,i)), abs(q(2,j,i)), &
                             real(a(2,j,i)), aimag(a(2,j,i)), abs(a(2,j,i))
          write(28,10) x(i), real(inprod(ndof, z1(:,j,i), dqdx(:,j,i) )),    &
                             aimag(inprod(ndof, z1(:,j,i), dqdx(:,j,i) )),   &
                             abs(inprod(ndof, z1(:,j,i), dqdx(:,j,i) ))
          write(29,10) x(i), real(inprod(ndof, z3(:,j,i), dqdxy(:,j,i) )),    &
                             aimag(inprod(ndof, z3(:,j,i), dqdxy(:,j,i) )),   &
                             abs(inprod(ndof, z3(:,j,i), dqdxy(:,j,i) ))
#endif
#undef DEBUG
        end do loop_h12

        close(20)  ! NPsigma
        close(21)  ! NPalpha
        
        deallocate( x, y, alpha, beta, omega, dalphadx )
        deallocate( q, dqdy, dqdx, dqdxy, a, ke, ind )
        deallocate( c1, c2, c3, z1, z2, z3, emax, demax )

        return
  10    format(20(1pe20.13,1x))
  20    format(i4,1x,8(1pe20.13,1x))

  100   write(*,"('ERROR opening output.dat')")
        call exit(1)

        end subroutine nonpar

!=============================================================================!
        subroutine fumax( x, g, d )
!=============================================================================!
!       Find the (local) maximum of a bsplined function
!=============================================================================!
        use bspline
        real :: x, f, g, d
!=============================================================================!
#ifdef IMSL
        f = BSDER( 0, x, kyord, yknot, nbs, bs(:) )
        g = BSDER( 1, x, kyord, yknot, nbs, bs(:) )
        d = BSDER( 2, x, kyord, yknot, nbs, bs(:) )
#else
        write(*,*) "ERROR:  IMSL BSDER is not available for fumax"
#endif
        return
        end subroutine fumax
