!=============================================================================!
      subroutine CONTE(nstep, tol, n, r, yo, yf, to, tf, eigfun, &
                       FHOMO, BC, y)
!=============================================================================!
!
!     First order linear boundary value problem solver using Conte's
!     method.  Fourth order Runge-Kutta is used for time advancement
!
!=============================================================================!
      implicit none

!.... Argument declarations

      integer :: nstep, n, r
      real    :: tol, to, tf
      complex :: yo(n,n), yf(n,n), bc(n), y(n,0:nstep)
      integer :: eigfun

!.... Other variables

      integer :: i, j, k, m, q, s, mi, mj, qend
      real    :: t, tq(0:nstep), h 
      complex :: B(n-r,0:nstep), Utemp(n), temp
      complex :: U(n,n-r,0:nstep), P(n-r,n-r,0:nstep), z(n,n-r)
      complex :: w(n-r), eta(n)
      complex :: ut(n,n-r)

!.... Temporary vars for automatic normalization

      real    :: aa, bb, cc, test
      logical :: norm

#ifdef USE_ODEINT

!.... externally defined routines

      complex, external :: inprod
      external FHOMO, RKQC

!.... Variables for ODEINT

      integer :: nbad, nok
#else
      complex, external :: inprod
      external FHOMO
#endif
!=============================================================================!

!.... initialize some variables

      U  = 0.0
      B  = 0.0
      P  = 0.0
      tq = 0.0
      y  = 0.0

!.... compute the step size

      h = (tf - to) / real(nstep)

!.... Begin the eigenvalue iteration loop

      q = 0
      tq(0) = to

!.... Set the initial conditions

      k = 0
      U = 0.0
      U(:,1:(n-r),k) = yo(:,1:(n-r))

!.... Gram-Schmidt

      w(1) = SQRT( inprod(n, U(1,1,k), U(1,1,k)) )
      do i = 1, n
        z(i,1) = U(i,1,k) / w(1)
      end do
      do mi = 2, (n-r)
        do i = 1, n
          eta(i) = U(i,mi,k)
        end do
        do mj = mi-1, 1, -1
          temp = inprod(n, U(1,mi,k), z(1,mj))
          do i = 1, n
            eta(i) = eta(i) - temp * z(i,mj)
          end do
        end do
        w(mi) = SQRT( inprod(n, eta, eta) )
        do i = 1, n
          z(i,mi) = eta(i) / w(mi)
        end do
      end do

!.... Now update the U matrix with the orthonormal values

      do i = 1, n
        do m = 1, n-r
          U(i,m,k) = z(i,m)
        end do
      end do

!.... Integrate the homogeneous equations

      do k = 1, nstep
        t = to + h*k

!.... Loop thru all homogeneous solutions

        do m = 1, n-r
#ifdef USE_ODEINT
          write(*,*) "Using ODEINT..."
          do i = 1, n
            Utemp(i) = U(i,m,k-1)
          end do
          call ODEINT(Utemp,n,t-h,t,1.E-5,h/2.0,1.e-20,nok,nbad, &
                      FHOMO,RKQC)
          write (*,*) k, nok, nbad
          do i = 1, n
            U(i,m,k) = Utemp(i)
          end do
#else
          call RK4(n, U(1,m,k-1), U(1,m,k), t-h, h, FHOMO)
#endif
        end do

!.... Test to see if normalization is required

        norm = .false.
        do mi = 1, n-r
          do mj = 1, n-r
            if (mi .ne. mj) then
              aa = ABS(inprod(n, U(1,mi,k), U(1,mi,k)))
              bb = ABS(inprod(n, U(1,mj,k), U(1,mj,k)))
              cc = ABS(inprod(n, U(1,mi,k), U(1,mj,k)))
              test = cc/SQRT(aa*bb)
              if (test .gt. tol) norm = .true.
            end if
          end do
        end do

!.... Perform normalization

        if ( norm .or. (k .eq. nstep) ) then
          q = q + 1
          tq(q) = t
          if (k .eq. nstep) then
            qend = q
          end if

!.... Gram-Schmidt

          w(1) = SQRT( inprod(n, U(1,1,k), U(1,1,k)) )
          do i = 1, n
            z(i,1) = U(i,1,k) / w(1)
          end do
          do mi = 2, (n-r)
            do i = 1, n
              eta(i) = U(i,mi,k)
            end do
            do mj = mi-1, 1, -1
              temp = inprod(n, U(1,mi,k), z(1,mj))
              do i = 1, n
                eta(i) = eta(i) - temp * z(i,mj)
              end do
            end do
            w(mi) = SQRT( inprod(n, eta, eta) )
            do i = 1, n
              z(i,mi) = eta(i) / w(mi)
            end do
          end do

!.... Now I have the orthonormal basis in z and 
!.... the norms in w so I can compute the P orthonormalization 
!.... matrix

          do j = 1, n-r
            do i = 1, j
              if (j .eq. i) then
                P(j,i,q) = 1.0 / w(j)
              else
                P(j,i,q) = 0.0
                do s = i, j-1
                  P(j,i,q) = P(j,i,q)-inprod(n,U(1,j,k),z(1,s)) / &
                             w(j) * P(s,i,q)
                end do
              end if
            end do
          end do

!.... Check the P matrix

          if (.false.) then
            do i = 1, n
              do m = 1, n-r
                ut(i,m) = 0.0
                do j = 1, n-r
                  ut(i,m) = ut(i,m) + U(i,j,k)*P(m,j,q)
                end do
              end do
            end do
            do i = 1,n
              write(*,*) i,(ut(i,m) - z(i,m), m = 1, n-r)
            end do
            write (*,*)
            write (*,*)
          end if

!.... Now update the U matrix with the orthonormal values

          do i = 1, n
            do m = 1, n-r
              U(i,m,k) = z(i,m)
            end do
          end do

        end if     ! norm
      end do       ! nstep

!.... return the solutions at the last node

#ifdef VERBOSE
      do k = 0, nstep
        write(*,*) "k = ",k," U(:,:,k) = ", U(:,:,k)
      enddo
#endif

      yf = 0.0
      yf(:,1:n-r) = U(:,1:n-r,nstep)

!.... If you would like to see the eigenfunction

      if (eigfun.eq.1) then

        B(:,q) = BC

        do i = 1, n
          y(i,nstep) = 0.0
          do m = 1, n-r
            y(i,nstep) = y(i,nstep) + U(i,m,nstep)*B(m,q)
          end do
        end do

        do m = 1, n-r
          B(m,q-1) = 0.0
          do j = 1, n-r
            B(m,q-1) = B(m,q-1) + P(j,m,q)*B(j,q) 
          end do
        end do

        do k = nstep-1, 0, -1
          t = to + h*k
          if ( t .gt. tq(q-1) ) then
            q = q - 1
            do m = 1, n-r
              B(m,q-1) = 0.0
              do j = 1, n-r
                B(m,q-1) = B(m,q-1) + P(j,m,q)*B(j,q)
              end do
            end do
          end if
          do i = 1, n
            y(i,k) = 0.0
            do m = 1, n-r
              y(i,k) = y(i,k) + U(i,m,k)*B(m,q-1)
            end do
          end do
        end do

      end if    ! eigfun

      return
      end

!=============================================================================!
      function INPROD(n, v1, v2)
!=============================================================================!
!
!.... Perform and inner product on two complex vectors, v1 and v2
!.... without taking the complex conjugate
!
!=============================================================================!
      implicit none

      integer :: n
      complex :: v1(n), v2(n)
      complex :: INPROD
      integer :: i
!=============================================================================!

!.... The analytic inner product should yield faster convergence
 
      INPROD = 0.0
      do i = 1, n
        INPROD = INPROD + v1(i) * v2(i)
      end do
      
      return
      end

!=============================================================================!
      function CINPROD(n, v1, v2)
!=============================================================================!
!
!.... Perform and inner product on two complex vectors, v1 and v2
!.... using the conjugate of the second vector
! 
!.... Note that this is a slightly weird definition of the 
!.... complex inner-produce
!=============================================================================!
      implicit none

      integer :: n
      complex :: v1(n), v2(n)
      complex :: CINPROD
      integer :: i
!=============================================================================!

      CINPROD = 0.0
      do i = 1, n
        CINPROD = CINPROD + v1(i) * conjg(v2(i))
      end do
      
      return
      end
