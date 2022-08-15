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
      complex :: yo(n,n), yf(n,n), BC(n), y(n,0:nstep)
      integer :: eigfun

!.... Other variables

      integer :: i, j, k, m, q, s, mi, mj, qend
      real    :: t, tq(0:nstep), h
      complex :: B(n-r,0:nstep), Utemp(n), temp
      !complex :: U(n,n-r,0:nstep), P(n-r,n-r,0:nstep), z(n,n-r)
      complex,allocatable :: U(:,:,:), P(:,:,:), z(:,:)
      complex :: w(n-r), eta(n)
      complex :: ut(n,n-r)

!.... Temporary vars for automatic normalization

      real    :: aa, bb, cc, test
      logical :: norm

!.... externally defined routines

      complex, external :: inprod

#ifdef USE_ODEINT
!.... ODEINT is for REAL -- not implemented for COMPLEX yet 
      external FHOMO, NR_CRKQC
      integer :: nbad, nok
#else
!.... RK4 works for COMPLEX 
      external FHOMO
#endif

!=============================================================================!

      allocate( U(n,n-r,0:nstep), P(n-r,n-r,0:nstep), z(n,n-r) )

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

      !write(*,*) "1:  I am here..."

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

      !write(*,*) "2: I am here..."

!.... Now update the U matrix with the orthonormal values

      do i = 1, n
        do m = 1, n-r
          U(i,m,k) = z(i,m)
        end do
      end do

!.... Integrate the homogeneous equations

      !write(*,*) "3: I am here..."

      do k = 1, nstep
#if VERBOSE>=2
        write(*,*) "Conte: k = ",k," with nstep = ",nstep
#endif
        t = to + h*k

!.... Loop thru all homogeneous solutions

        do m = 1, n-r
#ifdef USE_ODEINT
          write(*,*) "Using experimental NR_CODEINT..."
          do i = 1, n
            Utemp(i) = U(i,m,k-1)
          end do
          call NR_CODEINT(Utemp,n,t-h,t,1.E-5,h/2.0,1.e-20,nok,nbad, &
                          FHOMO,NR_CRKQC)
          write (*,*) k, nok, nbad
          do i = 1, n
            U(i,m,k) = Utemp(i)
          end do
#else
#if VERBOSE>=2
          write(*,*) "Starting RK4 with k = ",k," m = ",m
#endif
          call RK4(n, U(1,m,k-1), U(1,m,k), t-h, h, FHOMO)
#if VERBOSE>=2
          write(*,*) "Finished RK4 with k = ",k," m = ",m
#endif
#endif
        end do

        !write(*,*) "4: I am here..."

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

      !write(*,*) "5: I am here..."

!.... return the solutions at the last node

#if VERBOSE>=3
      do k = 0, nstep
        write(*,*) "k = ",k," U(:,:,k) = ", U(:,:,k)
      enddo
#endif

      yf = 0.0
      yf(:,1:n-r) = U(:,1:n-r,nstep)

!.... If you would like to see the eigenfunction

      if (eigfun.eq.1) then

        !B(:,q) = BC  ! SSC:  B(n-r,:) whereas BC(n), what to do?
        !B(:,q) = BC(r+1:n)
        B(:,q) = BC(1:n-r)

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

      !write(*,*) "6: I am here..."
     
      deallocate( U, P, z )
 
      return
      end

!=============================================================================!
      function inprod(n, v1, v2)
!=============================================================================!
!
!.... Pick which inner-product you wish to use 
!
!=============================================================================!
      implicit none
      integer :: n
      complex :: v1(n), v2(n)
      complex :: inprod
      complex, external :: cinprod, c1inprod, ninprod
!=============================================================================!
      !inprod =  cinprod(n, v1, v2) 
      !inprod = c1inprod(n, v1, v2) 
      inprod  =  ninprod(n, v1, v2) 
      return
      end
!=============================================================================!
      function ninprod(n, v1, v2)
!=============================================================================!
!
!.... Perform and inner product on two complex vectors, v1 and v2
!.... without taking the complex conjugate
!
!=============================================================================!
      implicit none
      integer :: n
      complex :: v1(n), v2(n)
      complex :: ninprod 
      integer :: i
!=============================================================================!
!.... This analytic inner product may yield faster convergence, but is not
!.... a real inner-product
      ninprod = 0.0
      do i = 1, n
        ninprod = ninprod + v1(i) * v2(i)
      end do
      return
      end
!=============================================================================!
      function c1inprod(n, v1, v2)
!=============================================================================!
!
!.... Perform and inner product on two complex vectors, v1 and v2
!.... using the conjugate of the first vector
!
!=============================================================================!
      implicit none
      integer :: n
      complex :: v1(n), v2(n)
      complex :: c1inprod
      integer :: i
!=============================================================================!
      c1inprod = 0.0
      do i = 1, n
        c1inprod = c1inprod + conjg(v1(i)) * v2(i)
      end do
      return
      end
!=============================================================================!
      function cinprod(n, v1, v2)
!=============================================================================!
!
!.... Perform and inner product on two complex vectors, v1 and v2
!.... using the conjugate of the second vector
!
!.... Note that this is a slightly weird definition of the
!.... complex inner-product
!=============================================================================!
      implicit none
      integer :: n
      complex :: v1(n), v2(n)
      complex :: cinprod
      integer :: i
!=============================================================================!
      cinprod = 0.0
      do i = 1, n
        cinprod = cinprod + v1(i) * conjg(v2(i))
      end do
      return
      end
