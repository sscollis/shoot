#include "port.h"
!=============================================================================!
       subroutine solve
!=============================================================================!
       use global
       implicit none

       integer :: ic, i, j
       complex :: Ui(neq,neq), Uf(neq,neq), BC(neq), norm

       complex, allocatable :: A(:,:)

       integer :: icount, mcount=20, ind

       real :: tol = 0.2, y, dy

!.... eigenvalue iteration variables

       complex :: ctemp, cm1, cm2, err, errm1, errm2
       complex :: At, Bt, Ct, qt, fd
       real    :: fdxr, fdxi, fdyr, fdyi

       external parder

!.... For Lapack eigensolver

       integer :: info
       integer :: lwork
       complex, allocatable :: work(:), eval(:), evec(:,:)
       real, allocatable    :: rwork(:)
!=============================================================================!
       mcount = maxIter

       BC = zero

!.... compute the initial condition

       call initial( Ui, ic )

#if VERBOSE>=3
       write(*,*) "Ui = ", Ui(:,:)
#endif

       lwork = 2 * ic
       allocate( A(ic,ic), work(lwork), rwork(lwork), eval(ic), &
                 evec(ic,ic) )

!.... Begin the eigenvalue iteration loop

       ievec  = 0          ! don't compute efun
       icount = 0          ! number of iteration
       err    = one

100    continue
         icount = icount + 1

!.... integrate the equations using orthonomalization
#if VERBOSE>=3
         write(*,*) "Starting Conte integration"
#endif
         call conte( ny-1, tol, neq, ic, Ui, Uf, ymax, zero, ievec, &
                     parder, BC, efun )
#if VERBOSE>=3
         write(*,*) "Uf = ", Uf(:,:)
#endif

!.... form the boundary condition matrix for an isothermal wall

         if (ic.ne.4) then
           write(*,*) "solve:  dimensions of A don't match Uf range"
           write(*,*) "ic = ",ic," which is not 4"
           call exit(1)
         endif
#if 1
         A(:,:) = Uf(2:5,1:ic)
#else
         A(:,:) = 0
         write(*,*) "ic = ", ic
         A(2:5,1:ic) = Uf(2:5,1:ic)
#endif

!.... compute the eigenvalues (only) of A and select the minimum eval

#if VERBOSE>=3
         write(*,*) "Calling [CZ]EEV: ic = ",ic
         write(*,*) "A = ", A(:,:)
         write(*,*) "CGEEV"
#endif
         call CGEEV('N', 'N', ic, A, ic, eval, evec, &
                    ic, evec, ic, work, lwork, rwork, info)
#if VERBOSE>=3
         write(*,*) "Finished [CZ]EEV..."
#endif
         err = eval(1)
         ind = 1
         do i = 2, ic
           if ( abs(eval(i)) .le. abs(err) ) then
             err = eval(i)
             ind = i
           end if
         end do

         write (*,30) icount, ind, real(alpha), aimag(alpha), &
                      real(err), aimag(err), abs(err)
30       format (1x,i2,1x,i2,2(2x,1pe13.6,1x,1pe13.6),2x,1pe13.6)

         if ( (abs(err) .ge. eps8) .and. (icount .lt. mcount) ) then

!.... compute a new guess for the eigenvalue
!....   1.  Just make a perturbation on the first iteration
!....   2.  A Newton correction
!....   3.  A quadratic correction

         ctemp = alpha
         if (icount .eq. 1) then
           cm2 = alpha
           errm2 = err
           alpha = alpha * 1.00000001
         else if (icount .eq. 2) then
           cm1 = alpha
           errm1 = err
           fd = (err-errm2)/(alpha-cm2)
           alpha = alpha - err/fd
         else
           qt = (alpha-cm1)/(cm1-cm2)
           At = qt*err-qt*(one+qt)*errm1+qt**2*errm2
           Bt = (two*qt+one)*err-(one+qt)**2*errm1+qt**2*errm2
           Ct = (one+qt)*err
           if ( ABS(Bt+SQRT(Bt**2-four*At*Ct)) .gt. &
                ABS(Bt-SQRT(Bt**2-four*At*Ct)) )  then
             alpha = ctemp-(ctemp-cm1)*two*Ct/(Bt+SQRT(Bt**2-four*At*Ct))
           else
             alpha = ctemp-(ctemp-cm1)*two*Ct/(Bt-SQRT(Bt**2-four*At*Ct))
           end if
           cm2 = cm1
           cm1 = ctemp
           errm2 = errm1
           errm1 = err
         end if
         goto 100
       endif

       write(*,*)
       write(*,"('alpha = ',1pe20.13,1x,1pe20.13)") real(alpha), &
                                                    aimag(alpha)
       write(*,*)

!.... if converged, determine the combination of the independent solutions
!.... that satisfies the boundary conditions by solving an eigensystem
!.... to determine the eigenvector cooresponding to the zero eigenvalue.

#if 1
       A(:,:) = Uf(2:5,1:ic)
#else
       A(:,:) = 0
       A(2:5,1:ic) = Uf(2:5,1:ic)
#endif

       !write(*,*) "2: CGEEV"
       call CGEEV('N', 'V', ic, A, ic, eval, evec, &
                  ic, evec, ic, work, lwork, rwork, info)

       i=1
       err = eval(1)
       BC(1:ic) = evec(:,i)
       do i = 2, ic
         if ( abs(eval(i)) .le. abs(err) ) then
           err = eval(i)
           BC(1:ic) = evec(:,i)
         end if
       end do
       !write(*,*) "Solve deallocate"
       deallocate( A, work, rwork, eval, evec )

       ievec = 1
       !write(*,*) "ny-1=",ny-1," neq=",neq," ic=",ic," ymax=",ymax
       call conte( ny-1, tol, neq, ic, Ui, Uf, ymax, zero, ievec, &
                   parder, BC, efun)

       !write(*,*) "Finished Conte"

!.... make the phase reference consistent

       j = ny/2
       efun = efun / exp( im * atan2(aimag(efun(2,j)),real(efun(2,j))))

!.... normalize and output the eigenfunction if desired

       if (efun_out) then

!.... determine the normalization factor

       norm = zero
       do i = 1, ndof
         do j = 1, ny
           if ( ABS(efun(i,j)) .gt. ABS(norm) ) norm = efun(i,j)
         end do
       end do

!.... normalize and output the eigenfunction

       efun = efun / norm
       open(15,file='efun.out')
       write(15,"('# alpha = ',1pe20.13,1x,1pe20.13)") alpha
       dy = ymax / real(ny-1)
       do j = 1, ny
         y = ymax - dy * real(j-1)
         write(15,"(17(1pe13.6,1x))") y, &
           (real(efun(i,j)), aimag(efun(i,j)), i = 1, neq )
       end do
       close(15)

       end if

       return
       end subroutine solve
