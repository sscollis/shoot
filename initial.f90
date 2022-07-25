#include "port.h"
!=============================================================================!
      subroutine initial( Uic, ic  )
!=============================================================================!
      use global
      implicit none

      complex :: Uic(neq,neq)
      integer :: ic

!.... local variables

      complex :: Eh(neq,neq), Fh(neq,neq), Ehinv(neq,neq)
      complex :: eval(neq), evec(neq,neq)
      integer :: ieq

!.... local variables for Lapack routines

      integer :: ipiv(neq), info
      integer :: lwork
      complex, allocatable :: work(:)
      real, allocatable    :: rwork(:)
!=============================================================================!
      lwork = 2 * neq
      allocate( work(lwork), rwork(lwork) )

!.... Compute the matrices in the farfield

      call parallel( ymax, Eh, Fh, 0 )

!.... Multiply Fh by Eh^{-1}

      call inverse( neq, Eh, Ehinv)
      Fh = matmul( Ehinv, Fh )

!     call CGETRF( neq, neq, Eh(:,:), neq, ipiv, info)
!     if (info.ne.0) write(*,*) 'CGETRF: ', info
!     call CGETRS('N', neq, neq, Eh(:,:), neq, ipiv, Fh(:,:), neq, info)
!     if (info.ne.0) write(*,*) 'CGETRS: ',info

!.... Negate the matrix to get the correct eigenvalues

      Fh = -Fh

!.... solve the eigenproblem

      call CGEEV('N', 'V', neq, Fh, neq, eval, evec, &
                 neq, evec, neq, work, lwork, rwork, info)

!.... use the solutions that are damped to infinity to form the initial

      Uic = zero
      ic = 0
      do ieq = 1, neq
        if ( real(eval(ieq)) .lt. zero ) then
          ic = ic + 1
          Uic(:,ic) = evec(:,ieq) * exp( eval(ieq) * ymax )
        end if
      end do
#if VERBOSE
      write(*,"('Completed initial')")
#endif
      write(*,*) "Deallocate initial"
      deallocate( work, rwork )
      return
      end
