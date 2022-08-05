#include "port.h"
!=============================================================================!
      subroutine initial( Uic, ic )
!=============================================================================!
      use global
      implicit none

      complex :: Uic(neq,neq)
      integer :: ic, j

!.... local variables

      complex :: Eh(neq,neq), Fh(neq,neq), Ehinv(neq,neq), Id(neq,neq)
      complex :: eval(neq), evec(neq,neq)
      integer :: ieq

!.... local variables for Lapack routines

      integer :: ipiv(neq), jpiv(neq), info
      integer :: lwork
      complex, allocatable :: work(:)
      real, allocatable    :: rwork(:)
!=============================================================================!
#define USE_LOCAL_INVERSE
!#define DEBUG_MATRICES

      lwork = 2 * neq
      allocate( work(lwork), rwork(lwork) )

!.... Initialize data-structures

      Eh=zero; Fh=zero; Ehinv=zero; eval=zero; evec=zero; ipiv=zero
      jpiv=zero; work=zero; rwork=zero

!.... Compute the matrices in the farfield

      call parallel( ymax, Eh, Fh, 0 )
#ifdef DEBUG_MATRICES
      write(*,*) "Eh = "
      do j=1,neq
       write(*,"(*('('sf8.2xspf8.2x'i)':x))") Eh(:,j)
      end do
#endif

!.... Multiply Fh by Eh^{-1}

#ifdef USE_LOCAL_INVERSE
      call inverse( neq, Eh, Ehinv)
#ifdef DEBUG_MATRICES
      Id = matmul(Ehinv, Eh) 
      write(*,*) "Ehinv * Eh = "
      do j=1,neq
       write(*,"(*('('sf8.2xspf8.2x'i)':x))") Id(:,j)
      end do
      !write(*,*) "Fh(1) = ", Fh
      !write(*,*) "Ehinv = ", Ehinv
#endif
      Fh = matmul( Ehinv, Fh )
#else  
#ifdef PARTIAL_PIVOTING
      call CGETRF( neq, neq, Eh, neq, ipiv, info)
#else
      call ZGETC2( neq, Eh, neq, ipiv, jpiv, info)
#endif
      if (info.ne.0) then
        write(*,*) 'CGETRF/ZGETC2: ', info
        call exit(info)
      endif
      call CGETRS('N', neq, neq, Eh, neq, ipiv, Fh, neq, info)
      if (info.ne.0) then
        write(*,*) 'CGETRS: ',info
        call exit(info)
      endif
#endif

!.... Negate the matrix to get the correct eigenvalues

      Fh = -Fh
#ifdef DEBUG_MATRICES
      write(*,*) "-Ehinv * Fh = "
      do j=1,neq
       write(*,"(*('('sf8.2xspf8.2x'i)':x))") Fh(:,j)
      end do
#endif

!.... solve the eigenproblem

      call CGEEV('N', 'V', neq, Fh, neq, eval, evec, &
                 neq, evec, neq, work, lwork, rwork, info)

!.... use the solutions that are damped to infinity to form the initial

      Uic = zero
      ic = 0
      do ieq = 1, neq
!        write(*,*) ieq, eval(ieq), real(eval(ieq))
        if ( real(eval(ieq)) .lt. zero ) then
          ic = ic + 1
          Uic(:,ic) = evec(:,ieq) * exp( eval(ieq) * ymax )
        end if
      end do

      if (ic.ne.4) then
        write(*,*) "initial:  there should be 4 damped modes"
        write(*,*) "          but we have ic = ", ic," dampled modes"
        call exit(1)
      endif
#ifdef VERBOSE
      write(*,"('Completed initial')")
#endif
      deallocate( work, rwork )
      return
      end
