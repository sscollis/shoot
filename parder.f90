!=============================================================================!
      subroutine parder( neq, Q, y, dQ )
!=============================================================================!
      implicit none
      
      integer :: neq
      real    :: y
      complex :: Q(neq), dQ(neq)
      
!.... local variables

      complex :: Eh(neq,neq), Fh(neq,neq), Ehinv(neq,neq)

!.... local variables for Lapack routines

      integer :: ipiv(neq), info

      real :: cpu
      real, external :: second
!=============================================================================!
      call parallel( y, Eh, Fh, 0 ) 

!.... Multiply Fh by Eh^{-1}

      call inverse( neq, Eh, Ehinv)
      Fh = matmul( Ehinv, Fh )

!     call CGETRF( neq, neq, Eh, neq, ipiv, info)
!     if (info.ne.0) write(*,*) 'CGETRF: ', info
!     call CGETRS('N', neq, neq, Eh, neq, ipiv, Fh, neq, info)
!     if (info.ne.0) write(*,*) 'CGETRS: ',info
      
!.... Form the derivative
      
      dQ = -matmul( Fh, Q )
      
      return
      end

!=============================================================================!
      subroutine inverse( neq, Eh, Ehinv )
!
!     Hardwired inverse for Eh.  Note that this takes advantage of the 
!     sparce structure of Eh.  If, for some wierd reason, the definition
!     of Eh changes then this routine WILL produce incorrect results.
!
!=============================================================================!
      implicit none
      integer :: neq
      complex :: Eh(neq,neq), Ehinv(neq,neq), fact
      real, parameter :: zero = 0.0, one = 1.0
!=============================================================================!
      fact = one / ( Eh(3,1) * Eh(4,7) )
      Ehinv(1,1) = ( Eh(3,7) * Eh(4,3) - Eh(3,3) * Eh(4,7) ) / &
                   ( Eh(1,3) * Eh(3,1) * Eh(4,7) )
      Ehinv(1,2) = zero
      Ehinv(1,3) = one / Eh(3,1)
      Ehinv(1,4) = -Eh(3,7) * fact
      Ehinv(1,5) = zero
      Ehinv(1,6) = -Eh(3,2) / Eh(3,1)
      Ehinv(1,7) = ( Eh(3,7) * Eh(4,4) - Eh(3,4) * Eh(4,7) ) * fact
      Ehinv(1,8) = ( Eh(3,7) * Eh(4,5) - Eh(3,5) * Eh(4,7) ) * fact

      Ehinv(2,1) = zero
      Ehinv(2,2) = zero
      Ehinv(2,3) = zero
      Ehinv(2,4) = zero
      Ehinv(2,5) = zero
      Ehinv(2,6) = one
      Ehinv(2,7) = zero
      Ehinv(2,8) = zero

      Ehinv(3,1) = one / Eh(1,3)
      Ehinv(3,2) = zero
      Ehinv(3,3) = zero
      Ehinv(3,4) = zero
      Ehinv(3,5) = zero
      Ehinv(3,6) = zero
      Ehinv(3,7) = zero
      Ehinv(3,8) = zero

      Ehinv(4,1) = zero
      Ehinv(4,2) = zero
      Ehinv(4,3) = zero
      Ehinv(4,4) = zero
      Ehinv(4,5) = zero
      Ehinv(4,6) = zero
      Ehinv(4,7) = one
      Ehinv(4,8) = zero

      Ehinv(5,1) = zero
      Ehinv(5,2) = zero
      Ehinv(5,3) = zero
      Ehinv(5,4) = zero
      Ehinv(5,5) = zero
      Ehinv(5,6) = zero
      Ehinv(5,7) = zero
      Ehinv(5,8) = one

      fact = one / Eh(2,6)
      Ehinv(6,1) = -Eh(2,3) / ( Eh(1,3) * Eh(2,6) )
      Ehinv(6,2) = fact
      Ehinv(6,3) = zero
      Ehinv(6,4) = zero
      Ehinv(6,5) = zero
      Ehinv(6,6) = -Eh(2,2) * fact
      Ehinv(6,7) = zero
      Ehinv(6,8) = -Eh(2,5) * fact

      fact = one / Eh(4,7)
      Ehinv(7,1) = -Eh(4,3) / ( Eh(1,3) * Eh(4,7) )
      Ehinv(7,2) = zero
      Ehinv(7,3) = zero
      Ehinv(7,4) = fact
      Ehinv(7,5) = zero
      Ehinv(7,6) = zero
      Ehinv(7,7) = -Eh(4,4) * fact
      Ehinv(7,8) = -Eh(4,5) * fact

      fact = one / Eh(5,8)
      Ehinv(8,1) = -Eh(5,3) / ( Eh(1,3) * Eh(5,8) )
      Ehinv(8,2) = zero
      Ehinv(8,3) = zero
      Ehinv(8,4) = zero
      Ehinv(8,5) = fact
      Ehinv(8,6) = -Eh(5,2) * fact
      Ehinv(8,7) = -Eh(5,4) * fact
      Ehinv(8,8) = -Eh(5,5) * fact

      return
      end
