!=============================================================================!
      subroutine adjder( neq, Q, y, dQ )
!=============================================================================!
      implicit none
      
      integer :: neq
      real    :: y
      complex :: Q(neq), dQ(neq)
      
!.... local variables

      complex :: Eh(neq,neq), Fh(neq,neq), Ehinv(neq,neq)

!.... local variables for Lapack routines

      integer :: ipiv(neq), info
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
      
      Fh = transpose( Fh )
      dQ = matmul( Fh, Q )
      
      return
      end
