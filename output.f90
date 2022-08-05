!=============================================================================!
       subroutine output( iout )
!=============================================================================!
       use global
       implicit none

       integer :: iout

!.... local variables

       complex :: Ah(neq,neq), Bh(neq,neq), Ch(neq,neq)
       complex :: Eh(neq,neq), Fh(neq,neq), Ahp(neq,neq)

       complex :: Ehp(neq,neq), Fhp(neq,neq), Ehpinv(neq,neq)

       complex :: Ehn(neq,neq), Fhn(neq,neq)

       complex :: vec(neq)
       complex, external :: inprod

       complex :: g2efun(neq,ny), c1(ny), c2(ny), c3(ny)
       complex :: Z(neq), Z1(neq,ny), Z2(neq,ny), Z3(neq,ny)

       real    :: y, dy
       integer :: j, k
!=============================================================================!

       dy = ymax / real(ny-1)

!.... loop over the normal index, go backwards in y since the eigenfunctions
!.... are computed and stored in that manner.

       loop_j: do j = 1, ny

         y = ymax - dy * real(j-1)

!.... compute "parallel" matrices

         call parallel( y, Ehp, Fhp, 1, Ahp)

#ifdef OUTPUT_MATRICES
         write(*,*) "Ehp = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Ehp(:,k)
         end do
         write(*,*) "Fhp = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Fhp(:,k)
         end do
         write(*,*) "Ahp = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Ahp(:,k)
         end do
#endif

!.... compute the "full" matrices

         call matrix( y, Ah, Bh, Ch, Eh, Fh )

#ifdef OUTPUT_MATRICES
         write(*,*) "Ah = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Ah(:,k)
         end do
         write(*,*) "Bh = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Bh(:,k)
         end do
         write(*,*) "Ch = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Ch(:,k)
         end do
         write(*,*) "Eh = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Eh(:,k)
         end do
         write(*,*) "Fh = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Fh(:,k)
         end do
#endif

!.... form the "nonparallel" matrices which are just the difference
!.... of the full matrices and the parallel matrices (i.e. they are the
!.... nonparallel mean-flow terms)

         Ehn = Eh - Ehp
         Fhn = Fh - Fhp
         
#ifdef OUTPUT_MATRICES
         write(*,*) "Ehn = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Ehn(:,k)
         end do
         write(*,*) "Fhn "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Fhn(:,k)
         end do
#endif

!.... now compute the inverse of Eh parallel and multiply through

         call inverse( neq, Ehp, Ehpinv)
         Fhp = matmul( Ehpinv, Fhp )

#ifdef OUTPUT_MATRICES
         write(*,*) "Ehpinv * Fh = "
         do k=1,neq
           write(*,"(*('('sf8.2xspf8.2x'i)':x))") Fhp(:,k)
         end do
#endif

!.... use the parallel equations to compute the derivative of efun wrt y

         g2efun(:,j) = -matmul( Fhp, efun(:,j) )

!.... Z is the adjoint times the inverse Eh parallel matrix

         Z = matmul( adj(:,j), Ehpinv )

!==============================================================================
!                      C 1    C o m p u t a t i o n
!==============================================================================

!.... c1 is the part of h1 that doesn't require g1efun
!.... note that the full matrices are used here

         vec   = matmul( Ah, efun(:,j) ) + matmul( Ch, g2efun(:,j) )
         c1(j) = inprod( neq, Z, vec )

!.... Z1 is the part of h1 that must be dotted with g1efun (zero)

         Z1(:,j) = zero
!==============================================================================
!                      C 2    C o m p u t a t i o n
!==============================================================================

!.... c2 is the part of h2 that doesn't require g1efun or g12efun but
!.... does require g1alpha

         vec   = -im * matmul( Bh, efun(:,j) )
         c2(j) = inprod( neq, Z, vec )

!.... c3 is the part of h2 that doesn't require g1efun, g12efun, or g1alpha
!.... These terms are due to the nonparallel meanflow:  -L1(\hat U_0)

         vec = -one * ( matmul( Ehn, g2efun(:,j) ) + matmul( Fhn, efun(:,j) ) )
         c3(j) = inprod( neq, Z, vec )

         Z2(:,j) = matmul( Z, Ah )

!.... z3 is the part of h2 that must be dotted with g12efun

         Z3(:,j) = matmul( Z, Ch )

       end do loop_j

!.... output the results to an unformatted file for processing later on

       write(iout) iver, sl, alpha, beta, omega, efun, g2efun, adj, &
                   c1, c2, c3, Z1, Z2, Z3

       return
     end subroutine output
