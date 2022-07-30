!=============================================================================!
      module global
!
!     Global variables for shoot
!
!     Revised:  9-25-96
!=============================================================================!

!.... problem dimensions

      integer :: ny
      integer :: ndof = 5, neq = 8, nsd = 3
      real    :: ymax
      integer :: maxIter = 50
      logical :: useParallel = .true.
      logical :: outputNonparallel = .true.

!.... flags

      integer :: ievec, iver, npcalc
      logical :: efun_out = .true.

!.... solution storage

      complex, allocatable :: efun(:,:), adj(:,:)

!.... flow parameters

      real    :: Ma, Re, Pr, T0
      complex :: alpha, beta, omega
      real    :: alphar, alphai, betar, betai, omegar, omegai
      real    :: sl
      integer :: curve=1

!.... flow constants

      real, parameter :: gamma  = 1.4
      real, parameter :: gamma1 = 0.4
      real, parameter :: cv     = 716.5
      real, parameter :: cp     = 1003.1
      real, parameter :: Rgas   = 286.6

      integer :: mattyp
      real    :: datmat(3)

!.... edge conditions

      real :: Te, rmue, dmue, d2mue, rlme, dlme, d2lme, cone
      real :: dcone, d2cone

!.... Useful constants

      real, parameter    :: zero    = 0.0000000000000000000d+00
      real, parameter    :: eps12   = 1.0000000000000000000d-12
      real, parameter    :: eps8    = 1.0000000000000000000d-08
      real, parameter    :: pt25    = 2.5000000000000000000d-01
      real, parameter    :: pt33    = 3.3333333333333333333d-01
      real, parameter    :: pt5     = 5.0000000000000000000d-01
      real, parameter    :: pt66    = 6.6666666666666666666d-01
      real, parameter    :: one     = 1.0000000000000000000d+00
      real, parameter    :: onept25 = 1.2500000000000000000d+00
      real, parameter    :: onept33 = 1.3333333333333333333d+00
      real, parameter    :: onept5  = 1.5000000000000000000d+00
      real, parameter    :: two     = 2.0000000000000000000d+00
      real, parameter    :: three   = 3.0000000000000000000d+00
      real, parameter    :: four    = 4.0000000000000000000d+00
      real, parameter    :: pi      = 3.1415926535897932385d+00
      complex, parameter :: im      = (0.0,1.0)

      end module global

!=============================================================================!
      module material
!
!     Interfaces for material routines
!
!     Revised:  4-16-96
!=============================================================================!
      interface getmat
        subroutine getmat(t, mu, lm, con, &
                          dmu, d2mu, dlm, d2lm, dcon, d2con)
          real t(:), mu(:), lm(:), con(:), dmu(:), d2mu(:)
          real dlm(:), d2lm(:), dcon(:), d2con(:)
        end subroutine getmat
        subroutine sgetmat(t, mu, lm, con, &
                           dmu, d2mu, dlm, d2lm, dcon, d2con)
          real t, mu, lm, con, dmu, d2mu, dlm, d2lm
          real dcon, d2con
        end subroutine sgetmat
      end interface

      end module material
