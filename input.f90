!==============================================================================
subroutine input
!==============================================================================
  use global
  use material
  implicit none

  namelist /inputparam/ mattyp, Ma, Re, Pr, ny, ymax, ievec, alphar, alphai, &
                        betar, betai, omegar, omegai, iver, sl, npcalc, &
                        efun_out
!==============================================================================

  write (*,"('SHOOT Compressible Flow Linear Stability Solver')")

!.... Constant mu or Sutherland's law

  write (*,"('(0) for constant Mu, (1) for Sutherland ==> ',$)")
  read (*,*) mattyp

  if (mattyp .eq. 1) then   ! for Sutherland's law (AIR)

!.... get the freestream stagnation temperature

    write (*,"('Enter the freestream T0 (K) ==> ',$)")
    read (*,*) T0
    te = t0 / ( one + pt5 * gamma1 * Ma**2 )

!.... set fluid properties

    datmat(1) = 1.715336725523065e-05
    datmat(2) = 273.0
    datmat(3) = 110.4
  else                                    ! constant viscosity
    datmat(1) = one
    datmat(2) = zero
    datmat(3) = zero
  end if

  write (*,"('Enter Ma, Re, Pr ==> ',$)")
  read (*,*) Ma, Re, Pr

!.... get the fluid properties at the reference state

  call getmat(te, rmue, rlme, cone, dmue, d2mue, &
              dlme, d2lme, dcone, d2cone)

  write (*,"('Enter ny, Ymax ==> ',$)")
  read (*,*) ny, ymax
  write (*,"('Compute eigenfuntions (0/1) ==> ',$)")
  read (*,*) ievec

  write (*,"('Enter alphar, alphai ==> ',$)")
  read (*,*) alphar, alphai
  write (*,"('Enter betar, betai ==> ',$)")
  read (*,*) betar, betai
  write (*,"('Enter omegar, omegai ==> ',$)")
  read (*,*) omegar, omegai

!.... set alpha and beta and omega

  alpha = alphar + im * alphai
  beta  = betar  + im * betai
  omega = omegar + im * omegai

!.... These two variables are now read from the parm.dat file
!.... along with alphar, alphai

  write (*,"('Enter profile index ==> ',$)")
  read (*,*) iver

  write(*,"('Enter s (-1 = no curvature) ==> ',$)")
  read(*,*) sl

!.... npcalc = 0   -->  polish eigenvalues, solve adjoint, compute nonparallel
!.... npcalc = 1   -->  compute nonparallel terms assuming output.dat is avail
!.... npcalc = 2   -->  polish, solve adjoint, write output.dat

  write(*,"('Enter npcalc ==> ',$)")
  read(*,*) npcalc

!.... output efunctions

  write(*,"('Enter efun_out ==> ',$)")
  read(*,*) efun_out 

!.... echo input using Namelist format

  write(*,*)
  write(*,NML=inputparam)
  write(*,*)

  return
  end subroutine input
