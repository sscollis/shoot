!==============================================================================
subroutine input
!==============================================================================
  use global
  use material
  implicit none

  namelist /shoot/ mattyp, Ma, Re, Pr, T0, ny, ymax, ievec, &
                   betar, betai, omegar, omegai, npcalc, &
                   efun_out, maxIter, curve, useParallel
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

  !write (*,"('Enter alphar, alphai ==> ',$)")
  !read (*,*) alphar, alphai
  alphar = 0
  alphai = 0
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

  iver = 0
  !write (*,"('Enter profile index ==> ',$)")
  !read (*,*) iver

  sl = 0
  !write(*,"('Enter s (-1 = no curvature) ==> ',$)")
  !read(*,*) sl
  
  write(*,"('Use curvative (0/1)  ==> ',$)")
  read(*,*) curve 
  
  write(*,"('Use parallel flow (T/F)  ==> ',$)")
  read(*,*) useParallel 

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
  write(*,NML=shoot)
  write(*,*)

  return
  end subroutine input
