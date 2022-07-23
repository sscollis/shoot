!=============================================================================!
        subroutine parallel( y, Eh, Fh, iflag, Ah)
!
!       Compute the parallel flow matrices
!
!       Note:  you can give this routine a profile that is not parallel and
!              it will still compute the parallel matrices...get it?
!
!       iflag controls what is returned by the routine
!
!       iflag = 0     return only Eh, Fh
!       iflag = 1     return Eh, Fh, Ah
!
!       The iflag = 1 is only used (right now) to compute the nonparallel terms
!
!       S. Scott Collis
!
!       Revised: 2-24-97
!=============================================================================!
        use global
        use material
        implicit none

        integer :: iflag
        real    :: y
        complex :: Eh(neq,neq), Fh(neq,neq), Ah(neq,neq)

        real    :: vm(ndof)
        real    :: A(ndof,ndof), C(ndof,ndof)
        real    :: D(ndof,ndof), G(ndof,ndof)
        complex :: B(ndof,ndof)

        real :: Vxx(ndof,ndof), Vxy(ndof,ndof), Vyy(ndof,ndof)
        real :: Vxz(ndof,ndof), Vyz(ndof,ndof), Vzz(ndof,ndof)

        real :: g1vm(ndof), g2vm(ndof), g3vm(ndof), g22vm(ndof)

        real :: g2rhom, g2u1m, g2u3m, g2tm

        real :: rhom, u1m, u3m, tm

        real :: rmu,   dmu,   d2mu
        real :: rlm,   dlm,   d2lm
        real :: con,   dcon,  d2con
        real :: g2mu,  g2lm,  g2con
        real :: g2dmu, g2dlm, g2dcon

        real :: S1jj, S2jj, S3jj, S12, S22, S23, Lapt

        real :: fact

!.... metrics for the curved wall

        real :: h, dhds, dhdr, dhdsr, dhdrr, hinv, hinv2, hinv3
!=============================================================================!

        call getmeanp( iver, y, vm, g2vm, g22vm )

!.... compute the curvature metrics

        call calch( sl, 1, y, h, dhds, dhdr, dhdsr, dhdrr )

        hinv  = one / h
        hinv2 = hinv**2
        hinv3 = hinv**3

!.... some abbreviations

        rhom = vm(1)
        u1m  = vm(2)
        u3m  = vm(4)
        tm   = vm(5)

!.... initialize derivatives

        g2rhom = g2vm(1)
        g2u1m  = g2vm(2)
        g2u3m  = g2vm(4)
        g2tm   = g2vm(5)

!.... compute some stuff that is useful for the viscous terms (curve*)

        S12 = pt5 * ( (-u1m * dhdr) * hinv + g2u1m )
        S22 = zero
        S23 = pt5 * g2u3m

        S1jj = -pt5 * ( dhdr**2 + dhdrr * h ) * hinv2 * u1m + &
               pt5 * dhdr * g2u1m * hinv + pt5 * g22vm(2)

        S2jj = pt5 * ( dhdr * dhds - h * dhdsr ) * hinv3 * u1m

        S3jj = pt5 * dhdr * g2u3m * hinv

!.... Laplacian of Tm (curve*)

        LapT = g22vm(5) + g2tm * dhdr * hinv

!.... compute mean material properties

        call getmat(tm*te, rmu, rlm, con, dmu, d2mu, dlm, d2lm, dcon, d2con)

!.... nondimensionalize

        rmu   = rmu / rmue
        dmu   = dmu * Te / rmue
        d2mu  = d2mu * Te**2 / rmue

        con   = con / cone
        dcon  = dcon * Te / cone
        d2con = d2con * Te**2 / cone

        rlm   = rlm / rlme
        dlm   = dlm * Te / rlme
        d2lm  = d2lm * Te**2 / rlme

!.... compute gradients of viscosity using chain-rule

        g2mu  = dmu * g2tm
        g2dmu = d2mu * g2tm

!.... compute gradients of conductivity using chain-rule

        g2con  = dcon * g2tm
        g2dcon = d2con * g2tm

!.... compute gradients of bulk viscosity using chain-rule

        g2lm  = dlm * g2tm
        g2dlm = d2lm * g2tm

!==============================================================================
        G   = zero
        A   = zero
        B   = zero
        C   = zero
        D   = zero
        Vxx = zero
        Vxy = zero
        Vyy = zero
        Vxz = zero
        Vyz = zero
        Vzz = zero

!.... Continuity equation (curve*)

        G(1,1) = one

        A(1,1) = u1m * hinv
        A(1,2) = rhom * hinv

        B(1,3) = rhom

        C(1,1) = u3m
        C(1,4) = rhom

        D(1,3) = g2rhom + rhom * dhdr * hinv

!.... Momentum equation -- x_1 (convective + pressure) (curve*)

        G(2,2) = rhom

        A(2,1) = tm/(h * gamma * Ma**2)
        A(2,2) = rhom * u1m * hinv
        A(2,5) = rhom/(h * gamma * Ma**2)

        C(2,2) = rhom * u3m

        D(2,2) = zero
        D(2,3) = rhom * ( g2u1m + u1m * dhdr * hinv )

!.... (viscous lambda) (curve*)

        fact = rlme / (rmue * Re)

        A(2,2) = A(2,2) - fact * ( -rlm * hinv3 * dhds )
        A(2,3) = A(2,3) - fact * rlm * hinv2 * dhdr

        D(2,3) = D(2,3) - fact * ( -rlm * hinv3 * dhds * dhdr + &
                                    rlm * hinv2 * dhdsr )

        Vxx(2,2) = fact * rlm * hinv2

        Vxy(2,3) = fact * rlm * hinv

        Vxz(2,4) = fact * rlm * hinv

!.... (viscous mu) (curve*)

        fact = one / Re

        A(2,2) = A(2,2) - fact * ( -two * rmu * dhds * hinv3 )
        A(2,3) = A(2,3) - fact * ( g2mu * hinv + &
                              rmu * 3.0 * dhdr * hinv2 )

        B(2,2) = B(2,2) - fact * ( g2mu + rmu * dhdr * hinv )
        B(2,5) = B(2,5) - fact * dmu * two * S12

        D(2,2) = D(2,2) - fact * ( g2mu * hinv * (-dhdr) - &
                              rmu * ( dhdr**2 + dhdrr * h ) * hinv2 )
        D(2,3) = D(2,3) - fact * ( &
                              two * rmu * ( dhdsr * h - dhds * dhdr ) * hinv3 )
        D(2,5) = D(2,5) - fact * two * ( &
                              g2dmu * S12 + &
                              dmu * S1jj )

        Vxx(2,2) = Vxx(2,2) + fact * two * rmu * hinv2
        Vxy(2,3) = Vxy(2,3) + fact * rmu * hinv
        Vyy(2,2) = Vyy(2,2) + fact * rmu
        Vxz(2,4) = Vxz(2,4) + fact * rmu * hinv
        Vzz(2,2) = Vzz(2,2) + fact * rmu

!.... Momentum equation -- x_2 (convective + pressure) (curve*)

        G(3,3) = rhom

        A(3,3) = rhom * u1m * hinv

        B(3,1) = tm/(gamma * Ma**2)
        B(3,3) = zero
        B(3,5) = rhom/(gamma * Ma**2)

        C(3,3) = rhom * u3m

        D(3,1) = u1m * hinv * ( -u1m * dhdr ) + &
                 g2tm / (gamma * Ma**2)
        D(3,2) = rhom * ( -two * u1m * dhdr ) * hinv
        D(3,5) = g2rhom / (gamma * Ma**2)

!.... (viscous lambda) (curve*)

        fact = rlme / (rmue * Re)

        A(3,2) = A(3,2) - fact * ( g2lm * hinv - rlm * dhdr * hinv2 )

        B(3,3) = B(3,3) - fact * ( g2lm + rlm * dhdr * hinv )

        C(3,4) = C(3,4) - fact * ( g2lm )

        D(3,3) = D(3,3) - fact * ( g2lm * hinv * dhdr - &
                              rlm * dhdr * hinv2 * dhdr + &
                              rlm * hinv * dhdrr )

        Vxy(3,2) = fact * rlm * hinv

        Vyy(3,3) = fact * rlm

        Vyz(3,4) = fact * rlm

!.... (viscous mu) (curve*)

        fact = one / Re

        A(3,2) = A(3,2) + fact * rmu * 3.0 * dhdr * hinv2
        A(3,3) = A(3,3) - fact * ( -rmu * dhds * hinv3 )
        A(3,5) = A(3,5) - fact * dmu * two * S12 * hinv

        B(3,3) = B(3,3) - fact * ( two * g2mu + two * rmu * dhdr * hinv )
        B(3,5) = B(3,5) - fact * dmu * two * S22

        C(3,5) = C(3,5) - fact * dmu * two * S23

        D(3,2) = D(3,2) - fact * ( &
                              rmu * (dhds * dhdr - h * dhdsr) * hinv3 )
        D(3,3) = D(3,3) + fact * two * rmu * dhdr**2 * hinv2
        D(3,5) = D(3,5) - fact * two * ( g2dmu * S22 + dmu * S2jj )

        Vxx(3,3) = Vxx(3,3) + fact * rmu * hinv2
        Vxy(3,2) = Vxy(3,2) + fact * rmu * hinv
        Vyy(3,3) = Vyy(3,3) + fact * two * rmu
        Vyz(3,4) = Vyz(3,4) + fact * rmu
        Vzz(3,3) = Vzz(3,3) + fact * rmu

!.... now use continuity to remove Vyy(:,3,3) term

        fact = Vyy(3,3) / rhom

        A(3,1) = A(3,1) + fact * ( g2u1m * hinv - u1m * dhdr * hinv2 )
        A(3,2) = A(3,2) + fact * ( g2rhom * hinv - rhom * dhdr * hinv2 )

        B(3,1) = B(3,1) + fact * ( -im * omega )
        B(3,3) = B(3,3) + fact * ( two * g2rhom + rhom * dhdr * hinv )

        C(3,1) = C(3,1) + fact * ( g2u3m )
        C(3,4) = C(3,4) + fact * ( g2rhom )

        D(3,3) = D(3,3) + fact * ( g2rhom * dhdr * hinv - &
                              rhom * dhdr**2 * hinv2 + rhom * dhdrr * hinv + &
                              g22vm(1) )

        Vxy(3,1) = Vxy(3,1) - fact * u1m * hinv
        Vxy(3,2) = Vxy(3,2) - fact * rhom * hinv
        Vyz(3,1) = Vyz(3,1) - fact * u3m
        Vyz(3,4) = Vyz(3,4) - fact * rhom

!.... zero the offending term

        Vyy(3,3) = zero

!.... Momentum equation -- x_3 (convective + pressure) (curve*)

        G(4,4) = rhom

        A(4,4) = rhom * u1m * hinv

        C(4,1) = tm/(gamma * Ma**2)
        C(4,4) = rhom * u3m
        C(4,5) = rhom/(gamma * Ma**2)

        D(4,3) = rhom * g2u3m

!.... (viscous lambda) (curve*)

        fact = rlme / (rmue * Re)

        C(4,3) = C(4,3) - fact * rlm * hinv * dhdr

        D(4,3) = D(4,3) - fact * ( g2lm * hinv * dhdr )

        Vxz(4,2) = fact * rlm * hinv

        Vyz(4,3) = fact * rlm

        Vzz(4,4) = fact * rlm

!.... (viscous mu) (curve*)

        fact = one / Re

        A(4,4) = A(4,4) - fact * ( -rmu * dhds * hinv3 )

        B(4,4) = B(4,4) - fact * ( g2mu + rmu * dhdr * hinv )
        B(4,5) = B(4,5) - fact * dmu * two * S23

        C(4,3) = C(4,3) - fact * ( g2mu + rmu * dhdr * hinv )

        D(4,5) = D(4,5) - fact * two * ( g2dmu * S23 + dmu * S3jj )

        Vxx(4,4) = Vxx(4,4) + fact * rmu * hinv2
        Vyy(4,4) = Vyy(4,4) + fact * rmu
        Vxz(4,2) = Vxz(4,2) + fact * rmu * hinv
        Vyz(4,3) = Vyz(4,3) + fact * rmu
        Vzz(4,4) = Vzz(4,4) + fact * two * rmu

!.... Energy equation (Advection + pressure) (curve*)

        G(5,5) = rhom

        A(5,2) = rhom * gamma1 * tm * hinv
        A(5,5) = rhom * u1m * hinv

        B(5,3) = rhom * gamma1 * tm

        C(5,4) = rhom * gamma1 * tm
        C(5,5) = rhom * u3m

        D(5,3) = rhom * g2tm + rhom * gamma1 * tm * dhdr * hinv

!.... diffusion (curve*)

        fact = gamma / (Pr * Re)

        A(5,5) = A(5,5) - fact * ( -con * dhds * hinv3 )
        B(5,5) = B(5,5) - fact * (g2con + dcon * g2tm + &
                              con * dhdr * hinv )
        D(5,5) = D(5,5) - fact * ( g2dcon * g2tm + dcon * Lapt )

        Vxx(5,5) = fact * con * hinv2
        Vyy(5,5) = fact * con
        Vzz(5,5) = fact * con

!.... dissipation (lambda) (curve*)

        fact = gamma * gamma1 * Ma**2 * rlme / (Re * rmue)

!.... dissipation (mu) (curve*)

        fact = gamma * gamma1 * Ma**2 / Re

        A(5,3) = A(5,3) - fact * four * rmu * S12 * hinv

        B(5,2) = B(5,2) - fact * four * rmu * S12
        B(5,3) = B(5,3) - fact * four * rmu * S22
        B(5,4) = B(5,4) - fact * four * rmu * S23

        C(5,3) = C(5,3) - fact * four * rmu * S23

        D(5,2) = D(5,2) + fact * four * rmu * S12 * dhdr * hinv
        D(5,5) = D(5,5) - fact * two * dmu * ( &
                       S12**2 + S12**2 + &
                       S22**2 + S23**2 + &
                       S23**2 )

!==============================================================================
!.... form the extended matrices
!==============================================================================
        Eh = zero
        Fh = zero

        Eh(1:5,1:5) = B - im * alpha * Vxy - im * beta * Vyz
        Eh(2,6)     = -Vyy(2,2)
        Eh(4,7)     = -Vyy(4,4)
        Eh(5,8)     = -Vyy(5,5)
        Eh(6,2)     = one
        Eh(7,4)     = one
        Eh(8,5)     = one

        Fh(1:5,1:5) = D + im * alpha * A + im * beta * C + &
                      alpha**2 * Vxx + alpha * beta * Vxz + &
                      beta**2 * Vzz - im * omega * G
        Fh(6,6)     = -one
        Fh(7,7)     = -one
        Fh(8,8)     = -one

!.... terms needed for the nonparallel correction

        if (iflag .eq. 1) then
          Ah = zero
          Ah(1:5,1:5) = two * im * alpha * Vxx + im * beta * Vxz - A
        end if

        return
        end
