!=============================================================================!
        subroutine matrix( y, Ah, Bh, Ch, Eh, Fh )
!
!       Compute the full flow matrices.
!
!       S. Scott Collis
!
!       Revised: 9-30-96
!=============================================================================!
        use global
        use material
        implicit none

        real :: y
        complex :: Ah(neq,neq), Bh(neq,neq), Ch(neq,neq), &
             Eh(neq,neq), Fh(neq,neq)

        real :: vm(ndof)
        real :: A(ndof,ndof), C(ndof,ndof)
        real :: D(ndof,ndof), G(ndof,ndof)
        complex :: B(ndof,ndof)

        real :: Vxx(ndof,ndof), Vxy(ndof,ndof), Vyy(ndof,ndof)
        real :: Vxz(ndof,ndof), Vyz(ndof,ndof), Vzz(ndof,ndof)

        real :: g1vm(ndof),  g2vm(ndof),  g3vm(ndof)
        real :: g11vm(ndof), g12vm(ndof), g13vm(ndof)
        real :: g22vm(ndof), g23vm(ndof), g33vm(ndof)

        real :: gum(nsd,nsd), grhom(nsd), gtm(nsd)
        real :: divum, g1divum, g2divum, g3divum

        real :: rhom, u1m, u2m, u3m, tm

        real :: rmu,    dmu,    d2mu
        real :: rlm,    dlm,    d2lm
        real :: con,    dcon,   d2con
        real :: g1mu,   g2mu,   g3mu
        real :: g1lm,   g2lm,   g3lm
        real :: g1con,  g2con,  g3con
        real :: g1dmu,  g2dmu,  g3dmu
        real :: g1dlm,  g2dlm,  g3dlm
        real :: g1dcon, g2dcon, g3dcon

        real :: S1jj, S2jj, S3jj, S(nsd,nsd), Lapt

        real :: fact, lsl

!.... metrics for the curved wall

        real :: h, dhds, dhdr, dhdsr, dhdrr, hinv, hinv2, hinv3
!=============================================================================!

        call getmean( iver, y, vm, g1vm, g2vm, g11vm, g12vm, g22vm )

!=============================================================================

!.... hack-in parallel flow temporarily

!       call getmeanp( iver, y, vm, g2vm, g22vm )

!       vm(2)    = zero
!       g1vm     = zero
!       g2vm(3)  = zero
!       g11vm    = zero
!       g12vm    = zero
!       g22vm(3) = zero

!=============================================================================

!.... compute the curvature metrics

       if (curve.eq.0) then
         lsl = -one
       else
         lsl = sl
       endif

        call calch( lsl, 1, y, h, dhds, dhdr, dhdsr, dhdrr )

        hinv  = one / h
        hinv2 = hinv**2
        hinv3 = hinv**3

!.... some abbreviations

        rhom = vm(1)
        u1m  = vm(2)
        u2m  = vm(3)
        u3m  = vm(4)
        tm   = vm(5)

!.... initialize derivatives in the spanwise direction (infinite span)

        g3vm  = zero
        g13vm = zero
        g23vm = zero
        g33vm = zero

!.... initialize gradient of mean velocity (not the velocity gradient tensor)

        gum(1,1) = g1vm(2)
        gum(1,2) = g2vm(2)
        gum(1,3) = g3vm(2)

        gum(2,1) = g1vm(3)
        gum(2,2) = g2vm(3)
        gum(2,3) = g3vm(3)

        gum(3,1) = g1vm(4)
        gum(3,2) = g2vm(4)
        gum(3,3) = g3vm(4)

!.... compute the mean divergence (curve*)

        divum = ( gum(1,1) + u2m * dhdr ) * hinv + gum(2,2) + gum(3,3)

!.... initialize gradient of rho and T in the mean
!.... not the gradient for curvilinear coordinates

        grhom(1) = g1vm(1)
        grhom(2) = g2vm(1)
        grhom(3) = g3vm(1)

        gtm(1) = g1vm(5)
        gtm(2) = g2vm(5)
        gtm(3) = g3vm(5)

!.... compute the gradient of the divergence of um (curve*)

        g1divum = -dhds * hinv3 * ( gum(1,1) + u2m * dhdr ) + &
                  hinv2 * ( g11vm(2) + gum(2,1) * dhdr + &
                  u2m * dhdsr ) + hinv * ( g12vm(3) + g13vm(4) )

        g2divum = -dhdr * hinv2 * ( gum(1,1) + u2m * dhdr ) + &
                  hinv * ( g12vm(2) + gum(2,2) * dhdr + &
                  u2m * dhdrr ) + ( g22vm(3) + g23vm(4) )

        g3divum = hinv * ( g13vm(2) + gum(2,3) * dhdr ) + g23vm(3) + g33vm(4)

!.... compute some stuff that is useful for the viscous terms (curve*)

        S(1,1) = ( gum(1,1) + u2m * dhdr ) * hinv
        S(1,2) = pt5 * ( ( gum(2,1) - u1m * dhdr ) * hinv + gum(1,2) )
        S(1,3) = pt5 * ( gum(3,1) * hinv + gum(1,3) )
        S(2,1) = S(1,2)
        S(2,2) = gum(2,2)
        S(2,3) = pt5 * ( gum(3,2) + gum(2,3) )
        S(3,1) = S(1,3)
        S(3,2) = S(2,3)
        S(3,3) = gum(3,3)

        S1jj = -pt5 * ( dhdr**2 + dhdrr * h ) * hinv2 * u1m + &
               pt5 * dhdr * gum(1,2) * hinv + pt5 * g22vm(2) - &
               dhds * gum(1,1) * hinv3 + g11vm(2) * hinv2 + &
               pt5 * g33vm(2) + ( h * dhdsr - dhdr * dhds ) * hinv3 * u2m + &
               3.0 * dhdr * gum(2,1) / (two * h**2) + &
               pt5 * g12vm(3) * hinv + pt5 * g13vm(4) * hinv

        S2jj = pt5 * ( dhdr * dhds - h * dhdsr ) * hinv3 * u1m - &
               3.0 * dhdr * gum(1,1) / (two * h**2) + &
               pt5 * g12vm(2) * hinv - dhdr**2 * u2m * hinv2 + &
               dhdr * gum(2,2) * hinv + g22vm(3) - &
               pt5 * dhds * gum(2,1) * hinv3 + pt5 * g11vm(3) * hinv2 + &
               pt5 * g33vm(3) + pt5 * g23vm(4)

        S3jj = pt5 * g13vm(2) * hinv + pt5 * g23vm(3) + &
               pt5 * dhdr * gum(2,3) * hinv + pt5 * dhdr * gum(3,2) * hinv + &
               pt5 * g22vm(4) + pt5 * g11vm(4) * hinv2 + g33vm(4) - &
               pt5 * dhds * gum(3,1) * hinv3

!.... Laplacian of Tm (curve*)

        LapT = hinv * ( -dhds/h**2 * gtm(1) + hinv * g11vm(5) + &
               h * g22vm(5) + gtm(2) * dhdr + h * g33vm(5) )

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

        g1mu = dmu * gtm(1)
        g2mu = dmu * gtm(2)
        g3mu = dmu * gtm(3)

        g1dmu = d2mu * gtm(1)
        g2dmu = d2mu * gtm(2)
        g3dmu = d2mu * gtm(3)

!.... compute gradients of conductivity using chain-rule

        g1con = dcon * gtm(1)
        g2con = dcon * gtm(2)
        g3con = dcon * gtm(3)

        g1dcon = d2con * gtm(1)
        g2dcon = d2con * gtm(2)
        g3dcon = d2con * gtm(3)

!.... compute gradients of bulk viscosity using chain-rule

        g1lm = dlm * gtm(1)
        g2lm = dlm * gtm(2)
        g3lm = dlm * gtm(3)

        g1dlm = d2lm * gtm(1)
        g2dlm = d2lm * gtm(2)
        g3dlm = d2lm * gtm(3)

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

        B(1,1) = u2m
        B(1,3) = rhom

        C(1,1) = u3m
        C(1,4) = rhom

        D(1,1) = divum
        D(1,2) = grhom(1) * hinv
        D(1,3) = grhom(2) + rhom * dhdr * hinv
        D(1,4) = grhom(3)

!.... Momentum equation -- x_1 (convective + pressure) (curve*)

        G(2,2) = rhom

        A(2,1) = tm/(h * gamma * Ma**2)
        A(2,2) = rhom * u1m * hinv
        A(2,5) = rhom/(h * gamma * Ma**2)

        B(2,2) = rhom * u2m

        C(2,2) = rhom * u3m

        D(2,1) = u1m * hinv * ( gum(1,1) + u2m * dhdr ) + &
                 u2m * gum(1,2) + u3m * gum(1,3) + &
                 gtm(1) / (h * gamma * Ma**2)
        D(2,2) = rhom * ( gum(1,1) + u2m * dhdr ) * hinv
        D(2,3) = rhom * ( gum(1,2) + u1m * dhdr * hinv )
        D(2,4) = rhom * gum(1,3)
        D(2,5) = grhom(1) / (h * gamma * Ma**2)

!.... (viscous lambda) (curve*)

        fact = rlme / (rmue * Re)

        A(2,2) = A(2,2) - fact * ( g1lm * hinv2 - rlm * hinv3 * dhds )
        A(2,3) = A(2,3) - fact * rlm * hinv2 * dhdr
        A(2,5) = A(2,5) - fact * dlm * divum * hinv

        B(2,3) = B(2,3) - fact * ( g1lm * hinv )

        C(2,4) = C(2,4) - fact * ( g1lm * hinv )

        D(2,3) = D(2,3) - fact * ( g1lm * dhdr * hinv2 - &
                                   rlm * hinv3 * dhds * dhdr + &
                                   rlm * hinv2 * dhdsr )
        D(2,5) = D(2,5) - fact * ( g1dlm * divum * hinv + dlm * g1divum )

        Vxx(2,2) = fact * rlm * hinv2

        Vxy(2,3) = fact * rlm * hinv

        Vxz(2,4) = fact * rlm * hinv

!.... (viscous mu) (curve*)

        fact = one / Re

        A(2,2) = A(2,2) - fact * ( two * g1mu * hinv2 - &
                              two * rmu * dhds * hinv3 )
        A(2,3) = A(2,3) - fact * ( g2mu * hinv + &
                              rmu * 3.0 * dhdr * hinv2 )
        A(2,4) = A(2,4) - fact * g3mu * hinv
        A(2,5) = A(2,5) - fact * dmu * two * S(1,1) * hinv

        B(2,2) = B(2,2) - fact * ( g2mu + rmu * dhdr * hinv )
        B(2,5) = B(2,5) - fact * dmu * two * S(1,2)

        C(2,2) = C(2,2) - fact * g3mu
        C(2,5) = C(2,5) - fact * dmu * two * S(1,3)

        D(2,2) = D(2,2) - fact * ( g2mu * hinv * (-dhdr) - &
                              rmu * ( dhdr**2 + dhdrr * h ) * hinv2 )
        D(2,3) = D(2,3) - fact * ( two * g1mu * hinv2 * dhdr + &
                              two * rmu * ( dhdsr * h - dhds * dhdr ) * hinv3 )
        D(2,5) = D(2,5) - fact * two * ( g1dmu * hinv * S(1,1) + &
                              g2dmu * S(1,2) + g3dmu * S(1,3) + &
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
        B(3,3) = rhom * u2m
        B(3,5) = rhom/(gamma * Ma**2)

        C(3,3) = rhom * u3m

        D(3,1) = u1m * hinv * ( gum(2,1) - u1m * dhdr ) + &
                   u2m * gum(2,2) + u3m * gum(2,3) + &
                   gtm(2) / (gamma * Ma**2)
        D(3,2) = rhom * ( gum(2,1) - two * u1m * dhdr ) * hinv
        D(3,3) = rhom * gum(2,2)
        D(3,4) = rhom * gum(2,3)
        D(3,5) = grhom(2) / (gamma * Ma**2)

!.... (viscous lambda) (curve*)

        fact = rlme / (rmue * Re)

        A(3,2) = A(3,2) - fact * ( g2lm * hinv - rlm * dhdr * hinv2 )

        B(3,3) = B(3,3) - fact * ( g2lm + rlm * dhdr * hinv )
        B(3,5) = B(3,5) - fact * dlm * divum

        C(3,4) = C(3,4) - fact * ( g2lm )

        D(3,3) = D(3,3) - fact * ( g2lm * hinv * dhdr - &
                              rlm * dhdr * hinv2 * dhdr + &
                              rlm * hinv * dhdrr )
        D(3,5) = D(3,5) - fact * ( g2dlm * divum + dlm * g2divum )

        Vxy(3,2) = fact * rlm * hinv

        Vyy(3,3) = fact * rlm

        Vyz(3,4) = fact * rlm

!.... (viscous mu) (curve*)

        fact = one / Re

        A(3,2) = A(3,2) + fact * rmu * 3.0 * dhdr * hinv2
        A(3,3) = A(3,3) - fact * ( g1mu * hinv2 - &
                              rmu * dhds * hinv3 )
        A(3,5) = A(3,5) - fact * dmu * two * S(2,1) * hinv

        B(3,2) = B(3,2) - fact * g1mu * hinv
        B(3,3) = B(3,3) - fact * ( two * g2mu + two * rmu * dhdr * hinv )
        B(3,4) = B(3,4) - fact * g3mu
        B(3,5) = B(3,5) - fact * dmu * two * S(2,2)

        C(3,3) = C(3,3) - fact * g3mu
        C(3,5) = C(3,5) - fact * dmu * two * S(2,3)

        D(3,2) = D(3,2) - fact * ( g1mu * hinv2 * (-dhdr) + &
                              rmu * (dhds * dhdr - h * dhdsr) * hinv3 )
        D(3,3) = D(3,3) + fact * two * rmu * dhdr**2 * hinv2
        D(3,5) = D(3,5) - fact * two * ( g1dmu * hinv * S(2,1) + &
                              g2dmu * S(2,2) + g3dmu * S(2,3) + &
                              dmu * S2jj )

        Vxx(3,3) = Vxx(3,3) + fact * rmu * hinv2
        Vxy(3,2) = Vxy(3,2) + fact * rmu * hinv
        Vyy(3,3) = Vyy(3,3) + fact * two * rmu
        Vyz(3,4) = Vyz(3,4) + fact * rmu
        Vzz(3,3) = Vzz(3,3) + fact * rmu

!.... now use continuity to remove Vyy(:,3,3) term

        fact = Vyy(3,3) / rhom

        A(3,1) = A(3,1) + fact * ( gum(1,2) * hinv - u1m * dhdr * hinv2 )
        A(3,2) = A(3,2) + fact * ( grhom(2) * hinv - rhom * dhdr * hinv2 )

        B(3,1) = B(3,1) + fact * ( divum + gum(2,2) - im * omega )
        B(3,2) = B(3,2) + fact * ( grhom(1) * hinv )
        B(3,3) = B(3,3) + fact * ( two * grhom(2) + rhom * dhdr * hinv )
        B(3,4) = B(3,4) + fact * ( grhom(3) )

        C(3,1) = C(3,1) + fact * ( gum(3,2) )
        C(3,4) = C(3,4) + fact * ( grhom(2) )

        D(3,1) = D(3,1) + fact * ( g2divum )
        D(3,2) = D(3,2) + fact * ( g12vm(1) * hinv - &
                                        grhom(1) * dhdr * hinv2 )
        D(3,3) = D(3,3) + fact * ( grhom(2) * dhdr * hinv - &
                              rhom * dhdr**2 * hinv2 + rhom * dhdrr * hinv + &
                              g22vm(1) )
        D(3,4) = D(3,4) +  fact * ( g23vm(1) )

        Vxy(3,1) = Vxy(3,1) - fact * u1m * hinv
        Vxy(3,2) = Vxy(3,2) - fact * rhom * hinv
        Vyy(3,1) = Vyy(3,1) - fact * u2m
        Vyz(3,1) = Vyz(3,1) - fact * u3m
        Vyz(3,4) = Vyz(3,4) - fact * rhom

!.... zero the offending term

        Vyy(3,3) = zero

!.... Momentum equation -- x_3 (convective + pressure) (curve*)

        G(4,4) = rhom

        A(4,4) = rhom * u1m * hinv

        B(4,4) = rhom * u2m

        C(4,1) = tm/(gamma * Ma**2)
        C(4,4) = rhom * u3m
        C(4,5) = rhom/(gamma * Ma**2)

        D(4,1) = u1m * gum(3,1) * hinv + u2m * gum(3,2) + &
                   u3m * gum(3,3) + gtm(3) / (gamma * Ma**2)
        D(4,2) = rhom * gum(3,1) * hinv
        D(4,3) = rhom * gum(3,2)
        D(4,4) = rhom * gum(3,3)
        D(4,5) = grhom(3) / (gamma * Ma**2)

!.... (viscous lambda) (curve*)

        fact = rlme / (rmue * Re)

        A(4,2) = A(4,2) - fact * g3lm * hinv

        B(4,3) = B(4,3) - fact * g3lm

        C(4,4) = C(4,4) - fact * g3lm
        C(4,3) = C(4,3) - fact * rlm * hinv * dhdr
        C(4,5) = C(4,5) - fact * dlm * divum

        D(4,3) = D(4,3) - fact * ( g2lm * hinv * dhdr )
        D(4,5) = D(4,5) - fact * ( g3dlm * divum + dlm * g3divum )

        Vxz(4,2) = fact * rlm * hinv

        Vyz(4,3) = fact * rlm

        Vzz(4,4) = fact * rlm

!.... (viscous mu) (curve*)

        fact = one / Re

        A(4,4) = A(4,4) - fact * ( g1mu * hinv2 - rmu * dhds * hinv3 )
        A(4,5) = A(4,5) - fact * dmu * two * S(3,1) * hinv

        B(4,4) = B(4,4) - fact * ( g2mu + rmu * dhdr * hinv )
        B(4,5) = B(4,5) - fact * dmu * two * S(3,2)

        C(4,2) = C(4,2) - fact * g1mu * hinv
        C(4,3) = C(4,3) - fact * ( g2mu + rmu * dhdr * hinv )
        C(4,4) = C(4,4) - fact * two * g3mu
        C(4,5) = C(4,5) - fact * dmu * two * S(3,3)

        D(4,5) = D(4,5) - fact * two * ( g1dmu * hinv * S(3,1) + &
                              g2dmu * S(3,2) + g3dmu * S(3,3) + &
                              dmu * S3jj )

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
        B(5,5) = rhom * u2m

        C(5,4) = rhom * gamma1 * tm
        C(5,5) = rhom * u3m

        D(5,1) = u1m * hinv * gtm(1) + u2m * gtm(2) + u3m * gtm(3) + &
                   gamma1 * tm * divum
        D(5,2) = rhom * gtm(1) * hinv
        D(5,3) = rhom * gtm(2) + rhom * gamma1 * tm * dhdr * hinv
        D(5,4) = rhom * gtm(3)
        D(5,5) = rhom * gamma1 * divum

!.... diffusion (curve*)

        fact = gamma / (Pr * Re)

        A(5,5) = A(5,5) - fact * (g1con * hinv2 + dcon * gtm(1) * hinv2 - &
                              con * dhds * hinv3 )
        B(5,5) = B(5,5) - fact * (g2con + dcon * gtm(2) + &
                              con * dhdr * hinv )
        C(5,5) = C(5,5) - fact * (g3con + dcon * gtm(3))
        D(5,5) = D(5,5) - fact * (g1dcon * gtm(1) * hinv2 + &
                              g2dcon * gtm(2) + g3dcon * gtm(3) + &
                              dcon * Lapt )

        Vxx(5,5) = fact * con * hinv2
        Vyy(5,5) = fact * con
        Vzz(5,5) = fact * con

!.... dissipation (lambda) (curve*)

        fact = gamma * gamma1 * Ma**2 * rlme / (Re * rmue)

        A(5,2) = A(5,2) - fact * two * rlm * divum * hinv
        B(5,3) = B(5,3) - fact * two * rlm * divum
        C(5,4) = C(5,4) - fact * two * rlm * divum
        D(5,3) = D(5,3) - fact * two * rlm * divum * dhdr * hinv
        D(5,5) = D(5,5) - fact * dlm * divum * divum

!.... dissipation (mu) (curve*)

        fact = gamma * gamma1 * Ma**2 / Re

        A(5,2) = A(5,2) - fact * four * rmu * S(1,1) * hinv
        A(5,3) = A(5,3) - fact * four * rmu * S(2,1) * hinv
        A(5,4) = A(5,4) - fact * four * rmu * S(3,1) * hinv

        B(5,2) = B(5,2) - fact * four * rmu * S(1,2)
        B(5,3) = B(5,3) - fact * four * rmu * S(2,2)
        B(5,4) = B(5,4) - fact * four * rmu * S(3,2)

        C(5,2) = C(5,2) - fact * four * rmu * S(1,3)
        C(5,3) = C(5,3) - fact * four * rmu * S(2,3)
        C(5,4) = C(5,4) - fact * four * rmu * S(3,3)

        D(5,2) = D(5,2) + fact * four * rmu * S(2,1) * dhdr * hinv
        D(5,3) = D(5,3) - fact * four * rmu * S(1,1) * dhdr * hinv
        D(5,5) = D(5,5) - fact * two * dmu * ( S(1,1)**2 + &
                       S(1,2)**2 + S(1,3)**2 + S(2,1)**2 + &
                       S(2,2)**2 + S(2,3)**2 + S(3,1)**2 + &
                       S(3,2)**2 + S(3,3)**2)

!==============================================================================
!.... form the extended matrices
!==============================================================================
        Ah = zero
        Bh = zero
        Ch = zero
        Eh = zero
        Fh = zero

        Ah(1:5,1:5) = two * im * alpha * Vxx + im * beta * Vxz - A
        Bh(1:5,1:5) = Vxx
        Ch(1:5,1:5) = Vxy

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

        return
        end
