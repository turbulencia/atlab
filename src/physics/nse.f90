!########################################################################
!#
!# Evolution equations per unit volume, nonlinear term in convective form and the
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!#
!########################################################################
subroutine NavierStokes_PerVolume(hq, hs, dte, remove_divergence)
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_NN
    use TLab_Memory, only: imax, jmax, kmax, isize_field, inb_flow, inb_scal
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3
    use TLab_Arrays, only: s
    use OPR_Partial
    use NSE_Burgers
    use OPR_Elliptic, only: OPR_Poisson
    implicit none

    real(wp), intent(inout) :: hq(isize_field, inb_flow)
    real(wp), intent(inout) :: hs(isize_field, inb_scal)
    real(wp), intent(in) :: dte
    logical, intent(in) :: remove_divergence

    ! -----------------------------------------------------------------------
    integer(wi) is

    ! #######################################################################
    ! Diffusion and advection terms
    ! #######################################################################
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, w, hq(:, 3), tmp1, tmp3)               ! store rho w in tmp3
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, u, hq(:, 1), tmp1, tmp2, rhou_in=tmp3)
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, v, hq(:, 2), tmp1, tmp2, rhou_in=tmp3)
    do is = 1, inb_scal
        call NSE_AddBurgers_PerVolume_Z(is, imax, jmax, kmax, s(:, is), hs(:, is), tmp1, tmp2, rhou_in=tmp3)
    end do

    call NSE_AddBurgers_PerVolume_X(0, imax, jmax, kmax, u, hq(:, 1), tmp1, tmp3)               ! store rho u transposed in tmp3
    call NSE_AddBurgers_PerVolume_X(0, imax, jmax, kmax, v, hq(:, 2), tmp1, tmp2, rhou_in=tmp3)
    call NSE_AddBurgers_PerVolume_X(0, imax, jmax, kmax, w, hq(:, 3), tmp1, tmp2, rhou_in=tmp3)
    do is = 1, inb_scal
        call NSE_AddBurgers_PerVolume_X(is, imax, jmax, kmax, s(:, is), hs(:, is), tmp1, tmp2, rhou_in=tmp3)
    end do

    call NSE_AddBurgers_PerVolume_Y(0, imax, jmax, kmax, v, hq(:, 2), tmp1, tmp3)               ! store rho v transposed in tmp3
    call NSE_AddBurgers_PerVolume_Y(0, imax, jmax, kmax, u, hq(:, 1), tmp1, tmp2, rhou_in=tmp3)
    call NSE_AddBurgers_PerVolume_Y(0, imax, jmax, kmax, w, hq(:, 3), tmp1, tmp2, rhou_in=tmp3)
    do is = 1, inb_scal
        call NSE_AddBurgers_PerVolume_Y(is, imax, jmax, kmax, s(:, is), hs(:, is), tmp1, tmp2, rhou_in=tmp3)
    end do

    ! #######################################################################
    ! Pressure term
    ! #######################################################################
    ! Forcing term
    if (remove_divergence) then ! remove residual divergence
        call Add_Residual_Divergence(w, hq(1, 3), tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp2, result=tmp1)     ! accumulate forcing term in tmp1
        call Add_Residual_Divergence(v, hq(1, 2), tmp2)
        call OPR_Partial_Y(OPR_P1_ADD, imax, jmax, kmax, tmp2, result=tmp1)
        call Add_Residual_Divergence(u, hq(1, 1), tmp2)
        call OPR_Partial_X(OPR_P1_ADD, imax, jmax, kmax, tmp2, result=tmp1)

    else
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, hq(:, 3), result=tmp1)
        call OPR_Partial_Y(OPR_P1_ADD, imax, jmax, kmax, hq(:, 2), result=tmp1)
        call OPR_Partial_X(OPR_P1_ADD, imax, jmax, kmax, hq(:, 1), result=tmp1)

    end if

    ! Solution of Poisson equation: pressure in tmp1
    call OPR_Poisson(imax, jmax, kmax, BCS_NN, tmp1, tmp2, tmp3, &
                     bcs_hb=hq(1:imax*jmax, 3), &                           ! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
                     bcs_ht=hq(isize_field - imax*jmax + 1:isize_field, 3))

    ! Add pressure gradient
    call OPR_Partial_X(OPR_P1_SUBTRACT, imax, jmax, kmax, tmp1, result=hq(:, 1))
    call OPR_Partial_Y(OPR_P1_SUBTRACT, imax, jmax, kmax, tmp1, result=hq(:, 2))
    call OPR_Partial_Z(OPR_P1_SUBTRACT, imax, jmax, kmax, tmp1, result=hq(:, 3))

    return

contains
    subroutine Add_Residual_Divergence(q, hq, result)
        use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC, DNS_EQNS_BOUSSINESQ
        use Thermo_Anelastic, only: rbackground
        real(wp), intent(in) :: q(imax*jmax, kmax)
        real(wp), intent(in) :: hq(imax*jmax, kmax)
        real(wp), intent(out) :: result(imax*jmax, kmax)

        integer(wi) k
        real(wp) dummy

        dummy = 1.0_wp/dte

        select case (nse_eqns)
        case (DNS_EQNS_ANELASTIC)
            do k = 1, kmax
                result(:, k) = hq(:, k) + q(:, k)*dummy*rbackground(k)
            end do

        case (DNS_EQNS_BOUSSINESQ)
            do k = 1, kmax
                result(:, k) = hq(:, k) + q(:, k)*dummy
            end do

        end select

        return
    end subroutine

end subroutine NavierStokes_PerVolume
