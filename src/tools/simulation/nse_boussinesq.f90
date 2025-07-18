!########################################################################
!#
!# Evolution equations, nonlinear term in convective form and the
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!#
!# It is written such that u and w transposes are calculated first for the
!# Ox and Oy momentum equations, stored in tmp5 and tmp6 and then used as needed.
!# This saves 2 MPI transpositions.
!# Includes the scalar to benefit from the same reduction
!#
!########################################################################
subroutine NSE_Boussinesq()
    use TLab_Constants, only: wp, wi, BCS_NN
    use TLab_Memory, only: imax, jmax, kmax, inb_flow, inb_scal
    use TLab_Arrays, only: s
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
    use DNS_Arrays
    use TimeMarching, only: dte, remove_divergence
    use BoundaryConditions
    use OPR_Partial
    use NSE_Burgers
    use OPR_Elliptic

    implicit none

    ! -----------------------------------------------------------------------
    integer(wi) iq, is
    integer ibc
    real(wp) dummy

    ! #######################################################################
    ! Preliminaries for Scalar BC
    ! (flow BCs initialized below as they are used for pressure in between)
    ! #######################################################################
    BcsScalKmin%ref(:, :, :) = 0.0_wp ! default is no-slip (dirichlet)
    BcsScalKmax%ref(:, :, :) = 0.0_wp

    ! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
    if (any(BcsScalKmin%SfcType(1:inb_scal) == DNS_SFC_LINEAR) .or. &
        any(BcsScalKmax%SfcType(1:inb_scal) == DNS_SFC_LINEAR)) then
        do is = 1, inb_scal
            if (BcsScalKmin%SfcType(is) == DNS_SFC_LINEAR) BcsScalKmin%ref(:, :, is) = p_hs(:, :, 1, is)
            if (BcsScalKmax%SfcType(is) == DNS_SFC_LINEAR) BcsScalKmax%ref(:, :, is) = p_hs(:, :, kmax, is)
        end do
    end if

    ! #######################################################################
    ! Diffusion and advection terms
    ! #######################################################################
    ! Diagonal terms and transposed velocity arrays
    call NSE_Burgers_X(0, imax, jmax, kmax, u, tmp1, tmp4) ! store u transposed in tmp4
    call NSE_Burgers_Y(0, imax, jmax, kmax, v, tmp2, tmp5) ! store v transposed in tmp5
    call NSE_Burgers_Z(0, imax, jmax, kmax, w, tmp3, w)

    ! Ox momentum equation
    call NSE_Burgers_Y(0, imax, jmax, kmax, u, tmp6, tmp8, tmp5) ! tmp5 contains v transposed
    call NSE_Burgers_Z(0, imax, jmax, kmax, u, tmp7, w)
    hq(:, 1) = hq(:, 1) + tmp1(:) + tmp6(:) + tmp7(:)

    ! Oy momentum equation
    call NSE_Burgers_X(0, imax, jmax, kmax, v, tmp6, tmp8, tmp4) ! tmp4 contains u transposed
    call NSE_Burgers_Z(0, imax, jmax, kmax, v, tmp7, w)
    hq(:, 2) = hq(:, 2) + tmp2(:) + tmp6(:) + tmp7(:)

    ! Oz momentum equation
    call NSE_Burgers_X(0, imax, jmax, kmax, w, tmp6, tmp8, tmp4) ! tmp4 contains u transposed
    call NSE_Burgers_Y(0, imax, jmax, kmax, w, tmp7, tmp8, tmp5) ! tmp5 contains v transposed
    hq(:, 3) = hq(:, 3) + tmp3(:) + tmp6(:) + tmp7(:)

    ! Scalar equations
    do is = 1, inb_scal
        call NSE_Burgers_X(is, imax, jmax, kmax, s(:, is), tmp1, tmp8, tmp4) ! tmp4 contains u transposed
        call NSE_Burgers_Y(is, imax, jmax, kmax, s(:, is), tmp2, tmp8, tmp5) ! tmp5 contains v transposed
        call NSE_Burgers_Z(is, imax, jmax, kmax, s(:, is), tmp3, w)
        hs(:, is) = hs(:, is) + tmp1(:) + tmp2(:) + tmp3(:)

    end do

    ! #######################################################################
    ! Pressure term
    ! #######################################################################
    ! Forcing term
    if (remove_divergence) then ! remove residual divergence
        dummy = 1.0_wp/dte
        tmp2(:) = hq(:, 1) + u(:)*dummy
        tmp3(:) = hq(:, 2) + v(:)*dummy
        tmp4(:) = hq(:, 3) + w(:)*dummy

        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, tmp2, tmp1)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, tmp3, tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp4, tmp3)

    else
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, hq(:, 1), tmp1)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, hq(:, 2), tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, hq(:, 3), tmp3)

    end if
    tmp1(:) = tmp1(:) + tmp2(:) + tmp3(:) ! forcing term in tmp1

    ! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    BcsFlowKmin%ref(:, :, 3) = p_hq(:, :, 1, 3)
    BcsFlowKmax%ref(:, :, 3) = p_hq(:, :, kmax, 3)

    ! Solution of Poisson equation: pressure in tmp1
    call OPR_Poisson(imax, jmax, kmax, BCS_NN, tmp1, tmp2, tmp3, BcsFlowKmin%ref(:, :, 3), BcsFlowKmax%ref(:, :, 3))

    ! Add pressure gradient
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, tmp1, tmp2)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, tmp1, tmp3)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp1, tmp4)
    hq(:, 1) = hq(:, 1) - tmp2(:)
    hq(:, 2) = hq(:, 2) - tmp3(:)
    hq(:, 3) = hq(:, 3) - tmp4(:)

    ! #######################################################################
    ! Boundary conditions
    ! #######################################################################
    BcsFlowKmin%ref = 0.0_wp ! default is no-slip (dirichlet)
    BcsFlowKmax%ref = 0.0_wp ! Scalar BCs initialized at start of routine

    do iq = 1, inb_flow
        ibc = 0
        if (BcsFlowKmin%type(iq) == DNS_BCS_Neumann) ibc = ibc + 1
        if (BcsFlowKmax%type(iq) == DNS_BCS_Neumann) ibc = ibc + 2
        if (ibc > 0) then
            call BCS_Neumann_Z(ibc, imax, jmax, kmax, hq(:, iq), &
                                        BcsFlowKmin%ref(:, :, iq), BcsFlowKmax%ref(:, :, iq), tmp1)
        end if

        p_hq(:, :, 1, iq) = BcsFlowKmin%ref(:, :, iq)
        p_hq(:, :, kmax, iq) = BcsFlowKmax%ref(:, :, iq)

    end do

    do is = 1, inb_scal
        ibc = 0
        if (BcsScalKmin%type(is) == DNS_BCS_Neumann) ibc = ibc + 1
        if (BcsScalKmax%type(is) == DNS_BCS_Neumann) ibc = ibc + 2
        if (ibc > 0) then
            call BCS_Neumann_Z(ibc, imax, jmax, kmax, hs(:, is), &
                                        BcsScalKmin%ref(:, :, is), BcsScalKmax%ref(:, :, is), tmp1)
        end if

        if (BcsScalJmin%type(is) /= DNS_SFC_STATIC .or. &
            BcsScalKmax%type(is) /= DNS_SFC_STATIC) then
            call BCS_SURFACE_Z(is, s, hs, tmp1, tmp2)
        end if

        p_hs(:, :, 1, is) = BcsScalKmin%ref(:, :, is)
        p_hs(:, :, kmax, is) = BcsScalKmax%ref(:, :, is)

    end do

    return
end subroutine NSE_Boussinesq
