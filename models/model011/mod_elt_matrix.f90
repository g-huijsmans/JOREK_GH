!> Analytic direct-harmonic electrostatic Poisson element matrix.
module mod_elt_matrix
  implicit none
contains

  subroutine element_matrix(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, &
                            R_xpoint, Z_xpoint, ELM, RHS, tid, i_tor_min, i_tor_max, aux_nodes)
    use constants,                  only: ATOMIC_MASS_UNIT, EL_CHG
    use mod_parameters,            only: n_vertex_max, n_degrees, var_psi, var_rho
    use data_structure,            only: type_element, type_node
    use gauss,                     only: n_gauss, wgauss
    use basis_at_gaussian,         only: H, H_s, H_t, H_ss, H_st, H_tt
    use phys_module,               only: F0, central_mass, filter_perp, filter_hyper, filter_par, &
                                         filter_perp_n0, filter_hyper_n0, filter_par_n0, mode, mode_type
    use mod_poisson_element_kernel, only: accumulate_poisson_n0_block, &
                                          accumulate_poisson_nonzero_blocks, scatter_poisson_harmonics

    implicit none

    type(type_element), intent(in)           :: element
    type(type_node), intent(in)              :: nodes(n_vertex_max)
    type(type_node), intent(in), optional    :: aux_nodes(n_vertex_max)
    logical, intent(in)                      :: xpoint2
    integer, intent(in)                      :: xcase2, tid, i_tor_min, i_tor_max
    real*8, intent(in)                       :: R_axis, Z_axis, psi_axis, psi_bnd
    real*8, intent(in)                       :: R_xpoint(2), Z_xpoint(2)
    real*8, dimension(:,:), allocatable      :: ELM
    real*8, dimension(:), allocatable        :: RHS

    integer, parameter :: basis_size = n_vertex_max*n_degrees
    integer            :: i, j, ms, mt, a, n_tor_local
    real*8             :: xjac, xjac_x, xjac_y, big_r, bb2, factor, weight
    real*8             :: psi_g, psi_s, psi_t, psi_x, psi_y, psi_norm, rho, filter_par_centre
    real*8             :: x_g, x_s, x_t, x_ss, x_st, x_tt
    real*8             :: y_g, y_s, y_t, y_ss, y_st, y_tt
    real*8             :: value(basis_size), deriv_x(basis_size), deriv_y(basis_size)
    real*8             :: deriv_s(basis_size), deriv_t(basis_size)
    real*8             :: deriv_ss(basis_size), deriv_st(basis_size), deriv_tt(basis_size)
    real*8             :: deriv_xx(basis_size), deriv_yy(basis_size), laplace_star(basis_size)
    real*8             :: block_n0(basis_size,basis_size), block_a(basis_size,basis_size)
    real*8             :: block_b(basis_size,basis_size), block_c(basis_size,basis_size)
    integer, allocatable :: component_modes(:)
    character(len=3), allocatable :: component_types(:)

    ELM = 0.d0
    RHS = 0.d0
    block_n0 = 0.d0; block_a = 0.d0; block_b = 0.d0; block_c = 0.d0

    ! Operator coefficients are equilibrium/background quantities and are
    ! deliberately reconstructed from physical storage component n=0 only.
    do ms = 1, n_gauss
      do mt = 1, n_gauss
        x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; x_ss = 0.d0; x_st = 0.d0; x_tt = 0.d0
        y_g = 0.d0; y_s = 0.d0; y_t = 0.d0; y_ss = 0.d0; y_st = 0.d0; y_tt = 0.d0
        psi_g = 0.d0; psi_s = 0.d0; psi_t = 0.d0; rho = 0.d0

        do i = 1, n_vertex_max
          do j = 1, n_degrees
            a = (i-1)*n_degrees + j
            value(a)    = element%size(i,j)*H(i,j,ms,mt)
            deriv_s(a)  = element%size(i,j)*H_s(i,j,ms,mt)
            deriv_t(a)  = element%size(i,j)*H_t(i,j,ms,mt)
            deriv_ss(a) = element%size(i,j)*H_ss(i,j,ms,mt)
            deriv_st(a) = element%size(i,j)*H_st(i,j,ms,mt)
            deriv_tt(a) = element%size(i,j)*H_tt(i,j,ms,mt)

            x_g  = x_g  + nodes(i)%x(1,j,1)*value(a)
            x_s  = x_s  + nodes(i)%x(1,j,1)*deriv_s(a)
            x_t  = x_t  + nodes(i)%x(1,j,1)*deriv_t(a)
            x_ss = x_ss + nodes(i)%x(1,j,1)*deriv_ss(a)
            x_st = x_st + nodes(i)%x(1,j,1)*deriv_st(a)
            x_tt = x_tt + nodes(i)%x(1,j,1)*deriv_tt(a)
            y_g  = y_g  + nodes(i)%x(1,j,2)*value(a)
            y_s  = y_s  + nodes(i)%x(1,j,2)*deriv_s(a)
            y_t  = y_t  + nodes(i)%x(1,j,2)*deriv_t(a)
            y_ss = y_ss + nodes(i)%x(1,j,2)*deriv_ss(a)
            y_st = y_st + nodes(i)%x(1,j,2)*deriv_st(a)
            y_tt = y_tt + nodes(i)%x(1,j,2)*deriv_tt(a)

            psi_g = psi_g + nodes(i)%values(1,j,var_psi)*value(a)
            psi_s = psi_s + nodes(i)%values(1,j,var_psi)*deriv_s(a)
            psi_t = psi_t + nodes(i)%values(1,j,var_psi)*deriv_t(a)
            rho   = rho   + nodes(i)%values(1,j,var_rho)*value(a)
          enddo
        enddo

        xjac = x_s*y_t - x_t*y_s
        xjac_x = (x_ss*y_t**2 - y_ss*x_t*y_t - 2.d0*x_st*y_s*y_t &
                 +y_st*(x_s*y_t+x_t*y_s) + x_tt*y_s**2-y_tt*x_s*y_s)/xjac
        xjac_y = (y_tt*x_s**2 - x_tt*y_s*x_s - 2.d0*y_st*x_t*x_s &
                 +x_st*(y_t*x_s+y_s*x_t) + y_ss*x_t**2-x_ss*y_t*x_t)/xjac

        big_r = x_g
        psi_x = ( y_t*psi_s-y_s*psi_t)/xjac
        psi_y = (-x_t*psi_s+x_s*psi_t)/xjac
        bb2 = (F0*F0 + psi_x*psi_x + psi_y*psi_y)/big_r**2
        factor = central_mass*ATOMIC_MASS_UNIT*rho/(EL_CHG*bb2)
        psi_norm = (psi_g-psi_axis)/(psi_bnd-psi_axis)
        filter_par_centre = 0.d0
        if (psi_norm.lt.0.64d0) filter_par_centre = 1.d0

        do a = 1, basis_size
          deriv_x(a) = ( y_t*deriv_s(a)-y_s*deriv_t(a))/xjac
          deriv_y(a) = (-x_t*deriv_s(a)+x_s*deriv_t(a))/xjac
          deriv_xx(a) = (deriv_ss(a)*y_t**2 - 2.d0*deriv_st(a)*y_s*y_t + deriv_tt(a)*y_s**2 &
                        +deriv_s(a)*(y_st*y_t-y_tt*y_s) + deriv_t(a)*(y_st*y_s-y_ss*y_t))/xjac**2 &
                        -xjac_x*(deriv_s(a)*y_t-deriv_t(a)*y_s)/xjac**2
          deriv_yy(a) = (deriv_ss(a)*x_t**2 - 2.d0*deriv_st(a)*x_s*x_t + deriv_tt(a)*x_s**2 &
                        +deriv_s(a)*(x_st*x_t-x_tt*x_s) + deriv_t(a)*(x_st*x_s-x_ss*x_t))/xjac**2 &
                        -xjac_y*(-deriv_s(a)*x_t+deriv_t(a)*x_s)/xjac**2
          laplace_star(a) = deriv_xx(a) + deriv_x(a)/big_r + deriv_yy(a)
        enddo

        weight = wgauss(ms)*wgauss(mt)
        call accumulate_poisson_n0_block(weight,big_r,xjac,factor,bb2,psi_x,psi_y, &
             filter_perp_n0,filter_hyper_n0,filter_par_n0+filter_par_centre, &
             value,deriv_x,deriv_y,laplace_star,block_n0)
        call accumulate_poisson_nonzero_blocks(weight,big_r,xjac,factor,bb2,F0,psi_x,psi_y, &
             filter_perp,filter_hyper,filter_par,value,deriv_x,deriv_y,laplace_star, &
             block_a,block_b,block_c)
      enddo
    enddo

    n_tor_local = i_tor_max-i_tor_min+1
    allocate(component_modes(n_tor_local))
    allocate(component_types(n_tor_local))
    component_modes = mode(i_tor_min:i_tor_max)
    component_types = mode_type(i_tor_min:i_tor_max)
    call scatter_poisson_harmonics(block_n0,block_a,block_b,block_c,component_modes,component_types,ELM)
    deallocate(component_modes, component_types)
  end subroutine element_matrix
end module mod_elt_matrix
