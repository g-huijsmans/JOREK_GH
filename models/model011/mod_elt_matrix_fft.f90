!> Compatibility entry point for generic assembly; model011 never performs FFTs.
module mod_elt_matrix_fft
  implicit none
contains
  subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, &
                                R_xpoint, Z_xpoint, ELM, RHS, tid, ELM_p, ELM_n, ELM_k, ELM_kn, &
                                RHS_p, RHS_k, eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, &
                                delta_g, delta_s, delta_t, i_tor_min, i_tor_max, aux_nodes, ELM_pnn)
    use mod_parameters, only: n_vertex_max, n_plane, n_var, n_degrees, n_eq_var
    use gauss, only: n_gauss
    use data_structure, only: type_element, type_node
    use mod_elt_matrix, only: element_matrix
    implicit none

    type(type_element), intent(in)        :: element
    type(type_node), intent(in)           :: nodes(n_vertex_max)
    type(type_node), intent(in), optional :: aux_nodes(n_vertex_max)
    logical, intent(in)                   :: xpoint2
    integer, intent(in)                   :: xcase2, tid, i_tor_min, i_tor_max
    real*8, intent(in)                    :: R_axis, Z_axis, psi_axis, psi_bnd
    real*8, intent(in)                    :: R_xpoint(2), Z_xpoint(2)
    real*8, allocatable                   :: ELM(:,:), RHS(:)
    real*8                                :: ELM_p(n_plane,n_vertex_max*n_var*n_degrees,n_vertex_max*n_var*n_degrees)
    real*8                                :: ELM_n(n_plane,n_vertex_max*n_var*n_degrees,n_vertex_max*n_var*n_degrees)
    real*8                                :: ELM_k(n_plane,n_vertex_max*n_var*n_degrees,n_vertex_max*n_var*n_degrees)
    real*8                                :: ELM_kn(n_plane,n_vertex_max*n_var*n_degrees,n_vertex_max*n_var*n_degrees)
    real*8                                :: ELM_pnn(n_plane,n_vertex_max*n_var*n_degrees,n_vertex_max*n_var*n_degrees)
    real*8                                :: RHS_p(n_plane,n_vertex_max*n_var*n_degrees)
    real*8                                :: RHS_k(n_plane,n_vertex_max*n_var*n_degrees)
    real*8                                :: eq_g(n_plane,n_eq_var,n_gauss,n_gauss)
    real*8                                :: eq_s(n_plane,n_eq_var,n_gauss,n_gauss)
    real*8                                :: eq_t(n_plane,n_eq_var,n_gauss,n_gauss)
    real*8                                :: eq_p(n_plane,n_eq_var,n_gauss,n_gauss)
    real*8                                :: eq_ss(n_plane,n_eq_var,n_gauss,n_gauss)
    real*8                                :: eq_st(n_plane,n_eq_var,n_gauss,n_gauss)
    real*8                                :: eq_tt(n_plane,n_eq_var,n_gauss,n_gauss)
    real*8                                :: delta_g(n_plane,n_eq_var,n_gauss,n_gauss)
    real*8                                :: delta_s(n_plane,n_eq_var,n_gauss,n_gauss)
    real*8                                :: delta_t(n_plane,n_eq_var,n_gauss,n_gauss)

    call element_matrix(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, &
                        R_xpoint, Z_xpoint, ELM, RHS, tid, i_tor_min, i_tor_max, aux_nodes)
  end subroutine element_matrix_fft
end module mod_elt_matrix_fft
