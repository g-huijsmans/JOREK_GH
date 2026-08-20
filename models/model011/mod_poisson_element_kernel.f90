!> Small, independently testable kernel for the model011 Poisson bilinear form.
module mod_poisson_element_kernel
  implicit none
  private
  public :: accumulate_poisson_blocks, scatter_poisson_harmonics
  public :: accumulate_poisson_load_mass, apply_poisson_load_harmonics
contains

  !> Add one poloidal Gaussian point to A and B in K_n = A + mode(n)^2 B.
  subroutine accumulate_poisson_blocks(weight, big_r, xjac, factor, bb2, psi_x, psi_y, &
                                       filter_hyper, filter_parallel, value, deriv_x, deriv_y, &
                                       laplace_star, block_a, block_b)
    real*8, intent(in)    :: weight, big_r, xjac, factor, bb2, psi_x, psi_y
    real*8, intent(in)    :: filter_hyper, filter_parallel
    real*8, intent(in)    :: value(:), deriv_x(:), deriv_y(:), laplace_star(:)
    real*8, intent(inout) :: block_a(:,:), block_b(:,:)
    integer               :: i, j
    real*8                 :: bgrad_i, bgrad_j

    do j = 1, size(value)
      bgrad_j = (deriv_x(j)*psi_y - deriv_y(j)*psi_x) / big_r
      do i = 1, size(value)
        bgrad_i = (deriv_x(i)*psi_y - deriv_y(i)*psi_x) / big_r
        block_a(i,j) = block_a(i,j) + weight*xjac*big_r * ( &
            factor*(deriv_x(i)*deriv_x(j) + deriv_y(i)*deriv_y(j)) &
          + filter_hyper*laplace_star(i)*laplace_star(j)           &
          + filter_parallel*bgrad_i*bgrad_j/bb2)
        block_b(i,j) = block_b(i,j) + weight*xjac/big_r * factor*value(i)*value(j)
      enddo
    enddo
  end subroutine accumulate_poisson_blocks


  !> Scatter only diagonal real-Fourier blocks for the selected components.
  subroutine scatter_poisson_harmonics(block_a, block_b, component_modes, component_types, n_plane, elm)
    real*8, intent(in)    :: block_a(:,:), block_b(:,:)
    integer, intent(in)   :: component_modes(:), n_plane
    character(len=*), intent(in) :: component_types(:)
    real*8, intent(inout) :: elm(:,:)
    integer               :: basis_size, component, i, j, row, col
    real*8                 :: toroidal_norm, mode_squared

    basis_size = size(block_a,1)
    if (size(component_types).ne.size(component_modes)) error stop 'Invalid Poisson harmonic metadata.'
    do component = 1, size(component_modes)
      if (component_modes(component).eq.0) then
        if (component_types(component).ne.'cos') error stop 'The n=0 Poisson component must be cosine.'
        toroidal_norm = real(n_plane,8)
      else
        if (component_types(component).ne.'cos' .and. component_types(component).ne.'sin') &
          error stop 'Unknown real-Fourier Poisson component type.'
        toroidal_norm = 0.5d0*real(n_plane,8)
      endif
      mode_squared = real(component_modes(component),8)**2
      do j = 1, basis_size
        col = (j-1)*size(component_modes) + component
        do i = 1, basis_size
          row = (i-1)*size(component_modes) + component
          elm(row,col) = toroidal_norm*(block_a(i,j) + mode_squared*block_b(i,j))
        enddo
      enddo
    enddo
  end subroutine scatter_poisson_harmonics


  !> Add one poloidal Gaussian point to the mass-like Poisson load operator.
  !! The legacy model011 source was v*aux_rhs*R*J/(n0*1e20).
  subroutine accumulate_poisson_load_mass(weight, big_r, xjac, value, load_mass)
    real*8, intent(in)    :: weight, big_r, xjac
    real*8, intent(in)    :: value(:)
    real*8, intent(inout) :: load_mass(:,:)
    integer               :: i, j

    do j = 1, size(value)
      do i = 1, size(value)
        load_mass(i,j) = load_mass(i,j) + weight*big_r*xjac*value(i)*value(j)
      enddo
    enddo
  end subroutine accumulate_poisson_load_mass


  !> Apply the legacy plane-summed real-Fourier load to charge coefficients.
  subroutine apply_poisson_load_harmonics(load_mass, charge, component_modes, &
                                          component_types, n_plane, density_norm, element_rhs)
    real*8, intent(in)    :: load_mass(:,:), charge(:,:), density_norm
    integer, intent(in)   :: component_modes(:), n_plane
    character(len=*), intent(in) :: component_types(:)
    real*8, intent(out)   :: element_rhs(:,:)
    integer               :: component
    real*8                 :: toroidal_norm

    if (size(charge,1).ne.size(load_mass,2)) error stop 'Invalid Poisson charge basis size.'
    if (size(charge,2).ne.size(component_modes)) error stop 'Invalid Poisson charge harmonic count.'
    if (size(element_rhs,1).ne.size(load_mass,1) .or. &
        size(element_rhs,2).ne.size(component_modes)) error stop 'Invalid Poisson element RHS size.'
    if (size(component_types).ne.size(component_modes)) error stop 'Invalid Poisson harmonic metadata.'
    if (density_norm.le.0.d0) error stop 'Invalid Poisson density normalization.'

    do component = 1, size(component_modes)
      if (component_modes(component).eq.0) then
        if (component_types(component).ne.'cos') error stop 'The n=0 Poisson component must be cosine.'
        toroidal_norm = real(n_plane,8)
      else
        if (component_types(component).ne.'cos' .and. component_types(component).ne.'sin') &
          error stop 'Unknown real-Fourier Poisson component type.'
        toroidal_norm = 0.5d0*real(n_plane,8)
      endif
      element_rhs(:,component) = toroidal_norm*matmul(load_mass,charge(:,component))/density_norm
    enddo
  end subroutine apply_poisson_load_harmonics
end module mod_poisson_element_kernel
