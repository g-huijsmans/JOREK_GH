!> Small, independently testable kernel for the model011 Poisson bilinear form.
module mod_poisson_element_kernel
  implicit none
  private
  public :: accumulate_poisson_n0_block, accumulate_poisson_nonzero_blocks
  public :: scatter_poisson_harmonics
  public :: accumulate_poisson_load_mass, apply_poisson_load_harmonics
contains

  !> Add one poloidal Gaussian point to the original n=0 Poisson block.
  subroutine accumulate_poisson_n0_block(weight,big_r,xjac,factor,bb2,psi_x,psi_y, &
       filter_perp,filter_hyper,filter_parallel,value,deriv_x,deriv_y,laplace_star,block_n0)
    real*8, intent(in) :: weight,big_r,xjac,factor,bb2,psi_x,psi_y
    real*8, intent(in) :: filter_perp,filter_hyper,filter_parallel
    real*8, intent(in) :: value(:),deriv_x(:),deriv_y(:),laplace_star(:)
    real*8, intent(inout) :: block_n0(:,:)
    integer :: i,j
    real*8 :: bgrad_i,bgrad_j

    do j=1,size(value)
      bgrad_j=(deriv_x(j)*psi_y-deriv_y(j)*psi_x)/big_r
      do i=1,size(value)
        bgrad_i=(deriv_x(i)*psi_y-deriv_y(i)*psi_x)/big_r
        block_n0(i,j)=block_n0(i,j)+weight*xjac*big_r*( &
             (factor+filter_perp)*(deriv_x(i)*deriv_x(j)+deriv_y(i)*deriv_y(j)) &
             +filter_hyper*laplace_star(i)*laplace_star(j) &
             +filter_parallel*bgrad_i*bgrad_j/bb2)
      enddo
    enddo
  end subroutine accumulate_poisson_n0_block

  !> Add one Gaussian point to A, B and C for a nonzero real Fourier pair.
  !! K_cc=K_ss=A+k^2 B, K_cs=k C and K_sc=-k C. C is skew-symmetric,
  !! so the complete real operator remains symmetric. The -factor part of
  !! parallel_coefficient is the field-parallel subtraction in grad_perp.
  subroutine accumulate_poisson_nonzero_blocks(weight,big_r,xjac,factor,bb2,f0,psi_x,psi_y, &
       filter_perp,filter_hyper,filter_parallel,value,deriv_x,deriv_y,laplace_star, &
       block_a,block_b,block_c)
    real*8, intent(in) :: weight,big_r,xjac,factor,bb2,f0,psi_x,psi_y
    real*8, intent(in) :: filter_perp,filter_hyper,filter_parallel
    real*8, intent(in) :: value(:),deriv_x(:),deriv_y(:),laplace_star(:)
    real*8, intent(inout) :: block_a(:,:),block_b(:,:),block_c(:,:)
    integer :: i,j
    real*8 :: bgrad_i,bgrad_j,toroidal_b,parallel_coefficient,prefactor

    toroidal_b=f0/big_r**2
    parallel_coefficient=(filter_parallel-factor)/bb2
    prefactor=weight*xjac*big_r
    do j=1,size(value)
      bgrad_j=(deriv_x(j)*psi_y-deriv_y(j)*psi_x)/big_r
      do i=1,size(value)
        bgrad_i=(deriv_x(i)*psi_y-deriv_y(i)*psi_x)/big_r
        block_a(i,j)=block_a(i,j)+prefactor*( &
             (factor+filter_perp)*(deriv_x(i)*deriv_x(j)+deriv_y(i)*deriv_y(j)) &
             +filter_hyper*laplace_star(i)*laplace_star(j) &
             +parallel_coefficient*bgrad_i*bgrad_j)
        block_b(i,j)=block_b(i,j)+prefactor*( &
             factor*value(i)*value(j)/big_r**2 &
             +parallel_coefficient*toroidal_b**2*value(i)*value(j))
        block_c(i,j)=block_c(i,j)+prefactor*parallel_coefficient*toroidal_b*( &
             bgrad_i*value(j)-value(i)*bgrad_j)
      enddo
    enddo
  end subroutine accumulate_poisson_nonzero_blocks

  !> Scatter the physical real-Fourier integral into selected components.
  subroutine scatter_poisson_harmonics(block_n0,block_a,block_b,block_c, &
       component_modes,component_types,elm)
    real*8, intent(in) :: block_n0(:,:),block_a(:,:),block_b(:,:),block_c(:,:)
    integer, intent(in) :: component_modes(:)
    character(len=*), intent(in) :: component_types(:)
    real*8, intent(inout) :: elm(:,:)
    integer :: basis_size,component,other,i,j,row,col
    real*8 :: toroidal_norm,k

    basis_size=size(block_a,1)
    if (size(component_types).ne.size(component_modes)) error stop 'Invalid Poisson harmonic metadata.'
    do component=1,size(component_modes)
      k=real(component_modes(component),8)
      if (component_modes(component).eq.0) then
        if (component_types(component).ne.'cos') error stop 'The n=0 Poisson component must be cosine.'
        toroidal_norm=2.d0*acos(-1.d0)
      else
        if (component_types(component).ne.'cos' .and. component_types(component).ne.'sin') &
          error stop 'Unknown real-Fourier Poisson component type.'
        toroidal_norm=acos(-1.d0)
      endif
      do j=1,basis_size
        col=(j-1)*size(component_modes)+component
        do i=1,basis_size
          row=(i-1)*size(component_modes)+component
          if (component_modes(component).eq.0) then
            elm(row,col)=toroidal_norm*block_n0(i,j)
          else
            elm(row,col)=toroidal_norm*(block_a(i,j)+k*k*block_b(i,j))
          endif
        enddo
      enddo

      if (component_modes(component).eq.0) cycle
      do other=1,size(component_modes)
        if (component_modes(other).ne.component_modes(component)) cycle
        if (component_types(other).eq.component_types(component)) cycle
        do j=1,basis_size
          col=(j-1)*size(component_modes)+other
          do i=1,basis_size
            row=(i-1)*size(component_modes)+component
            if (component_types(component).eq.'cos') then
              elm(row,col)=toroidal_norm*k*block_c(i,j)
            else
              elm(row,col)=-toroidal_norm*k*block_c(i,j)
            endif
          enddo
        enddo
      enddo
    enddo
  end subroutine scatter_poisson_harmonics

  subroutine accumulate_poisson_load_mass(weight,big_r,xjac,value,load_mass)
    real*8, intent(in) :: weight,big_r,xjac
    real*8, intent(in) :: value(:)
    real*8, intent(inout) :: load_mass(:,:)
    integer :: i,j
    do j=1,size(value)
      do i=1,size(value)
        load_mass(i,j)=load_mass(i,j)+weight*big_r*xjac*value(i)*value(j)
      enddo
    enddo
  end subroutine accumulate_poisson_load_mass

  !> Apply the physical real-Fourier mass integral to projected coefficients.
  subroutine apply_poisson_load_harmonics(load_mass,charge,component_modes, &
       component_types,n_plane,density_norm,element_rhs)
    real*8, intent(in) :: load_mass(:,:),charge(:,:),density_norm
    integer, intent(in) :: component_modes(:),n_plane
    character(len=*), intent(in) :: component_types(:)
    real*8, intent(out) :: element_rhs(:,:)
    integer :: component
    real*8 :: toroidal_norm

    if (n_plane.lt.1) error stop 'Invalid Poisson plane count.'
    if (size(charge,1).ne.size(load_mass,2)) error stop 'Invalid Poisson charge basis size.'
    if (size(charge,2).ne.size(component_modes)) error stop 'Invalid Poisson charge harmonic count.'
    if (size(element_rhs,1).ne.size(load_mass,1) .or. &
        size(element_rhs,2).ne.size(component_modes)) error stop 'Invalid Poisson element RHS size.'
    if (size(component_types).ne.size(component_modes)) error stop 'Invalid Poisson harmonic metadata.'
    if (density_norm.le.0.d0) error stop 'Invalid Poisson density normalization.'

    do component=1,size(component_modes)
      if (component_modes(component).eq.0) then
        if (component_types(component).ne.'cos') error stop 'The n=0 Poisson component must be cosine.'
        toroidal_norm=2.d0*acos(-1.d0)
      else
        if (component_types(component).ne.'cos' .and. component_types(component).ne.'sin') &
          error stop 'Unknown real-Fourier Poisson component type.'
        toroidal_norm=acos(-1.d0)
      endif
      element_rhs(:,component)=toroidal_norm*matmul(load_mass,charge(:,component))/density_norm
    enddo
  end subroutine apply_poisson_load_harmonics
end module mod_poisson_element_kernel
