!> Particle pusher module with the modified Qin variational scheme.
!> This module contains routines for pushing particles in the RZPhi
!>  (cylindrical) coordinate system. 
module mod_gc_variational
  use mod_fields, only: calc_gyro_average_E_common
  use mod_gk_gc_boundary, only: gc_boundary_geometry, boundary_winding
  use mod_particle_types
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT
  implicit none  
  public copy_particle_gc_vpar
  public copy_particle_gc_Qin
  public copy_particle_gc_Qin_to_vpar
  public initialise_gc_Qin
  public push_gc_Qin
  public push_gc_rk4
  integer, parameter :: GC_STEP_VALID=0, GC_WALL_LOSS=1, GC_NUMERICAL_FAILURE=2
  public convert_leapfrog_to_gc_vpar
  public convert_gc_to_gc_vpar
! Gyro-point and ring statuses are distinct from GC step/loss statuses above.
  integer, parameter :: GYRO_VALID=10, GYRO_FALLBACK=11, GYRO_OUTSIDE=12, GYRO_UNRESOLVED=13, GYRO_INVALID=14
  integer, parameter :: RING_VALID=20, RING_PARTIAL_OUTSIDE=21, RING_OUTSIDE=22, RING_FAILURE=23
  integer, parameter :: GYRO_LOCATION_FAILURE=-1002, GYRO_NOT_LOCATED=-1003
  type :: gyro_point_result
    integer :: classification=GYRO_INVALID, locator_status=GYRO_NOT_LOCATED, locator_element=0, element=0
    real*8 :: locator_st(2)=0.d0
    logical :: recovered=.false.
  end type
  type :: gyro_ring_result
    integer :: classification=RING_FAILURE, n_valid=0, n_outside=0, n_unresolved=0, n_invalid=0, n_recovered=0
  end type
contains

! Shared deposition decision; numerical/empty rings remain the caller's failure.
logical function ion_ring_deposits(ring) result(accept)
 use phys_module, only: gk_ion_boundary_policy
 type(gyro_ring_result),intent(in)::ring
 accept=ring%classification==RING_VALID
 if(gk_ion_boundary_policy=='truncate')accept=accept.or.ring%classification==RING_PARTIAL_OUTSIDE
end function ion_ring_deposits

subroutine copy_particle_gc_vpar(particle_in, particle_out)
implicit none
type(particle_gc_vpar), intent(in)  :: particle_in
type(particle_gc_vpar)              :: particle_out
    
particle_out%i_elm = particle_in%i_elm
particle_out%mu    = particle_in%mu
particle_out%x     = particle_in%x
particle_out%vpar  = particle_in%vpar
particle_out%st    = particle_in%st
particle_out%q     = particle_in%q
particle_out%weight= particle_in%weight
particle_out%B_norm= particle_in%B_norm

return
end

! Preserve raw locator evidence, and only recover an independently inside point.
! Geometry is built once by the caller outside particle/OpenMP loops.
subroutine classify_gyro_point(nodes,elements,point,raw_status,result,geometry)
use data_structure, only: type_node_list,type_element_list
use mod_interp, only: interp_RZ
use mod_find_rz_nearby, only: find_RZ_nearby
use phys_module, only: find_RZ_nearby_tol
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
type(type_node_list),intent(in)::nodes
type(type_element_list),intent(in)::elements
type(particle_kinetic_leapfrog),intent(inout)::point
integer,intent(in)::raw_status
type(gyro_point_result),intent(out)::result
type(gc_boundary_geometry),optional,intent(in)::geometry
integer::winding,e,status,e2
real*8::r,z,s,t,rs,rt,zs,zt,s2,t2
logical::ok
result=gyro_point_result()
result%locator_status=raw_status;result%locator_element=point%i_elm;result%locator_st=point%st
if(.not.all(ieee_is_finite(point%x)))then
  point%i_elm=0
  return
endif
result%classification=GYRO_UNRESOLVED
if(point%i_elm>0.and.point%i_elm<=elements%n_elements.and.all(ieee_is_finite(point%st)))then
  if(minval(point%st)>=0.d0.and.maxval(point%st)<=1.d0)then
    ! Ordinary successful locator retains its own convergence checks.
    if(raw_status==0)then
      result%classification=GYRO_VALID;result%element=point%i_elm
      return
    endif
    if(valid_location(point%i_elm,point%st))then
      result%classification=GYRO_FALLBACK;result%element=point%i_elm
      return
    endif
  endif
endif
if(present(geometry))then
  call boundary_winding(geometry,point%x(:2),winding,ok)
  if(ok)then
    if(winding==0)then
      result%classification=GYRO_OUTSIDE
      ! Negative metadata is evidence only; it is never a GC loss.
      point%i_elm=min(0,point%i_elm)
      result%element=point%i_elm
      return
    endif
    if(winding==1)then
      e=0
      call find_RZ(nodes,elements,point%x(1),point%x(2),r,z,e,s,t,status)
      if(e>0.and.e<=elements%n_elements.and.all(ieee_is_finite([s,t])))then
        if(min(s,t)>=0.d0.and.max(s,t)<=1.d0)then
          call interp_RZ(nodes,elements,e,s,t,r,rs,rt,z,zs,zt)
          if(all(ieee_is_finite([r,z,rs,rt,zs,zt])).and.rs*zt-rt*zs>0.d0)then
            call find_RZ_nearby(nodes,elements,r,z,s,t,e,point%x(1),point%x(2),s2,t2,e2,status)
            if(valid_location(e2,[s2,t2]))then
              point%i_elm=e2;point%st=[s2,t2]
              result%classification=GYRO_FALLBACK;result%element=e2;result%recovered=.true.
              return
            endif
          endif
        endif
      endif
    endif
  endif
endif
point%i_elm=0
contains
logical function valid_location(elm,st) result(valid)
integer,intent(in)::elm
real*8,intent(in)::st(2)
real*8::rr,zz,a,b,c,d
valid=.false.
if(elm<=0.or.elm>elements%n_elements.or..not.all(ieee_is_finite(st)))return
if(minval(st)<0.d0.or.maxval(st)>1.d0)return
call interp_RZ(nodes,elements,elm,st(1),st(2),rr,a,b,zz,c,d)
if(.not.all(ieee_is_finite([rr,zz,a,b,c,d])))return
valid=a*d-b*c>0.d0.and.(rr-point%x(1))**2+(zz-point%x(2))**2<=find_RZ_nearby_tol
end function
end subroutine

subroutine accumulate_gyro_status(ring,point)
type(gyro_ring_result),intent(inout)::ring
type(gyro_point_result),intent(in)::point
select case(point%classification)
case(GYRO_VALID,GYRO_FALLBACK);ring%n_valid=ring%n_valid+1
case(GYRO_OUTSIDE);ring%n_outside=ring%n_outside+1
case(GYRO_UNRESOLVED);ring%n_unresolved=ring%n_unresolved+1
case default;ring%n_invalid=ring%n_invalid+1
end select
if(point%recovered)ring%n_recovered=ring%n_recovered+1
end subroutine

subroutine finish_gyro_status(ring,status)
use mod_pusher_tools,only:GC_GEOMETRY_INVALID
type(gyro_ring_result),intent(inout)::ring
integer,intent(out)::status
ring%classification=RING_FAILURE
status=GYRO_LOCATION_FAILURE
if(ring%n_invalid>0)then
  status=GC_GEOMETRY_INVALID
else if(ring%n_unresolved>0)then
  return
else if(ring%n_outside>0)then
  status=-1
  ring%classification=RING_PARTIAL_OUTSIDE
  if(ring%n_valid==0)ring%classification=RING_OUTSIDE
else if(ring%n_valid>0)then
  status=0
  ring%classification=RING_VALID
endif
end subroutine

subroutine convert_gc_to_gc_vpar(particle_in, B_norm, mass, particle_out)
  implicit none

  type(particle_gc), intent(in) :: particle_in
  real*8, intent(in)            :: B_norm      !< Norm magnetic field at  guiding (gyro) center position [T]
  real*8, intent(in)            :: mass        !< Mass of the particle [amu]
  type(particle_gc_vpar),intent(out) :: particle_out

  real*8 :: v2, v_par

  particle_out%i_elm  = particle_in%i_elm
  particle_out%x      = particle_in%x
  particle_out%st     = particle_in%st
  particle_out%q      = particle_in%q
  particle_out%weight = particle_in%weight

  v2 = 2.d0 * particle_in%E * EL_CHG / (mass*ATOMIC_MASS_UNIT) ![m/s]

  v_par = sqrt(v2 - 2.d0 * abs(particle_in%mu) * B_norm * EL_CHG /  (mass * ATOMIC_MASS_UNIT))    

  particle_out%vpar   = sign(v_par, particle_in%mu)
  particle_out%mu     = abs(particle_in%mu) * EL_CHG / (mass*ATOMIC_MASS_UNIT)
  particle_out%B_norm = B_norm

end

subroutine convert_leapfrog_to_gc_vpar(node_list, element_list, particle_in, B, mass, particle_out, my_ifail)
  use data_structure
  use mod_pusher_tools, only: particle_position_to_gc, GC_GEOMETRY_INVALID
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mod_find_rz_nearby
  implicit none

  type(type_node_list), intent(in)            :: node_list
  type(type_element_list), intent(in)         :: element_list
  type(particle_kinetic_leapfrog), intent(in) :: particle_in
  real*8, dimension(3), intent(in)            :: B !< Magnetic field at kinetic position [T]
  real*8, intent(in)                          :: mass !< Mass of the particle [amu]
  type(particle_gc_vpar)                      :: particle_out
  real*8  :: B_hat(3), B_norm, v_par, v2
  integer :: ifail
  integer, optional, intent(out) :: my_ifail

  call copy_particle_base(particle_in, particle_out)

  particle_out%q      = particle_in%q
  particle_out%weight = particle_in%weight

  if (present(my_ifail)) my_ifail = GC_GEOMETRY_INVALID
  particle_out%vpar = 0.d0
  particle_out%mu = 0.d0
  particle_out%B_norm = 0.d0
  if (.not. all(ieee_is_finite(B)) .or. .not. all(ieee_is_finite(particle_in%v)) .or. &
      .not. ieee_is_finite(mass) .or. mass <= 0.d0) then
    particle_out%i_elm = 0
    return
  endif
  B_norm = norm2(B)
  if (.not. ieee_is_finite(B_norm) .or. B_norm <= 0.d0) then
    particle_out%i_elm = 0
    return
  endif
  B_hat  = B/B_norm
  v_par  = dot_product(particle_in%v,B_hat)
  v2     = dot_product(particle_in%v,particle_in%v)

  ! Inverse displacement uses B and velocity in the kinetic position's basis.
  call particle_position_to_gc(node_list,element_list,particle_in%x,particle_in%st, &
      particle_in%i_elm,mass*particle_in%v,particle_in%q,B_hat,B_norm, &
      particle_out%x,particle_out%st,particle_out%i_elm)

  ! Calculate velocity-related variables
  particle_out%vpar  = v_par
  ! Perhaps we could also use the magnetic field in the guiding center to calculate
  ! v_perp, but don't do it for now.
  particle_out%mu = 0.5d0 * sum((particle_in%v-v_par*B_hat)**2) / B_norm

  particle_out%B_norm = B_norm
  if (present(my_ifail)) then
    my_ifail = 0
    if (particle_out%i_elm <= 0) my_ifail = -1
  endif
end

subroutine copy_particle_gc_Qin(particle_in,particle_out)
implicit none
class(particle_gc_Qin), intent(in)  :: particle_in
class(particle_gc_Qin), intent(out) :: particle_out
  
particle_out%i_elm = particle_in%i_elm
particle_out%mu    = particle_in%mu
particle_out%x     = particle_in%x
particle_out%vpar  = particle_in%vpar
particle_out%st    = particle_in%st
particle_out%q     = particle_in%q
particle_out%weight= particle_in%weight

particle_out%x_m      = particle_in%x_m
particle_out%vpar_m   = particle_in%vpar_m
particle_out%Astar_m  = particle_in%Astar_m 
particle_out%Astar_k  = particle_in%Astar_k 
particle_out%dAstar_k = particle_in%dAstar_k
particle_out%Bn_k     = particle_in%Bn_k
particle_out%dBn_k    = particle_in%dBn_k
particle_out%Bnorm_k  = particle_in%Bnorm_k
particle_out%E_k      = particle_in%E_k

return
end

subroutine copy_particle_gc_Qin_to_vpar(particle_in,particle_out)
use mod_particle_types
implicit none
class(particle_gc_Qin), intent(in)   :: particle_in
class(particle_gc_vpar), intent(out) :: particle_out
  
particle_out%i_elm = particle_in%i_elm
particle_out%mu    = particle_in%mu
particle_out%x     = particle_in%x
particle_out%vpar  = particle_in%vpar
particle_out%st    = particle_in%st
particle_out%q     = particle_in%q
particle_out%weight= particle_in%weight

return
end

subroutine copy_particle_gc_Vpar_to_Qin(particle_in,particle_out)
use mod_particle_types
implicit none
class(particle_gc_Vpar), intent(in) :: particle_in
class(particle_gc_Qin), intent(out) :: particle_out
    
particle_out%i_elm = particle_in%i_elm
particle_out%mu    = particle_in%mu
particle_out%x     = particle_in%x
particle_out%vpar  = particle_in%vpar
particle_out%st    = particle_in%st
particle_out%q     = particle_in%q
particle_out%weight= particle_in%weight

return
end
   
subroutine convert_gc_vpar_to_kinetic(node_list, element_list, particle_in, B, mass, n_phases, particle_out, my_ifail, gyro_phase, point_results, ring_result, boundary_geometry)
  use constants
  use data_structure
  use mod_pusher_tools, only: get_orthonormals, larmor_displacement, cylindrical_displacement, &
      GC_GEOMETRY_INVALID
  use mod_coordinate_transforms, only: vector_rotation
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mod_find_rz_nearby
  implicit none

  type(type_node_list), intent(in)    :: node_list
  type(type_element_list), intent(in) :: element_list
  type(particle_gc_vpar), intent(in)  :: particle_in
  real*8, intent(in)   :: B(3)        !< Magnetic field at GC position [T]
  real*8, intent(in)   :: mass        !< Mass of the particle [amu]
  integer, intent(in)  :: n_phases    !< number of points of the gyro orbit
  type(particle_kinetic_leapfrog), intent(out)  :: particle_out(:)
  integer, intent(out) :: my_ifail
  ! Optional deterministic phase in radians; absent retains exactly one RNG draw.
  real*8, optional, intent(in) :: gyro_phase
  type(gyro_point_result), optional, intent(out) :: point_results(:)
  type(gyro_ring_result), optional, intent(out) :: ring_result
  type(gc_boundary_geometry), optional, intent(in) :: boundary_geometry
  type(gyro_point_result) :: point
  type(gyro_ring_result) :: ring

  real*8  :: B_norm, v_perp, v_par, B_hat(3), e1(3), e2(3), chi, chi_start
  real*8 :: rho(3), delta_phi, velocity_at_gc(3)
  integer :: i, ifail

  if(present(ring_result))ring_result=gyro_ring_result()
  if(present(point_results))then
    point_results=gyro_point_result()
    if(size(point_results)<max(1,n_phases))then
      my_ifail=GC_GEOMETRY_INVALID
      return
    endif
  endif
  ring=gyro_ring_result()
  ifail = 0
  my_ifail = 0

  if (n_phases < 0 .or. size(particle_out) < max(1,n_phases)) then
    my_ifail = GC_GEOMETRY_INVALID
    return
  endif

  if (n_phases .eq. 0) then             ! return a kinetic particle at gyro-centre
    particle_out(1)%x      = particle_in%x
    particle_out(1)%st     = particle_in%st
    particle_out(1)%i_elm  = particle_in%i_elm
    particle_out(1)%q      = particle_in%q
    particle_out(1)%weight = particle_in%weight
    if(present(ring_result))then
      ring_result%classification=RING_VALID
      ring_result%n_valid=1
    endif
    if(present(point_results))then
      point_results(1)%classification=GYRO_VALID
      point_results(1)%element=particle_in%i_elm
      point_results(1)%locator_element=particle_in%i_elm
      point_results(1)%locator_st=particle_in%st
    endif
    return
  endif


  ! Consume the same draw as the original normal path, before validating data.
  if (present(gyro_phase)) then
    chi_start = gyro_phase
  else
    call random_number(chi_start)
    chi_start = chi_start * TWOPI
  endif
  do i=1,n_phases
    particle_out(i)%x = particle_in%x
    particle_out(i)%st = particle_in%st
    particle_out(i)%i_elm = 0
    particle_out(i)%q = particle_in%q
    particle_out(i)%weight = particle_in%weight
    particle_out(i)%v = 0.d0
  enddo
  my_ifail = GC_GEOMETRY_INVALID
  if(present(ring_result))ring_result%n_invalid=n_phases
  if (.not. all(ieee_is_finite(B)) .or. .not. ieee_is_finite(mass)) return
  if (.not. ieee_is_finite(particle_in%mu) .or. .not. ieee_is_finite(particle_in%vpar)) return
  if (.not. ieee_is_finite(chi_start) .or. mass <= 0.d0 .or. particle_in%mu < 0.d0) return
  B_norm = sqrt(dot_product(B,B))
  if (.not. ieee_is_finite(B_norm) .or. B_norm <= 0.d0) return
  B_hat  = B/B_norm

  v_perp = sqrt(2.d0 * particle_in%mu * B_norm) ! [m/s]
  
  ! Define chi as the angle of the velocity vector with b x r
  call get_orthonormals(B_hat, e1, e2)

  my_ifail = 0

  do i = 1, n_phases
    
    ifail=0
    chi = chi_start + real(i-1,8) / real(n_phases,8) * TWOPI

    particle_out(i)%v  = particle_in%vpar * B_hat + v_perp * (cos(chi) * e1 + sin(chi) * e2)

    particle_out(i)%q      = particle_in%q
    particle_out(i)%weight = particle_in%weight

    velocity_at_gc = particle_out(i)%v
    call larmor_displacement(mass*velocity_at_gc,particle_in%q,B_hat,B_norm,rho,ifail)
    if (ifail == 0) call cylindrical_displacement(particle_in%x,rho,particle_out(i)%x,delta_phi,ifail)
    if (ifail /= 0) then
      my_ifail = GC_GEOMETRY_INVALID
      particle_out(i)%i_elm = 0
      ring%n_invalid=ring%n_invalid+1
      if(present(point_results))point_results(i)=gyro_point_result()
      cycle
    endif
    particle_out(i)%v = vector_rotation(velocity_at_gc,-delta_phi)
    if (particle_in%q /= 0) then
      call find_RZ_nearby(node_list, element_list, &
          particle_in%x(1),particle_in%x(2),particle_in%st(1),particle_in%st(2),particle_in%i_elm, &
          particle_out(i)%x(1),particle_out(i)%x(2),particle_out(i)%st(1),particle_out(i)%st(2), &
          particle_out(i)%i_elm,ifail)

    else
      particle_out(i)%i_elm = particle_in%i_elm
    endif
    if(ifail==0.and.particle_out(i)%i_elm>0)then
      ! Ordinary convergence is already checked by find_RZ_nearby. Avoid an
      ! exceptional-classifier call and unused record construction per point.
      ring%n_valid=ring%n_valid+1
      if(present(point_results))then
        point_results(i)%classification=GYRO_VALID
        point_results(i)%locator_status=ifail
        point_results(i)%locator_element=particle_out(i)%i_elm
        point_results(i)%locator_st=particle_out(i)%st
        point_results(i)%element=particle_out(i)%i_elm
      endif
    else
      call classify_gyro_point(node_list,element_list,particle_out(i),ifail,point,boundary_geometry)
      if(present(point_results))point_results(i)=point
      call accumulate_gyro_status(ring,point)
    endif
  enddo
  call finish_gyro_status(ring,my_ifail)
  if(present(ring_result))ring_result=ring

end subroutine convert_gc_vpar_to_kinetic


subroutine initialise_gc_Qin(fields, particle_Qin, mass, timestep)
! initialising backward in time
use mod_particle_types
use mod_fields, only: fields_base
implicit none
class(fields_base)      :: fields
type(particle_gc_Qin)  :: particle_Qin
real*8, intent(in)     :: timestep ! [s]
real*8, intent(in)     :: mass     ! [amu]
type(particle_gc_Vpar) :: particle_Vpar

real*8 :: A_m(3), dA_m(3,3), B_m(3), dB_m(3,3), E_m(3), bn_m, dbn_m(3), bnorm_m(3), dbnorm_m(3,3)
real*8 :: A_k(3), dA_k(3,3), B_k(3), dB_k(3,3), E_k(3), bn_k, dbn_k(3), bnorm_k(3), dbnorm_k(3,3)
real*8 :: qom, time_0

time_0 = 0.d0
qom = particle_Qin%q * EL_CHG / (mass * ATOMIC_MASS_UNIT) 

call fields%calc_Qin(time_0, particle_Qin%i_elm, particle_Qin%st, particle_Qin%x(3), &
                     A_k, dA_k, B_k, dB_k, Bnorm_k, dBnorm_k, bn_k, dbn_k, E_k)
!call fields%calc_Qin_analytic(particle_Qin%x(1),particle_Qin%x(2), particle_Qin%x(3), A_k, dA_k, B_k, dB_k, Bnorm_k, dBnorm_k, bn_k, dbn_k, E_k)

particle_Qin%Astar_k  = qom *  A_k + particle_Qin%vpar *  bnorm_k
particle_Qin%dAstar_k = qom * dA_k + particle_Qin%vpar * dbnorm_k 
particle_Qin%Bn_k     = Bn_k   
particle_Qin%dBn_k    = dBn_k   
particle_Qin%Bnorm_k  = Bnorm_k 
particle_Qin%E_k      = E_k 

call copy_particle_gc_Qin_to_Vpar(particle_Qin, particle_Vpar)

call push_gc_rk4(fields, particle_Vpar, mass, -timestep, 1, 0)      ! should be only the very first call

particle_Qin%x_m    = particle_Vpar%x
particle_Qin%vpar_m = particle_Vpar%vpar

call fields%calc_Qin(time_0, particle_Vpar%i_elm, particle_Vpar%st, particle_Vpar%x(3), &
                     A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dbn_m, E_m)
!call fields%calc_Qin_analytic(particle_Vpar%x(1), particle_Vpar%x(2), particle_Vpar%x(3), A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dbn_m, E_m)

particle_Qin%vpar_m = (particle_Qin%Astar_k(3) - qom * A_m(3)) / Bnorm_m(3)
particle_Qin%Astar_m = qom * A_m + particle_Qin%vpar_m *  Bnorm_m

!particle_Qin%vpar = (particle_Qin%Astar_m(3) - qom * A_k(3)) / Bnorm_k(3)
!particle_Qin%Astar_k  = qom *  A_k + particle_Qin%vpar *  bnorm_k
!particle_Qin%dAstar_k = qom * dA_k + particle_Qin%vpar * dbnorm_k 

return
end

subroutine push_gc_Qin(fields, particle_Qin, mass, timestep, n_steps)
use nodes_elements
use mod_find_rz_nearby
use mod_fields, only: fields_base
implicit none
class(fields_base)      :: fields
type(particle_gc_Qin)  :: particle_Qin
real*8, intent(in)     :: timestep ! [s]
real*8, intent(in)     :: mass     ! [amu]
integer, intent(in)    :: n_steps
real*8                 :: time_0   ! not yet implemented
real*8                 :: qom

real*8 :: x_k(3), vpar_k, Astar_k(3), st_k(2), x_p(3), st_p(2)
real*8 :: A_p(3), dA_p(3,3), B_p(3), dB_p(3,3), E_p(3), bn_p, dbn_p(3), bnorm_p(3), dbnorm_p(3,3)
real*8 :: Astar_p(3), dAstar_p(3,3)
real*8 :: newton_rhs(3), newton_matrix(3,3), residue, mp_k, wp_k, mp_p, wp_p, mp_m, wp_m
integer :: info, ifail, ipiv(3), iter, i, j, it, i_elm_k, i_elm_p
  
if (particle_Qin%i_elm .le. 0) return

qom = particle_Qin%q * EL_CHG / (mass * ATOMIC_MASS_UNIT) 

do it=1, n_steps

  x_k     = particle_Qin%x
  vpar_k  = particle_Qin%vpar
  Astar_k = particle_Qin%Astar_k
  st_k    = particle_Qin%st
  i_elm_k = particle_Qin%i_elm
  
  !----------------------------- initial guess from linearised equations
  do j=1,3
    newton_rhs(j) = - 2.d0 * particle_Qin%bnorm_k(j) * (particle_Qin%vpar_m + vpar_k) &
                    + 2.d0 * timestep * (particle_Qin%mu * particle_Qin%dBn_k(j) - qom * particle_Qin%E_k(j))
    do i=1,3
      newton_matrix(i,j) = particle_Qin%dAstar_k(i,j) - particle_Qin%dAstar_k(j,i) &
                         - 2.d0 * particle_Qin%bnorm_k(j) *  particle_Qin%bnorm_k(i) / timestep
    enddo
  enddo
  call dgesv(3, 1, newton_matrix, 3, ipiv, newton_rhs, 3, info)

  particle_Qin%x    =   particle_Qin%x_m + newton_rhs
  particle_Qin%vpar = - particle_Qin%vpar_m - 2.d0 * particle_Qin%vpar &
                    + 2.d0 * dot_product(particle_Qin%bnorm_k, particle_Qin%x - particle_Qin%x_m) / timestep

  call find_RZ_nearby(node_list, element_list, x_k(1), x_k(2), st_k(1), st_k(2), i_elm_k, &
                      particle_Qin%x(1), particle_Qin%x(2), particle_Qin%st(1), particle_Qin%st(2), particle_Qin%i_elm, ifail)

  if (particle_Qin%i_elm .le. 0) return

  call fields%calc_Qin(time_0, particle_Qin%i_elm, particle_Qin%st, particle_Qin%x(3), &
                       A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dbn_p, E_p)  
 !call fields%calc_Qin_analytic(particle_Qin%x(1), particle_Qin%x(2), particle_Qin%x(3), A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dbn_p, E_p)  

  Astar_p  = qom *  A_p + particle_Qin%vpar *  bnorm_p
  dAstar_p = qom * dA_p + particle_Qin%vpar * dbnorm_p

  !---------------------------- Newton iterations
  do iter = 1,3
 
     do i=1,3

      newton_rhs(i) = (Astar_p(i) - particle_Qin%Astar_m(i)) &
                    
                    + 2.d0 * timestep * (particle_Qin%mu * particle_Qin%dBn_k(i)  -  qom * particle_Qin%E_k(i)) &

                    - dot_product(particle_Qin%dAstar_k(:,i), particle_Qin%x - particle_Qin%x_m)
     !if (i .eq. 3) write(*,*) it, iter, norm2(newton_rhs)
     do j=1,3
        
        newton_matrix(j,i) = + particle_Qin%dAstar_k(i,j) - qom * dA_p(j,i)           &
                             
                             + dbnorm_p(j,i) * (particle_Qin%vpar_m + 2.d0 * vpar_k)  &
        
                             - 2.d0 * bnorm_p(j) * particle_Qin%bnorm_k(i) / timestep &

                             - 2.d0 * dbnorm_p(j,i) * dot_product(particle_Qin%bnorm_k, particle_Qin%x - particle_Qin%x_m) / timestep
      enddo
    enddo

    call dgesv(3, 1, newton_matrix, 3, ipiv, newton_rhs, 3, info)

    x_p     = particle_Qin%x
    st_p    = particle_Qin%st
    i_elm_p = particle_Qin%i_elm
    
    particle_Qin%x    = particle_Qin%x + newton_rhs
    particle_Qin%vpar = - particle_Qin%vpar_m - 2.d0 * vpar_k + 2.d0 * dot_product(particle_Qin%bnorm_k, particle_Qin%x - particle_Qin%x_m) / timestep
  
    call find_RZ_nearby(node_list, element_list, x_p(1), x_p(2), st_p(1), st_p(2), i_elm_p, &
                        particle_Qin%x(1), particle_Qin%x(2), particle_Qin%st(1), particle_Qin%st(2), particle_Qin%i_elm, ifail)
  
    if (particle_Qin%i_elm .le. 0) return

      call fields%calc_Qin(time_0, particle_Qin%i_elm, particle_Qin%st, particle_Qin%x(3), &
                           A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dbn_p, E_p)  
    !call fields%calc_Qin_analytic(particle_Qin%x(1), particle_Qin%x(2), particle_Qin%x(3), A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dbn_p, E_p)  
  
    Astar_p  = qom *  A_p + particle_Qin%vpar *  bnorm_p
    dAstar_p = qom * dA_p + particle_Qin%vpar * dbnorm_p

  enddo ! newton iterations


  do j=1,3    
    newton_rhs(j) = (Astar_p(j) - particle_Qin%Astar_m(j)) &
                 + 2.d0 * timestep * (particle_Qin%mu * particle_Qin%dBn_k(j) - qom * particle_Qin%E_k(j)) &
                 - dot_product(particle_Qin%dAstar_k(:,j), particle_Qin%x - particle_Qin%x_m)
  enddo
  residue = norm2(newton_rhs)
!  if (residue .gt. 1.d-3) write(*,'(A,i4,e14.6)') ' Qin : residue : ',it, residue

!  mp_p = Astar_p(3)
!  wp_p = 0.25d0 * (vpar_k**2 + particle_Qin%vpar**2) + 0.5d0 * particle_Qin%mu * (particle_Qin%bn_k + bn_p)

  particle_Qin%x_m      = x_k
  particle_Qin%vpar_m   = vpar_k
  particle_Qin%Astar_m  = Astar_k 
  particle_Qin%Astar_k  = Astar_p 
  particle_Qin%dAstar_k = dAstar_p 
  particle_Qin%Bn_k     = Bn_p
  particle_Qin%dBn_k    = dBn_p
  particle_Qin%Bnorm_k  = Bnorm_p
  particle_Qin%E_k      = E_p

enddo

return
end

subroutine push_gc_rk4(fields, particle_gc, mass, timestep, n_steps, n_gyro_phases, gyro_shift, lost_st, lost_phi, gyro_phase, wall_status, locator_retries, &
                       integration_substeps, start_time, failed_substep, failed_stage, stage_time, field_calls, &
                       max_stage_distance, max_cell_distance, failure_reason, gyro_boundary)
use nodes_elements
use mod_pusher_tools, only: GC_GEOMETRY_INVALID
use mod_find_rz_nearby
use mod_fields, only: fields_base
use mod_boundary, only: wall_normal_vector
use mod_interp, only: interp_RZ
use phys_module, only: find_RZ_nearby_iter
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
class(fields_base)     :: fields
type(particle_gc_vpar) :: particle_gc
real*8, intent(in)     :: timestep      ! [s]
real*8, intent(in)     :: mass          ! [amu]
integer, intent(in)    :: n_steps       ! number of time steps 
integer, intent(in)    :: n_gyro_phases ! number of gyro phases for gyro-averaging
real*8, intent(out), optional :: gyro_shift(3) ! shift of the gyro centre from guiding centre
! Optional wall metadata; leave accepted particle state and legacy callers unchanged.
real*8, intent(out), optional :: lost_st(2), lost_phi
! Status follows the existing RK stage-loss convention; it is not a located event time.
integer, optional, intent(out) :: wall_status, locator_retries
integer, intent(in), optional :: integration_substeps
real*8, intent(in), optional :: start_time
! Stage-loss diagnostics identify the RK trial, not a continuous first impact.
integer, intent(out), optional :: failed_substep, failed_stage, field_calls, failure_reason
real*8, intent(out), optional :: stage_time, max_stage_distance, max_cell_distance
type(gc_boundary_geometry), intent(in), optional :: gyro_boundary
type(gyro_ring_result) :: gyro_result
real*8, intent(in), optional :: gyro_phase ! Controlled quadrature offset; absent uses production RNG.

real*8 :: A_0(3), dA_0(3,3), B_0(3), dB_0(3,3), bn_0, dbn_0(3), Bnorm_0(3), dBnorm_0(3,3), E_0(3)
real*8 :: A_1(3), dA_1(3,3), B_1(3), dB_1(3,3), bn_1, dbn_1(3), Bnorm_1(3), dBnorm_1(3,3), E_1(3)
real*8 :: A_2(3), dA_2(3,3), B_2(3), dB_2(3,3), bn_2, dbn_2(3), Bnorm_2(3), dBnorm_2(3,3), E_2(3)
real*8 :: A_3(3), dA_3(3,3), B_3(3), dB_3(3,3), bn_3, dbn_3(3), Bnorm_3(3), dBnorm_3(3,3), E_3(3)
real*8 :: x_0(3), x_1(3), x_2(3), x_3(3), u_0, u_1, u_2, u_3, time_0, time_1, time_2, time_3
real*8 :: delta_x1(3), delta_x2(3), delta_x3(3), delta_x4(3)
real*8 :: delta_u1, delta_u2, delta_u3, delta_u4, qom, p_phi, energy

type(particle_gc_vpar) :: p_0, p_1, p_2, p_3, step_start
type(particle_kinetic_leapfrog), allocatable  :: p_orbit(:) 
integer :: i, ifail, n_phases, n_internal, total_steps, current_stage
real*8 :: dt_gc, initial_time, current_time

if(present(wall_status))wall_status=GC_STEP_VALID
if(present(locator_retries))locator_retries=0

if(present(failure_reason))failure_reason=0
if(present(failed_substep))failed_substep=0
if(present(failed_stage))failed_stage=0
if(present(stage_time))stage_time=0.d0
if(present(field_calls))field_calls=0
if(present(max_stage_distance))max_stage_distance=0.d0
if(present(max_cell_distance))max_cell_distance=0.d0
n_internal=1
if(present(integration_substeps))n_internal=integration_substeps
initial_time=0.d0
if(present(start_time))initial_time=start_time
if(n_internal<1.or.n_steps<1.or..not.ieee_is_finite(timestep).or.timestep<0.d0.or. &
   .not.ieee_is_finite(initial_time))then
  if(present(failure_reason))failure_reason=1
  if(present(wall_status))then
    wall_status=GC_NUMERICAL_FAILURE
    return
  endif
  error stop 'Invalid GC RK4 integration interval'
endif
if(n_steps>huge(total_steps)/n_internal)error stop 'GC RK4 substep count overflow'
total_steps=n_steps*n_internal
dt_gc=timestep/n_internal

qom = particle_gc%q * EL_CHG / (mass * ATOMIC_MASS_UNIT) 

if (present(gyro_shift)) gyro_shift = 0.d0

if (present(lost_st)) lost_st = particle_gc%st
if (present(lost_phi)) lost_phi = particle_gc%x(3)
if (particle_gc%i_elm .le. 0) return

n_phases = n_gyro_phases
if (n_gyro_phases .eq. 0) then     ! evolve guiding centre instead of gyro centre
  n_phases = 1
endif

allocate(p_orbit(n_phases))

call copy_particle_gc_vpar(particle_gc,p_0)
call copy_particle_gc_vpar(particle_gc,p_1)
call copy_particle_gc_vpar(particle_gc,p_2)
call copy_particle_gc_vpar(particle_gc,p_3)

do i =1, total_steps
  time_0=initial_time+(i-1)*dt_gc
  time_1=time_0+0.5d0*dt_gc
  time_2=time_1
  time_3=time_0+dt_gc
  step_start=p_0
  current_stage=0;current_time=time_0
call fields%calc_RK4(time_0, p_0%i_elm, p_0%st, p_0%x(3), A_0, dA_0, B_0, dB_0, Bnorm_0, dBnorm_0, bn_0, dbn_0, E_0)
if(present(field_calls))field_calls=field_calls+1
!call fields%calc_RK4_analytic(p_0%x(1), p_0%x(2), p_0%x(3), A_0, dA_0, B_0, dB_0, Bnorm_0, dBnorm_0, bn_0, dbn_0, E_0)

call convert_gc_vpar_to_kinetic(node_list, element_list, p_0, B_0, mass, n_gyro_phases, p_orbit, ifail, gyro_phase=gyro_phase, ring_result=gyro_result, boundary_geometry=gyro_boundary)
if(gyro_result%classification==RING_FAILURE.or.gyro_result%n_valid==0)then
  if(present(failure_reason))failure_reason=5
  if(present(failed_substep))failed_substep=i
  if(present(failed_stage))failed_stage=current_stage
  if(present(stage_time))stage_time=current_time
  if(present(wall_status))then
    wall_status=GC_NUMERICAL_FAILURE
    return
  endif
  error stop 'Unresolved or empty ion gyro-ring; accepted GC retained'
endif
call gather_ion_common_basis(fields,time_0, p_orbit, n_gyro_phases, E_0, p_0%x(3))


  
  call rk4_step(p_0%x, p_0%vpar, qom, p_0%mu, E_0, B_0, Bnorm_0, dBnorm_0, dBn_0, delta_x1, delta_u1)

  p_1%x    = p_0%x    + 0.5d0 * dt_gc * delta_x1
  p_1%vpar = p_0%vpar + 0.5d0 * dt_gc * delta_u1

  current_stage=1;current_time=time_1
  call locate_rk_stage( p_0%x(1), p_0%x(2), p_0%st(1), p_0%st(2), p_0%i_elm, &
                                               p_1%x(1), p_1%x(2), p_1%st(1), p_1%st(2), p_1%i_elm, ifail)
                                               
  if (p_1%i_elm .le. 0) then
    if(p_1%i_elm==0)then
      if(present(wall_status))then
        wall_status=GC_NUMERICAL_FAILURE
        return ! original accepted particle and loss metadata remain intact
      endif
      error stop 'Unresolved GC RK4 locator failure; particle retained'
    endif
    if(present(wall_status))wall_status=GC_WALL_LOSS
    particle_gc%i_elm = p_1%i_elm
    if (present(lost_st)) lost_st = p_1%st
    if (present(lost_phi)) lost_phi = p_1%x(3)
    return
  endif
    
  call fields%calc_RK4(time_1, p_1%i_elm, p_1%st, p_1%x(3), A_1, dA_1, B_1, dB_1, Bnorm_1, dBnorm_1, bn_1, dbn_1, E_1)
if(present(field_calls))field_calls=field_calls+1
 !call fields%calc_RK4_analytic(p_1%x(1), p_1%x(2), p_1%x(3), A_1, dA_1, B_1, dB_1, Bnorm_1, dBnorm_1, bn_1, dbn_1, E_1)

  call convert_gc_vpar_to_kinetic(node_list, element_list, p_1, B_1, mass, n_gyro_phases, p_orbit, ifail, gyro_phase=gyro_phase, ring_result=gyro_result, boundary_geometry=gyro_boundary)
if(gyro_result%classification==RING_FAILURE.or.gyro_result%n_valid==0)then
  if(present(failure_reason))failure_reason=5
  if(present(failed_substep))failed_substep=i
  if(present(failed_stage))failed_stage=current_stage
  if(present(stage_time))stage_time=current_time
  if(present(wall_status))then
    wall_status=GC_NUMERICAL_FAILURE
    return
  endif
  error stop 'Unresolved or empty ion gyro-ring; accepted GC retained'
endif
  call gather_ion_common_basis(fields,time_1, p_orbit, n_gyro_phases, E_1, p_1%x(3))
  
  call rk4_step(p_1%x, p_1%vpar, qom, p_1%mu, E_1, B_1, Bnorm_1, dBnorm_1, dBn_1, delta_x2, delta_u2)

  p_2%x    = p_0%x    + 0.5d0 * dt_gc * delta_x2
  p_2%vpar = p_0%vpar + 0.5d0 * dt_gc * delta_u2

  current_stage=2;current_time=time_2
  call locate_rk_stage( p_0%x(1), p_0%x(2), p_0%st(1), p_0%st(2), p_0%i_elm, &
                                               p_2%x(1), p_2%x(2), p_2%st(1), p_2%st(2), p_2%i_elm, ifail)
  if (p_2%i_elm .le. 0) then
    if(p_2%i_elm==0)then
      if(present(wall_status))then
        wall_status=GC_NUMERICAL_FAILURE
        return ! original accepted particle and loss metadata remain intact
      endif
      error stop 'Unresolved GC RK4 locator failure; particle retained'
    endif
    if(present(wall_status))wall_status=GC_WALL_LOSS
    particle_gc%i_elm = p_2%i_elm
    if (present(lost_st)) lost_st = p_2%st
    if (present(lost_phi)) lost_phi = p_2%x(3)
    return
  endif
  
  call fields%calc_RK4(time_2, p_2%i_elm, p_2%st, p_2%x(3), A_2, dA_2, B_2, dB_2, Bnorm_2, dBnorm_2, bn_2, dbn_2, E_2)
if(present(field_calls))field_calls=field_calls+1
 !call fields%calc_RK4_analytic(p_2%x(1), p_2%x(2), p_2%x(3), A_2, dA_2, B_2, dB_2, Bnorm_2, dBnorm_2, bn_2, dbn_2, E_2)
            
  call convert_gc_vpar_to_kinetic(node_list, element_list, p_2, B_2, mass, n_gyro_phases, p_orbit, ifail, gyro_phase=gyro_phase, ring_result=gyro_result, boundary_geometry=gyro_boundary)
if(gyro_result%classification==RING_FAILURE.or.gyro_result%n_valid==0)then
  if(present(failure_reason))failure_reason=5
  if(present(failed_substep))failed_substep=i
  if(present(failed_stage))failed_stage=current_stage
  if(present(stage_time))stage_time=current_time
  if(present(wall_status))then
    wall_status=GC_NUMERICAL_FAILURE
    return
  endif
  error stop 'Unresolved or empty ion gyro-ring; accepted GC retained'
endif
  call gather_ion_common_basis(fields,time_2, p_orbit, n_gyro_phases, E_2, p_2%x(3))

  call rk4_step(p_2%x, p_2%vpar, qom, p_2%mu, E_2, B_2, Bnorm_2, dBnorm_2, dBn_2, delta_x3, delta_u3)

  p_3%x    = p_0%x    + dt_gc * delta_x3
  p_3%vpar = p_0%vpar + dt_gc * delta_u3

  current_stage=3;current_time=time_3
  call locate_rk_stage( p_0%x(1), p_0%x(2), p_0%st(1), p_0%st(2), p_0%i_elm, &
                                               p_3%x(1), p_3%x(2), p_3%st(1), p_3%st(2), p_3%i_elm, ifail)

  if (p_3%i_elm .le. 0) then
    if(p_3%i_elm==0)then
      if(present(wall_status))then
        wall_status=GC_NUMERICAL_FAILURE
        return ! original accepted particle and loss metadata remain intact
      endif
      error stop 'Unresolved GC RK4 locator failure; particle retained'
    endif
    if(present(wall_status))wall_status=GC_WALL_LOSS
    particle_gc%i_elm = p_3%i_elm
    if (present(lost_st)) lost_st = p_3%st
    if (present(lost_phi)) lost_phi = p_3%x(3)
    return
  endif
  
  call fields%calc_RK4(time_3, p_3%i_elm, p_3%st, p_3%x(3), A_3, dA_3, B_3, dB_3, Bnorm_3, dBnorm_3, bn_3, dbn_3, E_3)
if(present(field_calls))field_calls=field_calls+1
 !call fields%calc_RK4_analytic(p_3%x(1), p_3%x(2), p_3%x(3), A_3, dA_3, B_3, dB_3, Bnorm_3, dBnorm_3, bn_3, dbn_3, E_3)
    
  call convert_gc_vpar_to_kinetic(node_list, element_list, p_3, B_3, mass, n_gyro_phases, p_orbit, ifail, gyro_phase=gyro_phase, ring_result=gyro_result, boundary_geometry=gyro_boundary)
if(gyro_result%classification==RING_FAILURE.or.gyro_result%n_valid==0)then
  if(present(failure_reason))failure_reason=5
  if(present(failed_substep))failed_substep=i
  if(present(failed_stage))failed_stage=current_stage
  if(present(stage_time))stage_time=current_time
  if(present(wall_status))then
    wall_status=GC_NUMERICAL_FAILURE
    return
  endif
  error stop 'Unresolved or empty ion gyro-ring; accepted GC retained'
endif
  call gather_ion_common_basis(fields,time_3, p_orbit, n_gyro_phases, E_3, p_3%x(3))

  call rk4_step(p_3%x, p_3%vpar, qom, p_3%mu, E_3, B_3, Bnorm_3, dBnorm_3, dBn_3, delta_x4, delta_u4)
                                
  p_0%x    = p_0%x    + dt_gc * (delta_x1 + 2.d0*delta_x2 + 2.d0*delta_x3 + delta_x4)/6.d0
  p_0%vpar = p_0%vpar + dt_gc * (delta_u1 + 2.d0*delta_u2 + 2.d0*delta_u3 + delta_u4)/6.d0

  current_stage=4;current_time=time_3
  call locate_rk_stage( &
                      step_start%x(1), step_start%x(2), step_start%st(1), step_start%st(2), step_start%i_elm, &
                      p_0%x(1),  p_0%x(2),  p_0%st(1),  p_0%st(2),  p_0%i_elm, ifail)

  if (p_0%i_elm .le. 0) then
    if(p_0%i_elm==0)then
      if(present(wall_status))then
        wall_status=GC_NUMERICAL_FAILURE
        return ! original accepted particle and loss metadata remain intact
      endif
      error stop 'Unresolved GC RK4 locator failure; particle retained'
    endif
    if(present(wall_status))wall_status=GC_WALL_LOSS
    particle_gc%i_elm = p_0%i_elm
    if (present(lost_st)) lost_st = p_0%st
    if (present(lost_phi)) lost_phi = p_0%x(3)
    return
  endif

enddo

time_0=initial_time+n_steps*timestep
call fields%calc_RK4(time_0, p_0%i_elm, p_0%st, p_0%x(3), A_0, dA_0, B_0, dB_0, Bnorm_0, dBnorm_0, bn_0, dbn_0, E_0)
if(present(field_calls))field_calls=field_calls+1

if (present(gyro_shift)) then
  call convert_gc_vpar_to_kinetic(node_list, element_list, p_0, B_0, mass, n_gyro_phases, p_orbit, ifail, gyro_phase=gyro_phase, ring_result=gyro_result, boundary_geometry=gyro_boundary)
if(gyro_result%classification==RING_FAILURE.or.gyro_result%n_valid==0)then
  if(present(failure_reason))failure_reason=5
  if(present(failed_substep))failed_substep=i
  if(present(failed_stage))failed_stage=current_stage
  if(present(stage_time))stage_time=current_time
  if(present(wall_status))then
    wall_status=GC_NUMERICAL_FAILURE
    return
  endif
  error stop 'Unresolved or empty ion gyro-ring; accepted GC retained'
endif
  call gather_ion_common_basis(fields,time_0, p_orbit, n_gyro_phases, E_0, p_0%x(3))
  gyro_shift = (E_0 - dot_product(E_0,B_0)*B_0/bn_0**2) / bn_0**2 * mass*ATOMIC_MASS_UNIT / EL_CHG
endif

p_0%B_norm = norm2(B_0)
if(.not.all(ieee_is_finite([p_0%x,p_0%vpar,p_0%mu,p_0%B_norm])).or.p_0%B_norm<=0.d0)then
  if(present(failure_reason))failure_reason=4
  if(present(wall_status))then
    wall_status=GC_NUMERICAL_FAILURE
    return
  endif
  error stop 'Invalid GC RK4 final state; particle retained'
endif

call copy_particle_gc_vpar(p_0,particle_gc)
if(present(failed_substep))failed_substep=0
if(present(failed_stage))failed_stage=0
if(present(stage_time))stage_time=initial_time+n_steps*timestep

return
contains
! Retry the identical RK trial, with the same origin and geometry, only when
! the ordinary locator is unresolved. No RK force, timestep or RNG is repeated.
function current_phi() result(phi)
real*8::phi
select case(current_stage)
case(1);phi=p_1%x(3)
case(2);phi=p_2%x(3)
case(3);phi=p_3%x(3)
case default;phi=p_0%x(3)
end select
end function
logical function finite_trial() result(valid)
select case(current_stage)
case(1);valid=all(ieee_is_finite([p_1%x,p_1%vpar]))
case(2);valid=all(ieee_is_finite([p_2%x,p_2%vpar]))
case(3);valid=all(ieee_is_finite([p_3%x,p_3%vpar]))
case default;valid=all(ieee_is_finite([p_0%x,p_0%vpar]))
end select
end function
subroutine locate_rk_stage(ro,zo,so,to,eo,r,z,s,t,e,status)
real*8,intent(in)::ro,zo,so,to,r,z
integer,intent(in)::eo
real*8,intent(out)::s,t
integer,intent(out)::e,status
integer::budget,side,next
real*8::wall(2),normal(3),tol, rr,zz,rs,rt,zs,zt,jac,ds,dt
logical::valid_wall
if(present(failed_substep))failed_substep=i
if(present(failed_stage))failed_stage=current_stage
if(present(stage_time))stage_time=current_time
if(.not.finite_trial())then
  if(present(failure_reason))failure_reason=2
  e=0;status=999;s=so;t=to
  return
endif
if(present(max_stage_distance))max_stage_distance=max(max_stage_distance,norm2([r*cos(current_phi())-ro*cos(step_start%x(3)), &
  r*sin(current_phi())-ro*sin(step_start%x(3)),z-zo]))
if(present(max_cell_distance))then
  call interp_RZ(node_list,element_list,eo,so,to,rr,rs,rt,zz,zs,zt)
  jac=rs*zt-rt*zs
  if(ieee_is_finite(jac).and.jac>0.d0)then
    ds=((r-ro)*zt-(z-zo)*rt)/jac
    dt=((z-zo)*rs-(r-ro)*zs)/jac
    max_cell_distance=max(max_cell_distance,max(abs(ds),abs(dt)))
  endif
endif
call find_RZ_nearby(node_list,element_list,ro,zo,so,to,eo,r,z,s,t,e,status)
if(e/=0)return
! 16 -> 32 -> 64 in COMPASS. Larger configured defaults are not reduced.
do budget=1,2
  if(present(locator_retries))locator_retries=locator_retries+1
  call find_RZ_nearby(node_list,element_list,ro,zo,so,to,eo,r,z,s,t,e,status, &
                      iteration_limit=max(16,find_RZ_nearby_iter)*2**budget)
  if(e>0)return
  if(e<0)then
    ! Verify that negative metadata names an actual wall face. This retains
    ! the existing RK stage-loss model, not an exact trajectory intersection.
    valid_wall=.false.;tol=128.d0*epsilon(1.d0)
    if(-e>element_list%n_elements.or..not.all(ieee_is_finite([s,t])))exit
    if(min(s,t)<-tol.or.max(s,t)>1.d0+tol)exit
    do side=1,4
      next=mod(side,4)+1
      if(element_list%element(-e)%neighbours(side)/=0)cycle
      if(node_list%node(element_list%element(-e)%vertex(side))%boundary<=0.or. &
         node_list%node(element_list%element(-e)%vertex(next))%boundary<=0)cycle
      select case(side)
      case(1);valid_wall=valid_wall.or.abs(t)<=tol
      case(2);valid_wall=valid_wall.or.abs(s-1.d0)<=tol
      case(3);valid_wall=valid_wall.or.abs(t-1.d0)<=tol
      case(4);valid_wall=valid_wall.or.abs(s)<=tol
      end select
    enddo
    if(.not.valid_wall)exit
    call interp_RZ(node_list,element_list,-e,s,t,wall(1),wall(2))
    normal=wall_normal_vector(node_list,element_list,-e,s,t)
    if(.not.all(ieee_is_finite([wall,normal])).or.norm2(normal)<=0.d0)exit
    if(dot_product([r,z]-wall,normal(:2))>tol*max(1.d0,norm2(wall)))exit
    return
  endif
enddo
if(present(failure_reason))failure_reason=3
e=0;status=999
end subroutine locate_rk_stage
end subroutine push_gc_rk4


function cross(A,B) result(AcrossB)
implicit none
real*8 :: A(3), B(3), AcrossB(3)
    AcrossB(1) = A(2)*B(3) - A(3)*B(2)
    AcrossB(2) = A(3)*B(1) - A(1)*B(3)
    AcrossB(3) = A(1)*B(2) - A(2)*B(1)
  return
end

function rot_tmp(x,A,dA) result(rotA)
implicit none
real*8, intent(in) :: x(3), A(3), dA(3,3)
real*8             :: rotA(3)
  rotA(1) = dA(3,2) - dA(2,3) / x(1)
  rotA(2) = dA(1,3)/x(1) - dA(3,1) - A(3) / x(1)
  rotA(3) = dA(2,1) - dA(1,2) 
return
end    

subroutine rk4_step(x, vpar, qom, zmu, E, B, Bnorm, dBnorm, dB, delta_x, delta_u)
implicit none
real*8 :: x(3), vpar, qom, zmu, A(3), dA(3,3), E(3), Bnorm(3), dBnorm(3,3), B(3), dB(3)
real*8 :: Bstar(3), Estar(3), Bpar_star, delta_x(3), delta_u
  
  Bstar     = B + vpar * rot_tmp(x,Bnorm,dBnorm) / qom
  Bpar_star = dot_product(Bstar,Bnorm)
  Estar(1:2) = E(1:2) - zmu * dB(1:2)    /qom
  Estar(3)   = E(3)   - zmu * dB(3)/x(1) /qom
        
  delta_x = (Bstar * vpar  - cross(Bnorm, Estar)) / Bpar_star
  delta_u = dot_product(Bstar,Estar) * qom        / Bpar_star

  delta_x(3) = delta_x(3) / x(1)

return
end

subroutine rk4_step2(x, vpar, qom, zmu, E, B, Bnorm, dBnorm, dB, delta_x, delta_u)
! identical to rk4_step (but can be useful to identify different velicity contributions)
implicit none
real*8 :: x(3), vpar, qom, zmu, A(3), dA(3,3), E(3), Bnorm(3), dBnorm(3,3), B(3), dB(3)
real*8 :: Bstar(3), Estar(3), Bpar_star, delta_x(3), delta_u
real*8 :: Bn, v_E(3), v_gradB(3), v_curvature(3)
  
  Bn        = norm2(B)
  Bstar     = B + vpar * rot_tmp(x,Bnorm,dBnorm) / qom
  Bpar_star = dot_product(Bstar,Bnorm)
  Estar(1:2) = E(1:2) - zmu * dB(1:2)    /qom
  Estar(3)   = E(3)   - zmu * dB(3)/x(1) /qom

  v_E         =           - cross(B,E)  / Bn**2
  v_gradB     = zmu / qom * cross(B,dB) / Bn**2
  v_curvature = vpar**2 / (Bn * qom) * (rot_tmp(x,Bnorm,dBnorm) - dot_product(rot_tmp(x,Bnorm,dBnorm),Bnorm)* Bnorm)
        
  delta_x = Bnorm * vpar  + Bn / Bpar_star * (v_gradB + v_E + v_curvature)
  delta_u = - dot_product(Bnorm + Bn/(vpar*Bpar_star) * (v_gradB + v_E + v_curvature) , zmu * dB - qom * E)

  delta_x(3) = delta_x(3) / x(1)

return
end


! Candidate scope: concrete JOREK linear FE fields only. Preserve custom callbacks.
subroutine gather_ion_common_basis(fields,time,points,n,E,phi)
 use phys_module, only: gk_ion_boundary_policy
 use mod_fields, only: fields_base,calc_gyro_average_E_common
 use mod_fields_linear, only: jorek_fields_interp_linear
 class(fields_base),intent(in)::fields
 real*8,intent(in)::time,phi
 integer,intent(in)::n
 type(particle_kinetic_leapfrog),intent(in)::points(n)
 real*8,intent(inout)::E(3)
 select type(fields)
 type is(jorek_fields_interp_linear)
  call calc_gyro_average_E_common(fields,time,points,n,E,phi,original_weights=gk_ion_boundary_policy=='truncate')
 class default
  call fields%calc_gyro_average_E(time,points,n,E)
  if(gk_ion_boundary_policy=='truncate'.and.n>0)then
    if(count(points%i_elm>0)>0.and.count(points%i_elm>0)<n) &
      E=E*(real(count(points%i_elm>0),8)/real(n,8))
  endif
 end select
end subroutine gather_ion_common_basis

end module
