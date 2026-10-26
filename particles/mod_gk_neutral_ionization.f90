!> Deterministic whole-marker H0 -> H+ + e- infrastructure; not called by a main.
module mod_gk_neutral_ionization
use mpi
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use mod_particle_sim, only: particle_sim
use mod_particle_slot_pool
use mod_particle_types, only: particle_gc_vpar, particle_kinetic_leapfrog, codify_particle_type
use mod_find_rz_nearby, only: find_RZ_nearby
use mod_math_operators, only: cross_product
use mod_interp, only: interp_RZP
use constants, only: ATOMIC_MASS_UNIT, EL_CHG
implicit none
private
public :: gk_neutral_ionization_state, initialize_gk_neutral_ionization
public :: ionize_selected_gk_neutrals, report_gk_neutral_ionization
real*8, parameter :: energy_tolerance=128.d0*epsilon(1.d0)
type :: gk_neutral_ionization_state
  integer :: ion_group=0, electron_group=0, neutral_group=0
  logical :: initialized=.false.
  ! Selected includes rejected selections; failed counts originating rank-local attempts.
  integer :: counts(3)=0, global_counts(3)=0 ! selected, converted, failed
  ! Wn, Wi, We, En[J], Ei[J], Ee[J], Q[C], Wi-Wn, Wi-We, Ei+Ee-En.
  real*8 :: totals(10)=0.d0, global_totals(10)=0.d0
  real*8 :: reaction_time=0.d0
end type
contains

subroutine reject(sim,g,j,message)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: g,j
  character(len=*), intent(in) :: message
  write(*,*) 'ERROR: GK ionization rank=',sim%my_id,' group=',g,' particle=',j,': ',trim(message)
  if(allocated(sim%groups)) then
    if(g>=1.and.g<=size(sim%groups)) write(*,*) ' group ID=',sim%groups(g)%id
  endif
end subroutine

! All public transaction calls are collective, including empty local selections.
subroutine vote(local_error,error,ierr)
  integer, intent(in) :: local_error
  integer, intent(out) :: error
  integer, optional, intent(out) :: ierr
  integer :: mpierr
  call MPI_ALLREDUCE(local_error,error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpierr)
  if(present(ierr)) then
    ierr=error
  else if(error/=0) then
    call MPI_ABORT(MPI_COMM_WORLD,error,mpierr)
  endif
end subroutine

subroutine initialize_gk_neutral_ionization(sim,ion_group,electron_group,neutral_group,state,ierr)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: ion_group,electron_group,neutral_group
  type(gk_neutral_ionization_state), intent(out) :: state
  integer, optional, intent(out) :: ierr
  integer :: local_error,error
  state%ion_group=ion_group; state%electron_group=electron_group; state%neutral_group=neutral_group
  call validate_groups(sim,state,local_error)
  call vote(local_error,error,ierr)
  state%initialized=error==0
end subroutine

subroutine validate_groups(sim,state,error)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_ionization_state), intent(in) :: state
  integer, intent(out) :: error
  integer :: groups(3),k,g,j,charge
  real*8 :: mi,mn
  character(len=220) :: message
  error=1
  groups=[state%ion_group,state%electron_group,state%neutral_group]
  if(.not.allocated(sim%groups)) then
    call reject(sim,0,0,'groups unallocated; adiabatic-electron ionisation requires a separate implementation')
    return
  endif
  if(groups(1)==groups(2).or.groups(1)==groups(3).or.groups(2)==groups(3)) then
    call reject(sim,0,0,'explicit group indices must be distinct')
    return
  endif
  do k=1,3
    g=groups(k)
    if(g<1.or.g>size(sim%groups)) then
      message='invalid explicit group index'
      if(k==2) message='absent kinetic electrons: adiabatic-electron ionisation requires a separate implementation'
      call reject(sim,g,0,message)
      return
    endif
    if(.not.allocated(sim%groups(g)%particles)) then
      message='particle array unallocated'
      if(k==2) message='absent kinetic electrons: adiabatic-electron ionisation requires a separate implementation'
      call reject(sim,g,0,message)
      return
    endif
    if(.not.ieee_is_finite(sim%groups(g)%mass).or.sim%groups(g)%mass<=0.d0) then
      call reject(sim,g,0,'mass must be finite and positive')
      return
    endif
    if(k<=2) then
      charge=3-2*k
      select type(p=>sim%groups(g)%particles)
      type is(particle_gc_vpar)
        do j=1,size(p)
          if(p(j)%i_elm<=0.or.p(j)%q==charge) cycle
          write(message,*) 'active charge=',p(j)%q,'; expected ',charge
          call reject(sim,g,j,message)
          return
        enddo
      class default
        write(message,*) 'dynamic type=',codify_particle_type(sim%groups(g)%particles), &
                         '; expected particle_gc_vpar'
        if(k==2) message=trim(message)//'; adiabatic-electron ionisation requires a separate implementation'
        call reject(sim,g,0,message)
        return
      end select
    else
      select type(p=>sim%groups(g)%particles)
      type is(particle_kinetic_leapfrog)
        do j=1,size(p)
          if(p(j)%i_elm<=0.or.p(j)%q==0) cycle
          call reject(sim,g,j,'active neutral charge must be zero')
          return
        enddo
      class default
        write(message,*) 'dynamic type=',codify_particle_type(sim%groups(g)%particles), &
                         '; expected particle_kinetic_leapfrog'
        call reject(sim,g,0,message)
        return
      end select
      if(sim%groups(g)%coupling_scheme/='non') then
        call reject(sim,g,0,"configured neutral group requires coupling_scheme='non'")
        return
      endif
    endif
  enddo
  mi=sim%groups(groups(1))%mass; mn=sim%groups(groups(3))%mass
  if(abs(mi-mn)>32.d0*epsilon(mi)*max(mi,mn)) then
    call reject(sim,groups(3),0,'ion and neutral masses must agree within 32 epsilon relative tolerance')
    return
  endif
  error=0
end subroutine

!> Validate, stage, vote globally, then commit. No allocation/growth or migration of
!> particle slots; lost (negative-index) destinations are deliberately protected.
subroutine ionize_selected_gk_neutrals(sim,state,registry,neutral_indices,reaction_time,ierr)
  class(particle_sim), intent(inout) :: sim
  type(gk_neutral_ionization_state), intent(inout) :: state
  type(particle_slot_registry), intent(inout) :: registry
  integer, intent(in) :: neutral_indices(:)
  real*8, intent(in) :: reaction_time
  integer, optional, intent(out) :: ierr
  type(particle_gc_vpar), allocatable :: ions(:),electrons(:)
  integer, allocatable :: ion_slots(:),electron_slots(:)
  real*8 :: totals(10)
  integer :: local_error,error,k,n
  n=size(neutral_indices)
  state%counts(1)=state%counts(1)+n
  call stage_transaction(sim,state,registry,neutral_indices,reaction_time,ions,electrons, &
                         ion_slots,electron_slots,totals,local_error)
  if(local_error/=0) state%counts(3)=state%counts(3)+1
  call vote(local_error,error,ierr)
  if(error/=0) return
  ! Every rank has successfully staged both products. No failure paths below here.
  select type(p=>sim%groups(state%ion_group)%particles)
  type is(particle_gc_vpar)
    do k=1,n
      p(ion_slots(k))=ions(k)
    enddo
  end select
  select type(p=>sim%groups(state%electron_group)%particles)
  type is(particle_gc_vpar)
    do k=1,n
      p(electron_slots(k))=electrons(k)
    enddo
  end select
  do k=1,n
    sim%groups(state%neutral_group)%particles(neutral_indices(k))%i_elm=0
  enddo
  ! Common commit order: products, source deactivation, consumption, release, diagnostics.
  call commit_particle_slot_consumption(registry%pool(state%ion_group),ion_slots)
  call commit_particle_slot_consumption(registry%pool(state%electron_group),electron_slots)
  call commit_particle_slot_release(registry%pool(state%neutral_group),neutral_indices)
  state%counts(2)=state%counts(2)+n
  state%totals=state%totals+totals
  state%reaction_time=reaction_time
end subroutine

subroutine stage_transaction(sim,state,registry,indices,time,ions,electrons,ion_slots,electron_slots,totals,error)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_ionization_state), intent(in) :: state
  type(particle_slot_registry), intent(in) :: registry
  integer, intent(in) :: indices(:)
  real*8, intent(in) :: time
  type(particle_gc_vpar), allocatable, intent(out) :: ions(:),electrons(:)
  integer, allocatable, intent(out) :: ion_slots(:),electron_slots(:)
  real*8, intent(out) :: totals(10)
  integer, intent(out) :: error
  integer :: n,k,j,g,ni,ne,conversion_error
  real*8 :: en,ei
  character(len=240) :: message
  totals=0.d0; error=1
  if(.not.state%initialized) then
    call reject(sim,0,0,'ionization state is not initialized')
    return
  endif
  call validate_groups(sim,state,error)
  if(error/=0) return
  error=1; n=size(indices); g=state%neutral_group
  if(.not.ieee_is_finite(time).or.abs(time)>real(huge(0.e0),8)) then
    call reject(sim,g,0,'reaction_time must be finite and representable as real*4 t_birth')
    return
  endif
  call check_particle_slot_registry_groups(sim,registry, &
       [state%ion_group,state%electron_group,state%neutral_group],error)
  if(error/=0) then
    call reject(sim,g,0,'slot registry is uninitialized or does not match configured groups')
    return
  endif
  error=1
  ni=registry%pool(state%ion_group)%n_free
  ne=registry%pool(state%electron_group)%n_free
  if(n>ni.or.n>ne) then
    write(message,*) 'insufficient local capacity: selected=',n,' free ions=',ni,' free electrons=',ne
    call reject(sim,g,0,message)
    write(*,*) ' ion group ID=',sim%groups(state%ion_group)%id, &
               ' electron group ID=',sim%groups(state%electron_group)%id
    return
  endif
  allocate(ions(n),electrons(n),ion_slots(n),electron_slots(n))
  call peek_particle_slots(registry%pool(state%ion_group),n,ion_slots,error)
  if(error==0) call peek_particle_slots(registry%pool(state%electron_group),n,electron_slots,error)
  if(error==0) call check_particle_slots_for_release(registry%pool(g),indices,error)
  if(error/=0) then
    call reject(sim,g,0,'invalid slot consumption or source release transaction')
    return
  endif
  error=1
  if(any(sim%groups(state%ion_group)%particles(ion_slots)%i_elm/=0).or. &
     any(sim%groups(state%electron_group)%particles(electron_slots)%i_elm/=0)) then
    call reject(sim,g,0,'pool destination is not physically reusable')
    return
  endif
  if(n==0) then
    error=0
    return
  endif
  if(.not.allocated(sim%fields)) then
    call reject(sim,g,0,'fields unallocated')
    return
  endif
  if(.not.associated(sim%fields%node_list).or..not.associated(sim%fields%element_list)) then
    call reject(sim,g,0,'field mesh unassociated')
    return
  endif
  select type(p=>sim%groups(g)%particles)
  type is(particle_kinetic_leapfrog)
    do k=1,n
      j=indices(k)
      if(j<1.or.j>size(p)) then
        call reject(sim,g,j,'selected index outside neutral array')
        return
      endif
      if(any(indices(:k-1)==j)) then
        call reject(sim,g,j,'duplicate selected neutral index')
        return
      endif
      if(p(j)%i_elm<=0.or.p(j)%q/=0) then
        call reject(sim,g,j,'selected source must be active with q=0; no double conversion')
        return
      endif
      if(.not.ieee_is_finite(p(j)%weight).or.p(j)%weight<=0.d0.or. &
         .not.all(ieee_is_finite(p(j)%v))) then
        call reject(sim,g,j,'source weight must be finite positive and velocity finite')
        return
      endif
      if(.not.valid_position(sim,p(j)%x,p(j)%st,p(j)%i_elm)) then
        call reject(sim,g,j,'invalid source position, st or element')
        return
      endif
    enddo
    do k=1,n
      j=indices(k)
      call stage_products(sim,state,p(j),time,ions(k),electrons(k),en,ei,conversion_error)
      if(conversion_error/=0) then
        call reject(sim,g,j,'product staging failed; global transaction cancelled')
        return
      endif
      if(sim%groups(state%ion_group)%particles(ion_slots(k))%i_life==huge(0).or. &
         sim%groups(state%electron_group)%particles(electron_slots(k))%i_life==huge(0)) then
        call reject(sim,g,j,'destination lifetime would overflow')
        return
      endif
      ions(k)%i_life=sim%groups(state%ion_group)%particles(ion_slots(k))%i_life+1
      electrons(k)%i_life=sim%groups(state%electron_group)%particles(electron_slots(k))%i_life+1
      totals(1:3)=totals(1:3)+p(j)%weight
      totals(4)=totals(4)+en; totals(5)=totals(5)+ei
    enddo
  end select
  totals(7)=EL_CHG*(totals(2)-totals(3))
  totals(8)=totals(2)-totals(1); totals(9)=totals(2)-totals(3)
  totals(10)=totals(5)+totals(6)-totals(4)
  if(.not.all(ieee_is_finite(totals)).or. &
     abs(totals(10))>energy_tolerance*max(totals(4),totals(5),tiny(1.d0))) then
    call reject(sim,g,0,'nonfinite accumulated diagnostics or kinetic-energy mismatch')
    return
  endif
  error=0
end subroutine

logical function valid_position(sim,x,st,element) result(valid)
  class(particle_sim), intent(in) :: sim
  real*8, intent(in) :: x(3),st(2)
  integer, intent(in) :: element
  real*8 :: r,z
  valid=.false.
  if(element<=0.or.element>sim%fields%element_list%n_elements) return
  if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(st))) return
  if(x(1)<=0.d0.or.any(st<0.d0).or.any(st>1.d0)) return
  call interp_RZP(sim%fields%node_list,sim%fields%element_list,element,st(1),st(2),x(3),r,z)
  ! Consistency tolerance follows the mesh locator, not the energy tolerance.
  valid=all(ieee_is_finite([r,z])).and.maxval(abs(x(:2)-[r,z]))<=1.d-8*max(1.d0,abs(r),abs(z))
end function

subroutine magnetic_field(sim,time,x,st,element,B,Bnorm,error)
  class(particle_sim), intent(in) :: sim
  real*8, intent(in) :: time,x(3),st(2)
  integer, intent(in) :: element
  real*8, intent(out) :: B(3),Bnorm
  integer, intent(out) :: error
  real*8 :: E(3),psi,U
  call sim%fields%calc_EBpsiU(time,element,st,x(3),E,B,psi,U)
  Bnorm=norm2(B); error=1
  if(all(ieee_is_finite(B)).and.ieee_is_finite(Bnorm).and.Bnorm>0.d0) error=0
  if(error/=0) call reject(sim,0,0,'magnetic field must be finite with strictly positive norm')
end subroutine

subroutine stage_products(sim,state,source,time,ion,electron,en,ei,error)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_ionization_state), intent(in) :: state
  type(particle_kinetic_leapfrog), intent(in) :: source
  real*8, intent(in) :: time
  type(particle_gc_vpar), intent(out) :: ion,electron
  real*8, intent(out) :: en,ei
  integer, intent(out) :: error
  integer :: ifail
  real*8 :: B(3),Bnorm,b0(3),v2,perp2,displacement(3),radial,mi,mn
  error=1
  call magnetic_field(sim,time,source%x,source%st,source%i_elm,B,Bnorm,ifail)
  if(ifail/=0) return
  mi=sim%groups(state%ion_group)%mass*ATOMIC_MASS_UNIT
  mn=sim%groups(state%neutral_group)%mass*ATOMIC_MASS_UNIT
  b0=B/Bnorm; v2=dot_product(source%v,source%v)
  ion%vpar=dot_product(source%v,b0); perp2=v2-ion%vpar**2
  if(.not.ieee_is_finite(v2).or..not.ieee_is_finite(perp2)) return
  if(perp2 < -32.d0*epsilon(v2)*max(v2,tiny(v2))) return
  perp2=max(0.d0,perp2)
  ! Preserve the shared JOREK kinetic-to-GC sign in its (R,Z,phi) component order.
  displacement=mi*cross_product(source%v,b0)/(EL_CHG*Bnorm)
  radial=source%x(1)+displacement(1)
  ion%x=[norm2([radial,displacement(3)]),source%x(2)+displacement(2), &
         source%x(3)+atan2(displacement(3),radial)]
  if(.not.all(ieee_is_finite(ion%x)).or.ion%x(1)<=0.d0) return
  call find_RZ_nearby(sim%fields%node_list,sim%fields%element_list,source%x(1),source%x(2), &
                     source%st(1),source%st(2),source%i_elm,ion%x(1),ion%x(2), &
                     ion%st(1),ion%st(2),ion%i_elm,ifail,phi=ion%x(3))
  ! Current finite-Larmor-radius boundary limitation: both source and GC must be inside.
  if(ifail/=0.or.ion%i_elm<=0) then
    call reject(sim,state%ion_group,0,'guiding centre outside mesh or locator failure (finite-Larmor-radius limit)')
    return
  endif
  if(.not.valid_position(sim,ion%x,ion%st,ion%i_elm)) return
  call magnetic_field(sim,time,ion%x,ion%st,ion%i_elm,B,ion%B_norm,ifail)
  if(ifail/=0) return
  ion%mu=perp2/(2.d0*ion%B_norm)
  ion%q=1; ion%weight=source%weight; ion%t_birth=real(time,4); ion%i_life=0
  ! Stage B2 must supply a physical secondary-electron energy model and remove
  ! the ionisation potential from the electron population. Neither is included here.
  electron%x=source%x; electron%st=source%st; electron%i_elm=source%i_elm
  electron%q=-1; electron%weight=source%weight; electron%vpar=0.d0; electron%mu=0.d0
  electron%B_norm=Bnorm; electron%t_birth=real(time,4); electron%i_life=0
  en=0.5d0*mn*v2*source%weight
  ei=mi*(0.5d0*ion%vpar**2+ion%mu*ion%B_norm)*source%weight
  if(.not.all(ieee_is_finite([ion%vpar,ion%mu,en,ei]))) return
  if(abs(ei-en)>energy_tolerance*max(en,ei,tiny(en))) then
    call reject(sim,state%ion_group,0,'ion/source kinetic-energy mismatch (128 epsilon relative tolerance)')
    return
  endif
  error=0
end subroutine

subroutine report_gk_neutral_ionization(sim,state)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_ionization_state), intent(inout) :: state
  integer :: mpierr
  ! Always reduce local cumulative buffers; never feed global results back into them.
  call MPI_ALLREDUCE(state%counts,state%global_counts,3,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
  call MPI_ALLREDUCE(state%totals,state%global_totals,10,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)
  if(.not.all(ieee_is_finite(state%global_totals)).or. &
     any(state%global_totals(7:9)/=0.d0).or. &
     abs(state%global_totals(10))>energy_tolerance*max(state%global_totals(4), &
                                                     state%global_totals(5),tiny(1.d0))) then
    call reject(sim,0,0,'global conservation diagnostic failed')
    call MPI_ABORT(MPI_COMM_WORLD,1,mpierr)
  endif
  if(sim%my_id==0) then
    write(*,*) 'GK ionization reaction_time=',state%reaction_time,' selected/converted/failed=',state%global_counts
    write(*,'(A,10ES24.15)') ' Wn Wi We En Ei Ee Q nuclear_res charge_res energy_res: ',state%global_totals
  endif
end subroutine
end module
