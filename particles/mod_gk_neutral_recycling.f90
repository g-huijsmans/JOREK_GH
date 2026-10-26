!> Prompt hydrogenic recycling on the existing finite-element boundary.
module mod_gk_neutral_recycling
use mpi
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use mod_particle_sim, only: particle_sim
use mod_particle_slot_pool
use mod_particle_types, only: particle_gc_vpar, particle_kinetic_leapfrog, codify_particle_type
use mod_pcg32_rng, only: pcg32_rng
use mod_random_seed, only: random_seed
use mod_sampling, only: sample_cosine
use mod_boundary, only: wall_normal_vector
use mod_interp, only: interp_RZ
use constants, only: ATOMIC_MASS_UNIT, EL_CHG
implicit none
private
public :: gk_neutral_recycling_state, initialize_gk_neutral_recycling
public :: recycle_lost_gk_ions, report_gk_neutral_recycling

type :: gk_neutral_recycling_state
  integer :: ion_group=0, neutral_group=0
  logical :: enabled=.false.
  type(pcg32_rng) :: rng
  ! Counts: incident, created, unresolved, multiple boundary sides, failed attempts.
  integer :: counts(5)=0, global_counts(5)=0
  ! Totals: incident weight/energy[J], created weight/energy[J], wall charge[C].
  real*8 :: totals(5)=0.d0, global_totals(5)=0.d0
end type

contains

subroutine initialize_gk_neutral_recycling(sim, ion_group, neutral_group, state)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: ion_group, neutral_group
  type(gk_neutral_recycling_state), intent(out) :: state
  integer :: error, seed, mpierr
  state%ion_group=ion_group
  state%neutral_group=neutral_group
  if (neutral_group > size(sim%groups)) return ! charged-only configuration
  call validate_groups(sim,state,error)
  if (error /= 0) call MPI_ABORT(MPI_COMM_WORLD,1,mpierr)
  seed=random_seed()
  if (seed==0) then
    call reject(sim,ion_group,0,'PCG seed is zero; expected nonzero manual/random seed')
    call MPI_ABORT(MPI_COMM_WORLD,1,mpierr)
  endif
  ! Separate rank-local stream. Each product consumes exactly two uniform doubles.
  call state%rng%initialize(2,seed,sim%n_mpi,sim%my_id+1)
  state%enabled=.true.
end subroutine

subroutine reject(sim,g,j,message)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: g,j
  character(len=*), intent(in) :: message
  write(*,*) 'ERROR: GK recycling rank=',sim%my_id,' group=',g,' particle=',j,': ',trim(message)
  if (g>=1.and.g<=size(sim%groups)) write(*,*) ' group ID=',sim%groups(g)%id
end subroutine

subroutine validate_groups(sim,state,error)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_recycling_state), intent(in) :: state
  integer, intent(out) :: error
  integer :: g,j,k
  real*8 :: mi,mn
  character(len=180) :: message
  error=1
  do k=1,2
    g=state%ion_group
    if(k==2) g=state%neutral_group
    if(g<1.or.g>size(sim%groups)) then
      call reject(sim,g,0,'invalid group index; expected allocated group')
      return
    endif
    if(.not.allocated(sim%groups(g)%particles)) then
      call reject(sim,g,0,'particle array unallocated')
      return
    endif
  enddo
  g=state%ion_group
  select type(p=>sim%groups(g)%particles)
  type is(particle_gc_vpar)
    do j=1,size(p)
      if(p(j)%i_elm==0) cycle
      if(p(j)%q==1) cycle
      write(message,*) 'detected ion charge=',p(j)%q,'; expected q=+1'
      call reject(sim,g,j,message)
      return
    enddo
  class default
    write(message,*) 'detected type code=',codify_particle_type(sim%groups(g)%particles), &
                     '; expected particle_gc_vpar (3)'
    call reject(sim,g,0,message)
    return
  end select
  g=state%neutral_group
  select type(p=>sim%groups(g)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,size(p)
      if(p(j)%i_elm<=0) cycle
      if(p(j)%q==0) cycle
      write(message,*) 'detected neutral charge=',p(j)%q,'; expected q=0'
      call reject(sim,g,j,message)
      return
    enddo
  class default
    write(message,*) 'detected type code=',codify_particle_type(sim%groups(g)%particles), &
                     '; expected particle_kinetic_leapfrog (6)'
    call reject(sim,g,0,message)
    return
  end select
  mi=sim%groups(state%ion_group)%mass
  mn=sim%groups(state%neutral_group)%mass
  if(.not.ieee_is_finite(mi).or..not.ieee_is_finite(mn).or.mi<=0.d0.or.mn<=0.d0) then
    write(message,*) 'detected ion/neutral masses=',mi,mn,'; expected finite positive masses'
    call reject(sim,state%ion_group,0,message)
    call reject(sim,state%neutral_group,0,message)
    return
  endif
  if(abs(mi-mn)>32.d0*epsilon(mi)*max(mi,mn)) then
    write(message,*) 'detected ion/neutral masses=',mi,mn,'; expected equal hydrogenic masses'
    call reject(sim,state%ion_group,0,message)
    call reject(sim,state%neutral_group,0,message)
    return
  endif
  error=0
end subroutine

!> Same adjacent-node definition as boundary_from_grid; no second wall mesh.
subroutine boundary_sides(sim,element,sides)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: element
  logical, intent(out) :: sides(4)
  integer :: side,n1,n2
  do side=1,4
    n1=sim%fields%element_list%element(element)%vertex(side)
    n2=sim%fields%element_list%element(element)%vertex(mod(side,4)+1)
    sides(side)=sim%fields%node_list%node(n1)%boundary>0.and.sim%fields%node_list%node(n2)%boundary>0
  enddo
end subroutine

!> Collective transaction, including ranks without local losses. Pool operations are
!> rank-local; the error vote precedes all particle, pool, RNG and diagnostic commits.
subroutine recycle_lost_gk_ions(sim,state,registry,birth_time,dt,ierr)
  class(particle_sim), intent(inout) :: sim
  type(gk_neutral_recycling_state), intent(inout) :: state
  type(particle_slot_registry), intent(inout) :: registry
  real*8, intent(in) :: birth_time,dt
  integer, optional, intent(out) :: ierr
  type(gk_neutral_recycling_state) :: staged
  type(particle_kinetic_leapfrog), allocatable :: products(:)
  integer, allocatable :: sources(:),slots(:)
  integer :: local_error,error,mpierr,k
  if(present(ierr)) ierr=0
  if(.not.state%enabled) return
  call recycle_transaction(sim,state,registry,birth_time,dt,sources,slots,products,staged,local_error)
  call MPI_ALLREDUCE(local_error,error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpierr)
  if(error/=0) then
    if(local_error/=0) state%counts(5)=state%counts(5)+1
    if(present(ierr)) then
      ierr=error
    else
      call MPI_ABORT(MPI_COMM_WORLD,1,mpierr)
    endif
    return
  endif
  ! Commit order: complete products, source deactivation, consumption, release, state.
  select type(neutrals=>sim%groups(state%neutral_group)%particles)
  type is(particle_kinetic_leapfrog)
    do k=1,size(sources)
      neutrals(slots(k))=products(k)
    enddo
  end select
  sim%groups(state%ion_group)%particles(sources)%i_elm=0
  call commit_particle_slot_consumption(registry%pool(state%neutral_group),slots)
  call commit_particle_slot_release(registry%pool(state%ion_group),sources)
  state=staged
end subroutine

subroutine recycle_transaction(sim,state,registry,birth_time,dt,sources,slots,products,staged,error)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_recycling_state), intent(in) :: state
  type(particle_slot_registry), intent(in) :: registry
  type(gk_neutral_recycling_state), intent(out) :: staged
  real*8, intent(in) :: birth_time,dt
  integer, intent(out) :: error
  integer :: g,h,j,k,n,available,e,side,nside,slot,nmultiple
  integer, allocatable, intent(out) :: sources(:),slots(:)
  real*8, allocatable :: normals(:,:),speed2(:),energies(:)
  type(particle_kinetic_leapfrog), allocatable, intent(out) :: products(:)
  type(pcg32_rng) :: rng
  real*8 :: distances(4),r,rs,rt,z,zs,zt,jac,normal(3),u(2),direction(3),scale,emitted
  logical :: sides(4),allowed(4)
  character(len=180) :: message
  staged=state
  call validate_groups(sim,state,error)
  if(error/=0) return
  g=state%ion_group; h=state%neutral_group
  call check_particle_slot_registry_groups(sim,registry,[g,h],error)
  if(error/=0) then
    call reject(sim,h,0,'slot registry is uninitialized or does not match configured groups')
    return
  endif
  error=1
  if(.not.ieee_is_finite(birth_time).or.abs(birth_time)>real(huge(0.e0),8)) then
    call reject(sim,h,0,'invalid birth time; expected finite representable physical end time')
    return
  endif
  select type(ions=>sim%groups(g)%particles)
  type is(particle_gc_vpar)
  select type(neutrals=>sim%groups(h)%particles)
  type is(particle_kinetic_leapfrog)
    n=count(ions%i_elm<0)
    available=registry%pool(h)%n_free
    if(n>available) then
      write(message,*) 'insufficient local capacity: required=',n,' available i_elm==0 slots=',available
      call reject(sim,h,0,message)
      return
    endif
    allocate(sources(n),slots(n),normals(3,n),speed2(n),energies(n),products(n))
    call peek_particle_slots(registry%pool(h),n,slots,error)
    if(error/=0) then
      call reject(sim,h,0,'invalid neutral pool consumption')
      return
    endif
    error=1
    if(any(neutrals(slots)%i_elm/=0)) then
      call reject(sim,h,0,'pool destination is not physically reusable')
      return
    endif
    if(n==0) then
      error=0
      return
    endif
    nmultiple=0
    ! Pass 1: validate every source and prepare geometry/energy without mutation.
    k=0
    do j=1,size(ions)
      if(ions(j)%i_elm>=0) cycle
      k=k+1; sources(k)=j
      e=-ions(j)%i_elm
      if(e<1.or.e>sim%fields%element_list%n_elements) then
        write(message,*) 'invalid i_last=',e,'; expected 1..',sim%fields%element_list%n_elements
        call reject(sim,g,j,message)
        return
      endif
      call boundary_sides(sim,e,sides)
      nside=count(sides)
      if(nside==0) then
        staged%counts(3)=staged%counts(3)+1
        call reject(sim,g,j,'unresolved non-boundary last element; reduce the particle timestep')
        write(*,*) ' rank=',sim%my_id,' particle=',j,' i_last=',e,' x=',ions(j)%x, &
                   ' st=',ions(j)%st,' vpar=',ions(j)%vpar,' mu=',ions(j)%mu,' particle timestep=',dt
        write(*,*) 'Trajectory backtracking for skipped boundary elements is deferred.'
        return
      endif
      if(nside>1) nmultiple=nmultiple+1
      if(.not.all(ieee_is_finite(ions(j)%st)).or..not.all(ieee_is_finite(ions(j)%x))) then
        call reject(sim,g,j,'nonfinite retained wall coordinates')
        return
      endif
      distances=[ions(j)%st(2),1.d0-ions(j)%st(1),1.d0-ions(j)%st(2),ions(j)%st(1)]
      if(minval(distances)<-1.d-10.or.count(sides.and.abs(distances)<=1.d-10)/=1) then
        call reject(sim,g,j,'retained st must identify exactly one boundary side; ambiguous corner/interior st')
        return
      endif
      side=minloc(abs(distances),mask=sides,dim=1)
      ! wall_normal_vector selects among zero-neighbour sides; verify agreement.
      allowed=sim%fields%element_list%element(e)%neighbours==0
      if(.not.any(allowed)) allowed=.true.
      if(minloc(distances,mask=allowed,dim=1)/=side) then
        call reject(sim,g,j,'boundary-node criterion and wall normal side disagree')
        return
      endif
      call interp_RZ(sim%fields%node_list,sim%fields%element_list,e,ions(j)%st(1),ions(j)%st(2), &
                     r,rs,rt,z,zs,zt)
      ! Reject degenerate geometry before wall_normal_vector's fallback can hide it.
      jac=rs*zt-rt*zs
      select case(side)
      case(1,3)
        scale=norm2([rs,zs])
      case(2,4)
        scale=norm2([rt,zt])
      end select
      if(.not.all(ieee_is_finite([r,z,jac,scale])).or.jac<=0.d0.or.scale<=0.d0) then
        call reject(sim,g,j,'invalid wall geometry: expected finite positive Jacobian and nonzero tangent')
        return
      endif
      normal=wall_normal_vector(sim%fields%node_list,sim%fields%element_list,e,ions(j)%st(1),ions(j)%st(2))
      if(.not.all(ieee_is_finite(normal)).or.norm2(normal)<=0.d0) then
        call reject(sim,g,j,'invalid inward wall normal')
        return
      endif
      normals(:,k)=normal/norm2(normal)
      if(.not.all(ieee_is_finite([ions(j)%B_norm,ions(j)%vpar,ions(j)%mu,ions(j)%weight])).or. &
         ions(j)%B_norm<=0.d0.or.ions(j)%weight<0.d0) then
        call reject(sim,g,j,'expected finite B_norm>0, vpar, mu and nonnegative weight')
        return
      endif
      ! The RK loss return preserves accepted vpar, mu and B_norm, as used by GK loss diagnostics.
      speed2(k)=ions(j)%vpar**2+2.d0*ions(j)%mu*ions(j)%B_norm
      scale=max(ions(j)%vpar**2,abs(2.d0*ions(j)%mu*ions(j)%B_norm))
      if(.not.ieee_is_finite(speed2(k)).or.speed2(k)<-32.d0*epsilon(scale)*scale) then
        call reject(sim,g,j,'nonfinite or materially negative incident speed_squared')
        return
      endif
      speed2(k)=max(0.d0,speed2(k))
      energies(k)=0.5d0*sim%groups(g)%mass*ATOMIC_MASS_UNIT*ions(j)%weight*speed2(k)
      if(.not.ieee_is_finite(energies(k))) then
        call reject(sim,g,j,'nonfinite incident kinetic energy')
        return
      endif
      products(k)%x=[r,z,ions(j)%x(3)]
      products(k)%st=ions(j)%st
      products(k)%i_elm=e
    enddo
    call check_particle_slots_for_release(registry%pool(g),sources,error)
    if(error/=0) then
      call reject(sim,g,0,'invalid source-ion release transaction')
      return
    endif
    error=1
    ! Pass 2: sample on a copy of the stream, stage complete products, then commit.
    rng=state%rng
    do k=1,n
      j=sources(k); slot=slots(k)
      call rng%next(u) ! exactly two uniform doubles per created marker, no rejection sampling
      direction=sample_cosine(u,normals(:,k)) ! cylindrical (R,Z,phi), inward hemisphere
      if(.not.all(ieee_is_finite(direction)).or.dot_product(direction,normals(:,k))<=0.d0) then
        call reject(sim,g,j,'cosine sample must have a finite, strictly inward unit direction')
        return
      endif
      products(k)%v=sqrt(speed2(k))*direction
      products(k)%q=0
      products(k)%weight=ions(j)%weight
      if(neutrals(slot)%i_life==huge(0)) then
        call reject(sim,h,slot,'destination lifetime would overflow')
        return
      endif
      products(k)%i_life=neutrals(slot)%i_life+1
      products(k)%t_birth=real(birth_time,4)
      emitted=0.5d0*sim%groups(h)%mass*ATOMIC_MASS_UNIT*products(k)%weight*sum(products(k)%v**2)
      if(.not.ieee_is_finite(emitted).or.abs(emitted-energies(k))> &
         128.d0*epsilon(emitted)*max(emitted,energies(k),tiny(emitted))) then
        call reject(sim,g,j,'emitted and incident kinetic energies disagree')
        return
      endif
    enddo
    do k=1,n
      j=sources(k); slot=slots(k)
      staged%counts(1:2)=staged%counts(1:2)+1
      staged%totals=staged%totals+[ions(j)%weight,energies(k),products(k)%weight, &
           0.5d0*sim%groups(h)%mass*ATOMIC_MASS_UNIT*products(k)%weight*sum(products(k)%v**2), &
           EL_CHG*ions(j)%weight]
    enddo
    staged%counts(4)=staged%counts(4)+nmultiple
    staged%rng=rng
  end select
  end select
  error=0
end subroutine

!> Accumulated since previous report. Never reduce a previously reduced buffer.
subroutine report_gk_neutral_recycling(sim,state)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_recycling_state), intent(inout) :: state
  integer :: ierr
  if(.not.state%enabled) return
  call MPI_REDUCE(state%counts,state%global_counts,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(state%totals,state%global_totals,5,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if(sim%my_id==0) write(*,'(A,I0,A,I0,5(A,ES24.16),3(A,I0))') &
       'GK recycling ions=',state%global_counts(1),' neutrals=',state%global_counts(2), &
       ' incident_weight=',state%global_totals(1),' incident_energy_J=',state%global_totals(2), &
       ' created_weight=',state%global_totals(3),' created_energy_J=',state%global_totals(4), &
       ' wall_charge_C=',state%global_totals(5),' unresolved=',state%global_counts(3), &
       ' multiple_sides=',state%global_counts(4),' failed=',state%global_counts(5)
  state%counts=0
  state%totals=0.d0
end subroutine
end module mod_gk_neutral_recycling
