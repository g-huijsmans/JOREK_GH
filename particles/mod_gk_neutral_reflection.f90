!> Elastic, unit-probability diffuse reflection of current-step neutral wall hits.
!> No slot ownership, wall temperature, remaining-substep motion, or volume physics.
module mod_gk_neutral_reflection
use mpi
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use mod_particle_sim, only: particle_sim
use mod_particle_types, only: particle_kinetic_leapfrog
use mod_gk_neutrals, only: gk_neutral_state, get_gk_neutral_events, consume_gk_neutral_events
use mod_pcg32_rng, only: pcg32_rng
use mod_random_seed, only: random_seed
use mod_sampling, only: sample_cosine
use mod_boundary, only: wall_normal_vector
use mod_interp, only: interp_RZ
use constants, only: ATOMIC_MASS_UNIT
implicit none
private
public :: gk_neutral_reflection_state, initialize_gk_neutral_reflection
public :: reflect_newly_lost_gk_neutrals, report_gk_neutral_reflection
real*8, parameter :: st_tolerance=1.d-10, roundoff_tolerance=128.d0*epsilon(1.d0)
type :: gk_neutral_reflection_state
  integer :: neutral_group=0
  logical :: initialized=.false., enabled=.false.
  type(pcg32_rng) :: rng
  ! Cumulative: incident, reflected, grazing, unresolved, multiple-side, corner, failed.
  ! Failed attempts are diagnosed, but their staged counters cannot be committed:
  ! optional-ierr rollback preserves ALL diagnostics, including failure counters.
  integer :: counts(7)=0, global_counts(7)=0
  ! Incident/reflected weight and incident/reflected translational energy [J].
  real*8 :: totals(4)=0.d0, global_totals(4)=0.d0
end type
contains

subroutine reject(sim,g,j,message)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: g,j
  character(len=*), intent(in) :: message
  write(*,*) 'ERROR: GK reflection rank=',sim%my_id,' group=',g,' marker=',j,': ',trim(message)
end subroutine

subroutine vote(local_error,error,ierr)
  integer, intent(in) :: local_error
  integer, intent(out) :: error
  integer, optional, intent(out) :: ierr
  integer :: mpierr
  call MPI_ALLREDUCE(local_error,error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpierr)
  if(present(ierr)) then
    ierr=error
  else if(error/=0) then
    call MPI_ABORT(MPI_COMM_WORLD,1,mpierr)
  endif
end subroutine

subroutine validate_group(sim,g,error)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: g
  integer, intent(out) :: error
  error=1
  if(.not.allocated(sim%groups)) then
    call reject(sim,g,0,'groups unallocated')
    return
  endif
  if(g<1.or.g>size(sim%groups)) then
    call reject(sim,g,0,'invalid configured neutral index')
    return
  endif
  if(.not.allocated(sim%groups(g)%particles)) then
    call reject(sim,g,0,'neutral particle array unallocated')
    return
  endif
  select type(p=>sim%groups(g)%particles)
  type is(particle_kinetic_leapfrog)
  class default
    call reject(sim,g,0,'expected particle_kinetic_leapfrog')
    return
  end select
  if(.not.ieee_is_finite(sim%groups(g)%mass).or.sim%groups(g)%mass<=0.d0) then
    call reject(sim,g,0,'expected finite positive neutral mass')
    return
  endif
  if(sim%groups(g)%coupling_scheme/='non') then
    call reject(sim,g,0,"configured neutral group requires coupling_scheme='non'")
    return
  endif
  error=0
end subroutine

subroutine initialize_gk_neutral_reflection(sim,neutral_group,enabled,state,ierr)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: neutral_group
  logical, intent(in) :: enabled
  type(gk_neutral_reflection_state), intent(out) :: state
  integer, optional, intent(out) :: ierr
  integer :: local_error,error,seed,j
  state%neutral_group=neutral_group
  if(present(ierr)) ierr=0
  if(.not.enabled) then
    state%initialized=.true.
    return
  endif
  if(allocated(sim%groups)) then
    if(neutral_group>size(sim%groups)) then
      state%initialized=.true.
      return
    endif
  endif
  call validate_group(sim,neutral_group,local_error)
  if(local_error==0) then
    select type(p=>sim%groups(neutral_group)%particles)
    type is(particle_kinetic_leapfrog)
    do j=1,size(p)
      if(sim%groups(neutral_group)%particles(j)%i_elm<=0) cycle
      if(p(j)%q==0) cycle
      call reject(sim,neutral_group,j,'active neutral requires q=0')
      local_error=1
      exit
    enddo
    end select
  endif
  seed=0
  if(local_error==0) then
    seed=random_seed()
    if(seed==0) then
      call reject(sim,neutral_group,0,'PCG requires a nonzero seed')
      local_error=1
    endif
  endif
  call vote(local_error,error,ierr)
  if(error/=0) return
  ! PCG accepts signed stream indices. Negative indices are disjoint from
  ! recycling's positive rank+1 streams; n_streams is unused by this PCG wrapper.
  call state%rng%initialize(2,seed,sim%n_mpi,-sim%my_id-1)
  state%initialized=.true.; state%enabled=.true.
end subroutine

!> Collective, even on ranks with no local events. No collective during commit.
subroutine reflect_newly_lost_gk_neutrals(sim,neutrals,state,impact_time,dt,ierr)
  class(particle_sim), intent(inout) :: sim
  type(gk_neutral_state), intent(inout) :: neutrals
  type(gk_neutral_reflection_state), intent(inout) :: state
  real*8, intent(in) :: impact_time,dt
  integer, optional, intent(out) :: ierr
  type(gk_neutral_reflection_state) :: staged
  type(particle_kinetic_leapfrog), allocatable :: products(:)
  integer, allocatable :: indices(:)
  integer :: group,local_error,error,k
  if(present(ierr)) ierr=0
  if(.not.state%enabled) return
  call get_gk_neutral_events(neutrals,indices,group,local_error)
  if(local_error/=0) then
    call reject(sim,state%neutral_group,0,'invalid passive event buffer/count')
  else if(group/=state%neutral_group.or..not.state%initialized) then
    call reject(sim,state%neutral_group,0,'uninitialized state or event group mismatch')
    local_error=1
  else
    call stage_reflection(sim,state,indices,impact_time,dt,products,staged,local_error)
  endif
  call vote(local_error,error,ierr)
  if(error/=0) return
  select type(p=>sim%groups(state%neutral_group)%particles)
  type is(particle_kinetic_leapfrog)
    do k=1,size(indices)
      p(indices(k))=products(k)
    enddo
  end select
  state=staged
  call consume_gk_neutral_events(neutrals)
end subroutine

! Adjacent boundary-node criterion shared in convention with boundary_from_grid.
function boundary_sides(sim,e) result(sides)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: e
  logical :: sides(4)
  integer :: side,n1,n2
  do side=1,4
    n1=sim%fields%element_list%element(e)%vertex(side)
    n2=sim%fields%element_list%element(e)%vertex(mod(side,4)+1)
    sides(side)=sim%fields%node_list%node(n1)%boundary>0.and.sim%fields%node_list%node(n2)%boundary>0
  enddo
end function

subroutine stage_reflection(sim,state,indices,impact_time,dt,products,staged,error)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_reflection_state), intent(in) :: state
  integer, intent(in) :: indices(:)
  real*8, intent(in) :: impact_time,dt
  type(particle_kinetic_leapfrog), allocatable, intent(out) :: products(:)
  type(gk_neutral_reflection_state), intent(out) :: staged
  integer, intent(out) :: error
  integer :: g,n,j,k,e,side,nsides
  real*8, allocatable :: normals(:,:),speed2(:),energy(:)
  real*8 :: distances(4),r,rs,rt,z,zs,zt,jac,tangent,normal(3),speed,vn,u(2),direction(3),emitted
  logical :: sides(4),allowed(4),ordered
  character(len=160) :: message
  staged=state; g=state%neutral_group; n=size(indices); j=0
  call validate_group(sim,g,error)
  if(error/=0) return
  error=1
  message='nonfinite time or nonpositive/nonfinite substep'
  if(.not.ieee_is_finite(impact_time).or..not.ieee_is_finite(dt).or.dt<=0.d0) goto 900
  allocate(products(n),normals(3,n),speed2(n),energy(n))
  ! The passive producer emits increasing indices. Only a malformed/unordered
  ! supplied list needs the fallback duplicate check; neither path scans slots.
  ordered=.true.
  if(n>1) ordered=all(indices(2:)>indices(:n-1))
  select type(p=>sim%groups(g)%particles)
  type is(particle_kinetic_leapfrog)
    ! Phase 1: validate ALL events and geometry before sampling any random values.
    do k=1,n
      j=indices(k)
      message='invalid event index'
      if(j<1.or.j>size(p)) goto 900
      message='duplicate event index'
      if(.not.ordered) then
        if(any(indices(:k-1)==j)) goto 900
      endif
      message='event source must have negative i_elm'
      if(p(j)%i_elm>=0) goto 900
      message='invalid last element index'
      ! Check before negating to avoid abs(minimum integer) overflow.
      if(p(j)%i_elm < -sim%fields%element_list%n_elements) goto 900
      e=-p(j)%i_elm
      message='expected q=0 and finite strictly positive marker weight'
      if(p(j)%q/=0.or..not.ieee_is_finite(p(j)%weight).or.p(j)%weight<=0.d0) goto 900
      message='nonfinite retained position, st, or incident velocity'
      if(.not.all(ieee_is_finite(p(j)%x)).or..not.all(ieee_is_finite(p(j)%st)).or. &
         .not.all(ieee_is_finite(p(j)%v))) goto 900
      sides=boundary_sides(sim,e)
      nsides=count(sides)
      message='unresolved non-boundary last element; reduce the particle timestep (no trajectory backtracking)'
      if(nsides==0) goto 900
      distances=[p(j)%st(2),1.d0-p(j)%st(1),1.d0-p(j)%st(2),p(j)%st(1)]
      message='invalid retained st or ambiguous corner: require exactly one boundary side within 1e-10'
      if(minval(distances)<-st_tolerance.or.count(sides.and.abs(distances)<=st_tolerance)/=1) goto 900
      side=minloc(abs(distances),mask=sides,dim=1)
      allowed=sim%fields%element_list%element(e)%neighbours==0
      if(.not.any(allowed)) allowed=.true.
      message='boundary-node side and wall_normal_vector side disagree'
      if(minloc(distances,mask=allowed,dim=1)/=side) goto 900
      call interp_RZ(sim%fields%node_list,sim%fields%element_list,e,p(j)%st(1),p(j)%st(2),r,rs,rt,z,zs,zt)
      jac=rs*zt-rt*zs
      tangent=norm2([rs,zs])
      if(side==2.or.side==4) tangent=norm2([rt,zt])
      message='invalid wall geometry: require positive Jacobian, radius and nonzero tangent'
      if(.not.all(ieee_is_finite([r,z,jac,tangent])).or.r<=0.d0.or.jac<=0.d0.or.tangent<=0.d0) goto 900
      normal=wall_normal_vector(sim%fields%node_list,sim%fields%element_list,e,p(j)%st(1),p(j)%st(2))
      message='nonfinite or zero inward normal'
      if(.not.all(ieee_is_finite(normal)).or.norm2(normal)<=0.d0) goto 900
      normals(:,k)=normal/norm2(normal)
      speed2(k)=dot_product(p(j)%v,p(j)%v)
      message='invalid/nonfinite incident speed_squared or kinetic energy'
      if(.not.ieee_is_finite(speed2(k)).or.speed2(k)<0.d0) goto 900
      energy(k)=0.5d0*sim%groups(g)%mass*ATOMIC_MASS_UNIT*p(j)%weight*speed2(k)
      if(.not.ieee_is_finite(energy(k))) goto 900
      speed=sqrt(speed2(k)); vn=dot_product(p(j)%v,normals(:,k))
      message='materially inward incident velocity at a negative wall crossing'
      if(vn>roundoff_tolerance*speed) goto 900
      message='cumulative reflection counter would overflow'
      if(any(staged%counts==huge(0))) goto 900
      if(abs(vn)<=roundoff_tolerance*speed) staged%counts(3)=staged%counts(3)+1
      if(nsides>1) staged%counts(5)=staged%counts(5)+1
      products(k)=p(j)
      products(k)%i_elm=e
      ! The locator has no exact impact phi: retain continuous end-of-substep phi.
      products(k)%x(1:2)=[r,z]
      staged%counts(1:2)=staged%counts(1:2)+1
    enddo
    ! Phase 2: sample on the copied stream, exactly two uniforms per nonzero speed.
    ! Zero-speed events consume none. q/weight/st/i_life/t_birth are copied exactly.
    do k=1,n
      j=indices(k)
      if(speed2(k)>0.d0) then
        call staged%rng%next(u)
        direction=sample_cosine(u,normals(:,k))
        message='nonfinite or non-inward cosine sample'
        if(.not.all(ieee_is_finite(direction)).or.dot_product(direction,normals(:,k))<=0.d0) goto 900
        products(k)%v=sqrt(speed2(k))*direction
      endif
      emitted=0.5d0*sim%groups(g)%mass*ATOMIC_MASS_UNIT*p(j)%weight*sum(products(k)%v**2)
      message='incident and reflected kinetic energies disagree'
      if(.not.ieee_is_finite(emitted).or.abs(emitted-energy(k))> &
         roundoff_tolerance*max(emitted,energy(k),tiny(emitted))) goto 900
      staged%totals=staged%totals+[p(j)%weight,p(j)%weight,energy(k),emitted]
      message='cumulative reflection totals overflow'
      if(.not.all(ieee_is_finite(staged%totals))) goto 900
    enddo
  end select
  error=0
  return
900 continue
  call reject(sim,g,j,message)
  write(*,*) ' impact/end time=',impact_time,' dt=',dt,' mass=',sim%groups(g)%mass
  if(j>=1.and.j<=size(sim%groups(g)%particles)) then
    select type(p=>sim%groups(g)%particles)
    type is(particle_kinetic_leapfrog)
      write(*,*) ' signed last element=',p(j)%i_elm,' x=',p(j)%x,' st=',p(j)%st, &
                 ' v=',p(j)%v,' q=',p(j)%q,' weight=',p(j)%weight
    end select
  endif
end subroutine

!> Cumulative local totals never reset or replaced by a reduced buffer.
subroutine report_gk_neutral_reflection(sim,state)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_reflection_state), intent(inout) :: state
  integer :: error
  if(.not.state%enabled) return
  call MPI_ALLREDUCE(state%counts,state%global_counts,7,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,error)
  call MPI_ALLREDUCE(state%totals,state%global_totals,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
  if(sim%my_id==0) then
    write(*,*) 'GK reflection cumulative incident/reflected/grazing/unresolved/multiple/corner/failed=',state%global_counts
    write(*,*) 'GK reflection incident/reflected weight, incident/reflected energy_J=',state%global_totals
    write(*,*) 'Passive newly_lost diagnostics count wall crossings BEFORE reflection; no energy deposited at wall.'
  endif
end subroutine
end module mod_gk_neutral_reflection
