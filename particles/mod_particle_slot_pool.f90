!> Rank-local derived free-slot cache. Particle arrays remain the source of truth.
module mod_particle_slot_pool
use mod_particle_types, only: particle_base
use mod_particle_sim, only: particle_sim
implicit none
private
public :: particle_slot_pool, particle_slot_registry
public :: initialize_particle_slot_registry, rebuild_particle_slot_pool, peek_particle_slots
public :: check_particle_slots_for_consumption, check_particle_slots_for_release
public :: commit_particle_slot_consumption, commit_particle_slot_release, is_particle_slot_free
public :: validate_particle_slot_pool, validate_particle_slot_registry, finalize_particle_slot_registry
public :: check_particle_slot_registry_groups

type :: particle_slot_pool
  integer :: group_num=0, n_slots=0, n_free=0
  logical :: initialized=.false.
  integer, allocatable :: free_indices(:),free_position(:)
end type

type :: particle_slot_registry
  type(particle_slot_pool), allocatable :: pool(:)
  logical :: initialized=.false.
end type
contains

subroutine finish(error,operation,ierr)
  integer, intent(in) :: error
  character(len=*), intent(in) :: operation
  integer, optional, intent(out) :: ierr
  if(present(ierr)) then
    ierr=error
  else if(error/=0) then
    write(*,*) 'ERROR: particle slot pool: ',operation,' status=',error
    error stop 1
  endif
end subroutine

! Constant-time structural check; does not scan entries. Full validation is explicit.
logical function ready(pool)
  type(particle_slot_pool), intent(in) :: pool
  ready=.false.
  if(.not.pool%initialized.or.pool%group_num<1.or.pool%n_slots<0) return
  if(pool%n_free<0.or.pool%n_free>pool%n_slots) return
  if(.not.allocated(pool%free_indices).or..not.allocated(pool%free_position)) return
  if(size(pool%free_indices)/=pool%n_slots.or.size(pool%free_position)/=pool%n_slots) return
  ready=.true.
end function

subroutine rebuild_particle_slot_pool(particles,group_num,pool,ierr)
  class(particle_base), intent(in) :: particles(:)
  integer, intent(in) :: group_num
  type(particle_slot_pool), intent(inout) :: pool
  integer, optional, intent(out) :: ierr
  type(particle_slot_pool) :: staged
  integer :: j,error
  error=1
  if(group_num>0) then
    staged%n_slots=size(particles); staged%group_num=group_num
    allocate(staged%free_indices(staged%n_slots),staged%free_position(staged%n_slots),stat=error)
    if(error==0) then
      staged%free_indices=0; staged%free_position=0
      ! Descending build: initial peeks return the lowest free indices first.
      ! Later releases are pushed in supplied order and reused last-in, first-out.
      do j=staged%n_slots,1,-1
        if(particles(j)%i_elm/=0) cycle
        staged%n_free=staged%n_free+1
        staged%free_indices(staged%n_free)=j
        staged%free_position(j)=staged%n_free
      enddo
      staged%initialized=.true.
      call validate_particle_slot_pool(particles,staged,error)
      if(error==0) then
        call move_alloc(staged%free_indices,pool%free_indices)
        call move_alloc(staged%free_position,pool%free_position)
        pool%group_num=group_num; pool%n_slots=staged%n_slots; pool%n_free=staged%n_free
        pool%initialized=.true.
      endif
    endif
  endif
  call finish(error,'rebuild: invalid group/allocation/consistency',ierr)
end subroutine

subroutine initialize_particle_slot_registry(sim,group_indices,registry,ierr)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: group_indices(:)
  type(particle_slot_registry), intent(inout) :: registry
  integer, optional, intent(out) :: ierr
  type(particle_slot_registry) :: staged
  integer :: j,g,error
  error=1
  if(allocated(sim%groups)) then
    allocate(staged%pool(size(sim%groups)),stat=error)
    if(error==0) then
      do j=1,size(group_indices)
        g=group_indices(j); error=1
        if(g<1.or.g>size(sim%groups)) exit
        if(any(group_indices(:j-1)==g)) exit
        if(.not.allocated(sim%groups(g)%particles)) exit
        call rebuild_particle_slot_pool(sim%groups(g)%particles,g,staged%pool(g),error)
        if(error/=0) exit
      enddo
      if(error==0) then
        call move_alloc(staged%pool,registry%pool)
        registry%initialized=.true.
      endif
    endif
  endif
  call finish(error,'initialize registry: invalid group/array/allocation',ierr)
end subroutine

!> Cheap binding check for transactions; deliberately no full membership scan.
subroutine check_particle_slot_registry_groups(sim,registry,group_indices,ierr)
  class(particle_sim), intent(in) :: sim
  type(particle_slot_registry), intent(in) :: registry
  integer, intent(in) :: group_indices(:)
  integer, optional, intent(out) :: ierr
  integer :: j,g,error
  error=1
  if(registry%initialized.and.allocated(registry%pool).and.allocated(sim%groups)) then
    if(size(registry%pool)==size(sim%groups)) then
      error=0
      do j=1,size(group_indices)
        g=group_indices(j); error=1
        if(g<1.or.g>size(registry%pool)) exit
        if(.not.ready(registry%pool(g))) exit
        if(registry%pool(g)%group_num/=g) exit
        if(.not.allocated(sim%groups(g)%particles)) exit
        if(registry%pool(g)%n_slots/=size(sim%groups(g)%particles)) exit
        error=0
      enddo
    endif
  endif
  call finish(error,'registry/group binding mismatch',ierr)
end subroutine

subroutine peek_particle_slots(pool,n_needed,indices,ierr)
  type(particle_slot_pool), intent(in) :: pool
  integer, intent(in) :: n_needed
  integer, intent(inout) :: indices(:)
  integer, optional, intent(out) :: ierr
  integer :: j,error
  error=1
  if(ready(pool).and.n_needed>=0.and.size(indices)==n_needed) then
    if(n_needed<=pool%n_free) then
      ! Validate the slice before touching even the caller's output buffer.
      call check_particle_slots_for_consumption(pool, &
           pool%free_indices(pool%n_free:pool%n_free-n_needed+1:-1),error)
      if(error==0) then
        do j=1,n_needed
          indices(j)=pool%free_indices(pool%n_free-j+1)
        enddo
      endif
    endif
  endif
  call finish(error,'peek: uninitialized/invalid request/capacity/membership',ierr)
end subroutine

subroutine check_particle_slots_for_consumption(pool,indices,ierr)
  type(particle_slot_pool), intent(in) :: pool
  integer, intent(in) :: indices(:)
  integer, optional, intent(out) :: ierr
  integer :: j,error
  error=1
  if(ready(pool)) then
    error=0
    do j=1,size(indices)
      error=1
      if(.not.is_particle_slot_free(pool,indices(j))) exit
      if(any(indices(:j-1)==indices(j))) exit
      error=0
    enddo
  endif
  call finish(error,'consume: invalid/duplicate/nonfree index',ierr)
end subroutine

subroutine check_particle_slots_for_release(pool,indices,ierr)
  type(particle_slot_pool), intent(in) :: pool
  integer, intent(in) :: indices(:)
  integer, optional, intent(out) :: ierr
  integer :: j,k,error
  error=1
  if(ready(pool)) then
    ! Subtraction avoids overflow when checking the new count.
    if(size(indices)<=pool%n_slots-pool%n_free) then
      error=0
      do j=1,size(indices)
        k=indices(j); error=1
        if(k<1.or.k>pool%n_slots) exit
        if(pool%free_position(k)/=0) exit
        if(any(indices(:j-1)==k)) exit
        error=0
      enddo
    endif
  endif
  ! Caller validates the physical source and sets i_elm=0 before release commit.
  call finish(error,'release: invalid/duplicate/already-free index or capacity',ierr)
end subroutine

logical function is_particle_slot_free(pool,particle_index) result(free)
  type(particle_slot_pool), intent(in) :: pool
  integer, intent(in) :: particle_index
  integer :: p
  free=.false.
  if(.not.ready(pool)) return
  if(particle_index<1.or.particle_index>pool%n_slots) return
  p=pool%free_position(particle_index)
  if(p<1.or.p>pool%n_free) return
  free=pool%free_indices(p)==particle_index
end function

!> Commit-only API: caller must have checked all indices before its success vote.
!> No allocation, validation failure path, particle modification or MPI operation.
subroutine commit_particle_slot_consumption(pool,indices)
  type(particle_slot_pool), intent(inout) :: pool
  integer, intent(in) :: indices(:)
  integer :: j,k,p,last
  do j=1,size(indices)
    k=indices(j); p=pool%free_position(k)
    last=pool%free_indices(pool%n_free)
    pool%free_indices(p)=last; pool%free_position(last)=p
    pool%free_position(k)=0; pool%n_free=pool%n_free-1
  enddo
end subroutine

subroutine commit_particle_slot_release(pool,indices)
  type(particle_slot_pool), intent(inout) :: pool
  integer, intent(in) :: indices(:)
  integer :: j,k
  do j=1,size(indices)
    k=indices(j); pool%n_free=pool%n_free+1
    pool%free_indices(pool%n_free)=k; pool%free_position(k)=pool%n_free
  enddo
end subroutine

subroutine validate_particle_slot_pool(particles,pool,ierr)
  class(particle_base), intent(in) :: particles(:)
  type(particle_slot_pool), intent(in) :: pool
  integer, optional, intent(out) :: ierr
  integer :: j,k,p,error
  error=1
  if(ready(pool)) then
    if(pool%n_slots==size(particles)) then
      error=0
      do p=1,pool%n_free
        k=pool%free_indices(p); error=1
        if(k<1.or.k>pool%n_slots) exit
        if(pool%free_position(k)/=p.or.particles(k)%i_elm/=0) exit
        error=0
      enddo
      if(error==0) then
        do j=1,pool%n_slots
          error=1
          if(particles(j)%i_elm==0) then
            if(.not.is_particle_slot_free(pool,j)) exit
          else
            if(pool%free_position(j)/=0) exit
          endif
          error=0
        enddo
      endif
    endif
  endif
  call finish(error,'full validation: particle membership/inverse map mismatch',ierr)
end subroutine

subroutine validate_particle_slot_registry(sim,registry,group_indices,ierr)
  class(particle_sim), intent(in) :: sim
  type(particle_slot_registry), intent(in) :: registry
  integer, intent(in) :: group_indices(:)
  integer, optional, intent(out) :: ierr
  integer :: j,g,error
  call check_particle_slot_registry_groups(sim,registry,group_indices,error)
  if(error==0) then
    do j=1,size(group_indices)
      g=group_indices(j)
      call validate_particle_slot_pool(sim%groups(g)%particles,registry%pool(g),error)
      if(error/=0) exit
    enddo
  endif
  call finish(error,'full registry validation failed',ierr)
end subroutine

subroutine finalize_particle_slot_registry(registry)
  type(particle_slot_registry), intent(inout) :: registry
  if(allocated(registry%pool)) deallocate(registry%pool)
  registry%initialized=.false.
end subroutine
end module
