module mod_particle_slot_pool_test
use fruit
use mod_particle_types, only: particle_gc_vpar
use mod_particle_sim, only: particle_sim
use mod_particle_slot_pool
implicit none
private
public :: run_slot_pool_tests, assert_slot_registry_equal, assert_slot_pool_equal, snapshot_slot_registry
contains

! Explicit component snapshot avoids ifx's global-derived-allocatable assignment
! descriptor mismatch; comparison still includes every entry and allocation flag.
subroutine snapshot_slot_registry(source,destination)
  type(particle_slot_registry), intent(in) :: source
  type(particle_slot_registry), intent(out) :: destination
  integer :: g
  destination%initialized=source%initialized
  if(.not.allocated(source%pool)) return
  allocate(destination%pool(size(source%pool)))
  do g=1,size(source%pool)
    destination%pool(g)%initialized=source%pool(g)%initialized
    destination%pool(g)%group_num=source%pool(g)%group_num
    destination%pool(g)%n_slots=source%pool(g)%n_slots
    destination%pool(g)%n_free=source%pool(g)%n_free
    if(allocated(source%pool(g)%free_indices)) then
      allocate(destination%pool(g)%free_indices(size(source%pool(g)%free_indices)))
      destination%pool(g)%free_indices=source%pool(g)%free_indices
    endif
    if(allocated(source%pool(g)%free_position)) then
      allocate(destination%pool(g)%free_position(size(source%pool(g)%free_position)))
      destination%pool(g)%free_position=source%pool(g)%free_position
    endif
  enddo
end subroutine

subroutine assert_slot_pool_equal(a,b,label)
  type(particle_slot_pool), intent(in) :: a,b
  character(len=*), intent(in) :: label
  call assert_true(a%initialized.eqv.b%initialized,label//' initialized')
  call assert_equals(a%group_num,b%group_num,label//' group')
  call assert_equals(a%n_slots,b%n_slots,label//' size')
  call assert_equals(a%n_free,b%n_free,label//' free count')
  call assert_true(allocated(a%free_indices).eqv.allocated(b%free_indices),label//' stack allocated')
  call assert_true(allocated(a%free_position).eqv.allocated(b%free_position),label//' inverse allocated')
  if(allocated(a%free_indices).and.allocated(b%free_indices)) &
    call assert_true(all(a%free_indices==b%free_indices),label//' complete stack (including inactive entries)')
  if(allocated(a%free_position).and.allocated(b%free_position)) &
    call assert_true(all(a%free_position==b%free_position),label//' complete inverse map')
end subroutine

subroutine assert_slot_registry_equal(a,b,label)
  type(particle_slot_registry), intent(in) :: a,b
  character(len=*), intent(in) :: label
  integer :: g
  call assert_true(a%initialized.eqv.b%initialized,label//' registry initialized')
  call assert_true(allocated(a%pool).eqv.allocated(b%pool),label//' registry allocated')
  if(.not.allocated(a%pool).or..not.allocated(b%pool)) return
  call assert_equals(size(a%pool),size(b%pool),label//' registry size')
  do g=1,min(size(a%pool),size(b%pool))
    call assert_slot_pool_equal(a%pool(g),b%pool(g),label)
  enddo
end subroutine

subroutine run_slot_pool_tests()
  call run_test_case(pool_operations,'slot pool ordering, transactions and rejection')
  call run_test_case(corruption_checks,'slot pool full consistency validation and repair')
  call run_test_case(registry_checks,'registry lifecycle, selected groups and charged-only')
end subroutine

subroutine mixed(p,pool)
  type(particle_gc_vpar), intent(out) :: p(9)
  type(particle_slot_pool), intent(inout) :: pool
  integer :: error
  p%i_elm=1; p(3)%i_elm=-7; p([2,5,9])%i_elm=0
  call rebuild_particle_slot_pool(p,1,pool,error)
  call assert_equals(0,error,'build mixed states')
end subroutine

subroutine pool_operations()
  type(particle_gc_vpar) :: p(9),empty(0)
  type(particle_slot_pool) :: pool,saved
  integer :: error,slots(3),again(3),one(1),none(0),four(4),k
  call peek_particle_slots(pool,0,none,error)
  call assert_true(error/=0,'uninitialized pool rejected')
  call rebuild_particle_slot_pool(empty,1,pool,error)
  call assert_equals(0,error,'zero-length build')
  call validate_particle_slot_pool(empty,pool,error)
  call assert_equals(0,error,'zero-length validation')
  call peek_particle_slots(pool,0,none,error)
  call assert_equals(0,error,'empty no-op peek')
  call commit_particle_slot_consumption(pool,none)
  call commit_particle_slot_release(pool,none)
  call validate_particle_slot_pool(empty,pool,error)
  call assert_equals(0,error,'zero-length commits')
  call mixed(p,pool)
  saved=pool
  call assert_equals(3,pool%n_free,'only exact zero counted')
  call assert_true(.not.is_particle_slot_free(pool,3),'negative not free')
  call assert_true(.not.is_particle_slot_free(pool,1),'active not free')
  call peek_particle_slots(pool,3,slots,error)
  call assert_equals(0,error,'peek three')
  call assert_true(all(slots==[2,5,9]),'initial low-index order')
  call peek_particle_slots(pool,3,again,error)
  call assert_true(all(slots==again),'repeat peek identical')
  call assert_slot_pool_equal(pool,saved,'peek immutable')
  call peek_particle_slots(pool,0,none,error)
  call assert_equals(0,error,'zero request')
  call peek_particle_slots(pool,-1,none,error)
  call assert_true(error/=0,'negative request')
  call peek_particle_slots(pool,4,slots,error)
  call assert_true(error/=0,'capacity/output shape rejected')
  four=-99
  call peek_particle_slots(pool,4,four,error)
  call assert_true(error/=0,'insufficient capacity with correctly sized output')
  call assert_true(all(four==-99),'failed capacity peek preserves output')
  call rebuild_particle_slot_pool(p,0,pool,error)
  call assert_true(error/=0,'invalid rebuild group rejected')
  call assert_slot_pool_equal(pool,saved,'failed peeks immutable')
  call check_particle_slots_for_consumption(pool,[2,2],error)
  call assert_true(error/=0,'duplicate consume')
  call check_particle_slots_for_consumption(pool,[1],error)
  call assert_true(error/=0,'nonfree consume')
  call check_particle_slots_for_consumption(pool,[0],error)
  call assert_true(error/=0,'consume below range')
  call check_particle_slots_for_consumption(pool,[10],error)
  call assert_true(error/=0,'consume above range')
  call check_particle_slots_for_release(pool,[1,1],error)
  call assert_true(error/=0,'duplicate release')
  call check_particle_slots_for_release(pool,[2],error)
  call assert_true(error/=0,'already free release')
  call check_particle_slots_for_release(pool,[0],error)
  call assert_true(error/=0,'release below range')
  call check_particle_slots_for_release(pool,[10],error)
  call assert_true(error/=0,'release above range')
  call assert_slot_pool_equal(pool,saved,'all rejected checks immutable')
  do k=1,3
    call peek_particle_slots(pool,1,one,error)
    call assert_equals(slots(k),one(1),'sequential deterministic consume')
    call check_particle_slots_for_consumption(pool,one,error)
    call assert_equals(0,error,'consume preflight')
    p(one)%i_elm=1
    call commit_particle_slot_consumption(pool,one)
    call validate_particle_slot_pool(p,pool,error)
    call assert_equals(0,error,'consumed pool agrees with particles')
    call assert_equals(3-k,pool%n_free,'consume count')
  enddo
  call check_particle_slots_for_release(pool,[2,5],error)
  call assert_equals(0,error,'release preflight before physical deactivation')
  p([2,5])%i_elm=0
  call commit_particle_slot_release(pool,[2,5])
  call validate_particle_slot_pool(p,pool,error)
  call assert_equals(0,error,'release consistency')
  call peek_particle_slots(pool,1,one,error)
  call assert_equals(5,one(1),'last released is first reused')
  call assert_equals(2,pool%n_free,'release count')
  ! Arbitrary consumption uses inverse-map swap removal, not a scan.
  call check_particle_slots_for_consumption(pool,[2],error)
  p(2)%i_elm=1
  call commit_particle_slot_consumption(pool,[2])
  call validate_particle_slot_pool(p,pool,error)
  call assert_equals(0,error,'non-top consume mapping')
  call assert_true(is_particle_slot_free(pool,5),'other free slot preserved')
end subroutine

subroutine corruption_checks()
  type(particle_gc_vpar) :: p(9)
  type(particle_slot_pool) :: pool
  integer :: k,error
  do k=1,9
    call mixed(p,pool)
    select case(k)
    case(1); pool%n_free=10
    case(2); pool%free_indices(1)=pool%free_indices(2)
    case(3); pool%free_position(2)=1
    case(4); p(2)%i_elm=1
    case(5); p(4)%i_elm=0
    case(6); p(2)%i_elm=-1
    case(7); pool%free_indices(1)=10
    case(8); pool%n_slots=8
    case(9); pool%n_free=-1
    end select
    call validate_particle_slot_pool(p,pool,error)
    call assert_true(error/=0,'corrupted pool detected')
    call rebuild_particle_slot_pool(p,1,pool,error)
    call assert_equals(0,error,'rebuild repairs corruption')
    call validate_particle_slot_pool(p,pool,error)
    call assert_equals(0,error,'repaired pool validates')
  enddo
end subroutine

subroutine registry_checks()
  type(particle_sim) :: sim
  type(particle_slot_registry) :: registry,saved
  integer :: g,error,k
  do k=2,3
    if(allocated(sim%groups)) deallocate(sim%groups)
    allocate(sim%groups(k))
    do g=1,k
      allocate(particle_gc_vpar::sim%groups(g)%particles(4))
      sim%groups(g)%particles%i_elm=[1,0,-2,0]
    enddo
    call initialize_particle_slot_registry(sim,[(g,g=1,k)],registry,error)
    call assert_equals(0,error,'registry with existing charged/neutral groups')
    call validate_particle_slot_registry(sim,registry,[(g,g=1,k)],error)
    call assert_equals(0,error,'registry validates')
    saved=registry
    call initialize_particle_slot_registry(sim,[k+1],registry,error)
    call assert_true(error/=0,'missing group rejected without mutation')
    call assert_slot_registry_equal(registry,saved,'failed registry init')
    call initialize_particle_slot_registry(sim,[1,1],registry,error)
    call assert_true(error/=0,'duplicate registry group rejected')
    call assert_slot_registry_equal(registry,saved,'duplicate registry init')
    call finalize_particle_slot_registry(registry)
    call assert_true(.not.registry%initialized.and..not.allocated(registry%pool),'finalization')
    call finalize_particle_slot_registry(registry)
    call initialize_particle_slot_registry(sim,[1],registry,error)
    call assert_equals(0,error,'safe reinitialize selected subset')
    call assert_true(.not.registry%pool(2)%initialized,'unselected pool remains uninitialized')
    call validate_particle_slot_registry(sim,registry,[2],error)
    call assert_true(error/=0,'uninitialized selected pool rejected')
  enddo
end subroutine
end module
