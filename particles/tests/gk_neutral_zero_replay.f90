! Replay a saved binary FE mesh and text particle state; optional rollback/roundoff mode.
! Usage: mpirun -n 16 ./gk_neutral_zero_replay zero_mesh_rank0008.bin zero_event_rank0008.txt reflect
! The saved case-C state is exact. Roundoff mode uses the separately recorded,
! printed rank-2/marker-7 state; its decimal precision is validated against the FE map.
! Modes long / long_rollback accept the rank-15/marker-82 full-precision abort capture
! directly and test its 65-visit crossing through passive advance and reflection.
program gk_neutral_zero_replay
use mpi
use data_structure
use mod_particle_types,only:particle_kinetic_leapfrog
use mod_gk_neutrals
use mod_gk_neutral_reflection
use mod_particle_sim,only:particle_sim
use mod_fields_linear,only:jorek_fields_interp_linear
use mod_particle_types,only:particle_gc_vpar
use mod_particle_slot_pool
use mod_pcg32_rng,only:pcg32_rng
use mod_sampling,only:sample_cosine
use constants,only:ATOMIC_MASS_UNIT
use mod_boundary,only:wall_normal_vector
use phys_module,only:use_manual_random_seed,manual_seed
use mod_gk_neutral_trajectory
use mod_find_rz_nearby,only:find_RZ_nearby
use mod_element_rtree,only:populate_element_rtree
use mod_interp,only:interp_RZ
use phys_module,only:find_RZ_nearby_iter,find_RZ_nearby_tol,i_plane_rtree
implicit none
type(particle_sim)::sim
type(gk_neutral_state)::passive
type(gk_neutral_reflection_state)::reflection,reflection_before
type(particle_slot_registry)::pool,pool_before
type(pcg32_rng)::expected_rng
type(particle_kinetic_leapfrog),allocatable::snapshot(:)
type(gk_neutral_event),allocatable::events_before(:)
integer,allocatable::indices_before(:)
integer::nloss,g,status,bad,marker,n_events
real*8::wl,el,draw(2),want_v(3),normal(3),energy0
logical::rollback,active,long_case,refine_case,vertex_case
character(len=64)::mode
type(type_node_list),target::nodes
type(type_element_list),target::elements
type(particle_kinetic_leapfrog)::old,p,before,endpoint
type(side_correction)::hit,subhit
integer::rank,np,ierr,u,nn,ne,dims(5),i,j,e,owner,identity(6),expected(2),raw_error
real*8::dt,time,expected_x(3),expected_v(3),returned_st(2),mwt(3),xpost(3),vpost(3)
real*8::r,z,rs,rt,zs,zt
character(len=512)::mesh_file,state_file,line
character(len=64)::label
logical::ok,candidate_sides(4)
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr);call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
call init_threads()
call get_command_argument(1,mesh_file);call get_command_argument(2,state_file)
call get_command_argument(3,mode)
rollback=index(mode,'rollback')>0
vertex_case=index(mode,'vertex')==1
long_case=index(mode,'long')==1
refine_case=index(mode,'refine')==1
if(rollback.and.np<2) call MPI_ABORT(MPI_COMM_WORLD,90,ierr)
open(newunit=u,file=trim(state_file),status='old')
do
 read(u,'(A)',iostat=ierr)line
 if(ierr/=0)exit
 read(line,*)label
 select case(trim(label))
 case('identity');read(line,*)label,identity
 case('time_dt');read(line,*)label,time,dt
 case('old_element');read(line,*)label,old%i_elm
 case('old_st');read(line,*)label,old%st
 case('old_x');read(line,*)label,old%x
 case('old_v');read(line,*)label,old%v
 case('ballistic_x');read(line,*)label,expected_x;xpost=expected_x
 case('ballistic_v');read(line,*)label,expected_v;vpost=expected_v
 case('locator');read(line,*)label,expected
 case('locator_st');read(line,*)label,returned_st
 case('post_locator_x');read(line,*)label,xpost
 case('post_locator_v');read(line,*)label,vpost
 case('mass_weight_q');read(line,*)label,mwt
 case('weight_charge_life_birth');read(line,*)label,old%weight,old%q,old%i_life,old%t_birth
 case('life');read(line,*)label,old%i_life
 case('birth');read(line,*)label,old%t_birth
 case('locator_iterations');read(line,*)label,find_RZ_nearby_iter
 case('locator_tolerance');read(line,*)label,find_RZ_nearby_tol
 end select
enddo
close(u);old%weight=mwt(2);old%q=mwt(3)
owner=identity(1);if(owner>=np)owner=0
 open(newunit=u,file=trim(mesh_file),status='old',access='stream',form='unformatted')
 read(u)nn,ne;read(u)dims
 call init_node_list(nodes,nn,1,n_var);elements%n_elements=ne
 if(any(dims(:3)/=shape(nodes%node(1)%x)).or.any(dims(4:)/=shape(elements%element(1)%size))) &
  call MPI_ABORT(MPI_COMM_WORLD,10,ierr)
 do i=1,nn
  read(u)nodes%node(i)%x,nodes%node(i)%boundary
 enddo
 do e=1,ne
  read(u)elements%element(e)%size,elements%element(e)%vertex,elements%element(e)%neighbours
 enddo
 close(u);i_plane_rtree=1
 call populate_element_rtree(nodes,elements)
 if(trim(mode)=='roundoff') then
  owner=2;if(owner>=np)owner=0
  old%i_elm=31916;old%st=[0.524271007149745d0,-3.469446951953614d-18]
  old%x=[0.490710283162818d0,-0.375695212594015d0,0.379663992403519d0]
  old%v=[-124583.817190361d0,39312.6126229714d0,15017.8850049275d0]
  p=old;call push_gk_neutral_ballistic(p,dt)
  expected_x=p%x;expected_v=p%v;xpost=p%x;vpost=p%v
 endif
 bad=mod(owner+1,np);marker=identity(4);if(trim(mode)=='roundoff')marker=7
 p=old;call push_gk_neutral_ballistic(p,dt);endpoint=p
 ok=all(p%x==expected_x).and.all(p%v==expected_v).and.all(xpost==expected_x).and.all(vpost==expected_v)
 if(.not.ok)call MPI_ABORT(MPI_COMM_WORLD,11,ierr)
 call locate(old,p,raw_error)
 write(*,'(A,4I10,2ES26.17)')'RAW rank/marker/element/ifail/st=',rank,marker,p%i_elm,raw_error,p%st
 ok=.true.
 if(trim(mode)/='roundoff')ok=p%i_elm==0.and.raw_error==99
 if(refine_case)ok=p%i_elm==-31722.and.raw_error==-1
 if(vertex_case)ok=p%i_elm==expected(1).and.raw_error==expected(2)
 ok=ok.and.all(p%x==expected_x).and.all(p%v==expected_v)
 ok=ok.and.p%weight==old%weight.and.p%q==old%q.and.p%i_life==old%i_life.and.p%t_birth==old%t_birth
 if(.not.ok)call MPI_ABORT(MPI_COMM_WORLD,12,ierr)
 write(*,*)'EXACT CAPTURED ENDPOINT, RAW LOCATOR RESULT AND PARTICLE IDENTITY PASS'
 call verify_neutral_trajectory(nodes,elements,old%i_elm,old%st,old%x,old%v,dt,hit,endpoint%x,endpoint%v)
 write(*,'(A,5I10,3ES26.17)')'VERIFY status/element/side/visits/crossings/lambda/st=', &
  hit%status,hit%element,hit%side,hit%n_visits,hit%n_trace,hit%lambda,hit%st
 write(*,*)'VERIFY reason=',trim(hit%reason)
 write(*,'(A,6ES26.17)')'VERIFY wall/normal/vn=',hit%wall,hit%normal,hit%vn
 do i=1,hit%n_trace
  e=hit%trace_elements(i);j=hit%trace_faces(i)
  write(*,'(A,4I10,ES26.17)')'TRACE visit/element/face/neighbour/fraction=', &
   i,e,j,elements%element(e)%neighbours(j),hit%trace_fractions(i)
 enddo
 call interp_RZ(nodes,elements,old%i_elm,old%st(1),old%st(2),r,rs,rt,z,zs,zt)
 write(*,'(A,2ES26.17)')'PRE residual/J=',norm2([r,z]-old%x(:2)),rs*zt-rt*zs
 call require(hit%status==WALL_CROSSING,'shared verifier classifies a physical wall')
 if(vertex_case)then
  call require(hit%local_clusters>0.and..not.hit%fallback_attempted,'local tolerance cluster avoids global search')
  write(*,'(A,I8,2ES26.17)')'LOCAL_CLUSTER ',hit%local_clusters,hit%cluster_tolerance,hit%cluster_extent
  endpoint=old;call push_gk_neutral_ballistic(endpoint,hit%lambda*dt)
  call require(norm2(endpoint%x-hit%impact_x)<1.d-12,'full physical contact position')
  call require(all(endpoint%v==hit%impact_v),'incident velocity at contact angle')
  endpoint=old;call push_gk_neutral_ballistic(endpoint,dt)
  if(old%i_elm==31127)then
   call require(hit%element==31638.and.hit%side==2,'first capture independently verified wall')
   call require(abs(hit%lambda-.550622293075919637d0)<1.d-12,'first capture high-precision impact time')
  else
   call require(hit%element==32036.and.hit%side==1,'second capture independently verified wall')
   call require(abs(hit%lambda-.832045198834392804d0)<1.d-12,'second capture high-precision impact time')
  endif
 else if(refine_case) then
  call require(hit%element==31722.and.hit%side==2,'marker 104 first wall 31722 side 2')
  call require(hit%n_visits==10.and.hit%n_trace==10,'marker 104 complete ten-element traversal')
  call require(all(hit%trace_elements(:10)==[29683,29938,30193,30192,30447,30702,30957,31212,31467,31722]), &
    'marker 104 independent chronological element sequence')
  call require(all(hit%trace_faces(:10)==[2,2,1,2,2,2,2,2,2,2]),'marker 104 exit side sequence')
  call require(all(hit%trace_fractions(2:10)>hit%trace_fractions(1:9)),'marker 104 strict fraction ordering')
  call require(abs(hit%trace_fractions(4)-0.22907180112516258d0)<5.d-13,'refined visit-four fraction')
  call require(abs(hit%trace_u(4)-0.94146538684537465d0)<5.d-13,'refined visit-four side coordinate')
  call interp_RZ(nodes,elements,30192,1.d0,hit%trace_u(4),r,z)
  before=old;call push_gk_neutral_ballistic(before,dt*hit%trace_fractions(4))
  call require(norm2(before%x(:2)-[r,z])<2.d-15,'refined visit-four residual near roundoff')
  call require(abs(hit%lambda-0.8226199226895325d0)<2.d-10,'marker 104 wall fraction')
  call require(norm2(hit%st-[1.d0,0.1777579235981458d0])<2.d-11,'marker 104 wall coordinates')
  call require(hit%refinement_counts(1)>0.and.hit%refinement_counts(2)==0,'successful bounded refinement recorded')
  write(*,'(A,4I8,2ES26.17)') 'REFINEMENT counts/residuals=',hit%refinement_counts,hit%refinement_residuals
  write(*,'(A,4ES26.17)') 'REFINED visit-four u/lambda/residual/wall-lambda=', &
    hit%trace_u(4),hit%trace_fractions(4),norm2(before%x(:2)-[r,z]),hit%lambda
 else if(long_case) then
  call require(hit%element==31693.and.hit%side==2,'rank-15 first wall 31693 side 2')
  call require(hit%n_visits==65.and.hit%n_trace==65,'64 interior crossings and wall on visit 65')
  call require(abs(hit%lambda-0.8442583996100523d0)<1.d-12,'rank-15 wall fraction')
  call require(norm2(hit%st-[1.d0,0.001166304336192911d0])<1.d-11,'rank-15 wall coordinates')
  call require(all(hit%trace_fractions(2:65)>hit%trace_fractions(1:64)),'strict chronological fraction ordering')
  call require(hit%trace_elements(64)==31438.and.hit%trace_elements(65)==31693,'former limit crosses into wall element')
 else if(trim(mode)/='roundoff') then
  call require(hit%element==31621.and.hit%side==2,'captured first wall 31621 side 2')
  call require(hit%n_visits==18.and.hit%n_trace==18,'captured 17 interior interfaces, 18 visits')
  call require(abs(hit%lambda-0.40001300851850724d0)<1.d-12,'captured wall fraction')
  call require(norm2(hit%st-[1.d0,0.88469327940539966d0])<1.d-11,'captured wall coordinates')
 else
  call require(hit%start_clamped,'recorded roundoff start independently validated and clipped locally')
  call require(residual_start()<128*epsilon(1.d0),'recorded position matches boundary geometry')
  before=old;before%st(2)=-1.d-10
  call verify_neutral_trajectory(nodes,elements,before%i_elm,before%st,before%x,before%v,dt,subhit)
  call require(subhit%status==UNRESOLVED,'nearby material excursion is rejected')
  before=old;before%x(2)=before%x(2)+1.d-8
  call verify_neutral_trajectory(nodes,elements,before%i_elm,before%st,before%x,before%v,dt,subhit)
  call require(subhit%status==UNRESOLVED,'inconsistent physical start is rejected')
 endif
 sim%my_id=rank;sim%n_mpi=np;sim%time=time
 allocate(jorek_fields_interp_linear::sim%fields)
 sim%fields%node_list=>nodes;sim%fields%element_list=>elements
 allocate(sim%groups(3));sim%groups%coupling_scheme='non';sim%groups%mass=mwt(1)
 sim%groups(1)%id='ion';sim%groups(2)%id='ele';sim%groups(3)%id='neu'
 allocate(particle_gc_vpar::sim%groups(1)%particles(1),sim%groups(2)%particles(1))
 allocate(particle_kinetic_leapfrog::sim%groups(3)%particles(marker))
 active=rank==owner.or.(rollback.and.rank==bad)
 select type(q=>sim%groups(3)%particles)
 type is(particle_kinetic_leapfrog)
  if(active) q(marker)=old
 end select
 use_manual_random_seed=.true.;manual_seed=12345
 call initialize_gk_neutrals(sim,3,passive)
 call initialize_gk_neutral_reflection(sim,3,.true.,reflection)
 call initialize_particle_slot_registry(sim,[1,2,3],pool)
 pool_before=pool;reflection_before=reflection
 call advance_gk_neutrals(sim,passive,dt,1)
 call get_gk_neutral_loss_totals(passive,nloss,wl,el)
 call require(nloss==merge(1,0,active).and.passive%n_newly_lost_indices==nloss,'exactly one local wall event')
 if(vertex_case)then
  call require(passive%zero_count==0.and.passive%verified_wall_count==nloss,'one confirmed negative-locator wall event')
 else if(refine_case) then
  call require(passive%zero_count==0.and.passive%verified_wall_count==nloss,'negative-locator physical wall counted')
  if(active)call require(passive%refinement_counts(1)>0.and.passive%refinement_counts(2)==0,'passive refinement counters')
 else if(trim(mode)/='roundoff') then
  call require(passive%zero_count==nloss.and.passive%zero_wall_count==nloss,'real passive zero return counted')
 else
  call require(passive%roundoff_count==nloss,'real passive roundoff recovery counted')
 endif
 call require(passive%max_successful_traversed==merge(hit%n_visits,0,active),'maximum successful traversal')
 call require(sum(passive%visit_histogram(:,1))==int(nloss,8).and. &
   sum(passive%visit_histogram(:,2))==0_8,'one successful histogram entry per verified event')
 i=findloc(hit%n_visits<=neutral_visit_upper,.true.,dim=1)
 call require(passive%visit_histogram(i,1)==int(nloss,8),'correct visit histogram bin')
 call require(passive%recovered_count==0.and.passive%unresolved_count==0,'no false recovery or unresolved event')
 call require(all(transfer(reflection%rng,[0_1])==transfer(reflection_before%rng,[0_1])),'geometry consumes no RNG')
 select type(q=>sim%groups(3)%particles)
 type is(particle_kinetic_leapfrog)
  if(active) then
   call require(q(marker)%i_elm==-hit%element.and.all(q(marker)%st==hit%st),'passive verified wall metadata')
   call require(all(q(marker)%x==endpoint%x).and.all(q(marker)%v==endpoint%v),'passive exact overshoot x/v')
   call require(wl==old%weight,'exact lost weight')
   energy0=0.5d0*old%weight*mwt(1)*ATOMIC_MASS_UNIT*sum(endpoint%v**2)
   call require(abs(el/energy0-1.d0)<128*epsilon(1.d0),'lost energy exactly once')
   call require(all(passive%newly_lost_events(1)%x_old==old%x).and. &
     all(passive%newly_lost_events(1)%v_old==old%v).and.all(passive%newly_lost_events(1)%st_old==old%st), &
     'unaltered original event provenance')
   call require(passive%newly_lost_events(1)%element==old%i_elm.and. &
     passive%newly_lost_events(1)%marker==marker.and.passive%newly_lost_events(1)%dt==dt,'event identity')
  endif
  if(rollback.and.rank==bad) then
   if(trim(mode)=='refine_failure_rollback') then
    ! Same saved mesh/start, with a much smaller vZ: the real polynomial candidate
    ! cannot converge within the fixed refinement budget. This is not fault
    ! injection into the solver or a weakened tolerance.
    before=old;before%v(2)=1.d-4
    passive%newly_lost_events(1)%v_old=before%v
    q(marker)=before;call push_gk_neutral_ballistic(q(marker),dt)
    call verify_neutral_trajectory(nodes,elements,before%i_elm,before%st,before%x,before%v,dt,subhit)
    call require(subhit%status==UNRESOLVED.and.subhit%refinement_counts(2)>0,'actual failed refinement fixture')
    write(*,'(A,I6,4I8,1X,A)')'EXPECTED REFINEMENT FAILURE rank=',rank,subhit%refinement_counts,trim(subhit%reason)
    ! A wrong but geometrically valid locator wall with inward endpoint velocity
    ! invokes the existing reflection verifier fallback on this rank.
    ok=.false.
    search_wall: do e=1,elements%n_elements
     candidate_sides=boundary_sides(nodes,elements,e)
     if(count(candidate_sides)/=1)cycle
     do i=1,4
      if(.not.candidate_sides(i))cycle
      select case(i)
      case(1);q(marker)%st=[0.5d0,0.d0]
      case(2);q(marker)%st=[1.d0,0.5d0]
      case(3);q(marker)%st=[0.5d0,1.d0]
      case(4);q(marker)%st=[0.d0,0.5d0]
      end select
      normal=wall_normal_vector(nodes,elements,e,q(marker)%st(1),q(marker)%st(2))
      if(dot_product(q(marker)%v,normal)>0.d0) then
       q(marker)%i_elm=-e;ok=.true.;exit search_wall
      endif
     enddo
    enddo search_wall
    call require(ok,'valid wrong-wall metadata invokes reflection verifier')
   else
    q(marker)%q=1
   endif
  endif
  snapshot=q
 end select
 n_events=passive%n_newly_lost_indices
 events_before=passive%newly_lost_events;indices_before=passive%newly_lost_indices
 expected_rng=reflection%rng
 normal=wall_normal_vector(nodes,elements,hit%element,hit%st(1),hit%st(2));normal=normal/norm2(normal)
 if(active) then
  call expected_rng%next(draw);want_v=sqrt(sum(endpoint%v**2))*sample_cosine(draw,normal)
 endif
 call reflect_newly_lost_gk_neutrals(sim,passive,reflection,sim%time,dt,status)
 call require((status/=0).eqv.rollback,'MPI reflection transaction status')
 select type(q=>sim%groups(3)%particles)
 type is(particle_kinetic_leapfrog)
  if(rollback) then
   call require(all(transfer(q,[0_1])==transfer(snapshot,[0_1])),'rollback preserves all particle bytes')
   call require(all(transfer(reflection%rng,[0_1])==transfer(reflection_before%rng,[0_1])),'rollback preserves RNG')
   call require(all(reflection%counts==reflection_before%counts).and. &
     all(reflection%totals==reflection_before%totals),'rollback preserves diagnostics')
   call require(n_events==passive%n_newly_lost_indices.and.all(indices_before==passive%newly_lost_indices), &
     'rollback preserves event list')
   call require(all(transfer(events_before,[0_1])==transfer(passive%newly_lost_events,[0_1])), &
     'rollback preserves provenance bytes')
  else
   call require(passive%n_newly_lost_indices==0,'reflection consumes events once')
   call require(all(transfer(expected_rng,[0_1])==transfer(reflection%rng,[0_1])),'exactly two RNG draws per event')
   if(active) then
    call require(q(marker)%i_elm==hit%element.and.all(q(marker)%st==hit%st),'reflected active wall metadata')
    call require(norm2(q(marker)%x(:2)-hit%wall)<1.d-13,'reflected correct FE wall position')
    call require(dot_product(q(marker)%v,normal)>0.d0,'diffuse velocity points inward')
    call require(norm2(q(marker)%v-want_v)<1.d-8,'usual cosine sample on unchanged RNG stream')
    call require(abs(sum(q(marker)%v**2)/sum(endpoint%v**2)-1.d0)<128*epsilon(1.d0),'roundoff energy conservation')
    call require(q(marker)%weight==old%weight.and.q(marker)%q==old%q.and.q(marker)%i_life==old%i_life.and. &
      q(marker)%t_birth==old%t_birth,'exact weight charge life birth')
    write(*,'(A,3I10,8ES26.17)') 'REFLECT rank/marker/element/st/wall/v/energy_error=', &
      rank,marker,q(marker)%i_elm,q(marker)%st,q(marker)%x(:2),q(marker)%v, &
      sum(q(marker)%v**2)/sum(endpoint%v**2)-1.d0
   endif
  endif
 end select
 do g=1,3
  call require(pool%pool(g)%n_free==pool_before%pool(g)%n_free.and. &
    all(pool%pool(g)%free_indices==pool_before%pool(g)%free_indices).and. &
    all(pool%pool(g)%free_position==pool_before%pool(g)%free_position),'no free-slot transition')
 enddo
 call validate_particle_slot_registry(sim,pool,[1,2,3],status)
 call require(status==0,'pool physical membership remains valid')
 if(.not.rollback)then
  reflection_before=reflection
  call reflect_newly_lost_gk_neutrals(sim,passive,reflection,sim%time,dt,status)
  call require(status==0.and.all(reflection%counts==reflection_before%counts).and. &
    all(transfer(reflection%rng,[0_1])==transfer(reflection_before%rng,[0_1])),'no duplicate reflection')
 endif
 call report_gk_neutral_diagnostics(sim,passive)
 call report_gk_neutral_reflection(sim,reflection)
 if(rank==0)write(*,*)'CAPTURED PASSIVE / REFLECTION / CONSERVATION / RNG / POOL PASS mode=',trim(mode)
 call MPI_FINALIZE(ierr)
contains
subroutine require(condition,message)
 logical,intent(in)::condition
 character(len=*),intent(in)::message
 if(condition)return
 write(*,*)'FAIL rank=',rank,': ',message
 call MPI_ABORT(MPI_COMM_WORLD,20,ierr)
end subroutine
function residual_start() result(residual)
 real*8::residual,y(2),rs,rt,zs,zt
 call interp_RZ(nodes,elements,old%i_elm,old%st(1),0.d0,y(1),rs,rt,y(2),zs,zt)
 residual=norm2(y-old%x(:2))
 write(*,'(A,2ES26.17)')'ROUNDOFF boundary residual/J=',residual,rs*zt-rt*zs
 call require(rs*zt-rt*zs>0.d0,'positive Jacobian at clipped start')
end function
subroutine locate(start,finish,error)
 type(particle_kinetic_leapfrog),intent(in)::start
 type(particle_kinetic_leapfrog),intent(inout)::finish
 integer,intent(out)::error
 call find_RZ_nearby(nodes,elements,start%x(1),start%x(2),start%st(1),start%st(2),start%i_elm, &
  finish%x(1),finish%x(2),finish%st(1),finish%st(2),finish%i_elm,error,phi=finish%x(3))
end subroutine
end program
