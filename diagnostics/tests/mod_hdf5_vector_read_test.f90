! Standalone deterministic regression: make hdf5_vector_read_test, then run with
! OMP_NUM_THREADS=1 mpirun -n 1 (or 2/3) ./hdf5_vector_read_test.
! "positive_only" also runs against the historical reader to demonstrate failure
! without invoking its unchecked invalid/oversized allocation paths.
module mod_hdf5_vector_read_test
  use mpi
  use iso_fortran_env, only: int64
  use fruit
  use fruit_mpi, only: fruit_summary_mpi
  use hdf5_io_module
  implicit none
  private
  public :: run_hdf5_vector_read_test
  integer, parameter :: n=64
  integer :: rank,n_mpi,error,i,j,cycle,unit,scalar_i
  integer(int64) :: stamp
  integer(HID_T) :: file_id,space_id,dataset_id,plist
  integer(SIZE_T) :: objects_before,objects_after
  real*8 :: vector3(3,n),vector2(2,n),scalar_r
  real*8 :: source3(2,3,8),source4(2,3,2,8),source5(2,3,2,2,8)
  real*8, allocatable :: r1(:),r2(:,:),r3(:,:,:),r4(:,:,:,:),r5(:,:,:,:,:)
  real*4, allocatable :: s1(:)
  integer, allocatable :: i1(:),i2(:,:)
  character(len=160) :: filename
  character(len=32) :: mode
  logical :: ok,all_ok

contains
  subroutine run_hdf5_vector_read_test()
  call MPI_INIT(error)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,error)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_mpi,error)
  call init_fruit
  call fruit_hide_dots
  call get_command_argument(1,mode)
  do j=1,n
    do i=1,3
      vector3(i,j)=100*j+i
    enddo
  enddo
  vector2=vector3(1:2,:)
  source3=reshape([(real(i,8),i=1,size(source3))],shape(source3))
  source4=reshape([(real(i,8),i=1,size(source4))],shape(source4))
  source5=reshape([(real(i,8),i=1,size(source5))],shape(source5))
  call H5open_f(error)
  if (rank==0) then
    call system_clock(count=stamp)
    write(filename,'(A,I0,A)') 'hdf5_vector_read_',stamp,'.h5'
  endif
  call MPI_BCAST(filename,len(filename),MPI_CHARACTER,0,MPI_COMM_WORLD,error)
  if (rank==0) then
    call H5Fcreate_f(trim(filename),H5F_ACC_EXCL_F,file_id,error)
    if (error/=0) call MPI_ABORT(MPI_COMM_WORLD,1,error)
    call HDF5_integer_saving(file_id,17,'scalar_int')
    call HDF5_real_saving(file_id,1.25_8,'scalar_real')
    call HDF5_array1D_saving_int(file_id,int(vector3(1,:)),n,'integer1')
    call HDF5_array1D_saving_r4(file_id,real(vector3(1,:),4),n,'single1')
    call HDF5_array1D_saving(file_id,vector3(1,:),n,'double1')
    call HDF5_array2D_saving_int(file_id,int(vector3),3,n,'integer2')
    call HDF5_array2D_saving(file_id,vector3,3,n,'vector3')
    call HDF5_array2D_saving(file_id,vector2,2,n,'vector2')
    call HDF5_array3D_saving(file_id,source3,2,3,8,'double3')
    call HDF5_array4D_saving(file_id,source4,2,3,2,8,'double4')
    call HDF5_array5D_saving(file_id,source5,2,3,2,2,8,'double5')
    ! No payload is needed for a zero-extent dataset.
    call H5Screate_simple_f(2,[3_HSIZE_T,0_HSIZE_T],space_id,error)
    call H5Dcreate_f(file_id,'empty',H5T_NATIVE_DOUBLE,space_id,dataset_id,error)
    call H5Dclose_f(dataset_id,error)
    call H5Sclose_f(space_id,error)
    ! Huge logical extents with no allocated chunks: exercise byte-overflow
    ! rejection without creating a large fixture or attempting an allocation.
    call H5Screate_simple_f(2,[4294967296_HSIZE_T,536870912_HSIZE_T],space_id,error)
    call H5Pcreate_f(H5P_DATASET_CREATE_F,plist,error)
    call H5Pset_chunk_f(plist,2,[1_HSIZE_T,1_HSIZE_T],error)
    call H5Dcreate_f(file_id,'huge',H5T_STD_I8LE,space_id,dataset_id,error,dcpl_id=plist)
    call H5Dclose_f(dataset_id,error)
    call H5Pclose_f(plist,error)
    call H5Sclose_f(space_id,error)
    call H5Fclose_f(file_id,error)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,error)
  ! Two independent opens verify close/reopen behavior as well as repeated reads.
  do cycle=1,2
    call H5Fopen_f(trim(filename),H5F_ACC_RDONLY_F,file_id,error)
    call H5Fget_obj_count_f(file_id,H5F_OBJ_ALL_F,objects_before,error)
    call positive_reads()
    if (trim(mode)/='positive_only') call rejected_reads()
    call H5Fget_obj_count_f(file_id,H5F_OBJ_ALL_F,objects_after,error)
    call assert_true(objects_before==objects_after,'dataset/file identifiers do not leak')
    call H5Fclose_f(file_id,error)
    call assert_equals(0,error,'close after all reads')
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,error)
  if (rank==0) then
    open(newunit=unit,file=trim(filename),status='old')
    close(unit,status='delete')
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,error)
  call H5close_f(error)
  call fruit_summary_mpi(n_mpi,rank)
  call is_all_successful(ok)
  call MPI_ALLREDUCE(ok,all_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,error)
  if (.not.all_ok) call MPI_ABORT(MPI_COMM_WORLD,1,error)
  call MPI_FINALIZE(error)
  end subroutine run_hdf5_vector_read_test

  subroutine check2(expected,label)
    real*8, intent(in) :: expected(:,:)
    character(len=*), intent(in) :: label
    logical :: matches
    matches=.false.
    if (allocated(r2)) then
      if (all(shape(r2)==shape(expected))) matches=all(r2==expected)
    endif
    call assert_true(matches,label)
  end subroutine check2

  subroutine positive_reads()
    integer :: first,last,local_n
    integer(HSIZE_T) :: offset,count
    call HDF5_integer_reading(file_id,scalar_i,'scalar_int')
    call assert_equals(17,scalar_i,'integer scalar')
    call HDF5_real_reading(file_id,scalar_r,'scalar_real')
    call assert_true(scalar_r==1.25_8,'real scalar')
    call HDF5_allocatable_array1D_reading_int(file_id,i1,'integer1')
    call assert_true(all(i1==int(vector3(1,:))),'integer 1D full')
    call HDF5_allocatable_array1D_reading_r4(file_id,s1,'single1')
    call assert_true(all(s1==real(vector3(1,:),4)),'single 1D full')
    call HDF5_allocatable_array1D_reading(file_id,r1,'double1')
    call assert_true(all(r1==vector3(1,:)),'double 1D full')
    call HDF5_allocatable_array1D_reading_int(file_id,i1,'integer1',[5_HSIZE_T],[11_HSIZE_T])
    call assert_true(all(i1==int(vector3(1,12:16))),'integer 1D partial')
    call HDF5_allocatable_array1D_reading_r4(file_id,s1,'single1',[5_HSIZE_T],[11_HSIZE_T])
    call assert_true(all(s1==real(vector3(1,12:16),4)),'single 1D partial')
    call HDF5_allocatable_array1D_reading(file_id,r1,'double1',[5_HSIZE_T],[11_HSIZE_T])
    call assert_true(all(r1==vector3(1,12:16)),'double 1D partial')
    call HDF5_allocatable_array1D_reading_int(file_id,i1,'integer1',[0_HSIZE_T],[64_HSIZE_T])
    call assert_equals(0,size(i1),'integer 1D empty at end')
    call HDF5_allocatable_array1D_reading_r4(file_id,s1,'single1',[0_HSIZE_T],[64_HSIZE_T])
    call assert_equals(0,size(s1),'single 1D empty at end')
    call HDF5_allocatable_array1D_reading(file_id,r1,'double1',[0_HSIZE_T],[64_HSIZE_T])
    call assert_equals(0,size(r1),'double 1D empty at end')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3')
    call check2(vector3,'3-component full')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[3_HSIZE_T,64_HSIZE_T],[0_HSIZE_T,0_HSIZE_T])
    call check2(vector3,'explicit [3,64]')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[-1_HSIZE_T,64_HSIZE_T],[0_HSIZE_T,0_HSIZE_T])
    call check2(vector3,'sentinel [-1,64] must restore exact [3,64]')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector2',[-1_HSIZE_T,64_HSIZE_T],[0_HSIZE_T,0_HSIZE_T])
    call check2(vector2,'2-component sentinel')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[2_HSIZE_T,5_HSIZE_T],[1_HSIZE_T,11_HSIZE_T])
    call check2(vector3(2:3,12:16),'partial components and markers')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[-1_HSIZE_T,5_HSIZE_T],[0_HSIZE_T,11_HSIZE_T])
    call check2(vector3(:,12:16),'sentinel components, partial markers')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[2_HSIZE_T,-1_HSIZE_T],[1_HSIZE_T,0_HSIZE_T])
    call check2(vector3(2:3,:),'sentinel markers, partial components')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[0_HSIZE_T,-7_HSIZE_T],[0_HSIZE_T,0_HSIZE_T])
    call check2(vector3,'legacy nonpositive counts mean full dimensions')
    ! Either optional argument alone retains the historical whole-dataset read.
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',reqdims_in=[1_HSIZE_T,1_HSIZE_T])
    call check2(vector3,'request without start retains full read')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',start=[1_HSIZE_T,1_HSIZE_T])
    call check2(vector3,'start without request retains full read')
    call HDF5_allocatable_array2D_reading(file_id,r2,'absent')
    call check2(vector3,'missing optional dataset leaves caller array unchanged')
    ! Physical partitions differ between MPI sizes, including unequal blocks.
    first=rank*n/n_mpi+1; last=(rank+1)*n/n_mpi
    local_n=last-first+1; offset=int(first-1,HSIZE_T); count=int(local_n,HSIZE_T)
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[-1_HSIZE_T,count],[0_HSIZE_T,offset])
    call check2(vector3(:,first:last),'rank-local sentinel partition')
    call HDF5_allocatable_array2D_reading_int(file_id,i2,'integer2',[-1_HSIZE_T,count],[0_HSIZE_T,offset])
    ok=.false.
    if (allocated(i2)) then
      if (all(shape(i2)==[3,local_n])) ok=all(i2==int(vector3(:,first:last)))
    endif
    call assert_true(ok,'integer 2D sentinel partition')
    call HDF5_allocatable_array3D_reading(file_id,r3,'double3', &
         [-1_HSIZE_T,-1_HSIZE_T,3_HSIZE_T],[0_HSIZE_T,0_HSIZE_T,2_HSIZE_T])
    ok=.false.
    if (allocated(r3)) then
      if (all(shape(r3)==[2,3,3])) ok=all(r3==source3(:,:,3:5))
    endif
    call assert_true(ok,'3D two sentinels and marker offset')
    call HDF5_allocatable_array4D_reading(file_id,r4,'double4', &
         [-1_HSIZE_T,-1_HSIZE_T,-1_HSIZE_T,3_HSIZE_T],[0_HSIZE_T,0_HSIZE_T,0_HSIZE_T,2_HSIZE_T])
    ok=.false.
    if (allocated(r4)) then
      if (all(shape(r4)==[2,3,2,3])) ok=all(r4==source4(:,:,:,3:5))
    endif
    call assert_true(ok,'4D three sentinels and marker offset')
    call HDF5_allocatable_array5D_reading(file_id,r5,'double5', &
         [-1_HSIZE_T,-1_HSIZE_T,-1_HSIZE_T,-1_HSIZE_T,3_HSIZE_T], &
         [0_HSIZE_T,0_HSIZE_T,0_HSIZE_T,0_HSIZE_T,2_HSIZE_T])
    ok=.false.
    if (allocated(r5)) then
      if (all(shape(r5)==[2,3,2,2,3])) ok=all(r5==source5(:,:,:,:,3:5))
    endif
    call assert_true(ok,'5D four sentinels and marker offset')
    call HDF5_allocatable_array2D_reading(file_id,r2,'empty',[-1_HSIZE_T,-1_HSIZE_T],[0_HSIZE_T,0_HSIZE_T])
    call assert_true(all(shape(r2)==[3,0]),'zero-extent dataset retains [3,0]')
  end subroutine positive_reads

  subroutine rejected_reads()
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[4_HSIZE_T,1_HSIZE_T],[0_HSIZE_T,0_HSIZE_T])
    call assert_true(.not.allocated(r2),'component count beyond file rejected')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[1_HSIZE_T,1_HSIZE_T],[-1_HSIZE_T,0_HSIZE_T])
    call assert_true(.not.allocated(r2),'negative offset rejected')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[1_HSIZE_T,1_HSIZE_T],[0_HSIZE_T,65_HSIZE_T])
    call assert_true(.not.allocated(r2),'offset beyond file rejected')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[-1_HSIZE_T,1_HSIZE_T],[1_HSIZE_T,0_HSIZE_T])
    call assert_true(.not.allocated(r2),'sentinel means full extent, not remaining extent')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[1_HSIZE_T,huge(1_HSIZE_T)], &
                                       [0_HSIZE_T,1_HSIZE_T])
    call assert_true(.not.allocated(r2),'huge count rejected without count+offset overflow')
    call HDF5_allocatable_array2D_reading(file_id,r2,'vector3',[1_HSIZE_T,1_HSIZE_T], &
                                       [0_HSIZE_T,huge(1_HSIZE_T)])
    call assert_true(.not.allocated(r2),'huge offset rejected')
    call HDF5_allocatable_array2D_reading(file_id,r2,'double1')
    call assert_true(.not.allocated(r2),'rank mismatch rejected')
    call HDF5_allocatable_array2D_reading(file_id,r2,'huge')
    call assert_true(.not.allocated(r2),'allocation byte overflow rejected before allocate')
    call HDF5_allocatable_array1D_reading_int(file_id,i1,'integer1',[-1_HSIZE_T],[0_HSIZE_T])
    call assert_true(.not.allocated(i1),'rank-1 integer negative count rejected')
    call HDF5_allocatable_array1D_reading_r4(file_id,s1,'single1',[-1_HSIZE_T],[0_HSIZE_T])
    call assert_true(.not.allocated(s1),'rank-1 single negative count rejected')
    call HDF5_allocatable_array1D_reading(file_id,r1,'double1',[-1_HSIZE_T],[0_HSIZE_T])
    call assert_true(.not.allocated(r1),'rank-1 double negative count rejected')
    call HDF5_allocatable_array1D_reading(file_id,r1,'double1',[0_HSIZE_T],[65_HSIZE_T])
    call assert_true(.not.allocated(r1),'even an empty selection must have a valid offset')
  end subroutine rejected_reads
end module mod_hdf5_vector_read_test
