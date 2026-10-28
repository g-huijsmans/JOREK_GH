!---------------------------------------------
! file : HDF5_io.f90
! date : 25/01/2006
!  array saving and reading in HDF5 format
!
! Adapted from GYSELA routines
! Documentation concerning the HDF5 file handling
! routines implemented in the present module is
! reported in the JOREK wiki at the link below:
! https://www.jorek.eu/wiki/doku.php?id=hdf5tools
!---------------------------------------------
module hdf5_io_module
#ifdef USE_HDF5  
  use HDF5
  use iso_fortran_env, only: int64
  implicit none

  ! Parameters to be used within the module
  integer,parameter :: master_task=0
  private :: prepare_array_read, close_array_read

  !******************************
  contains
  !******************************

  !---------------------------------------- 
  ! create HDF5 file 
  !----------------------------------------
  ! Create an HDF5 file. For enabling Parallel I/O, the mpi_comm_in  and mpi_info
  ! arguments must be defined and the create_access_plist_in must be set to .true. 
  ! inputs:
  !   filename:               (character)(*) name of the file to be created
  !   create_access_plist_in: (logical)(optional) File access properties list (plist)
  !                           used to control different methods of performing
  !                           I/O in HDF5 file. Available data:
  !                           true: create the H5P_FILE_ACCESS_F property list required
  !                           for storing the MPI communicator information for parallel I/O
  !                           default: HDF5 default file access property list H5P_DEFAULT_F
  !   mpi_comm_in:            (integer)(optional) mpi communicator information
  !   mpi_info:               (integer)(optional) MPI_INFO object used for storing a set of
  !                           key,value table used for defining the behaviour of the 
  !                           MPI communicator (e.g. for performance tuning)
  ! outputs:
  !   file_id: (HID_T) identifier of the file
  !   ierr:    (integer) error code, if success =0
  !----------------------------------------
  subroutine HDF5_create(filename,file_id,ierr,create_access_plist_in,mpi_comm_in,mpi_info)
    implicit none
    character(LEN=*) , intent(in)  :: filename
    integer(HID_T)   , intent(out) :: file_id
    integer, optional, intent(out) :: ierr
    integer, optional, intent(in)  :: mpi_comm_in,mpi_info
    logical, optional, intent(in)  :: create_access_plist_in

    integer        :: ierr_HDF5
    integer(HID_T) :: plist
    logical        :: create_access_plist

    !*** Initialize fortran interface ***
    call H5open_f(ierr_HDF5)

    !*** define tha access properties
    create_access_plist = .false.;
    if(present(create_access_plist_in)) create_access_plist = create_access_plist_in
    if(create_access_plist) then
      call H5Pcreate_f(H5P_FILE_ACCESS_F,plist,ierr_HDF5)
    else
      plist = H5P_DEFAULT_F
    endif
    if(present(mpi_comm_in).and.present(mpi_info).and.create_access_plist) then
      call H5Pset_fapl_mpio_f(plist,mpi_comm_in,mpi_info,ierr_HDF5)
    endif

    !*** Create a new file using default properties ***
    call H5Fcreate_f(trim(filename)//char(0), &
      H5F_ACC_TRUNC_F,file_id,ierr_HDF5,access_prp=plist)
    if(plist.ne.H5P_DEFAULT_F) call H5Pclose_f(plist,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_create

  !---------------------------------------- 
  ! open HDF5 file 
  !----------------------------------------
  ! Open an HDF5 file. For enabling Parallel I/O, the mpi_comm_in  and mpi_info
  ! arguments must be defined and the create_access_plist_in must be set to .true. 
  ! inputs:
  !   filename:               (character)(*) name of the file to be created
  !   create_access_plist_in: (logical)(optional) File access properties list (plist)
  !                           used to control different methods of performing
  !                           I/O in HDF5 file. Available data:
  !                           true: create the H5P_FILE_ACCESS_F property list required
  !                           for storing the MPI communicator information for parallel I/O
  !                           default: HDF5 default file access property list H5P_DEFAULT_F
  !   mpi_comm_in:            (integer)(optional) mpi communicator information
  !   mpi_info:               (integer)(optional) MPI_INFO object used for storing a set of
  !                           key,value table used for defining the behaviour of the 
  !                           MPI communicator (e.g. for performance tuning)
  ! outputs:
  !   file_id: (HID_T) identifier of the file
  !   ierr:    (integer) error code, if success =0
  !----------------------------------------
  subroutine HDF5_open(filename,file_id,ierr,create_access_plist_in,mpi_comm_in,mpi_info)
    implicit none
    character(LEN=*) , intent(in)  :: filename  ! file name
    integer(HID_T)   , intent(out) :: file_id   ! file identifier
    integer, optional, intent(out) :: ierr
    integer, optional, intent(in)  :: mpi_comm_in,mpi_info
    logical, optional, intent(in)  :: create_access_plist_in

    integer        :: ierr_HDF5
    integer(HID_T) :: plist
    logical        :: create_access_plist

    !*** Initialize fortran interface ***
    call H5open_f(ierr_HDF5)

    !*** define tha access properties
    create_access_plist = .false.;
    if(present(create_access_plist_in)) create_access_plist = create_access_plist_in
    if(create_access_plist) then
      call H5Pcreate_f(H5P_FILE_ACCESS_F,plist,ierr_HDF5)
    else
      plist = H5P_DEFAULT_F
    endif
    if(present(mpi_comm_in).and.present(mpi_info).and.create_access_plist) then
      call H5Pset_fapl_mpio_f(plist,mpi_comm_in,mpi_info,ierr_HDF5)
    endif

    !*** open the HDF5 file ***
    call H5Fopen_f(trim(filename)//char(0), &
      H5F_ACC_RDONLY_F,file_id,ierr_HDF5,access_prp=plist)
    if(plist.ne.H5P_DEFAULT_F) call H5Pclose_f(plist,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_open

  !---------------------------------------- 
  ! close HDF5 file 
  !----------------------------------------
  ! inputs:
  !   file_id: (HID_T) identifier of the file
  !   ierr:    (integer) error code, if success =0
  ! outputs:
  !   ierr:    (integer) error code, if success =0
  !----------------------------------------
  subroutine HDF5_close(file_id)
    implicit none
    integer(HID_T), intent(in) :: file_id
    integer :: error
    call H5Fclose_f(file_id,error)
    call H5close_f(error)
  end subroutine HDF5_close

  !---------------------------------------- 
  ! create and set HDF5 parallel IO transfer properties 
  !----------------------------------------
  ! Create a HDF5 - MPI transfer property list and set the behaviour of the
  ! paralle I/O. Presently, available parallel I/O behaviour are:
  ! COLLECTIVE: all MPI tasks must performed the same operation
  ! INDEPENDENT: not all MPI tasks must performed the same operation
  ! actions modifying the structure of a file or of the metadata must be done collectively
  ! inputs: 
  !   mpio_collective_in: (logical)(optional) if true (default) MPIO operations are collective
  ! outputs:
  !   transfer_property: (HID_T) transfer property list defining the MPIO behaviour
  !----------------------------------------
  subroutine HDF5_set_parallel_io_properties(transfer_property,mpio_collective_in)
    implicit none
    integer(HID_T),intent(out)  :: transfer_property
    logical,intent(in),optional :: mpio_collective_in
    integer :: ierr
    logical :: mpio_collective
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    call H5Pcreate_f(H5P_DATASET_XFER_F,transfer_property,ierr)
    if(mpio_collective) then
      call H5Pset_dxpl_mpio_f(transfer_property,H5FD_MPIO_COLLECTIVE_F,ierr) 
    else
      call H5Pset_dxpl_mpio_f(transfer_property,H5FD_MPIO_INDEPENDENT_F,ierr)
    endif
  end subroutine HDF5_set_parallel_io_properties

  !----------------------------------------
  ! check if dataset is in HDF5 file
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset to be found
  ! outputs:
  !   in_file: (logical) true file contains the dataset
  !----------------------------------------
  subroutine HDF5_is_dataset_in_file(file_id,dsetname,in_file)
    use H5LT, only: h5ltfind_dataset_f
    !> inputs:
    integer(HID_T),intent(in)   :: file_id
    character(len=*),intent(in) :: dsetname
    !> outputs:
    logical,intent(out) :: in_file
    in_file = h5ltfind_dataset_f(file_id,trim(dsetname)).ge.1
  end subroutine HDF5_is_dataset_in_file

  !---------------------------------------- 
  ! extract dataset rank and shape
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset to be found
  ! outputs:
  !   rank: (integer) dimensionality of the dataset (number of indexes)
  !   dims: (HSIZE_T) dimensions (shape) of the dataset
  !----------------------------------------
  subroutine HDF5_extract_dataset_rank_shape(file_id,rank,dims,dsetname)
    !> inputs:
    character(len=*),intent(in) ::dsetname
    integer(HID_T), intent(in)  :: file_id
    !> outputs:
    integer,intent(out) :: rank !< rank (number of indexes) of the dataset
    integer(HSIZE_T),dimension(:),allocatable,intent(out) :: dims !< shape of the dataset
    !> variables:
    integer(HID_T) :: dataset_id,dataspace_id !< datatset,dataspace identifier
    integer        :: error !< error flag
    integer(HSIZE_T),dimension(:),allocatable :: maxdims !< maximum dataset shape
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    !> allocate shapes
    allocate(dims(rank)); allocate(maxdims(rank));
    !*** get space size ***
    call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
    !*** Closing ***
    call H5Dclose_f(dataset_id,error)
  end subroutine HDF5_extract_dataset_rank_shape

  !----------------------------------------
  ! Open or create HDF5 file
  ! Takes a property list as optional argument, which can be used to set parameters
  ! for the file opening or access, such as block size, offsets, buffering,
  ! parallel IO and data alignment.
  ! See https://support.hdfgroup.org/HDF5/doc1.6/Files.html for more information.
  !----------------------------------------
  ! inputs:
  !   filename:               (character)(*) name of the file to be created
  !   file_access:            (integer)(optional) integer identifying the property
  !                           specifying how an HDF5 is created and opened, default:
  !                           H5F_ACC_EXCL_F: exclusive access to the file
  !                           (fail if file exists, create otherwise)
  !   create_access_plist_in: (logical)(optional) File access properties list (plist)
  !                           used to control different methods of performing
  !                           I/O in HDF5 file. Available data:
  !                           true: create the H5P_FILE_ACCESS_F property list required
  !                           for storing the MPI communicator information for parallel I/O
  !                           default: HDF5 default file access property list H5P_DEFAULT_F
  !   mpi_comm_in:            (integer)(optional) mpi communicator information
  !   mpi_info:               (integer)(optional) MPI_INFO object used for storing a set of
  !                           key,value table used for defining the behaviour of the 
  !                           MPI communicator (e.g. for performance tuning)
  ! outputs:
  !   file_id: (HID_T) identifier of the file
  !   ierr:    (integer) error code, if success =0
  !----------------------------------------
  subroutine HDF5_open_or_create(filename,file_id,ierr,&
  file_access,create_access_plist_in,mpi_comm_in,mpi_info)
    use mpi, only: MPI_Abort
    implicit none
    character(LEN=*) , intent(in)  :: filename  !< file name
    integer(HID_T)   , intent(out) :: file_id   !< file identifier
    integer, optional, intent(out) :: ierr
    integer          , intent(in), optional :: file_access !< type of access for opening/creating an HDF5 file
    logical          , intent(in), optional :: create_access_plist_in !< Which features to use when opening
    integer          , intent(in), optional :: mpi_comm_in,mpi_info !< mpi communicator and info
    
    integer        :: ierr_HDF5,access_f
    logical        :: file_exists, is_hdf5, create_access_plist
    integer(HID_T) :: plist

    !*** Initialize fortran interface ***
    access_f = H5F_ACC_EXCL_F; if(present(file_access)) access_f = file_access;
    call H5open_f(ierr_HDF5)

    !*** define tha access properties
    create_access_plist = .false.;
    if(present(create_access_plist_in)) create_access_plist = create_access_plist_in
    if(create_access_plist) then
      call H5Pcreate_f(H5P_FILE_ACCESS_F,plist,ierr_HDF5)
    else
      plist = H5P_DEFAULT_F
    endif
    if(present(mpi_comm_in).and.present(mpi_info).and.create_access_plist) then
      call H5Pset_fapl_mpio_f(plist,mpi_comm_in,mpi_info,ierr_HDF5)
    endif

    !*** Test if the file exists ***
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists.and.(access_f.ne.H5F_ACC_TRUNC_F)) then
      !*** Test if it is an HDF5 file ***!
      call H5Fis_hdf5_f(trim(filename)//char(0), is_hdf5, ierr_HDF5)
      if (is_hdf5) then
        !*** Try to open the HDF5 file ***!
        call H5Fopen_f(trim(filename)//char(0), &
          H5F_ACC_RDWR_F,file_id,ierr_HDF5, access_prp=plist)
        if (present(ierr)) ierr = ierr_HDF5
      else
        !*** Present an error ***!
        write(*,*) "ERROR: Invalid HDF5 file, exiting"
        call MPI_Abort(mpi_comm_in, -1, ierr)
      end if
    else
      !*** Try to create an HDF5 file ***
      !*** Create a new file using default properties ***
      call H5Fcreate_f(trim(filename)//char(0), &
        access_f, file_id, ierr_HDF5, access_prp=plist)
      if (present(ierr)) ierr = ierr_HDF5
    end if
    if(plist.ne.H5P_DEFAULT_F) call H5Pclose_f(plist,ierr_HDF5)
  end subroutine HDF5_open_or_create

  !*************************************************
  !  HDF5 WRITING
  !*************************************************
  !---------------------------------------- 
  ! HDF5 saving for a character string. Parallel HDF5 is enabled
  ! if the argument varibales: mpi_rank, n_mpi_tasks, mpi_comm_in are
  ! are defined
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   charvar:              (character)(*) character to be written in file
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   mpi_rank:             (integer)(optional) identifier of the MPI task
  !   n_mpi_tasks:          (integer)(optional) number of MPI tasks
  !   mpi_comm_in:          (integer)(optional) identifier of the MPI communicator
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_char_saving(file_id,charvar,dsetname,mpi_rank,n_mpi_tasks,&
  mpi_comm_in,use_hdf5_parallel_in,mpio_collective_in)
#ifdef __GFORTRAN__
    use mpi_mod
#else
    use mpi
#endif
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    character(LEN=*), intent(in) :: charvar
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer,optional, intent(in) :: mpi_rank,n_mpi_tasks,mpi_comm_in
    logical,optional, intent(in) :: use_hdf5_parallel_in,mpio_collective_in

    integer              :: ierr_HDF5  ! error flag
    integer              :: rank       ! dataset rank
    integer              :: len_char   ! maximum length of strings locally
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: filespace  ! filespace identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
    integer(HID_T)       :: type_id    ! datatype identifier
    integer(HID_T)       :: transfer_property ! property for dset transfer
    logical              :: use_hdf5_parallel ! use parallel dset transfer
    logical              :: mpio_collective   ! MPIO collective/individual toggle

   !*** Find the maximum string length and create the hdf5 string type
   len_char = len(charvar)
   if(present(mpi_comm_in)) then
     call MPI_Allreduce(MPI_IN_PLACE,len_char,1,MPI_INTEGER,&
     MPI_MAX,mpi_comm_in,ierr_HDF5);
   endif

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,ierr_HDF5)
    call h5tset_size_f(type_id,int(len_char,kind=HSIZE_T),ierr_HDF5)

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

   !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1; if(present(n_mpi_tasks).and.present(mpi_rank)) dim(1) = n_mpi_tasks;
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,ierr_HDF5)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),type_id,&
    filespace,dataset,ierr_HDF5)

    !*** Write the character array data to the dataset using ***
    !***  default transfer properties                     ***
    if (present(mpi_rank).and.present(n_mpi_tasks).and.present(mpi_comm_in)) then
      call H5Screate_simple_f(rank,int([1],kind=HSIZE_T),dataspace,ierr_HDF5)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=[int(mpi_rank,kind=HSIZE_T)],count=[int(1,kind=HSIZE_T)],hdferr=ierr_HDF5)
      call H5Dwrite_f(dataset,type_id,[charvar],[int(1,kind=HSIZE_T)],ierr_HDF5, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,ierr_HDF5)
    else
      call h5dwrite_f(dataset,type_id,charvar,dim,ierr_HDF5)
    endif

    !*** Closing ***
    call H5Sclose_f(filespace,ierr_HDF5)
    call H5Dclose_f(dataset,ierr_HDF5)
    call H5Pclose_f(transfer_property,ierr_HDF5)
  end subroutine HDF5_char_saving

  !---------------------------------------- 
  ! HDF5 saving for a 1D array of characters
  ! For parallel applications the variable: mpi_comm_in and
  ! start must be defined.
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array1D:              (character)(*)(:) array of characters to be written in file
  !   dim1:                 (integer) size of the array
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (integer)(1)(optional) starting index of the input data chunk 
  !                         in the global dataset
  !   mpi_comm_in:          (integer)(optional) identifier of the MPI communicator
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array1D_saving_char(file_id,array1D,dim1,dsetname,start,&
  mpi_comm_in,use_hdf5_parallel_in,mpio_collective_in)
#ifdef __GFORTRAN__
    use mpi_mod
#else
    use mpi
#endif
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)                , intent(in) :: file_id   ! file identifier
    character(LEN=*), dimension(:), intent(in) :: array1D
    integer                       , intent(in) :: dim1
    character(LEN=*)              , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Begin position of data
    integer                       , intent(in), optional :: mpi_comm_in
    logical                       , intent(in), optional :: use_hdf5_parallel_in
    logical                       , intent(in), optional :: mpio_collective_in

    integer              :: ii
    integer              :: max_len_char   ! maximum length of strings locally
    integer              :: ierr_HDF5  ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim               ! dataset dimensions
    integer(HID_T)       :: dataset           ! dataset identifier
    integer(HID_T)       :: dataspace         ! dataspace identifier
    integer(HID_T)       :: type_id           ! char datatype identifier
    integer(HID_T)       :: filespace         ! filespace identifier
    integer(HID_T)       :: transfer_property ! Transfer property list identifier
    logical              :: use_hdf5_parallel,mpio_collective

    !*** compute and max reduce the lenght of all strings ***
    !*** for each MPI rank ***
    max_len_char = 0
    do ii=1,dim1
      max_len_char = max(max_len_char,len(array1D(ii)))
    enddo
    if(present(mpi_comm_in)) then
      call MPI_Allreduce(MPI_IN_PLACE,max_len_char,1,MPI_INTEGER,&
      MPI_MAX,mpi_comm_in,ierr_HDF5);
    endif

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Set the maximum char length of HDF5 IO **
    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,ierr_HDF5)
    call h5tset_size_f(type_id,int(max_len_char,kind=SIZE_T),ierr_HDF5)

    !*** Create and initialize filespaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,ierr_HDF5)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),type_id,&
    filespace,dataset,ierr_HDF5)

    !*** Write the character array data to the dataset using ***
    !***  default transfer properties                     ***
    if (present(start)) then
      call H5Screate_simple_f(1,shape(array1D,kind=HSIZE_T),dataspace,ierr_HDF5)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=ierr_HDF5)
      call H5Dwrite_f(dataset,type_id,array1D,shape(array1D,kind=HSIZE_T),ierr_HDF5, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,ierr_HDF5)
    else
      call H5Dwrite_f(dataset,type_id,array1D,dim,ierr_HDF5)
    endif

    !*** Closing ***
    call H5Sclose_f(filespace,ierr_HDF5)
    call H5Dclose_f(dataset,ierr_HDF5)
    call H5Pclose_f(transfer_property,ierr_HDF5)
  end subroutine HDF5_array1D_saving_char

  !---------------------------------------- 
  ! HDF5 saving for an integer. Parallel MPIO enabled 
  ! if the variable mpi_rank and n_mpi_tasks are enabled 
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   intv:                 (integer) integer to be written in HDF5 file
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   mpi_rank:             (integer)(optional) identifier of the MPI task
  !   n_mpi_tasks:          (integer)(optional) number of MPI tasks
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_integer_saving(file_id,intv,dsetname,mpi_rank,&
  n_mpi_tasks,use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    integer         , intent(in) :: intv
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer,optional, intent(in) :: mpi_rank,n_mpi_tasks
    logical,optional, intent(in) :: use_hdf5_parallel_in,mpio_collective_in

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: filespace ! filespace identifier
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
    integer(HID_T)       :: transfer_property ! Transfer property list identifier
    logical              :: use_hdf5_parallel,mpio_collective

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1; if(present(n_mpi_tasks).and.present(mpi_rank)) dim(1) = n_mpi_tasks;
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Create integer dataset ***
    call H5Dcreate_f(file_id,trim(dsetname), &
      H5T_NATIVE_INTEGER,filespace,dataset,error)

    !*** Write the integer data to the dataset ***
    if(present(mpi_rank).and.present(n_mpi_tasks)) then
      call H5Screate_simple_f(rank,int([1],kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=[int(mpi_rank,kind=HSIZE_T)],count=[int(1,kind=HSIZE_T)],hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,[intv],[int(1,kind=HSIZE_T)],error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,intv,dim,error)
    endif

    !*** Closing ***
    call H5Dclose_f(dataset,error)
    call H5Sclose_f(filespace,error)
  end subroutine HDF5_integer_saving

  !---------------------------------------- 
  ! HDF5 saving for a real double. Parallel MPIO enabled 
  ! if the variable mpi_rank and n_mpi_tasks are enabled 
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   rd:                   (real8) double float to be written in HDF5 file
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   mpi_rank:             (integer)(optional) identifier of the MPI task
  !   n_mpi_tasks:          (integer)(optional) number of MPI tasks
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_real_saving(file_id,rd,dsetname,mpi_rank,n_mpi_tasks,&
  use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    real*8          , intent(in) :: rd
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer,optional, intent(in) :: mpi_rank,n_mpi_tasks
    logical,optional, intent(in) :: use_hdf5_parallel_in,mpio_collective_in

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim               ! dataset dimensions
    integer(HID_T)       :: filespace
    integer(HID_T)       :: dataset           ! dataset identifier
    integer(HID_T)       :: dataspace         ! dataspace identifier
    integer(HID_T)       :: transfer_property ! Transfer property list identifier
    logical              :: mpio_collective,use_hdf5_parallel

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1; if(present(n_mpi_tasks).and.present(mpi_rank)) dim(1) = n_mpi_tasks;
    rank   = 1;
 
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Create double dataset ***
    call H5Dcreate_f(file_id,trim(dsetname), &
      H5T_NATIVE_DOUBLE,filespace,dataset,error)

    !*** Write the double data to the dataset ***
    if(present(mpi_rank).and.present(n_mpi_tasks)) then
      call H5Screate_simple_f(rank,int([1],kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=[int(mpi_rank,kind=HSIZE_T)],count=[int(1,kind=HSIZE_T)],hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,[rd],[int(1,kind=HSIZE_T)],error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error)
    endif

    !*** Closing ***
    call H5Dclose_f(dataset,error)
    call H5Sclose_f(filespace,error)
  end subroutine HDF5_real_saving

  !---------------------------------------- 
  ! HDF5 saving for a 1D array of integer. Parallel applications
  ! require the definition of the data chunk starting indexes
  ! in the global dataset table. Data compression is applied 
  ! if data chunking is not used ("start" is not defined)
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array1D:              (integer)(:) array of integers to be written in file
  !   dim1:                 (integer) size of the array
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array1D_saving_int(file_id,array1D,dim1,dsetname,start,compress_level,&
    use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)       , intent(in) :: file_id   ! file identifier
    integer, dimension(:), intent(in) :: array1D
    integer              , intent(in) :: dim1
    character(LEN=*)     , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Begin position of data
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in
    logical                       , intent(in), optional :: mpio_collective_in !< hdf5 dataset transfer properties

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim               ! dataset dimensions
    integer(HID_T)      :: dataset           ! dataset identifier
    integer(HID_T)      :: dataspace         ! dataspace identifier
    integer(HID_T)      :: filespace         ! dataspace identifier
    integer(HID_T)      :: property          ! Property list identifier 
    integer(HID_T)      :: transfer_property ! Transfer property list identifier
    logical             :: mpio_collective,use_hdf5_parallel

    !*** Check property for parallel IO
    mpio_collective   = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in 
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER,&
      filespace,dataset,error,property)

    !*** Write the integer array data to the dataset using ***
    !***  default transfer properties                     ***
    if (present(start)) then
      call H5Screate_simple_f(1,shape(array1D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array1D,shape(array1D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
    call H5Pclose_f(transfer_property,error)
  end subroutine HDF5_array1D_saving_int

  !---------------------------------------- 
  ! HDF5 saving for a 1D array of integer. if use_gatherv 
  ! is true and dim1_all_tasks, displs, mpi_rank, n_mpi,
  ! and mpi_comm_loc are defined,  parallelization
  ! based on the MPI gather of all data by the master task
  ! the data are written by the master task in HDF5 file.
  ! Otherwise, native HDF5 implementation is used
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array1D:              (integer)(:) array of integers of each MPI task
  !   dim1_tot:             (integer) total size of all arrays sum(dim1_all_tasks)
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   use_gatherv:          (logical) if true use gatherv parallelization, HDF5-IO is used is false
  !   dim1_all_tasks:       (integer)(n_mpi)(optional) size of the array of each task
  !   displs:               (integer)(n_mpi)(optional) each element specifies the displacement relative to
  !                         the receive MPI buffer at which to place the incoming data from processes
  !   mpi_rank:             (integer)(optional) identifier of the current MPI task
  !   n_mpi:                (integer)(optional) number of MPI tasks
  !   mpi_comm_loc:         (integer)(optional) MPI communicator identifier
  !   start:                (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer)(optional) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) use native HDF5 parallel if true and if
  !                         use_gatherv is false
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true (default)
  !----------------------------------------
  subroutine HDF5_array1D_saving_int_native_or_gatherv(file_id,array1D,dim1_tot,dsetname,&
    use_gatherv,dim1_all_tasks,displs,mpi_rank,n_mpi,mpi_comm_loc,start,compress_level,&
    use_hdf5_parallel_in,mpio_collective_in)
    use mpi
    implicit none
    integer(HID_T)            , intent(in) :: file_id   ! file identifier
    integer                   , intent(in) :: dim1_tot
    integer, dimension(:)     , intent(in) :: array1D
    character(LEN=*)          , intent(in) :: dsetname  ! dataset name
    logical                   , intent(in) :: use_gatherv
    integer                       , intent(in), optional :: mpi_rank,n_mpi,mpi_comm_loc
    integer, dimension(:)         , intent(in), optional :: dim1_all_tasks,displs
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Begin position of data
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in,mpio_collective_in
    integer                      :: ierr
    integer, dimension(dim1_tot) :: array1D_tot
    logical                      :: use_gatherv_loc,use_hdf5_parallel,mpio_collective

    ! check preset
    mpio_collective=.true.; if(present(mpio_collective_in)) mpio_collective=mpio_collective_in;
    ! check whether gatherv can/should be used default false
    use_gatherv_loc=use_gatherv.and.present(mpi_rank).and.present(n_mpi).and.&
    present(mpi_comm_loc).and.present(dim1_all_tasks).and.present(displs)
    use_hdf5_parallel=.false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel=use_hdf5_parallel_in
    if(use_gatherv_loc) use_hdf5_parallel=.not.use_gatherv_loc

    ! Gather all arrays in one
    if(use_gatherv_loc) then
      call MPI_Gatherv(array1D,dim1_all_tasks(mpi_rank+1),MPI_INTEGER,array1D_tot,&
      dim1_all_tasks,displs,MPI_INTEGER,master_task,mpi_comm_loc,ierr)
      if(mpi_rank.eq.master_task) then
        if(present(start).and.present(compress_level)) then
          call HDF5_array1D_saving_int(file_id,array1D_tot,dim1_tot,dsetname,start=start,&
          compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else if(present(start)) then
          call HDF5_array1D_saving_int(file_id,array1D_tot,dim1_tot,dsetname,start=start,&
          use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        else if(present(compress_level)) then
          call HDF5_array1D_saving_int(file_id,array1D_tot,dim1_tot,dsetname,compress_level=compress_level,&
          use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        else
          call HDF5_array1D_saving_int(file_id,array1D_tot,dim1_tot,dsetname,&
          use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        endif
      endif
    else
      if(present(start).and.present(compress_level)) then
        call HDF5_array1D_saving_int(file_id,array1D,dim1_tot,dsetname,&
        start=start,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else if(present(start)) then
        call HDF5_array1D_saving_int(file_id,array1D,dim1_tot,dsetname,&
        start=start,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else if(present(compress_level)) then
        call HDF5_array1D_saving_int(file_id,array1D,dim1_tot,dsetname,&
        compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else
        call HDF5_array1D_saving_int(file_id,array1D,dim1_tot,dsetname,&
        use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      endif
    endif
  end subroutine HDF5_array1D_saving_int_native_or_gatherv

  !---------------------------------------- 
  ! HDF5 saving for a 2D array of integer. Parallel applications
  ! require the definition of the data chunk starting indexes
  ! in the global dataset table. Data compression is applied 
  ! if data chunking is not used ("start" is not defined)
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array2D:              (integer)(:,:) array of integers to be written in file
  !   dim1:                 (integer) size of the first array dimension
  !   dim2:                 (integer) size of the second array dimension
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (HSIZE_T)(2)(optional) starting indexes of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array2D_saving_int(file_id,array2D,dim1,dim2,dsetname,start,&
  compress_level,use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)           , intent(in) :: file_id   ! file identifier
    integer, dimension(:,:)  , intent(in) :: array2D
    integer                  , intent(in) :: dim1, dim2
    character(LEN=*)         , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T)         , intent(in), optional :: start(2) !< Begin position of data
    integer                  , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                  , intent(in), optional :: use_hdf5_parallel_in !< use parallel hdf5 io
    logical                  , intent(in), optional :: mpio_collective_in !< HDF5 dataset MPI transfer property

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    integer(HID_T)      :: transfer_property ! Transfer property list identifier
    logical             :: use_hdf5_parallel,mpio_collective

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    rank   = 2
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER,&
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset using ***
    !***  default transfer properties                     ***
    if (present(start)) then
      call H5Screate_simple_f(2,shape(array2D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array2D,shape(array2D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array2D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
    call H5Pclose_f(transfer_property,error)
  end subroutine HDF5_array2D_saving_int

  !---------------------------------------- 
  ! HDF5 saving for a 2D array of integer. if use_gatherv 
  ! is true and dim3_all_tasks, displs, mpi_rank, n_mpi,
  ! and mpi_comm_loc are defined,  parallelization
  ! based on the MPI gather of all data by the master task
  ! the data are written by the master task in HDF5 file.
  ! Otherwise, native HDF5 implementation is used
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array2D:              (integer)(:,:) array of integers of each MPI task
  !   dim1_tot:             (integer) total size of first dimension of all arrays 
  !   dim2_tot              (integer) total size of second dimension of all arrays sum(dim1_all_tasks)
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   use_gatherv:          (logical) if true use gatherv parallelization, HDF5-IO is used is false
  !   dim2_all_tasks:       (integer)(n_mpi)(optional) size of the array of each task
  !   displs:               (integer)(n_mpi)(optional) each element specifies the displacement relative to
  !                         the receive MPI buffer at which to place the incoming data from processes
  !   mpi_rank:             (integer)(optional) identifier of the current MPI task
  !   n_mpi:                (integer)(optional) number of MPI tasks
  !   mpi_comm_loc:         (integer)(optional) MPI communicator identifier
  !   start:                (HSIZE_T)(2)(optional) starting index of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer)(optional) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) use native HDF5 parallel if true and if
  !                         use_gatherv is false
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true (default)
  !----------------------------------------
  subroutine HDF5_array2D_saving_int_native_or_gatherv(file_id,array2D,dim1_tot,dim2_tot,&
    dsetname,use_gatherv,dim2_all_tasks,displs,mpi_rank,n_mpi,mpi_comm_loc,start,&
    compress_level,use_hdf5_parallel_in,mpio_collective_in)
    use mpi
    implicit none
    integer(HID_T)            , intent(in) :: file_id   ! file identifier
    integer                   , intent(in) :: dim1_tot,dim2_tot
    integer, dimension(:,:)   , intent(in) :: array2D
    character(LEN=*)          , intent(in) :: dsetname  ! dataset name
    logical                   , intent(in) :: use_gatherv
    integer                       , intent(in), optional :: mpi_rank,n_mpi,mpi_comm_loc
    integer, dimension(:)         , intent(in), optional :: dim2_all_tasks,displs
    integer(HSIZE_T), dimension(2), intent(in), optional :: start !< Begin position of data
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in,mpio_collective_in
    integer                               :: ii,ierr
    integer, dimension(dim1_tot,dim2_tot) :: array2D_tot
    logical                               :: use_gatherv_loc,use_hdf5_parallel,mpio_collective

    ! check preset
    mpio_collective=.true.; if(present(mpio_collective_in)) mpio_collective=mpio_collective_in;
    ! check whether gatherv can/should be used default false
    use_gatherv_loc=use_gatherv.and.present(mpi_rank).and.present(n_mpi).and.&
    present(mpi_comm_loc).and.present(dim2_all_tasks).and.present(displs)
    use_hdf5_parallel=.false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel=use_hdf5_parallel_in
    if(use_gatherv_loc) use_hdf5_parallel=.not.use_gatherv_loc

    ! Gather all arrays in one
    if(use_gatherv_loc) then
      do ii=1,dim1_tot
        call MPI_Gatherv(array2D(ii,:),dim2_all_tasks(mpi_rank+1),MPI_INTEGER,&
        array2D_tot(ii,:),dim2_all_tasks,displs,MPI_INTEGER,&
        master_task,mpi_comm_loc,ierr)
      enddo
      if(mpi_rank.eq.master_task) then
        if(present(start).and.present(compress_level)) then
          call HDF5_array2D_saving_int(file_id,array2D_tot,dim1_tot,dim2_tot,dsetname,&
          start=start,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else if(present(start)) then
          call HDF5_array2D_saving_int(file_id,array2D_tot,dim1_tot,dim2_tot,dsetname,&
          start=start,use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        else if(present(compress_level)) then
          call HDF5_array2D_saving_int(file_id,array2D_tot,dim1_tot,dim2_tot,&
          dsetname,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else
          call HDF5_array2D_saving_int(file_id,array2D_tot,dim1_tot,dim2_tot,dsetname,&
          use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        endif
      endif
    else
      if(present(start).and.present(compress_level)) then
        call HDF5_array2D_saving_int(file_id,array2D,dim1_tot,dim2_tot,dsetname,start=start,&
        compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      else if(present(start)) then
        call HDF5_array2D_saving_int(file_id,array2D,dim1_tot,dim2_tot,dsetname,start=start,&
        use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      else if(present(compress_level)) then
        call HDF5_array2D_saving_int(file_id,array2D,dim1_tot,dim2_tot,dsetname,&
        compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else
        call HDF5_array2D_saving_int(file_id,array2D,dim1_tot,dim2_tot,dsetname,&
        use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      endif
    endif
  end subroutine HDF5_array2D_saving_int_native_or_gatherv

  !---------------------------------------- 
  ! gzip HDF5 saving for a 3D array integer. Parallel applications
  ! require the definition of the data chunk starting indexes
  ! in the global dataset table. Data compression is applied 
  ! if data chunking is not used ("start" is not defined)
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array3D:              (integer)(:,:,:) array of integers to be written in file
  !   dim1:                 (integer) size of the first array dimension
  !   dim2:                 (integer) size of the second array dimension
  !   dim3:                 (integer) size of the third array dimension
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (HSIZE_T)(3)(optional) starting indexes of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array3D_saving_int(file_id,array3D,dim1,dim2,dim3,dsetname,start,&
  compress_level,use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)           , intent(in) :: file_id   ! file identifier
    integer, dimension(:,:,:), intent(in) :: array3D
    integer                  , intent(in) :: dim1, dim2, dim3
    character(LEN=*)         , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(3), intent(in), optional :: start !< Offset of array to write
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in !< parallel hdf5 io 
    logical                       , intent(in), optional :: mpio_collective_in !< HDF5 dataset MPI transfer property

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(3)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    integer(HID_T)      :: transfer_property ! Transfer property list identifier
    logical             :: use_hdf5_parallel,mpio_collective

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    rank   = 3
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    if (present(start)) then
      call H5Screate_simple_f(3,shape(array3D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array3D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array3D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array3D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_saving_int

  !---------------------------------------- 
  ! gzip HDF5 saving for a 1D array of real*4. Parallel applications
  ! require the definition of the data chunk starting indexes
  ! in the global dataset table. Data compression is applied 
  ! if data chunking is not used ("start" is not defined)
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array1D:              (real4)(:) array of single precision floats to be written in file
  !   dim1:                 (integer) size of the first array dimension
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array1D_saving_r4(file_id,array1D,dim1,dsetname,start,compress_level,&
   use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)       , intent(in) :: file_id   ! file identifier
    real(4), dimension(:), intent(in) :: array1D
    integer              , intent(in) :: dim1
    character(LEN=*)     , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T)     , intent(in), optional :: start(1) !< Begin position of data
    integer              , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical              , intent(in), optional :: use_hdf5_parallel_in !< use parallel hdf5 io 
    logical              , intent(in), optional :: mpio_collective_in !< HDF5 dataset MPI transfer property

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    integer(HID_T)      :: transfer_property ! Transfer property list identifier
    logical             :: use_hdf5_parallel,mpio_collective

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_REAL,&
      filespace,dataset,error,property)

    !*** Write the real*4 array data to the dataset using ***
    !***  default transfer properties                     ***
    if (present(start)) then
      call H5Screate_simple_f(1,shape(array1D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_REAL,array1D,shape(array1D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_REAL,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
    call H5Pclose_f(transfer_property,error)
  end subroutine HDF5_array1D_saving_r4

  !---------------------------------------- 
  ! HDF5 saving for a 1D array of real*4. if use_gatherv 
  ! is true and dim1_all_tasks, displs, mpi_rank, n_mpi,
  ! and mpi_comm_loc are defined,  parallelization
  ! based on the MPI gather of all data by the master task
  ! the data are written by the master task in HDF5 file.
  ! Otherwise, native HDF5 implementation is used
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array1D:              (real4)(:) array of floats of each MPI task
  !   dim1_tot:             (integer) total size of all arrays sum(dim1_all_tasks)
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   use_gatherv:          (logical) if true use gatherv parallelization, HDF5-IO is used is false
  !   dim1_all_tasks:       (integer)(n_mpi)(optional) size of the array of each task
  !   displs:               (integer)(n_mpi)(optional) each element specifies the displacement relative to
  !                         the receive MPI buffer at which to place the incoming data from processes
  !   mpi_rank:             (integer)(optional) identifier of the current MPI task
  !   n_mpi:                (integer)(optional) number of MPI tasks
  !   mpi_comm_loc:         (integer)(optional) MPI communicator identifier
  !   start:                (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer)(optional) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) use native HDF5 parallel if true and if
  !                         use_gatherv is false
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true (default)
  !----------------------------------------
  subroutine HDF5_array1D_saving_r4_native_or_gatherv(file_id,array1D,dim1_tot,dsetname,&
    use_gatherv,dim1_all_tasks,displs,mpi_rank,n_mpi,mpi_comm_loc,start,compress_level,&
    use_hdf5_parallel_in,mpio_collective_in)
    use mpi
    implicit none
    integer(HID_T)            , intent(in) :: file_id   ! file identifier
    integer                   , intent(in) :: dim1_tot
    real*4, dimension(:)      , intent(in) :: array1D
    character(LEN=*)          , intent(in) :: dsetname  ! dataset name
    logical                   , intent(in) :: use_gatherv
    integer                       , intent(in), optional :: mpi_rank,n_mpi,mpi_comm_loc
    integer, dimension(:)         , intent(in), optional :: dim1_all_tasks,displs
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Begin position of data
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in,mpio_collective_in
    integer                      :: ierr
    real*4,  dimension(dim1_tot) :: array1D_tot
    logical                      :: use_gatherv_loc,use_hdf5_parallel,mpio_collective

    ! check preset
    mpio_collective=.true.; if(present(mpio_collective_in)) mpio_collective=mpio_collective_in;
    ! check whether gatherv can/should be used default false
    use_gatherv_loc=use_gatherv.and.present(mpi_rank).and.present(n_mpi).and.&
    present(mpi_comm_loc).and.present(dim1_all_tasks).and.present(displs)
    use_hdf5_parallel=.false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel=use_hdf5_parallel_in
    if(use_gatherv_loc) use_hdf5_parallel=.not.use_gatherv_loc

    ! Gather all arrays in one
    if(use_gatherv_loc) then
      call MPI_Gatherv(array1D,dim1_all_tasks(mpi_rank+1),MPI_REAL4,array1D_tot,&
      dim1_all_tasks,displs,MPI_REAL4,master_task,mpi_comm_loc,ierr)
      if(mpi_rank.eq.master_task) then
        if(present(start).and.present(compress_level)) then
          call HDF5_array1D_saving_r4(file_id,array1D_tot,dim1_tot,dsetname,&
          start=start,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else if(present(start)) then
          call HDF5_array1D_saving_r4(file_id,array1D_tot,dim1_tot,dsetname,start=start,&
          use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        else if(present(compress_level)) then
          call HDF5_array1D_saving_r4(file_id,array1D_tot,dim1_tot,dsetname,&
          compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else
          call HDF5_array1D_saving_r4(file_id,array1D_tot,dim1_tot,dsetname,&
          use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        endif
      endif
    else
      if(present(start).and.present(compress_level)) then
        call HDF5_array1D_saving_r4(file_id,array1D,dim1_tot,dsetname,start=start,&
        compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else if(present(start)) then
        call HDF5_array1D_saving_r4(file_id,array1D,dim1_tot,dsetname,start=start,&
        use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      else if(present(compress_level)) then
        call HDF5_array1D_saving_r4(file_id,array1D,dim1_tot,dsetname,&
        compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else
        call HDF5_array1D_saving_r4(file_id,array1D,dim1_tot,dsetname,&
        use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      endif
    endif
  end subroutine HDF5_array1D_saving_r4_native_or_gatherv

  !-------------------------------------------- 
  ! gzip HDF5 saving for a 1D array of real*8. Parallel applications
  ! require the definition of the data chunk starting indexes
  ! in the global dataset table. Data compression is applied 
  ! if data chunking is not used ("start" is not defined)
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array1D:              (real8)(:) array of double precision floats to be written in file
  !   dim1:                 (integer) size of the first array dimension
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array1D_saving(file_id,array1D,dim1,dsetname,start,compress_level,&
  use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)      , intent(in) :: file_id   ! file identifier
    real*8, dimension(:), intent(in) :: array1D
    integer             , intent(in) :: dim1
    character(LEN=*)    , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T)    , intent(in), optional :: start(1) !< Begin position of data
    integer             , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical             , intent(in), optional :: use_hdf5_parallel_in !< use parallel hdf5 io
    logical             , intent(in), optional :: mpio_collective_in !< HDF5 dataset MPI transfer property

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    integer(HID_T)      :: transfer_property ! Transfer property list identifier
    logical             :: use_hdf5_parallel,mpio_collective

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE,&
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    if (present(start)) then
      call H5Screate_simple_f(1,shape(array1D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array1D,shape(array1D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Sclose_f(filespace,error)
    call H5Pclose_f(property,error)
    call H5Dclose_f(dataset,error)
    call H5Pclose_f(transfer_property,error)
  end subroutine HDF5_array1D_saving

  !---------------------------------------- 
  ! HDF5 saving for a 1D array of real*8. if use_gatherv 
  ! is true and dim1_all_tasks, displs, mpi_rank, n_mpi,
  ! and mpi_comm_loc are defined,  parallelization
  ! based on the MPI gather of all data by the master task
  ! the data are written by the master task in HDF5 file.
  ! Otherwise, native HDF5 implementation is used
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array1D:              (real8)(:) array of floats of each MPI task
  !   dim1_tot:             (integer) total size of all arrays sum(dim1_all_tasks)
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   use_gatherv:          (logical) if true use gatherv parallelization, HDF5-IO is used is false
  !   dim1_all_tasks:       (integer)(n_mpi)(optional) size of the array of each task
  !   displs:               (integer)(n_mpi)(optional) each element specifies the displacement relative to
  !                         the receive MPI buffer at which to place the incoming data from processes
  !   mpi_rank:             (integer)(optional) identifier of the current MPI task
  !   n_mpi:                (integer)(optional) number of MPI tasks
  !   mpi_comm_loc:         (integer)(optional) MPI communicator identifier
  !   start:                (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer)(optional) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) use native HDF5 parallel if true and if
  !                         use_gatherv is false
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true (default)
  !----------------------------------------
  subroutine HDF5_array1D_saving_native_or_gatherv(file_id,array1D,dim1_tot,dsetname,&
    use_gatherv,dim1_all_tasks,displs,mpi_rank,n_mpi,mpi_comm_loc,start,compress_level,&
    use_hdf5_parallel_in,mpio_collective_in)
    use mpi
    implicit none
    integer(HID_T)            , intent(in) :: file_id   ! file identifier
    integer                   , intent(in) :: dim1_tot
    real*8, dimension(:)      , intent(in) :: array1D
    character(LEN=*)          , intent(in) :: dsetname  ! dataset name
    logical                   , intent(in) :: use_gatherv
    integer                       , intent(in), optional :: mpi_rank,n_mpi,mpi_comm_loc
    integer, dimension(:)         , intent(in), optional :: dim1_all_tasks,displs
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Begin position of data
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in,mpio_collective_in
    integer                      :: ierr
    real*8,  dimension(dim1_tot) :: array1D_tot
    logical                      :: use_gatherv_loc,use_hdf5_parallel,mpio_collective

    ! check preset
    mpio_collective=.true.; if(present(mpio_collective_in)) mpio_collective=mpio_collective_in;
    ! check whether gatherv can/should be used default false
    use_gatherv_loc=use_gatherv.and.present(mpi_rank).and.present(n_mpi).and.&
    present(mpi_comm_loc).and.present(dim1_all_tasks).and.present(displs)
    use_hdf5_parallel=.false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel=use_hdf5_parallel_in
    if(use_gatherv_loc) use_hdf5_parallel=.not.use_gatherv_loc

    ! Gather all arrays in one
    if(use_gatherv_loc) then
      call MPI_Gatherv(array1D,dim1_all_tasks(mpi_rank+1),MPI_REAL8,array1D_tot,&
      dim1_all_tasks,displs,MPI_REAL8,master_task,mpi_comm_loc,ierr)
      if(mpi_rank.eq.master_task) then
        if(present(start).and.present(compress_level)) then
          call HDF5_array1D_saving(file_id,array1D_tot,dim1_tot,dsetname,start=start,&
          compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else if(present(start)) then
          call HDF5_array1D_saving(file_id,array1D_tot,dim1_tot,dsetname,start=start,&
          use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        else if(present(compress_level)) then
          call HDF5_array1D_saving(file_id,array1D_tot,dim1_tot,dsetname,&
          compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else
          call HDF5_array1D_saving(file_id,array1D_tot,dim1_tot,dsetname,&
          use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        endif
      endif
    else
      if(present(start).and.present(compress_level)) then
        call HDF5_array1D_saving(file_id,array1D,dim1_tot,dsetname,start=start,&
        compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else if(present(start)) then
        call HDF5_array1D_saving(file_id,array1D,dim1_tot,dsetname,&
        start=start,mpio_collective_in=mpio_collective)
      else if(present(compress_level)) then
        call HDF5_array1D_saving(file_id,array1D,dim1_tot,dsetname,&
        compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else
        call HDF5_array1D_saving(file_id,array1D,dim1_tot,dsetname,&
        use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      endif
    endif
  end subroutine HDF5_array1D_saving_native_or_gatherv

  !---------------------------------------- 
  ! gzip HDF5 saving for a 2D array. Parallel applications
  ! require the definition of the data chunk starting indexes
  ! in the global dataset table. Data compression is applied 
  ! if data chunking is not used ("start" is not defined)
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array2D:              (real8)(:,:) array of double precision floats to be written in file
  !   dim1:                 (integer) size of the first array dimension
  !   dim2:                 (integer) size of the second array dimension
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (HSIZE_T)(2)(optional) starting indexes of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array2D_saving(file_id,array2D,dim1,dim2,dsetname,start,compress_level,&
  use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)        , intent(in) :: file_id   ! file identifier
    real*8, dimension(:,:), intent(in) :: array2D
    integer               , intent(in) :: dim1, dim2
    character(LEN=*)      , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(2), intent(in), optional :: start !< Offset of array to write
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in !< use parallel hdf5
    logical                       , intent(in), optional :: mpio_collective_in !< HDF5 dataset MPI transfer property

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
    integer(HID_T)       :: filespace  ! dataspace identifier
    integer(HID_T)       :: property   ! Property list identifier 
    integer(HID_T)       :: transfer_property ! Transfer property list identifier
    logical              :: use_hdf5_parallel,mpio_collective

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    rank   = 2
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    if (present(start)) then
      call H5Screate_simple_f(2,shape(array2D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array2D,shape(array2D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
    call H5Pclose_f(transfer_property,error)
  end subroutine HDF5_array2D_saving

  !---------------------------------------- 
  ! HDF5 saving for a 2D array of real8. if use_gatherv 
  ! is true and dim2_all_tasks, displs, mpi_rank, n_mpi,
  ! and mpi_comm_loc are defined,  parallelization
  ! based on the MPI gather of all data by the master task
  ! the data are written by the master task in HDF5 file.
  ! Otherwise, native HDF5 implementation is used
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array2D:              (real8)(:,:) array of integers of each MPI task
  !   dim1_tot:             (integer) total size of first dimension of all arrays 
  !   dim2_tot              (integer) total size of second dimension of all arrays sum(dim2_all_tasks)
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   use_gatherv:          (logical) if true use gatherv parallelization, HDF5-IO is used is false
  !   dim2_all_tasks:       (integer)(n_mpi)(optional) size of the array of each task
  !   displs:               (integer)(n_mpi)(optional) each element specifies the displacement relative to
  !                         the receive MPI buffer at which to place the incoming data from processes
  !   mpi_rank:             (integer)(optional) identifier of the current MPI task
  !   n_mpi:                (integer)(optional) number of MPI tasks
  !   mpi_comm_loc:         (integer)(optional) MPI communicator identifier
  !   start:                (HSIZE_T)(2)(optional) starting index of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer)(optional) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) use native HDF5 parallel if true and if
  !                         use_gatherv is false
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true (default)
  !----------------------------------------
  subroutine HDF5_array2D_saving_native_or_gatherv(file_id,array2D,dim1_tot,dim2_tot,&
    dsetname,use_gatherv,dim2_all_tasks,displs,mpi_rank,n_mpi,mpi_comm_loc,start,&
    compress_level,use_hdf5_parallel_in,mpio_collective_in)
    use mpi
    implicit none
    integer(HID_T)            , intent(in) :: file_id   ! file identifier
    integer                   , intent(in) :: dim1_tot,dim2_tot
    real*8, dimension(:,:)    , intent(in) :: array2D
    character(LEN=*)          , intent(in) :: dsetname  ! dataset name
    logical                   , intent(in) :: use_gatherv
    integer                       , intent(in), optional :: mpi_rank,n_mpi,mpi_comm_loc
    integer, dimension(:)         , intent(in), optional :: dim2_all_tasks,displs
    integer(HSIZE_T), dimension(2), intent(in), optional :: start !< Begin position of data
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in,mpio_collective_in
    integer                              :: ii,ierr
    real*8, dimension(dim1_tot,dim2_tot) :: array2D_tot
    logical                              :: use_gatherv_loc,use_hdf5_parallel,mpio_collective

    ! check preset
    mpio_collective=.true.; if(present(mpio_collective_in)) mpio_collective=mpio_collective_in;
    ! check whether gatherv can/should be used default false
    use_gatherv_loc=use_gatherv.and.present(mpi_rank).and.present(n_mpi).and.&
    present(mpi_comm_loc).and.present(dim2_all_tasks).and.present(displs)
    use_hdf5_parallel=.false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel=use_hdf5_parallel_in
    if(use_gatherv_loc) use_hdf5_parallel=.not.use_gatherv_loc

    ! Gather all arrays in one
    if(use_gatherv_loc) then
      do ii=1,dim1_tot
        call MPI_Gatherv(array2D(ii,:),dim2_all_tasks(mpi_rank+1),MPI_REAL8,&
        array2D_tot(ii,:),dim2_all_tasks,displs,MPI_REAL8,&
        master_task,mpi_comm_loc,ierr)
      enddo
      if(mpi_rank.eq.master_task) then
        if(present(start).and.present(compress_level)) then
          call HDF5_array2D_saving(file_id,array2D_tot,dim1_tot,dim2_tot,dsetname,&
          start=start,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else if(present(start)) then
          call HDF5_array2D_saving(file_id,array2D_tot,dim1_tot,dim2_tot,dsetname,&
          start=start,mpio_collective_in=mpio_collective,use_hdf5_parallel_in=use_hdf5_parallel)
        else if(present(compress_level)) then
          call HDF5_array2D_saving(file_id,array2D_tot,dim1_tot,dim2_tot,&
          dsetname,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else
          call HDF5_array2D_saving(file_id,array2D_tot,dim1_tot,dim2_tot,dsetname,&
          use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
        endif
      endif
    else
      if(present(start).and.present(compress_level)) then
        call HDF5_array2D_saving(file_id,array2D,dim1_tot,dim2_tot,dsetname,&
        start=start,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else if(present(start)) then
        call HDF5_array2D_saving(file_id,array2D,dim1_tot,dim2_tot,dsetname,start=start,&
        use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      else if(present(compress_level)) then
        call HDF5_array2D_saving(file_id,array2D,dim1_tot,dim2_tot,dsetname,&
        compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else
        call HDF5_array2D_saving(file_id,array2D,dim1_tot,dim2_tot,dsetname,&
        use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      endif
    endif
  end subroutine HDF5_array2D_saving_native_or_gatherv

  !---------------------------------------- 
  ! gzip HDF5 saving for a 3D array. Parallel applications
  ! require the definition of the data chunk starting indexes
  ! in the global dataset table. Data compression is applied 
  ! if data chunking is not used ("start" is not defined)
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array3D:              (real8)(:,:,:) array of double precision floats to be written in file
  !   dim1:                 (integer) size of the first array dimension
  !   dim2:                 (integer) size of the second array dimension
  !   dim3:                 (integer) size of the third array dimension
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (HSIZE_T)(3)(optional) starting indexes of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array3D_saving(file_id,array3D,dim1,dim2,dim3,dsetname,start,&
  compress_level,use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)          , intent(in) :: file_id   ! file identifier
    real*8, dimension(:,:,:), intent(in) :: array3D
    integer                 , intent(in) :: dim1, dim2, dim3
    character(LEN=*)        , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(3), intent(in), optional :: start !< Offset of array to write
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in !< use parallel hdf5
    logical                       , intent(in), optional :: mpio_collective_in !< HDF5 dataset MPI transfer property

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(3)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    integer(HID_T)      :: transfer_property ! Transfer property list identifier
    logical             :: use_hdf5_parallel,mpio_collective

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    rank   = 3
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    if (present(start)) then
      call H5Screate_simple_f(3,shape(array3D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array3D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
    call H5Pclose_f(transfer_property,error)
  end subroutine HDF5_array3D_saving

  !---------------------------------------- 
  ! HDF5 saving for a 3D array of real8. if use_gatherv 
  ! is true and dim3_all_tasks, displs, mpi_rank, n_mpi,
  ! and mpi_comm_loc are defined,  parallelization
  ! based on the MPI gather of all data by the master task
  ! the data are written by the master task in HDF5 file.
  ! Otherwise, native HDF5 implementation is used
  !----------------------------------------
  ! inputs:
  !   file_id:            (HID_T) file identifier
  !   array3D:            (real8)(:,:,:) array of integers of each MPI task
  !   dim1_tot:           (integer) total size of first dimension all arrays 
  !   dim2_tot            (integer) total size of second dimension of all arrays
  !   dim3_tot            (integer) total size of third dimension of all arrays sum(dim3_all_tasks)
  !   dsetname:           (character)(*) name of the dataset in which the data are written
  !   use_gatherv:        (logical) if true use gatherv parallelization, HDF5-IO is used is false
  !   dim3_all_tasks:     (integer)(n_mpi)(optional) size of the array of each task
  !   displs:             (integer)(n_mpi)(optional) each element specifies the displacement relative to
  !                       the receive MPI buffer at which to place the incoming data from processes
  !   mpi_rank:           (integer)(optional) identifier of the current MPI task
  !   n_mpi:              (integer)(optional) number of MPI tasks
  !   mpi_comm_loc:       (integer)(optional) MPI communicator identifier
  !   start:              (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !                       in the global dataset
  !   compress_level:     (integer)(optional) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) use native HDF5 parallel if true and if
  !                         use_gatherv is false
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true (default)
  !----------------------------------------
  subroutine HDF5_array3D_saving_native_or_gatherv(file_id,array3D,dim1_tot,dim2_tot,&
    dim3_tot,dsetname,use_gatherv,dim3_all_tasks,displs,mpi_rank,n_mpi,mpi_comm_loc,start,&
    compress_level,use_hdf5_parallel_in,mpio_collective_in)
    use mpi
    implicit none
    integer(HID_T)            , intent(in) :: file_id   ! file identifier
    integer                   , intent(in) :: dim1_tot,dim2_tot,dim3_tot
    real*8, dimension(:,:,:)  , intent(in) :: array3D
    character(LEN=*)          , intent(in) :: dsetname  ! dataset name
    logical                   , intent(in) :: use_gatherv
    integer                       , intent(in), optional :: mpi_rank,n_mpi,mpi_comm_loc
    integer, dimension(:)         , intent(in), optional :: dim3_all_tasks,displs
    integer(HSIZE_T), dimension(3), intent(in), optional :: start !< Begin position of data
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in,mpio_collective_in
    integer                                       :: ii,jj,ierr
    real*8, dimension(dim1_tot,dim2_tot,dim3_tot) :: array3D_tot
    logical                                       :: use_gatherv_loc,use_hdf5_parallel,mpio_collective

    ! check preset
    mpio_collective=.true.; if(present(mpio_collective_in)) mpio_collective=mpio_collective_in;
    ! check whether gatherv can/should be used default false
    use_gatherv_loc=use_gatherv.and.present(mpi_rank).and.present(n_mpi).and.&
    present(mpi_comm_loc).and.present(dim3_all_tasks).and.present(displs)
    use_hdf5_parallel=.false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel=use_hdf5_parallel_in
    if(use_gatherv_loc) use_hdf5_parallel=.not.use_gatherv_loc

    ! Gather all arrays in one
    if(use_gatherv_loc) then
    do jj=1,dim2_tot
      do ii=1,dim1_tot
        call MPI_Gatherv(array3D(ii,jj,:),dim3_all_tasks(mpi_rank+1),MPI_REAL8,&
        array3D_tot(ii,jj,:),dim3_all_tasks,displs,MPI_REAL8,&
        master_task,mpi_comm_loc,ierr)
      enddo
    enddo
      if(mpi_rank.eq.master_task) then
        if(present(start).and.present(compress_level)) then
          call HDF5_array3D_saving(file_id,array3D_tot,dim1_tot,dim2_tot,dim3_tot,dsetname,&
          start=start,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel)
        else if(present(start)) then
          call HDF5_array3D_saving(file_id,array3D_tot,dim1_tot,dim2_tot,dim3_tot,&
          dsetname,start=start,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else if(present(compress_level)) then
          call HDF5_array3D_saving(file_id,array3D_tot,dim1_tot,dim2_tot,dim3_tot,&
          dsetname,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        else
          call HDF5_array3D_saving(file_id,array3D_tot,dim1_tot,dim2_tot,&
          dim3_tot,dsetname,use_hdf5_parallel_in=use_hdf5_parallel,&
          mpio_collective_in=mpio_collective)
        endif
      endif
    else
      if(present(start).and.present(compress_level)) then
        call HDF5_array3D_saving(file_id,array3D,dim1_tot,dim2_tot,dim3_tot,dsetname,&
        start=start,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else if(present(start)) then
        call HDF5_array3D_saving(file_id,array3D,dim1_tot,dim2_tot,dim3_tot,dsetname,&
        start=start,use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      else if(present(compress_level)) then
        call HDF5_array3D_saving(file_id,array3D,dim1_tot,dim2_tot,dim3_tot,&
        dsetname,compress_level=compress_level,use_hdf5_parallel_in=use_hdf5_parallel,&
        mpio_collective_in=mpio_collective)
      else
        call HDF5_array3D_saving(file_id,array3D,dim1_tot,dim2_tot,dim3_tot,dsetname,&
        use_hdf5_parallel_in=use_hdf5_parallel,mpio_collective_in=mpio_collective)
      endif
    endif
  end subroutine HDF5_array3D_saving_native_or_gatherv

  !---------------------------------------- 
  ! gzip HDF5 saving for a 4D array. Parallel applications
  ! require the definition of the data chunk starting indexes
  ! in the global dataset table. Data compression is applied 
  ! if data chunking is not used ("start" is not defined)
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array4D:              (real8)(:,:,:,:) array of double precision floats to be written in file
  !   dim1:                 (integer) size of the first array dimension
  !   dim2:                 (integer) size of the second array dimension
  !   dim3:                 (integer) size of the third array dimension
  !   dim4:                 (integer) size of the fourth array dimension
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (HSIZE_T)(4)(optional) starting indexes of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array4D_saving(file_id,array4d,dim1,dim2,dim3,dim4,&
  dsetname,start,compress_level,use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)            , intent(in) :: file_id   ! file identifier
    real*8, dimension(:,:,:,:), intent(in) :: array4d
    integer                   , intent(in) :: dim1, dim2, dim3, dim4
    character(LEN=*)          , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(4), intent(in), optional :: start !< Offset of array to write
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in !< parallel hdf5 io
    logical                       , intent(in), optional :: mpio_collective_in !< HDF5 dataset MPI transfer property

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(4)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier
    integer(HID_T)      :: transfer_property ! Transfer property list identifier
    logical             :: use_hdf5_parallel,mpio_collective

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    dim(4) = dim4
    rank   = 4
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***  using default transfer properties ***
    if (present(start)) then
      call H5Screate_simple_f(4,shape(array4D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array4D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array4D,shape(array4D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
    call H5Pclose_f(transfer_property,error)
  end subroutine HDF5_array4D_saving

  !---------------------------------------- 
  ! gzip HDF5 saving for a 5D array. Parallel applications
  ! require the definition of the data chunk starting indexes
  ! in the global dataset table. Data compression is applied 
  ! if data chunking is not used ("start" is not defined)
  !----------------------------------------
  ! inputs:
  !   file_id:              (HID_T) file identifier
  !   array5D:              (real8)(:,:,:,:,:) array of double precision floats to be written in file
  !   dim1:                 (integer) size of the first array dimension
  !   dim2:                 (integer) size of the second array dimension
  !   dim3:                 (integer) size of the third array dimension
  !   dim4:                 (integer) size of the fourth array dimension
  !   dim5:                 (integer) size of the fifth array dimension
  !   dsetname:             (character)(*) name of the dataset in which the data are written
  !   start:                (HSIZE_T)(5)(optional) starting indexes of the input data chunk 
  !                         in the global dataset
  !   compress_level:       (integer) level of data compression to be used
  !   use_hdf5_parallel_in: (logical)(optional) if true, dataset transfer property is
  !                         set to parallel MPIO, H5P_DEFAULT is used otherwise
  !   mpio_collective_in:   (logical)(optional) toggle MPIO collective actions if true
  !----------------------------------------
  subroutine HDF5_array5D_saving(file_id,array5d,dim1,dim2,dim3,dim4,dim5,&
  dsetname,start,compress_level,use_hdf5_parallel_in,mpio_collective_in)
    use hdf5, only: H5P_DEFAULT_F
    implicit none
    integer(HID_T)              , intent(in) :: file_id  ! file identifier
    real*8, dimension(:,:,:,:,:), intent(in) :: array5d
    integer                     , intent(in) :: dim1, dim2
    integer                     , intent(in) :: dim3, dim4, dim5
    character(LEN=*)            , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(5), intent(in), optional :: start !< Offset of array to write
    integer                       , intent(in), optional :: compress_level !< if set and start is not provided compress with this level
    logical                       , intent(in), optional :: use_hdf5_parallel_in !< use parallel hdf5 io
    logical                       , intent(in), optional :: mpio_collective_in !< HDF5 dataset MPI transfer property

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(5)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    integer(HID_T)      :: transfer_property ! Transfer property list identifier
    logical             :: use_hdf5_parallel,mpio_collective

    !*** Check property for parallel IO
    mpio_collective = .true.
    if(present(mpio_collective_in)) mpio_collective = mpio_collective_in
    transfer_property = H5P_DEFAULT_F; use_hdf5_parallel = .false.
    if(present(use_hdf5_parallel_in)) use_hdf5_parallel = use_hdf5_parallel_in
    if(use_hdf5_parallel) call HDF5_set_parallel_io_properties(transfer_property,mpio_collective)

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    dim(4) = dim4
    dim(5) = dim5
    rank   = 5
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***  using default transfer properties         ***
    if (present(start)) then
      call H5Screate_simple_f(5,shape(array5D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array5D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array5D,shape(array5D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace,xfer_prp=transfer_property)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
    call H5Pclose_f(transfer_property,error)
  end subroutine HDF5_array5D_saving

  !************************************************
  !  HDF5 READING
  !************************************************
  !---------------------------------------- 
  ! HDF5 reading for a character. Parallel compatible
  ! reading is enabled defining the argumets mpi_rank
  ! and n_mpi_tasks. 
  !----------------------------------------
  ! inputs:
  !   file_id:     (HID_T) file identifier
  !   dsetname:    (character)(*) name of the dataset from which data are read
  !   mpi_rank:    (integer)(optional) identifier of the MPI task
  !   n_mpi_tasks: (integer)(optional) number of MPI tasks
  ! outputs:
  !   charvar:     (character)(*) character to be read from file
  !----------------------------------------
  subroutine HDF5_char_reading(file_id,charvar,dsetname,mpi_rank,n_mpi_tasks)
    implicit none
    integer(HID_T)  , intent(in)  :: file_id   ! file identifier
    character(LEN=*), intent(out) :: charvar
    character(LEN=*), intent(in)  :: dsetname  ! dataset name
    integer,intent(in),optional   :: mpi_rank,n_mpi_tasks
    integer                       :: error     ! error flag
    integer                       :: rank      ! dataset rank
    integer(HID_T)                :: dataset   ! dataset identifier
    integer(HID_T)                :: filespace ! filespace identifier
    integer(HID_T)                :: dataspace ! dataspace identifier
    integer(HID_T)                :: type_id   ! string type identifier
    integer(HSIZE_T),dimension(1) :: dim       ! dimensions
    integer(HSIZE_T),dimension(1) :: dims,maxdims ! sized of the datsets
    logical                       :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    !*** set the string datatype
    call H5Dget_type_f(dataset,type_id,error)
    !*** get size of the dataspace
    call H5dget_space_f(dataset,dataspace,error)
    call H5Sget_simple_extent_ndims_f(dataspace,rank,error)
    if(rank.eq.1) call H5Sget_simple_extent_dims_f(dataspace,dims,maxdims,error)
    !*** read the string data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = int(1,kind=HSIZE_T)
    if(present(mpi_rank).and.present(n_mpi_tasks).and.(&
    (rank.eq.1).and.(dims(1).gt.1))) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=int([mpi_rank],kind=HSIZE_T),count=int([1],kind=HSIZE_T),hdferr=error)
      call H5Dread_f(dataset,type_id,charvar,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,type_id,charvar,dim,error)
    endif
    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_char_reading

  !---------------------------------------- 
  ! HDF5 reading for a allocatable character. Parallel compatible
  ! reading is enabled defining the argumets mpi_rank
  ! and n_mpi_tasks. 
  !----------------------------------------
  ! inputs:
  !   file_id:     (HID_T) file identifier
  !   dsetname:    (character)(*) name of the dataset from which data are read
  !   mpi_rank:    (integer)(optional) identifier of the MPI task
  !   n_mpi_tasks: (integer)(optional) number of MPI tasks
  ! outputs:
  !   charvar:     (character)(:)(allocatable) character to be read from file
  !----------------------------------------
  subroutine HDF5_allocatable_char_reading(file_id,charvar,dsetname,mpi_rank,n_mpi_tasks)
    implicit none
    integer(HID_T)  , intent(in)                        :: file_id  ! file identifier
    character(LEN=:), allocatable, intent(out)          :: charvar
    character(LEN=*)             , intent(in)           :: dsetname ! dataset name
    integer                      , intent(in), optional :: mpi_rank,n_mpi_tasks
    
    integer                       :: error        ! error flag
    integer                       :: rank         ! rank of the array
    integer(HID_T)                :: dataset      ! dataset identifier
    integer(HID_T)                :: filespace    ! filespace identifier
    integer(HID_T)                :: dataspace    ! dataspace identifier
    integer(HID_T)                :: type_id      ! string type identifier
    integer(HSIZE_T)              :: len_char     ! character length
    integer(HSIZE_T),dimension(1) :: dim          ! dimensions
    integer(HSIZE_T),dimension(1) :: dims,maxdims ! dimensions of the array
    logical                       :: exists       ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    !*** set the string datatype
    call H5Dget_type_f(dataset,type_id,error)
    call H5Tget_size_f(type_id,len_char,error)
    !*** allocate character 
    if(allocated(charvar)) deallocate(charvar)
    allocate(character(len=int(len_char))::charvar)
    !*** get size of the dataspace
    call H5dget_space_f(dataset,dataspace,error)
    call H5Sget_simple_extent_ndims_f(dataspace,rank,error)
    if(rank.eq.1) call H5Sget_simple_extent_dims_f(dataspace,dims,maxdims,error)
    !*** read the string data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = int(1,kind=HSIZE_T)
    if(present(mpi_rank).and.present(n_mpi_tasks).and.&
    (rank.eq.1).and.(dims(1).gt.int(1,kind=HSIZE_T))) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=int([mpi_rank],kind=HSIZE_T),count=int([1],kind=HSIZE_T),hdferr=error)
      call H5Dread_f(dataset,type_id,charvar,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,type_id,charvar,dim,error)
    endif
    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_allocatable_char_reading

  !---------------------------------------- 
  ! HDF5 reading for a character array 1D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset in which the data are written
  !   start:    (HSIZE_T)(1)(optional) starting index of 
  !             the input data chunk in the global dataset
  ! outputs:
  !   array1D:  (character)(*)(:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array1D_reading_char(file_id,array1D,dsetname,start)
    implicit none
    integer(HID_T),   intent(in) :: file_id   ! file identifier
    character(len=*), dimension(:), intent(out) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Offset of array to read
    integer             :: len_char
    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim       ! dataset dimension
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: type_id   ! string datatype identifier
    integer(HID_T)      :: filespace ! filespace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    !*** set the string datatype
    call H5Dget_type_f(dataset,type_id,error)
    !*** read the character data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,type_id,array1D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,type_id,array1D,dim,error)
    endif

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading_char

  subroutine HDF5_array1D_reading_char_len(file_id, array1D, char_len, dsetname)
    use HDF5
    implicit none

    ! Arguments
    integer(HID_T), intent(in)                    :: file_id       ! File identifier
    character(LEN=*), dimension(:), intent(out)   :: array1D       ! Output character array
    integer, intent(in)                           :: char_len      ! Length of array entries
    character(LEN=*), intent(in)                  :: dsetname      ! Dataset name

    ! Local variables
    integer :: ierr_HDF5                          ! Error flag
    integer(HID_T) :: dataset                     ! Dataset identifier
    integer(HID_T) :: dataspace                   ! Dataspace identifier
    integer(HID_T) :: string_type                 ! Custom string datatype
    integer(HSIZE_T), dimension(1) :: dim        ! Dataset dimensions

    ! Open the dataset
    call H5Dopen_f(file_id, trim(dsetname), dataset, ierr_HDF5)

    ! Get the dataspace
    call H5Dget_space_f(dataset, dataspace, ierr_HDF5)

    dim(1) = size(array1D,1)

    ! Define a custom string datatype matching the length of each character
    call H5Tcopy_f(H5T_NATIVE_CHARACTER, string_type, ierr_HDF5)
    call H5Tset_size_f(string_type, int(char_len,kind=SIZE_T), ierr_HDF5)

    ! Read the dataset into the Fortran array
    call H5Dread_f(dataset, string_type, array1D, dim, ierr_HDF5)

    ! Close resources
    call H5Tclose_f(string_type, ierr_HDF5)
    call H5Sclose_f(dataspace, ierr_HDF5)
    call H5Dclose_f(dataset, ierr_HDF5)

end subroutine HDF5_array1D_reading_char_len


  !---------------------------------------- 
  ! HDF5 reading for an integer. Parallel compatible
  ! reading is enabled defining the argumets mpi_rank
  ! and n_mpi_tasks. 
  !----------------------------------------
  ! inputs:
  !   file_id:     (HID_T) file identifier
  !   dsetname:    (character)(*) name of the dataset from which data are read
  !   mpi_rank:    (integer)(optional) identifier of the MPI task
  !   n_mpi_tasks: (integer)(optional) number of MPI tasks
  ! outputs:
  !   intv:        (integer) character to be read from file
  !---------------------------------------- 
  subroutine HDF5_integer_reading(file_id,intv,dsetname,mpi_rank,&
  n_mpi_tasks)
    implicit none
    integer(HID_T)  , intent(in)          :: file_id  ! file identifier
    integer         , intent(out)         :: intv
    character(LEN=*), intent(in)          :: dsetname ! dataset name
    integer,          intent(in),optional :: mpi_rank,n_mpi_tasks 
    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim       ! desired dataset dimensions
    integer(HSIZE_T), &
      dimension(1)      :: dims,maxdims ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: data_type
    integer(HID_T)      :: filespace ! filespace identifier
    logical             :: exists    ! true if dataset exists

    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    !*** get size of the dataspace
    call H5dget_space_f(dataset,dataspace,error)
    call H5Sget_simple_extent_ndims_f(dataspace,rank,error)
    if(rank.eq.1) call H5Sget_simple_extent_dims_f(dataspace,dims,maxdims,error)
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = int(1,kind=HSIZE_T)
    if(present(mpi_rank).and.present(n_mpi_tasks).and.&
    (rank.eq.1).and.(dims(1).gt.1)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=int([mpi_rank],kind=HSIZE_T),count=int([1],kind=HSIZE_T),hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,intv,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,intv,dim,error)
    endif

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_integer_reading

  !---------------------------------------- 
  ! HDF5 reading for an integer array 1D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset in which the data are written
  !   start:    (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !             in the global dataset
  ! outputs:
  !   array1D:  (integer)(:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array1D_reading_int(file_id,array1D,dsetname,start)
    implicit none
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    integer         , intent(out), dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start ! Offset of array to read

    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim       ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: filespace ! filespace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading_int

  !---------------------------------------- 
  ! HDF5 reading for an integer allocatable array 1D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task and the size of the MPI task array in which the
  ! data chuck is stored.
  !----------------------------------------
  ! inputs:
  !   file_id:            (HID_T) file identifier
  !   dsetname:           (character)(*) name of the dataset in which the data are written
  !   reqdims_in:         (HSIZE_T)(1) size of the data chunck to be loaded
  !   start:              (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !                       in the global dataset
  ! outputs:
  !   array1D:            (integer)(:)(allocatable) array of characters to be read from file
  !---------------------------------------- 
  ! Prepare numeric allocatable reads without copying signed sentinels into
  ! HDF5 counts. In particular, keep allocation extents separate from counts:
  ! ifx 2023.2 with -O2 -axCOMMON-AVX512 can reuse a pre-normalization extent
  ! when a copied request array is conditionally overwritten in place.
  ! Legacy semantics: ranks >= 2 interpret ANY nonpositive request as the full
  ! file dimension (not the remainder after start). Rank 1 has no sentinel;
  ! zero is an empty selection. A slab requires BOTH optional arguments.
  subroutine prepare_array_read(file_id,dsetname,expected_rank,element_bytes,reqdims_in,start, &
                                dataset_id,file_space,mem_space,counts,extents,exists,ready)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname
    integer, intent(in) :: expected_rank,element_bytes
    integer(HSIZE_T), optional, intent(in) :: reqdims_in(expected_rank),start(expected_rank)
    integer(HID_T), intent(out) :: dataset_id,file_space,mem_space
    integer(HSIZE_T), intent(out) :: counts(expected_rank)
    integer(int64), intent(out) :: extents(expected_rank)
    logical, intent(out) :: exists,ready
    integer(HSIZE_T) :: file_dims(expected_rank),max_dims(expected_rank),offset(expected_rank)
    integer(int64) :: requested,extent,points
    integer(HSSIZE_T) :: file_points,mem_points
    integer :: rank,error,i,bad_dim
    logical :: slab
    character(len=80) :: reason

    dataset_id=-1; file_space=-1; mem_space=-1
    counts=0; extents=0; file_dims=0; offset=0
    ready=.false.; exists=.false.; rank=-1; bad_dim=0
    reason='dataset existence query failed'
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if (error<0) goto 900
    if (.not.exists) then
      write(*,*) 'WARNING: dataset ',dsetname,' does not exists!'
      return
    endif
    reason='dataset open failed'
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    if (error<0) goto 900
    reason='dataspace query failed'
    call H5Dget_space_f(dataset_id,file_space,error)
    if (error<0) goto 900
    reason='dataset rank query failed'
    call H5Sget_simple_extent_ndims_f(file_space,rank,error)
    if (error<0) goto 900
    reason='dataset rank does not match reader rank'
    if (rank/=expected_rank) goto 900
    reason='dataset extent query failed'
    call H5Sget_simple_extent_dims_f(file_space,file_dims,max_dims,error)
    if (error<0) goto 900
    slab=present(reqdims_in).and.present(start)
    if (slab) offset=start
    do i=1,expected_rank
      bad_dim=i
      reason='file extent is negative or exceeds signed allocation range'
      if (file_dims(i)<0.or.file_dims(i)>huge(extent)) goto 900
      extent=int(file_dims(i),int64)
      if (slab) then
        ! Convert only after the signed request is known to fit. Negative
        ! requests at ranks >= 2 are sentinels, never HDF5 extents/counts.
        reason='requested count exceeds signed allocation range'
        if (reqdims_in(i)>huge(requested)) goto 900
        if (expected_rank>=2.and.reqdims_in(i)<=0) then
          extent=int(file_dims(i),int64)
        else
          requested=int(reqdims_in(i),int64)
          reason='negative count is not supported by the rank-1 reader'
          if (requested<0) goto 900
          extent=requested
        endif
      endif
      reason='offset is negative or beyond the file extent'
      if (offset(i)<0.or.offset(i)>file_dims(i)) goto 900
      reason='requested count exceeds the file extent after offset'
      ! Subtract after checking offset; do not form the overflowing sum
      ! offset+count. A full-dimension sentinel therefore requires start=0.
      if (extent>file_dims(i)-offset(i)) goto 900
      extents(i)=extent
      counts(i)=int(extent,HSIZE_T)
    enddo
    bad_dim=0
    ! Check the allocation byte count and signed HDF5 selected-point range.
    ! Empty arrays still have their individual extents validated above.
    points=0
    if (all(extents>0)) then
      points=1
      do i=1,expected_rank
        bad_dim=i
        reason='selected element count overflows signed integer range'
        if (points>huge(points)/extents(i)) goto 900
        points=points*extents(i)
      enddo
      reason='allocation byte count or HDF5 point count overflows'
      if (points>huge(points)/element_bytes.or.points>huge(file_points)) goto 900
    endif
    bad_dim=0
    reason='memory dataspace creation failed'
    call H5Screate_simple_f(expected_rank,counts,mem_space,error)
    if (error<0) goto 900
    reason='file hyperslab selection failed'
    if (points==0) then
      call H5Sselect_none_f(file_space,error)
      if (error<0) goto 900
      call H5Sselect_none_f(mem_space,error)
    else
      call H5Sselect_hyperslab_f(file_space,H5S_SELECT_SET_F,start=offset,count=counts,hdferr=error)
    endif
    if (error<0) goto 900
    reason='selected file/memory point count query failed'
    call H5Sget_select_npoints_f(file_space,file_points,error)
    if (error<0) goto 900
    call H5Sget_select_npoints_f(mem_space,mem_points,error)
    if (error<0) goto 900
    reason='selected file and memory point counts do not match allocation'
    if (file_points/=points.or.mem_points/=points) goto 900
    ready=.true.
    return
900 continue
    write(*,*) 'WARNING: HDF5 array read rejected: ',trim(dsetname),': ',trim(reason)
    write(*,*) '  dataset rank / reader rank / failing dimension: ',rank,expected_rank,bad_dim
    if (rank==expected_rank) write(*,*) '  file extents: ',file_dims
    if (present(reqdims_in)) write(*,*) '  requested counts: ',reqdims_in
    if (present(start)) write(*,*) '  requested offsets: ',start
    write(*,*) '  normalized counts (zero for dimensions not yet validated): ',counts
    call close_array_read(dataset_id,file_space,mem_space)
  end subroutine prepare_array_read

  subroutine close_array_read(dataset_id,file_space,mem_space)
    integer(HID_T), intent(inout) :: dataset_id,file_space,mem_space
    integer :: error
    if (mem_space>=0) call H5Sclose_f(mem_space,error)
    if (file_space>=0) call H5Sclose_f(file_space,error)
    if (dataset_id>=0) call H5Dclose_f(dataset_id,error)
    dataset_id=-1; file_space=-1; mem_space=-1
  end subroutine close_array_read

  subroutine HDF5_allocatable_array1D_reading_int(file_id,array1D,dsetname,reqdims_in,start)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname
    integer(HSIZE_T), optional, intent(in) :: reqdims_in(1),start(1)
    integer, allocatable, intent(inout) :: array1D(:)
    integer(HID_T) :: dataset_id,file_space,mem_space
    integer(HSIZE_T) :: counts(1)
    integer(int64) :: extents(1)
    integer :: error
    logical :: exists,ready

    call prepare_array_read(file_id,dsetname,1,storage_size(0)/8,reqdims_in,start, &
                            dataset_id,file_space,mem_space,counts,extents,exists,ready)
    ! Preserve the existing missing-dataset behavior, including the caller's array.
    if (.not.exists) return
    if (allocated(array1D)) deallocate(array1D)
    if (.not.ready) return
    allocate(array1D(extents(1)),stat=error)
    if (error==0) then
      ! An empty independent read needs no transfer; retain its allocated shape.
      if (all(extents>0)) then
        array1D=0
        call H5Dread_f(dataset_id,H5T_NATIVE_INTEGER,array1D,counts,error, &
                       file_space_id=file_space,mem_space_id=mem_space)
      endif
    endif
    if (error/=0) then
      write(*,*) 'WARNING: HDF5 array allocation/read failed: ',dsetname,' extents=',extents,' status=',error
      if (allocated(array1D)) deallocate(array1D)
    endif
    call close_array_read(dataset_id,file_space,mem_space)
  end subroutine HDF5_allocatable_array1D_reading_int


  !---------------------------------------- 
  ! HDF5 reading for an integer array 2D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset in which the data are written
  !   start:    (HSIZE_T)(2)(optional) starting indexes of 
  !             the input data chunk in the global dataset
  ! outputs:
  !   array2D:  (integer)(:,:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array2D_reading_int(file_id,array2D,dsetname,start)
    implicit none
    integer(HID_T)  , intent(in)     :: file_id   ! file identifier
    integer         , dimension(:,:) :: array2D
    character(LEN=*), intent(in)     :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(2), intent(in), optional :: start ! Offset of array to read

    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim       ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: filespace ! dataspace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
     if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif  
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim = shape(array2D)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(2,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array2D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array2D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_reading_int

  !---------------------------------------- 
  ! HDF5 reading for an allocatable integer array 2D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task and the size of the MPI task array in which the
  ! data chuck is stored.
  !----------------------------------------
  ! inputs:
  !   file_id:    (HID_T) file identifier
  !   dsetname:   (character)(*) name of the dataset in which the data are written
  !   reqdims_in: (HSIZE_T)(2) size of the data chunck to be loaded
  !   start:      (HSIZE_T)(2)(optional) starting indexes of 
  !               the input data chunk in the global dataset
  ! outputs:
  !   array2D:            (integer)(:,:)(allocatable) array of characters to be read from file
  !----------------------------------------  
  subroutine HDF5_allocatable_array2D_reading_int(file_id,array2D,dsetname,reqdims_in,start)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname
    integer(HSIZE_T), optional, intent(in) :: reqdims_in(2),start(2)
    integer, allocatable, intent(inout) :: array2D(:,:)
    integer(HID_T) :: dataset_id,file_space,mem_space
    integer(HSIZE_T) :: counts(2)
    integer(int64) :: extents(2)
    integer :: error
    logical :: exists,ready

    call prepare_array_read(file_id,dsetname,2,storage_size(0)/8,reqdims_in,start, &
                            dataset_id,file_space,mem_space,counts,extents,exists,ready)
    ! Preserve the existing missing-dataset behavior, including the caller's array.
    if (.not.exists) return
    if (allocated(array2D)) deallocate(array2D)
    if (.not.ready) return
    allocate(array2D(extents(1),extents(2)),stat=error)
    if (error==0) then
      ! An empty independent read needs no transfer; retain its allocated shape.
      if (all(extents>0)) then
        array2D=0
        call H5Dread_f(dataset_id,H5T_NATIVE_INTEGER,array2D,counts,error, &
                       file_space_id=file_space,mem_space_id=mem_space)
      endif
    endif
    if (error/=0) then
      write(*,*) 'WARNING: HDF5 array allocation/read failed: ',dsetname,' extents=',extents,' status=',error
      if (allocated(array2D)) deallocate(array2D)
    endif
    call close_array_read(dataset_id,file_space,mem_space)
  end subroutine HDF5_allocatable_array2D_reading_int

  !---------------------------------------- 
  ! HDF5 reading for an integer array 3D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:            (HID_T) file identifier
  !   dsetname:           (character)(*) name of the dataset in which the data are written
  !   start:              (HSIZE_T)(3)(optional) starting indexes of 
  !                       the input data chunk in the global dataset
  ! outputs:
  !   array1D:            (integer)(:,:,:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array3D_reading_int(file_id,array3D,dsetname,&
       in1, in2, in3,start)
    implicit none
    integer(HID_T)     , intent(in)       :: file_id   ! file identifier
    integer            , dimension(:,:,:) :: array3D
    character(LEN=*)   , intent(in)       :: dsetname  ! dataset name
    integer, intent(in), optional  :: in1, in2, in3
    integer(HSIZE_T), dimension(3), intent(in), optional :: start ! Offset of array to read
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(3)      :: dim       ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: filespace ! dataspace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
     if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif  
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    if (present(in1) .and. present(in2) .and. present(in3)) then
      dim = [in1, in2, in3]
    else
      dim = shape(array3D)
    end if
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(3,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array3D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array3D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array3D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_reading_int

  !---------------------------------------- 
  ! HDF5 reading for a real double. Parallel compatible
  ! reading is enabled defining the argumets mpi_rank
  ! and n_mpi_tasks. 
  !----------------------------------------
  ! inputs:
  !   file_id:     (HID_T) file identifier
  !   dsetname:    (character)(*) name of the dataset from which data are read
  !   mpi_rank:    (integer)(optional) identifier of the MPI task
  !   n_mpi_tasks: (integer)(optional) number of MPI tasks
  ! outputs:
  !   rd:          (real8) character to be read from file
  !---------------------------------------- 
  subroutine HDF5_real_reading(file_id,rd,dsetname,mpi_rank,n_mpi_tasks)
    implicit none
    integer(HID_T)  , intent(in)          :: file_id  ! file identifier
    real*8          , intent(out)         :: rd
    character(LEN=*), intent(in)          :: dsetname ! dataset name
    integer,          intent(in),optional :: mpi_rank,n_mpi_tasks     
    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim       ! dataset dimensions
    integer(HSIZE_T),dimension(1) :: dims,maxdims ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: data_type
    integer(HID_T)      :: filespace
    logical             :: exists    ! true if dataset exists

    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    !*** get size of the dataspace
    call H5dget_space_f(dataset,dataspace,error)
    call H5Sget_simple_extent_ndims_f(dataspace,rank,error)
    if(rank.eq.1) call H5Sget_simple_extent_dims_f(dataspace,dims,maxdims,error)
    !*** read the double data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = int(1,kind=HSIZE_T)
    if(present(mpi_rank).and.present(n_mpi_tasks).and.&
    (rank.eq.1).and.(dims(1).gt.int(1,kind=HSIZE_T))) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=int([mpi_rank],kind=HSIZE_T),count=int([1],kind=HSIZE_T),hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error)
    endif

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_real_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 1D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:            (HID_T) file identifier
  !   dsetname:           (character)(*) name of the dataset in which the data are written
  !   start:              (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !                       in the global dataset
  ! outputs:
  !   array1D:            (real8)(:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array1D_reading(file_id,array1D,dsetname,start)
    implicit none
    integer(HID_T), intent(in)   :: file_id   ! file identifier
    real*8        , dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start ! Offset of array to read

    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim       ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: filespace ! dataspace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 1D (real*4). Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset in which the data are written
  !   start:    (HSIZE_T)(1)(optional) starting index of the input data chunk 
  !             in the global dataset
  ! outputs:
  !   array1D:  (real4)(:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array1D_reading_r4(file_id,array1D,dsetname,start)
    implicit none
    integer(HID_T), intent(in)   :: file_id   ! file identifier
    real*4        , dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start ! Offset of array to read
    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim       ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: filespace ! dataspace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=int(start,kind=HID_T), count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_REAL,array1D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_REAL,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading_r4

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 1D (real*4). Parallel application 
  ! requires the setting of the starting indexes "start" of the data chunk 
  ! within the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task and the size of the MPI task array in which the
  ! data chuck is stored.
  !----------------------------------------
  ! inputs:
  !   file_id:    (HID_T) file identifier
  !   dsetname:   (character)(*) name of the dataset in which the data are written
  !   reqdims_in: (HSIZE_T)(1) size of the data chunck to be loaded
  !   start:      (HSIZE_T)(1)(optional) starting index of 
  !               the input data chunk in the global dataset
  ! outputs:
  !   array1D:    (real4)(:)(allocatable) array of characters to be read from file
  !---------------------------------------- 
  subroutine HDF5_allocatable_array1D_reading_r4(file_id,array1D,dsetname,reqdims_in,start)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname
    integer(HSIZE_T), optional, intent(in) :: reqdims_in(1),start(1)
    real*4, allocatable, intent(inout) :: array1D(:)
    integer(HID_T) :: dataset_id,file_space,mem_space
    integer(HSIZE_T) :: counts(1)
    integer(int64) :: extents(1)
    integer :: error
    logical :: exists,ready

    call prepare_array_read(file_id,dsetname,1,storage_size(0.0_4)/8,reqdims_in,start, &
                            dataset_id,file_space,mem_space,counts,extents,exists,ready)
    ! Preserve the existing missing-dataset behavior, including the caller's array.
    if (.not.exists) return
    if (allocated(array1D)) deallocate(array1D)
    if (.not.ready) return
    allocate(array1D(extents(1)),stat=error)
    if (error==0) then
      ! An empty independent read needs no transfer; retain its allocated shape.
      if (all(extents>0)) then
        array1D=0.0_4
        call H5Dread_f(dataset_id,H5T_NATIVE_REAL,array1D,counts,error, &
                       file_space_id=file_space,mem_space_id=mem_space)
      endif
    endif
    if (error/=0) then
      write(*,*) 'WARNING: HDF5 array allocation/read failed: ',dsetname,' extents=',extents,' status=',error
      if (allocated(array1D)) deallocate(array1D)
    endif
    call close_array_read(dataset_id,file_space,mem_space)
  end subroutine HDF5_allocatable_array1D_reading_r4

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 1D. Parallel application 
  ! requires the setting of the starting indexes "start" of the data chunk 
  ! within the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task and the size of the MPI task array in which the
  ! data chuck is stored.
  !----------------------------------------
  ! inputs:
  !   file_id:    (HID_T) file identifier
  !   dsetname:   (character)(*) name of the dataset in which the data are written
  !   reqdims_in: (HSIZE_T)(1) size of the data chunck to be loaded
  !   start:      (HSIZE_T)(1)(optional) starting index of 
  !               the input data chunk in the global dataset
  ! outputs:
  !   array1D:    (real8)(:)(allocatable) array of characters to be read from file
  !----------------------------------------  
  subroutine HDF5_allocatable_array1D_reading(file_id,array1D,dsetname,reqdims_in,start)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname
    integer(HSIZE_T), optional, intent(in) :: reqdims_in(1),start(1)
    real*8, allocatable, intent(inout) :: array1D(:)
    integer(HID_T) :: dataset_id,file_space,mem_space
    integer(HSIZE_T) :: counts(1)
    integer(int64) :: extents(1)
    integer :: error
    logical :: exists,ready

    call prepare_array_read(file_id,dsetname,1,storage_size(0.0_8)/8,reqdims_in,start, &
                            dataset_id,file_space,mem_space,counts,extents,exists,ready)
    ! Preserve the existing missing-dataset behavior, including the caller's array.
    if (.not.exists) return
    if (allocated(array1D)) deallocate(array1D)
    if (.not.ready) return
    allocate(array1D(extents(1)),stat=error)
    if (error==0) then
      ! An empty independent read needs no transfer; retain its allocated shape.
      if (all(extents>0)) then
        array1D=0.0_8
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array1D,counts,error, &
                       file_space_id=file_space,mem_space_id=mem_space)
      endif
    endif
    if (error/=0) then
      write(*,*) 'WARNING: HDF5 array allocation/read failed: ',dsetname,' extents=',extents,' status=',error
      if (allocated(array1D)) deallocate(array1D)
    endif
    call close_array_read(dataset_id,file_space,mem_space)
  end subroutine HDF5_allocatable_array1D_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 2D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset in which the data are written
  !   start:    (HSIZE_T)(2)(optional) starting indexes of the input data chunk 
  !             in the global dataset
  ! outputs:
  !   array2D:  (real8)(:,:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array2D_reading(file_id,array2D,dsetname,start)
    implicit none
    integer(HID_T)  , intent(in)     :: file_id   ! file identifier
    real*8          , dimension(:,:) :: array2D
    character(LEN=*), intent(in)     :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(2), intent(in), optional :: start ! Offset of array to read

    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim       ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: filespace ! dataspace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
     if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif  
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim = shape(array2D)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(2,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_reading

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 2D. Parallel application 
  ! requires the setting of the starting indexes "start" of the data chunk 
  ! within the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task and the size of the MPI task array in which the
  ! data chuck is stored.
  !----------------------------------------
  ! inputs:
  !   file_id:    (HID_T) file identifier
  !   dsetname:   (character)(*) name of the dataset in which the data are written
  !   reqdims_in: (HSIZE_T)(2) size of the data chunck to be loaded
  !   start:      (HSIZE_T)(2)(optional) starting indexes of the input data chunk 
  !               in the global dataset
  ! outputs:
  !   array2D:    (real4)(:,:)(allocatable) array of characters to be read from file
  !---------------------------------------- 
  subroutine HDF5_allocatable_array2D_reading(file_id,array2D,dsetname,reqdims_in,start)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname
    integer(HSIZE_T), optional, intent(in) :: reqdims_in(2),start(2)
    real*8, allocatable, intent(inout) :: array2D(:,:)
    integer(HID_T) :: dataset_id,file_space,mem_space
    integer(HSIZE_T) :: counts(2)
    integer(int64) :: extents(2)
    integer :: error
    logical :: exists,ready

    call prepare_array_read(file_id,dsetname,2,storage_size(0.0_8)/8,reqdims_in,start, &
                            dataset_id,file_space,mem_space,counts,extents,exists,ready)
    ! Preserve the existing missing-dataset behavior, including the caller's array.
    if (.not.exists) return
    if (allocated(array2D)) deallocate(array2D)
    if (.not.ready) return
    allocate(array2D(extents(1),extents(2)),stat=error)
    if (error==0) then
      ! An empty independent read needs no transfer; retain its allocated shape.
      if (all(extents>0)) then
        array2D=0.0_8
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array2D,counts,error, &
                       file_space_id=file_space,mem_space_id=mem_space)
      endif
    endif
    if (error/=0) then
      write(*,*) 'WARNING: HDF5 array allocation/read failed: ',dsetname,' extents=',extents,' status=',error
      if (allocated(array2D)) deallocate(array2D)
    endif
    call close_array_read(dataset_id,file_space,mem_space)
  end subroutine HDF5_allocatable_array2D_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 3D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset in which the data are written
  !   start:    (HSIZE_T)(3)(optional) starting indexes of the input 
  !             data chunk in the global dataset
  ! outputs:
  !   array3D:  (real8)(:,:,:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array3D_reading(file_id,array3D,dsetname,&
       in1, in2, in3,start)
    implicit none
    integer(HID_T)     , intent(in)       :: file_id   ! file identifier
    real*8             , dimension(:,:,:) :: array3D
    character(LEN=*)   , intent(in)       :: dsetname  ! dataset name
    integer, intent(in), optional  :: in1, in2, in3
    integer(HSIZE_T), dimension(3), intent(in), optional :: start ! Offset of array to read
    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(3)      :: dim       ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: filespace ! dataspace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
     if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif  
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    if (present(in1) .and. present(in2) .and. present(in3)) then
      dim = [in1, in2, in3]
    else
      dim = shape(array3D)
    end if
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(3,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array3D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_reading

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 3D. Parallel application 
  ! requires the setting of the starting indexes "start" of the data chunk 
  ! within the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task and the size of the MPI task array in which the
  ! data chuck is stored.
  !----------------------------------------
  ! inputs:
  !   file_id:    (HID_T) file identifier
  !   dsetname:   (character)(*) name of the dataset in which the data are written
  !   reqdims_in: (HSIZE_T)(3) size of the data chunck to be loaded
  !   start:      (HSIZE_T)(3)(optional) starting indexes of the 
  !               input data chunk in the global dataset
  ! outputs:
  !   array3D:    (real8)(:,:,:)(allocatable) array of characters to be read from file
  !---------------------------------------- 
  subroutine HDF5_allocatable_array3D_reading(file_id,array3D,dsetname,reqdims_in,start)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname
    integer(HSIZE_T), optional, intent(in) :: reqdims_in(3),start(3)
    real*8, allocatable, intent(inout) :: array3D(:,:,:)
    integer(HID_T) :: dataset_id,file_space,mem_space
    integer(HSIZE_T) :: counts(3)
    integer(int64) :: extents(3)
    integer :: error
    logical :: exists,ready

    call prepare_array_read(file_id,dsetname,3,storage_size(0.0_8)/8,reqdims_in,start, &
                            dataset_id,file_space,mem_space,counts,extents,exists,ready)
    ! Preserve the existing missing-dataset behavior, including the caller's array.
    if (.not.exists) return
    if (allocated(array3D)) deallocate(array3D)
    if (.not.ready) return
    allocate(array3D(extents(1),extents(2),extents(3)),stat=error)
    if (error==0) then
      ! An empty independent read needs no transfer; retain its allocated shape.
      if (all(extents>0)) then
        array3D=0.0_8
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array3D,counts,error, &
                       file_space_id=file_space,mem_space_id=mem_space)
      endif
    endif
    if (error/=0) then
      write(*,*) 'WARNING: HDF5 array allocation/read failed: ',dsetname,' extents=',extents,' status=',error
      if (allocated(array3D)) deallocate(array3D)
    endif
    call close_array_read(dataset_id,file_space,mem_space)
  end subroutine HDF5_allocatable_array3D_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 4D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset in which the data are written
  !   start:    (HSIZE_T)(4)(optional) starting indexes of the input 
  !             data chunk in the global dataset
  ! outputs:
  !   array4D:  (real8)(:,:,:,:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array4D_reading(file_id,array4D,dsetname,ierr,start)
    implicit none
    integer(HID_T), intent(in)     :: file_id   ! file identifier
    real*8, dimension(:,:,:,:)     :: array4D
    character(LEN=*), intent(in)   :: dsetname  ! dataset name
    integer, optional, intent(out) :: ierr
    integer(HSIZE_T), dimension(4), intent(in), optional :: start ! Offset of array to read

    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
      dimension(4)      :: dim       ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: filespace ! dataspace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
     if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif  
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim = shape(array4D)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(4,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array4D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
    if (present(ierr)) ierr = error
  end subroutine HDF5_array4D_reading

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 4D. Parallel application 
  ! requires the setting of the starting indexes "start" of the data chunk 
  ! within the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task and the size of the MPI task array in which the
  ! data chuck is stored.
  !----------------------------------------
  ! inputs:
  !   file_id:    (HID_T) file identifier
  !   dsetname:   (character)(*) name of the dataset in which the data are written
  !   reqdims_in: (HSIZE_T)(4) size of the data chunck to be loaded
  !   start:      (HSIZE_T)(4)(optional) starting indexes of the 
  !               input data chunk in the global dataset
  ! outputs:
  !   array4D:    (real8)(:,:,:,:)(allocatable) array of characters to be read from file
  !---------------------------------------- 
  subroutine HDF5_allocatable_array4D_reading(file_id,array4D,dsetname,reqdims_in,start)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname
    integer(HSIZE_T), optional, intent(in) :: reqdims_in(4),start(4)
    real*8, allocatable, intent(inout) :: array4D(:,:,:,:)
    integer(HID_T) :: dataset_id,file_space,mem_space
    integer(HSIZE_T) :: counts(4)
    integer(int64) :: extents(4)
    integer :: error
    logical :: exists,ready

    call prepare_array_read(file_id,dsetname,4,storage_size(0.0_8)/8,reqdims_in,start, &
                            dataset_id,file_space,mem_space,counts,extents,exists,ready)
    ! Preserve the existing missing-dataset behavior, including the caller's array.
    if (.not.exists) return
    if (allocated(array4D)) deallocate(array4D)
    if (.not.ready) return
    allocate(array4D(extents(1),extents(2),extents(3),extents(4)),stat=error)
    if (error==0) then
      ! An empty independent read needs no transfer; retain its allocated shape.
      if (all(extents>0)) then
        array4D=0.0_8
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array4D,counts,error, &
                       file_space_id=file_space,mem_space_id=mem_space)
      endif
    endif
    if (error/=0) then
      write(*,*) 'WARNING: HDF5 array allocation/read failed: ',dsetname,' extents=',extents,' status=',error
      if (allocated(array4D)) deallocate(array4D)
    endif
    call close_array_read(dataset_id,file_space,mem_space)
  end subroutine HDF5_allocatable_array4D_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 5D. Parallel application requires
  ! the setting of the starting indexes "start" of the data chunk within 
  ! the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task.
  !----------------------------------------
  ! inputs:
  !   file_id:  (HID_T) file identifier
  !   dsetname: (character)(*) name of the dataset in which the data are written
  !   start:    (HSIZE_T)(5)(optional) starting indexes of 
  !             the input data chunk in the global dataset
  ! outputs:
  !   array5D:  (real8)(:,:,:,:,:) array of characters to be read from file
  !----------------------------------------
  subroutine HDF5_array5D_reading(file_id,array5D,dsetname,start)
    implicit none
    integer(HID_T), intent(in)   :: file_id   
    real*8, dimension(:,:,:,:,:) :: array5D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(4), intent(in), optional :: start ! Offset of array to read

    integer             :: error     ! error flag
    integer             :: rank      ! dataset rank
    integer(HSIZE_T), &
           dimension(5) :: dim       ! dataset dimensions
    integer(HID_T)      :: dataset   ! dataset identifier
    integer(HID_T)      :: dataspace ! dataspace identifier
    integer(HID_T)      :: filespace ! dataspace identifier
    logical             :: exists    ! true if dataset exists
    !*** check if dataset exists otherwise return ***
    call H5Lexists_f(file_id,trim(dsetname),exists,error)
    if(.not.exists) then
      write(*,*) "WARNING: dataset ",dsetname," does not exists!"
      return
    endif
    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim = shape(array5D)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(5,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array5D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array5D_reading

  !---------------------------------------- 
  ! HDF5 get dataset rank and dimensions
  !----------------------------------------
  ! inputs:
  !   file_id:    (HID_T) file identifier
  !   dsetname:   (character)(*) name of the dataset in which the data are written
  ! outputs:
  !   rank:    (integer) number of dimensions of the dataset
  !   dims:    (HSIZE_T)(rank)(allocatable) actual size of each dimensions
  !            of the dataset
  !   maxdims: (HSIZE_T)(rank)(allocatable) maximum size of each dimensions
  !            of the dataset
  !----------------------------------------
  subroutine HDF5_get_dataset_rank_dims(file_id,dsetname,rank,dims,maxdims)
    !> inputs:
    integer(HID_T),intent(in)   :: file_id
    character(len=*),intent(in) :: dsetname
    !> outputs:
    integer,intent(out)                                   :: rank
    integer(HSIZE_T),dimension(:),allocatable,intent(out) :: dims,maxdims
    !> variables:
    integer        :: error
    integer(HID_T) :: dataset_id,dataspace_id
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    if(error.ne.0) then
      write(*,*) "WARNING: error in opening the dataset ",dsetname,"!"
      return
    endif
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    !*** datasapce dimensions
    allocate(dims(rank)); allocate(maxdims(rank));
    call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine HDF5_get_dataset_rank_dims

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 5D. Parallel application 
  ! requires the setting of the starting indexes "start" of the data chunk 
  ! within the global HDF5 dataset corresponding to the data to be loaded by 
  ! each specific MPI task and the size of the MPI task array in which the
  ! data chuck is stored.
  !----------------------------------------
  ! inputs:
  !   file_id:    (HID_T) file identifier
  !   dsetname:   (character)(*) name of the dataset in which the data are written
  !   reqdims_in: (HSIZE_T)(5) size of the data chunck to be loaded
  !   start:      (HSIZE_T)(5)(optional) starting indexes of the 
  !               input data chunk in the global dataset
  ! outputs:
  !   array5D:    (real8)(:,:,:,:,:)(allocatable) array of characters to be read from file
  !----------------------------------------  
  subroutine HDF5_allocatable_array5D_reading(file_id,array5D,dsetname,reqdims_in,start)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dsetname
    integer(HSIZE_T), optional, intent(in) :: reqdims_in(5),start(5)
    real*8, allocatable, intent(inout) :: array5D(:,:,:,:,:)
    integer(HID_T) :: dataset_id,file_space,mem_space
    integer(HSIZE_T) :: counts(5)
    integer(int64) :: extents(5)
    integer :: error
    logical :: exists,ready

    call prepare_array_read(file_id,dsetname,5,storage_size(0.0_8)/8,reqdims_in,start, &
                            dataset_id,file_space,mem_space,counts,extents,exists,ready)
    ! Preserve the existing missing-dataset behavior, including the caller's array.
    if (.not.exists) return
    if (allocated(array5D)) deallocate(array5D)
    if (.not.ready) return
    allocate(array5D(extents(1),extents(2),extents(3),extents(4),extents(5)),stat=error)
    if (error==0) then
      ! An empty independent read needs no transfer; retain its allocated shape.
      if (all(extents>0)) then
        array5D=0.0_8
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array5D,counts,error, &
                       file_space_id=file_space,mem_space_id=mem_space)
      endif
    endif
    if (error/=0) then
      write(*,*) 'WARNING: HDF5 array allocation/read failed: ',dsetname,' extents=',extents,' status=',error
      if (allocated(array5D)) deallocate(array5D)
    endif
    call close_array_read(dataset_id,file_space,mem_space)
  end subroutine HDF5_allocatable_array5D_reading

  !----------------------------------------
  ! Get HDF5 property list. Create setand returns a HDF5 property list
  !----------------------------------------
  ! inputs:
  !   rank:
  !   chunksize: (HSIZE_T)(rank) size of the data chunk for
  !              each dimension of the dataset
  !   gzip:      (logical) if true, set a data compression level
  !              different the default one of 6
  !   level:     (integer)(optional) data compression level
  ! outputs:
  !   property:  (HID_T) indentifier of the HDF5 property list
  !----------------------------------------
  function get_HDF5_plist(rank, chunksize, gzip, level) result(property)
    implicit none
    integer, intent(in) :: rank
    integer(HSIZE_T), dimension(rank), intent(in) :: chunksize
    logical, intent(in) :: gzip
    integer, intent(in), optional :: level
    integer(HID_T) :: property

    integer        :: error, i
    integer(HSIZE_T), dimension(rank) :: chk
    integer, parameter :: cmpr_default = 6
    integer :: cmpr_level

    ! Check chunksize, if too large make it smaller.
    ! Aim for 10000 elements per chunk, by reducing the last dimension first
    chk = chunksize
    i   = rank
    do while (product(chk) .gt. 100000)
      if (chk(i) .gt. 1) then
        chk(i) = chk(i) / 2
      else
        i = i-1
      end if
    end do

    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)

    ! If we are not doing parallel IO we can enable compression
    if (gzip) then
      if (present(level)) then
        cmpr_level = level
      else
        cmpr_level = cmpr_default
      end if
      if (cmpr_level .gt. 0) then
        call H5Pset_chunk_f(property,rank,chk,error)
        call H5Pset_deflate_f(property,cmpr_level,error)
      end if
    end if
  end function get_HDF5_plist

#endif
end module hdf5_io_module
