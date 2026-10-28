program gk_neutral_reflection_test
use mpi
use fruit
use fruit_mpi, only: fruit_summary_mpi
use data_structure, only: init_threads
use mod_gk_neutral_reflection_test
implicit none
integer :: error,rank,tasks
logical :: ok,all_ok
character(len=40) :: mode
call MPI_INIT(error)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,error)
call MPI_COMM_SIZE(MPI_COMM_WORLD,tasks,error)
call init_threads()
call get_command_argument(1,mode)
if(len_trim(mode)>0) call run_reflection_invalid(trim(mode))
call init_fruit
call fruit_hide_dots
call run_reflection_tests()
call fruit_summary_mpi(tasks,rank)
call is_all_successful(ok)
call MPI_ALLREDUCE(ok,all_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,error)
if(.not.all_ok) call MPI_ABORT(MPI_COMM_WORLD,1,error)
call MPI_FINALIZE(error)
end program
