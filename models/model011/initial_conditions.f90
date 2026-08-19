!> Required model hook for electrostatic gyrokinetic startup.
subroutine initial_conditions(my_id, node_list, element_list, bnd_node_list, &
                              bnd_elm_list, xpoint2, xcase2)
  use data_structure, only: type_node_list, type_element_list, &
                            type_bnd_node_list, type_bnd_element_list

  implicit none

  integer, intent(in)                      :: my_id, xcase2
  logical, intent(in)                      :: xpoint2
  type(type_node_list), intent(inout)      :: node_list
  type(type_element_list), intent(inout)   :: element_list
  type(type_bnd_node_list), intent(inout)  :: bnd_node_list
  type(type_bnd_element_list), intent(inout) :: bnd_elm_list

  ! The model011 fields are imported by read_jorek_fields_interp_linear.
  ! There is no independent MHD initial-condition model to apply here, and
  ! changing node values would destroy the imported equilibrium/background.
end subroutine initial_conditions
