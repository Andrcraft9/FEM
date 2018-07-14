program main
    use fem

    implicit none

    integer*4, parameter :: N = 300
    real*8, parameter :: time = 50.0d0
    real*8, parameter :: dtime = 1.0d0
    
    call generate_mesh(N, 0.0d0, -100.0d0, 200.0d0, 100.0d0)
    call init_fem(time, dtime, 1.0d0, 0.0d0, 1.0d0, 0.3d0)
    call allocate_fem()
    call init_quad()

!   call generate_fem_matrix()
!   call show_fem_matrix()
 
    call fem_time_integrate()

!   print *, sp_cell_gbf_bf(1, 0, 2, 0, 1, 0, 1.0d0, 1.0d0)
!   print *, sp_gbf_bf(1, 0, 2, 0, 1.0d0, 1.0d0)
    
end program
