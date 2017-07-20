program main
	use fem

	implicit none

	integer*4, parameter :: N = 256
	real*8, parameter :: time = 50.0d0
	real*8, parameter :: dtime = 1.0d0
	
	call generate_mesh(N, 0.0d0, -100.0d0, 200.0d0, 100.0d0)
	call init_fem(time, dtime, 1.0d0, 0.0d0, 1.0d-3, 1.0d-4)
	call allocate_fem()
	call init_quad()

	call fem_time_integrate()

	call print_analytic_solution(time, 512)

end program
