! version 3.
! Optimazed v2.
module fem
    use solver
    implicit none

    integer*4, private :: fem_N
    integer*4, private :: fem_T
    real*8, private :: fem_h
    real*8, private :: fem_dt
    real*8, private :: fem_x_init
    real*8, private :: fem_y_init
    real*8, private :: fem_x_end
    real*8, private :: fem_y_end

    real*8, private :: fem_bx
    real*8, private :: fem_by
    real*8, private :: fem_dx
    real*8, private :: fem_dy

    ! Arrays for the matrix sparse row (CSR) format
    integer*4, allocatable :: fem_IA(:), fem_JA(:)
    real*8, allocatable :: fem_A(:)

    real*8, allocatable :: fem_C(:) 
    real*8, allocatable :: fem_C0(:)
    real*8, allocatable :: fem_RHS(:)

    real*8, allocatable :: fem_sp_bf_bf(:, :)
    real*8, allocatable :: fem_sp_gbf_bf(:, :)
    real*8, allocatable :: fem_sp_gbf_gbf(:, :)

    ! Parameters for integrate
    real*8, private :: quad_A(3, 9)
    real*8, private :: quad_w1, quad_w2

    contains

    subroutine generate_mesh(N, x_init, y_init, x_end, y_end)
        implicit none

        integer*4 :: N
        real*8 :: x_init, y_init, x_end, y_end

        fem_N = N
        fem_x_init = x_init
          fem_y_init = y_init
         fem_x_end = x_end
         fem_y_end = y_end

        fem_h = (x_end - x_init) / N

        return
    end subroutine

    subroutine init_fem(T, dt, bx, by, dx, dy)
        implicit none
        
        real*8 :: T
        real*8 :: dt
        real*8 :: bx
        real*8 :: by
        real*8 :: dx
        real*8 :: dy
        real*8 :: c0

        fem_T = T
        fem_dt = dt
        fem_bx = bx
        fem_by = by
        fem_dx = dx
        fem_dy = dy

        return
    end subroutine

    subroutine allocate_fem()
        implicit none

        allocate(fem_C((fem_N - 1)*(fem_N - 1)))
        fem_C = 0.0d0
        allocate(fem_C0((fem_N - 1)*(fem_N - 1)))
        fem_C0 = 0.0d0
        allocate(fem_RHS((fem_N - 1)*(fem_N - 1)))
        fem_RHS = 0.0d0

        allocate(fem_IA((fem_N - 1)*(fem_N - 1) + 1))
        allocate(fem_JA(10*(fem_N - 1)*(fem_N - 1)))
        allocate(fem_A(10*(fem_N - 1)*(fem_N - 1)))

        allocate(fem_sp_bf_bf((fem_N - 1)*(fem_N - 1), 9))
        allocate(fem_sp_gbf_bf((fem_N - 1)*(fem_N - 1), 9))
        allocate(fem_sp_gbf_gbf((fem_N - 1)*(fem_N - 1), 9))

        return
    end subroutine

    subroutine init_quad()
        implicit none

        real*8 :: a11, a12, a21, a22, a23
            
        a11 = 0.12495d0
        a12 = 0.43752d0
        a21 = 0.79711d0

        a22 = 0.16541d0
        a23 = 0.03748d0

        quad_w1 = 0.20595d0
        quad_w2 = 0.06369d0

        ! Init quad_A
        quad_A(1,1) = a11 
        quad_A(2,1) = a12
        quad_A(3,1) = a12

        quad_A(1,2) = a12 
        quad_A(2,2) = a11
        quad_A(3,2) = a12

        quad_A(1,3) = a12 
        quad_A(2,3) = a12
        quad_A(3,3) = a11

        quad_A(1,4) = a21 
        quad_A(2,4) = a22
        quad_A(3,4) = a23

        quad_A(1,5) = a21 
        quad_A(2,5) = a23
        quad_A(3,5) = a22

        quad_A(1,6) = a22 
        quad_A(2,6) = a21
        quad_A(3,6) = a23

        quad_A(1,7) = a22 
        quad_A(2,7) = a23
        quad_A(3,7) = a21

        quad_A(1,8) = a23 
        quad_A(2,8) = a21
        quad_A(3,8) = a22

        quad_A(1,9) = a23 
        quad_A(2,9) = a22
        quad_A(3,9) = a21

        return
    end subroutine

    ! Mesh
    real*8 function meshx(i)
        implicit none
        integer*4 :: i
        meshx = fem_x_init + fem_h*i
        return
    end function

    real*8 function meshy(j)
        implicit none
        integer*4 :: j
        meshy = fem_y_init + fem_h*j
        return
    end function
        
    ! Check point (x, y) in (i, j) cell
    logical function check_cell(x, y, i_cell, j_cell)
        implicit none
        
        integer*4 :: N
        real*8 :: x, y
        integer*4 :: i_cell, j_cell

        check_cell = .false.
        if (meshx(i_cell) <= x .and. x <= meshx(i_cell + 1)) then
            if (meshy(j_cell) <= y .and. y <= meshy(j_cell + 1)) then
                check_cell = .true.
            endif
        endif

        return
    end function

    ! Base functions
    real*8 function base_f(i, j, x, y)
        implicit none

        integer*4 :: i, j
        real*8 :: x, y

        base_f = 0.0d0

        ! Function in (i ,j) cell
        if (check_cell(x, y, i, j)) then
            base_f = (1.0 / (fem_h**2)) * (meshx(i + 1) - x) * (meshy(j + 1) - y)
            return
        endif

        ! Function in (i ,j - 1) cell
        if (check_cell(x, y, i, j - 1)) then
            base_f = (1.0 / (fem_h**2)) * (meshx(i + 1) - x) * (y - meshy(j - 1))
            return
        endif

        ! Function in (i - 1 ,j - 1) cell
        if (check_cell(x, y, i - 1, j - 1)) then
            base_f = (1.0 / (fem_h**2)) * (x - meshx(i - 1)) * (y - meshy(j - 1))
            return
        endif

        ! Function in (i - 1 ,j) cell
        if (check_cell(x, y, i - 1, j)) then
            base_f = (1.0 / (fem_h**2)) * (x - meshx(i - 1)) * (meshy(j + 1) - y)
            return
        endif

!		print *, i, j, x, y, base_f

        return
    end function

    ! Gradient of base functions and scalar product by vector b = (bx, by)
    real*8 function grad_base_f(i, j, x, y, bx, by)
        implicit none

        integer*4 :: i, j
        real*8 :: x, y
        real*8 :: bx, by

        grad_base_f = 0.0d0

        ! Function in (i ,j) cell
        if (check_cell(x, y, i, j)) then
            grad_base_f = (1.0 / (fem_h**2)) * ((y - meshy(j + 1))*bx + (x - meshx(i + 1))*by) 
            return
        endif

        ! Function in (i ,j - 1) cell
        if (check_cell(x, y, i, j - 1)) then
            grad_base_f = (1.0 / (fem_h**2)) * ((meshy(j - 1) - y)*bx + (meshx(i + 1) - x)*by) 
            return
        endif

        ! Function in (i - 1 ,j - 1) cell
        if (check_cell(x, y, i - 1, j - 1)) then
            grad_base_f = (1.0 / (fem_h**2)) * ((y - meshy(j - 1))*bx + (x - meshx(i - 1))*by) 
            return
        endif

        ! Function in (i - 1 ,j) cell
        if (check_cell(x, y, i - 1, j)) then
            grad_base_f = (1.0 / (fem_h**2)) * ((meshy(j + 1) - y)*bx + (meshx(i - 1) - x)*by) 
            return
        endif

        return
    end function


    ! Boundary condition
    real*8 function g_D(x, y, t)
        implicit none
        
        real*8 :: x, y, t
        real*8, parameter :: a = 10.0d0
        real*8, parameter :: c0 = 1.0d0

        g_D = 0.0d0
        if (x .eq. fem_x_init) then
            if (dabs(y) < a) then
                g_D = c0
            else
                g_D = 0.0d0
            endif
!		else
!			if (t > 0.0d0) then
!				g_D = analytic_solution(x, y, t, 256)		
!			endif
        endif

        return
    end function

    ! L2 scalar product in cell: (base function1, base function2)
    real*8 function sp_cell_bf_bf(i1, j1, i2, j2, i_cell, j_cell)
        implicit none
        
        integer*4 :: i1, j1, i2, j2
        integer*4 :: i_cell, j_cell
        
        integer*4 :: m
        real*8 :: S, T
        real*8 :: x1, x2, x3, y1, y2, y3, x, y

        ! Program
        S = 0.0d0

        ! Check boundary
        if (i_cell < 0 .or. j_cell < 0) then
            sp_cell_bf_bf = 0.0d0
            return
        endif

        if (i_cell >= fem_N .or. j_cell >= fem_N) then
            sp_cell_bf_bf = 0.0d0
            return
        endif

        ! Compute integral
        T = 0.5d0 * fem_h * fem_h												
    
        x1 = meshx(i_cell)
        y1 = meshy(j_cell)
        x3 = meshx(i_cell + 1)
        y3 = meshy(j_cell + 1)			

        ! Down-triangle
        x2 = meshx(i_cell + 1)
        y2 = meshy(j_cell)
        do m = 1, 3
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)

            S = S + T*quad_w1*( base_f(i1, j1, x, y) * base_f(i2, j2, x, y) )				
        enddo

        do m = 4, 9
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)
            
            S = S + T*quad_w2*( base_f(i1, j1, x, y) * base_f(i2, j2, x, y) )
        enddo

        ! Up-triangle
        x2 = meshx(i_cell)
        y2 = meshy(j_cell + 1)
        do m = 1, 3
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)

            S = S + T*quad_w1*( base_f(i1, j1, x, y) * base_f(i2, j2, x, y) )				
        enddo
        do m = 4, 9
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)

            S = S + T*quad_w2*( base_f(i1, j1, x, y) * base_f(i2, j2, x, y) )
        enddo
            
        sp_cell_bf_bf = S

        return
    end function

    ! L2 scalar product in cell: ( (vec(b), grad of base function1), base function2 )
    real*8 function sp_cell_gbf_bf(i1, j1, i2, j2, i_cell, j_cell, bx, by)
        implicit none
        
        integer*4 :: i1, j1, i2, j2
        integer*4 :: i_cell, j_cell
        real*8 :: bx, by
        
        integer*4 :: m
        real*8 :: S, T
        real*8 :: x1, x2, x3, y1, y2, y3, x, y

        ! Program
        S = 0.0d0

        ! Check boundary
        if (i_cell < 0 .or. j_cell < 0) then
            sp_cell_gbf_bf = 0.0d0
            return
        endif

        if (i_cell >= fem_N .or. j_cell >= fem_N) then
            sp_cell_gbf_bf = 0.0d0
            return
        endif

        ! Compute integral
        T = 0.5d0 * fem_h * fem_h												
    
        x1 = meshx(i_cell)
        y1 = meshy(j_cell)
        x3 = meshx(i_cell + 1)
        y3 = meshy(j_cell + 1)			

        ! Down-triangle
        x2 = meshx(i_cell + 1)
        y2 = meshy(j_cell)
        do m = 1, 3
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)

            S = S + T*quad_w1*( grad_base_f(i1, j1, x, y, bx, by) * base_f(i2, j2, x, y) )				
        enddo

        do m = 4, 9
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)
            
            S = S + T*quad_w2*( grad_base_f(i1, j1, x, y, bx, by) * base_f(i2, j2, x, y) )
        enddo

        ! Up-triangle
        x2 = meshx(i_cell)
        y2 = meshy(j_cell + 1)
        do m = 1, 3
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)

            S = S + T*quad_w1*( grad_base_f(i1, j1, x, y, bx, by) * base_f(i2, j2, x, y) )				
        enddo
        do m = 4, 9
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)

            S = S + T*quad_w2*( grad_base_f(i1, j1, x, y, bx, by) * base_f(i2, j2, x, y) )
        enddo
            
        sp_cell_gbf_bf = S

        return
    end function

    ! L2 scalar product in cell: ( diag(b) * grad of base function1, grad of base function2 )
    real*8 function sp_cell_gbf_gbf(i1, j1, i2, j2, i_cell, j_cell, bx, by)
        implicit none
        
        integer*4 :: i1, j1, i2, j2
        integer*4 :: i_cell, j_cell
        real*8 :: bx, by
        
        integer*4 :: m
        real*8 :: S, T
        real*8 :: S1, S2
        real*8 :: x1, x2, x3, y1, y2, y3, x, y

        ! Program
        S = 0.0d0

        ! Check boundary
        if (i_cell < 0 .or. j_cell < 0) then
            sp_cell_gbf_gbf = 0.0d0
            return
        endif

        if (i_cell >= fem_N .or. j_cell >= fem_N) then
            sp_cell_gbf_gbf = 0.0d0
            return
        endif

        ! Compute integral
        T = 0.5d0 * fem_h * fem_h												
    
        x1 = meshx(i_cell)
        y1 = meshy(j_cell)
        x3 = meshx(i_cell + 1)
        y3 = meshy(j_cell + 1)			

        ! Down-triangle
        x2 = meshx(i_cell + 1)
        y2 = meshy(j_cell)
        do m = 1, 3
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)

!			S = S + T*quad_w1*( grad_base_f(i1, j1, x, y, b1x, b1y) * grad_base_f(i2, j2, x, y, b2x, b2y) )							
            S1 = bx * grad_base_f(i1, j1, x, y, 1.0d0, 0.0d0) * grad_base_f(i2, j2, x, y, 1.0d0, 0.0d0)
            S2 = by * grad_base_f(i1, j1, x, y, 0.0d0, 1.0d0) * grad_base_f(i2, j2, x, y, 0.0d0, 1.0d0)
            S = S + T*quad_w1*(S1 + S2)								
        enddo

        do m = 4, 9
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)
            
!			S = S + T*quad_w2*( grad_base_f(i1, j1, x, y, b1x, b1y) * grad_base_f(i2, j2, x, y, b2x, b2y) )
            S1 = bx * grad_base_f(i1, j1, x, y, 1.0d0, 0.0d0) * grad_base_f(i2, j2, x, y, 1.0d0, 0.0d0)
            S2 = by * grad_base_f(i1, j1, x, y, 0.0d0, 1.0d0) * grad_base_f(i2, j2, x, y, 0.0d0, 1.0d0)
            S = S + T*quad_w2*(S1 + S2)				
        enddo

        ! Up-triangle
        x2 = meshx(i_cell)
        y2 = meshy(j_cell + 1)
        do m = 1, 3
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)

!			S = S + T*quad_w1*( grad_base_f(i1, j1, x, y, b1x, b1y) * grad_base_f(i2, j2, x, y, b2x, b2y) )				
            S1 = bx * grad_base_f(i1, j1, x, y, 1.0d0, 0.0d0) * grad_base_f(i2, j2, x, y, 1.0d0, 0.0d0)
            S2 = by * grad_base_f(i1, j1, x, y, 0.0d0, 1.0d0) * grad_base_f(i2, j2, x, y, 0.0d0, 1.0d0)
            S = S + T*quad_w1*(S1 + S2)				
        enddo
        do m = 4, 9
            x = x1*quad_A(1, m) + x2*quad_A(2, m) + x3*quad_A(3, m)					
            y = y1*quad_A(1, m) + y2*quad_A(2, m) + y3*quad_A(3, m)

!			S = S + T*quad_w2*( grad_base_f(i1, j1, x, y, b1x, b1y) * grad_base_f(i2, j2, x, y, b2x, b2y) )
            S1 = bx * grad_base_f(i1, j1, x, y, 1.0d0, 0.0d0) * grad_base_f(i2, j2, x, y, 1.0d0, 0.0d0)
            S2 = by * grad_base_f(i1, j1, x, y, 0.0d0, 1.0d0) * grad_base_f(i2, j2, x, y, 0.0d0, 1.0d0)
            S = S + T*quad_w2*(S1 + S2)				
        enddo
            
        sp_cell_gbf_gbf = S

        return
    end function

    ! L2 scalar product: (base function1, base function2)
    real*8 function sp_bf_bf(i1, j1, i2, j2)
        implicit none
        
        integer*4 :: i1, j1, i2, j2
        integer*4 :: i_cell, j_cell
        
        real*8 :: S1, S2, S3, S4

        S1 = 0.0d0
        S2 = 0.0d0
        S3 = 0.0d0
        S4 = 0.0d0

        S1 = sp_cell_bf_bf(i1, j1, i2, j2, i1, j1)
        S2 = sp_cell_bf_bf(i1, j1, i2, j2, i1, j1 - 1)
        S3 = sp_cell_bf_bf(i1, j1, i2, j2, i1 - 1, j1 - 1)
        S4 = sp_cell_bf_bf(i1, j1, i2, j2, i1 - 1, j1)

        sp_bf_bf = S1 + S2 + S3 + S4		

        return
    end function

    ! L2 scalar product: ( (vec(b), grad of base function1), base function2 )
    real*8 function sp_gbf_bf(i1, j1, i2, j2, bx, by)
        implicit none
        
        integer*4 :: i1, j1, i2, j2
        real*8 :: bx, by
        integer*4 :: i_cell, j_cell
        
        real*8 :: S1, S2, S3, S4

        S1 = 0.0d0
        S2 = 0.0d0
        S3 = 0.0d0
        S4 = 0.0d0

        S1 = sp_cell_gbf_bf(i1, j1, i2, j2, i1, j1, bx, by)
        S2 = sp_cell_gbf_bf(i1, j1, i2, j2, i1, j1 - 1, bx, by)
        S3 = sp_cell_gbf_bf(i1, j1, i2, j2, i1 - 1, j1 - 1, bx, by)
        S4 = sp_cell_gbf_bf(i1, j1, i2, j2, i1 - 1, j1, bx, by)

        sp_gbf_bf = S1 + S2 + S3 + S4		

        return
    end function

    ! L2 scalar product in cell: ( diag(b) * grad of base function1, grad of base function2 )
    real*8 function sp_gbf_gbf(i1, j1, i2, j2, bx, by)
        implicit none
        
        integer*4 :: i1, j1, i2, j2
        real*8 :: bx, by
        integer*4 :: i_cell, j_cell
        
        real*8 :: S1, S2, S3, S4

        S1 = 0.0d0
        S2 = 0.0d0
        S3 = 0.0d0
        S4 = 0.0d0

        S1 = sp_cell_gbf_gbf(i1, j1, i2, j2, i1, j1, bx, by)		
        S2 = sp_cell_gbf_gbf(i1, j1, i2, j2, i1, j1 - 1, bx, by)		
        S3 = sp_cell_gbf_gbf(i1, j1, i2, j2, i1 - 1, j1 - 1, bx, by)
        S4 = sp_cell_gbf_gbf(i1, j1, i2, j2, i1 - 1, j1, bx, by)

        sp_gbf_gbf = S1 + S2 + S3 + S4		

        return
    end function


    subroutine generate_sp()
        implicit none

        integer*4 :: i, j, k, l
        integer*4 :: indx, sp_indx

        do j = 1, fem_N-1
            do i = 1, fem_N-1		  
                   indx = (j-1)*(fem_N-1) + i
                do l = j - 1, j + 1
                    do k = i - 1, i + 1
                        sp_indx = (l - (j-1))*3 + k - (i-1) + 1
                        fem_sp_bf_bf(indx, sp_indx) = sp_bf_bf(k, l, i, j)
                        fem_sp_gbf_bf(indx, sp_indx) = sp_gbf_bf(k, l, i, j, fem_bx, fem_by)
                        fem_sp_gbf_gbf(indx, sp_indx) = sp_gbf_gbf(k, l, i, j, fem_dx, fem_dy)
                    enddo
                enddo
            enddo
        enddo

        return
    end subroutine

    subroutine generate_fem_matrix()
        implicit none

        integer*4 :: i, j, k, l
        integer*4 :: indx_col, indx_row
        integer*4 :: sp_indx
        integer*4 :: num_col
        real*8 :: S

        fem_IA(1) = 1
        do j = 1, fem_N-1
            do i = 1, fem_N-1		

                ! Generate a matrix row
                   indx_row = (j-1)*(fem_N-1) + i
                num_col = 0
                do l = max(j - 1, 1), min(j + 1, fem_N - 1)
                    do k = max(i - 1, 1), min(i + 1, fem_N - 1)

                        sp_indx = (l - (j-1))*3 + k - (i-1) + 1
                        S = (1.0d0/fem_dt) * fem_sp_bf_bf(indx_row, sp_indx)
                        S = S + 0.5*fem_sp_gbf_bf(indx_row, sp_indx)
                        S = S + 0.5*fem_sp_gbf_gbf(indx_row, sp_indx)
                        
                        if (S .ne. 0.0d0) then
                              indx_col = (l-1)*(fem_N-1) + k
                            fem_JA(fem_IA(indx_row) + num_col) = indx_col
                            fem_A(fem_IA(indx_row) + num_col) = S
                            num_col = num_col + 1 
                        endif 

                    enddo
                enddo
                fem_IA(indx_row + 1) = fem_IA(indx_row) + num_col

            enddo
        enddo

        return
    end subroutine

    subroutine show_fem_matrix()
        implicit none

!		print *, fem_N		
!		print *, fem_IA
!		print *, fem_JA
        call draw_matrix((fem_N-1)*(fem_N-1), fem_IA, fem_JA, "matrix.ps")

        return
    end subroutine

    subroutine generate_fem_rhs(i_T)
        implicit none

        real*8 :: i_T
        integer*4 :: i, j, indx, indx_in
        integer*4 :: sp_indx
        integer*4 :: k, l
        real*8 :: S, S1, S2
    
        do j = 1, fem_N-1
            do i = 1, fem_N-1		  

                    indx = (j-1)*(fem_N-1) + i
                fem_RHS(indx) = 0.0d0
                do l = j - 1, j + 1
                    do k = i - 1, i + 1

                        sp_indx = (l - (j-1))*3 + k - (i-1) + 1
                        indx_in = (l-1)*(fem_N-1) + k
                        ! Check (k, l) in boundary or not
                        if ( (l == 0 .or. l == fem_N) .or. (k == 0 .or. k == fem_N) ) then
                            S1 = -(1.0d0/fem_dt) * fem_sp_bf_bf(indx, sp_indx)
                            S1 = S1 - 0.5d0*fem_sp_gbf_bf(indx, sp_indx)
                            S1 = S1 - 0.5d0*fem_sp_gbf_gbf(indx, sp_indx)
                            S1 = g_D(meshx(k), meshy(l), i_T) * S1
                            
                            S2 = (1.0d0/fem_dt) * fem_sp_bf_bf(indx, sp_indx)
                            S2 = S2 - 0.5d0*fem_sp_gbf_bf(indx, sp_indx)
                            S2 = S2 - 0.5d0*fem_sp_gbf_gbf(indx, sp_indx)
                            S2 = g_D(meshx(k), meshy(l), i_T) * S2

                            S = S1 + S2
                        else
                            S = (1.0d0/fem_dt) * fem_sp_bf_bf(indx, sp_indx)
                            S = S - 0.5d0*fem_sp_gbf_bf(indx, sp_indx)
                            S = S - 0.5d0*fem_sp_gbf_gbf(indx, sp_indx)
                            S = fem_C0(indx_in) * S
                        endif
                        fem_RHS(indx) = fem_RHS(indx) + S

                    enddo
                enddo

            enddo
        enddo					

        return
    end subroutine

    subroutine fem_time_integrate()
        implicit none

        integer*4 :: i, j		
        real*8 :: i_T
        real*8 :: Db, Peh

!---------Pekle number-----------------------------------------------------------------------------
        Db = (fem_dx*fem_bx**2 + fem_dy*fem_by**2) / (fem_bx**2 + fem_by**2)
        Peh = (fem_h * sqrt(fem_bx**2 + fem_by**2)) / Db
        print *, 'Pekle number is ', Peh

        call generate_sp()
        print *, 'SPs were generated'

        call generate_fem_matrix()
        print *, size(fem_IA), size(fem_JA), size(fem_A)

!---------Init precondition for matrix-------------------------------------------------------------
        call set_precond((fem_N-1)*(fem_N-1), fem_IA, fem_JA, fem_A)

!---------Output file------------------------------------------------------------------------------
        open(12, FILE='result')
        write(12, *) fem_N-1
        write(12, *) int(fem_T/fem_dt) + 1

!---------Local Output file------------------------------------------------------------------------
        open(13, FILE='vresult')
        write(13, *) fem_N-1
        write(13, *) int(fem_T/fem_dt) + 1

        i_T = 0.0d0
        do while (i_T <= fem_T)

            call generate_fem_rhs(i_T)

            call solve_gmres((fem_N-1)*(fem_N-1), fem_IA, fem_JA, fem_A, fem_RHS, fem_C)
            fem_C0 = fem_C		

!--------------Print local solve
!			do j = 1, fem_N-1
!				do i = 1, fem_N-1
!					write(13, *) fem_C((j-1)*(fem_N-1) + i)
!				enddo
!			enddo

            i_T = i_T + fem_dt
        enddo

!---------Print solve
        do j = 1, fem_N-1
            do i = 1, fem_N-1
                write(12, *) fem_C((j-1)*(fem_N-1) + i)
            enddo	
        enddo
    
        return
    end subroutine

    real*8 function integral_func(x, y, tau)
        implicit none

        real*8 :: x, y, tau, a
        real*8 :: I

        a = 10.0d0
        
        I = tau**(-1.5d0) * (erf( (a + y) / ((4.0d0 * fem_dy * tau)**0.5d0)) +   &
                             erf( (a - y) / ((4.0d0 * fem_dy * tau)**0.5d0)))

        I = I * exp(-((x - tau) / ((4.0 * fem_dx * tau)**0.5d0))**2)

        integral_func = I
        return

    end function

    real*8 function analytic_solution(x, y, t, n_I)
        implicit none

        real*8 :: x, y, t
        integer :: n_I, j
        real*8 :: S, zeta
        real*8 :: pi

        pi = 4.0d0 * datan(1.0d0)
        
        S = 0.0d0
        do j = 1, n_I
            zeta = cos((2*j - 1)*pi / (2.0*n_I))

            S = S + ((1.0 - zeta**2)**0.5d0) * integral_func(x, y, t*0.5d0*(zeta + 1))
        enddo

        S = S * ((pi*t) / (2.0*n_I))
        S = S * (x / ((16.0d0 * pi * fem_dx)**0.5d0))
        
        analytic_solution = S
        return

    end function

    subroutine print_analytic_solution(t, n_I)
        implicit none

        real*8 :: t
        integer :: n_I, i, j		


        open(14, FILE='aresult')
        write(14, *) fem_N-1
        write(14, *) int(t)

!---------Print analytic solve
        do j = 1, fem_N-1
            do i = 1, fem_N-1
                write(14, *) analytic_solution(meshx(i), meshy(j), t, n_I)
            enddo	
        enddo

        open(15, FILE='erresult')
        write(15, *) fem_N-1
        write(15, *) int(t)

!---------Print error 
        do j = 1, fem_N-1
            do i = 1, fem_N-1
                write(15, *) analytic_solution(meshx(i), meshy(j), t, n_I) - fem_C((j-1)*(fem_N-1) + i)
            enddo	
        enddo

        return
    end subroutine	


end module

