program main

use func

implicit none

real(8) a, b
real(8) t_0, t_b, dt, dx, dy
real(8), allocatable :: u(:,:,:)
integer i, n, ier



call Initialization(t_0, t_b, a, b) !Task parameter initialization

dx = a/r
dy = b/r


allocate(u(r+2, r+2, 3), stat = ier); if (ier /= 0) stop 3

call Set_IC(u, t_0, t_b) !Setting initial conditions


dt = 2 !Time step definition


do i = 1, 999

	call Step(u, dt, dx, dy) !Calculation of the desired function at a new time step
	
	call Update_IC(u) !Updating Initial Conditions
	
	call Save_Data(u, i) !Saving calculation results to a file
		
enddo



deallocate(u, stat = ier); if (ier /= 0) stop 10



end program main

