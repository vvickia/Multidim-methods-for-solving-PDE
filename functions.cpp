module func

implicit none


real(8), parameter :: eps = 1e-8
integer, parameter :: r = 100
real(8), parameter :: D = 1.11 * 1e-4


contains

	subroutine Initialization(t_0, t_b, a, b) !Task parameter initialization.
	
	real(8) t_0, t_b, a, b
	
	
	open(1, file = 'Init_file.dat') !Reading model parameters from a file

		read(1, *) a, b
		read(1, *) t_0
		read(1, *) t_b
	
	close(1)
	
	end subroutine Initialization
	
	
	
	
	subroutine Thomas(a, b, res)

	real(8), dimension(0:, :) :: a
	real(8), dimension(0:) :: b, res
	real(8), allocatable :: coef(:,:)
	real(8) eps
	integer i, lenght, ier
	
	lenght = size(b)
	
	! Memory allocation for array

	allocate(coef(0:lenght - 2, 0:1), stat = ier); if(ier /= 0) stop 4

	! X_N-2 and Y_N-2 calculation

	coef(lenght - 2, 0) = -a(lenght - 1, 1)/a(lenght - 1, 2)
	coef(lenght - 2, 1) = b(lenght - 1)/a(lenght - 1, 2)

	! Straight sweep

	do i = lenght - 3, 0, -1

		if (abs(a(i + 1, 3) * coef(i + 1, 0) + a(i + 1, 2)) < eps) then
		
			write(*,*) "Have not any results, or infinitely many results"; stop 10
			
		endif
		
		coef(i,0) = -a(i + 1, 1)/(a(i + 1, 3) * coef(i + 1, 0) + a(i + 1, 2))
		coef(i,1) = (b(i + 1) - (a(i + 1, 3) * coef(i + 1, 1)))/(a(i + 1, 3) * coef(i + 1, 0) + a(i + 1, 2))

	enddo
	
	! U_0 calculation
	if (abs(a(0,2) + a(0,3) * coef(0,0)) < eps) then
		
			write(*,*) "Have not any results, or infinitely many results"; stop 20
			
	endif
	res(0) = (b(0) - a(0,3) * coef(0,1))/ (a(0,2) + a(0,3) * coef(0,0))
	
	

	! Reverse sweep

	do i = 1, lenght - 1

		res(i) = coef(i - 1, 0) * res(i - 1) + coef(i - 1, 1)

	enddo

	deallocate(coef, stat = ier); if(ier /= 0) stop 5

	end subroutine Thomas
	
	
	


	subroutine Set_IC(u, t_0, t_b)
	
	real(8), dimension(:,:,:) :: u
	real(8) t_0, t_b
	integer i, j
	
	do i = 1, r + 2 
	
		u(i,1,1) = t_b; u(i,1,2) = t_b; u(i,1,3) = t_b
		u(i,r + 2,1) = t_b; u(i,r + 2,2) = t_b; u(i,r + 2,3) = t_b !Border conditions
		
		do j = 2, r + 1
		
			u(i,j,1) = t_0; u(i,j,2) = t_0; u(i,j,3) = t_0 !Initial
		
		enddo
	
	enddo
	
	
	
	end subroutine Set_IC
	
	
	
	
	subroutine Step(u, dt, dx, dy)

	
	
	real(8), allocatable :: a(:,:), b(:)
	real(8), dimension(:,:,:) :: u
	real(8) dt, dx, dy
	integer i, j, ier
	
	
	allocate(a(0:r-1,3), b(0:r-1), stat = ier); if(ier /= 0) stop 16
	
	!Stage 1
	
	do i = 1, r - 2
	
		a(i,1) = -D/dy**2; a(i,3) = -D/dy**2
		a(i,2) = 2*D/dy**2 + 2/dt
	
	enddo
	
	a(0,1) = 0; a(0,2) = -1; a(0,3) = 1
	a(r-1,1) = -1; a(r-1,2) = 1; a(r-1,3) = 0
	
	!a(0,1) = 0; a(0,2) = 2*D/dy**2 + 2/dt; a(0,3) = -D/dy**2
	!a(r-1,1) = -D/dy**2; a(r-1,2) = 2*D/dy**2 + 2/dt; a(r-1,3) = 0
	
	
	do j = 2, r + 1
	
		do i = 2, r + 1
		
			b(i - 2) = D/dx**2 * (u(i+1,j,1) - 2*u(i,j,1) + u(i-1,j,1)) + 2/dt*u(i,j,1)
			
		enddo
		
		b(0) = 0
		b(r-1) = 0
		
		call Thomas(a, b, u(2:r+1,j,2))
	
	enddo
	
	call Set_BC(u(:,:,2)) !Refresh BC
	
	
	! Stage 2
	
	do i = 1, r - 2
	
		a(i,1) = -D/dx**2; a(i,3) = -D/dx**2
		a(i,2) = 2*D/dx**2 + 2/dt
	
	enddo
	
	a(0,1) = 0; a(0,2) = 1; a(0,3) = 0
	a(r-1,1) = 0; a(r-1,2) = 1; a(r-1,3) = 0
	
	!a(0,1) = 0; a(0,2) = 2*D/dx**2 + 2/dt; a(0,3) = -D/dx**2
	!a(r-1,1) = -D/dx**2; a(r-1,2) = 2*D/dx**2 + 2/dt; a(r-1,3) = 0
	
	do i = 2, r + 1
	
		do j = 2, r + 1
		
			b(j - 2) = D/dy**2 * (u(i,j+1,2) - 2*u(i,j,2) + u(i,j-1,2)) + 2/dt*u(i,j,2)
			
		enddo
		
		b(0) = 25.0
		b(r-1) = 25.0
		
		call Thomas(a, b, u(i,2:r+1,3))
	
	enddo
	
	call Set_BC(u(:,:,3)) !Refresh BC
	
	
	
	deallocate(a, b, stat = ier); if(ier /= 0) stop 23
	
	end subroutine Step
	
	
	
	
	
	
	subroutine Set_BC(u)
	
	
	real(8), dimension(:,:) :: u
	
	!u(1,:) = 0; u(r+2,:) = 0
	
	u(1,:) = u(2,:); u(r + 2, :) = u(r + 1, :)
	
	end subroutine Set_BC	
	
	
	
	
	subroutine Update_IC(u)
	
	
	real(8), dimension(:,:,:) :: u
	
	u(:,:,1) = u(:,:,3)
	
	end subroutine Update_IC
	
	
	
	subroutine Save_Data(u, iter)
	
	!Saving calculation results to a file.
	!Store the 2D matrix of results for each iteration step. 
	
	real(8), dimension(:,:,:) :: u
	integer iter, i, j
	character(3) str
	
	if (iter < 10) then
		write(str,"(I1.1)") iter
	elseif (iter < 100) then
		write(str,"(I2.2)") iter
	else
		write(str,"(I3.3)") iter
	endif
	
	open(1, file = 'result'//trim(str)//'.dat')
	
		
		do i = 2, r + 1
		
			write(1, *) u(i,:,3)
			
		enddo
	
	close(1)
	
	
	end subroutine Save_Data
	
	
	



end module func
