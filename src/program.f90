program sph

	use gnufor2
	use module1

	implicit none
	integer, parameter	:: N1=50
	real(kind=8)		:: x1(N1), f1(N1)
	integer :: i
! generate data for 2D plots
	do  i=1,N1
		x1(i)=5.0*i/N1
	end do
	f1=sin(2*x1)


	print *,'********************************************************************'
	print *,'Example 1: simple 2D graph'
	print *,'call plot(x1,f1)'
	call plot(x1,f1)
	print *,'press ENTER to go to the next example'
	read  *


end program sph
