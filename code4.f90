! This code is free & open-source, 
! distributed under the terms of GNU General Public 
! License (Version 3, 29 June 2007), 
! https://www.gnu.org/licenses/gpl-3.0.txt .
!   
!   Programmer: Ashkan Shekaari
! 
!      Contact: shekaari.theory@gmail.com
!               shekaari@email.kntu.ac.ir
! 
! Social Links: https://rt.academia.edu/AshkanShekaari
!               https://orcid.org/0000-0002-7434-467X
!  
! >>> PLEASE, cite the main paper in case of using this code, as follows:
!     A. Shekaari, M. Jafari (2021). "Theory and simulation of the Ising model". arXiv:2105.00841v1 [cond-mat.stat-mech]. doi:10.48550/arXiv.2105.00841.
! 
! PURPOSE: Ising model in 3 dimensions.
! 
module shared_data
implicit none
save
integer i, j, k ! Loop counters
integer x, y, z ! Spins' position coordinates
integer up, down, right, left, above, below
integer, parameter :: L = 4	        ! Size of the square lattice 
integer, parameter :: D = 3         ! Lattice dimension
integer, parameter :: N = L**D      ! Number of spins
integer, parameter :: NN = 6        ! Number of nearest-neighbors
integer, parameter :: mcstp = 100   ! Number of Monte Carlo sweeps 
integer, parameter :: eqstp = 100   ! Number of Monte Carlo sweeps for equilibration
real, parameter :: T_0 = 0.5        ! Initial temperature
real, parameter :: T_f = 10.0       ! Final temperature
real, parameter :: norm = 1./(mcstp*N*NN)
real, parameter :: nr2 = norm/mcstp
real, parameter :: dT = 0.1         ! Temperature step
end module
!--------------------------------
subroutine METROPOLIS(spin, beta)
! The Metropolis algorithm
use shared_data
implicit none
integer, intent(inout) :: spin(L, L, L)	! Spins
real, intent(in) :: beta	! Beta = 1/kB*T
integer s	! Box for spin
real R
real dE		! Energy difference
real x0, y0, z0, rand	! Random numbers
do i = 1, L 
	do j = 1, L 
		do k = 1, L
			call random_number(x0)
			call random_number(y0)
			call random_number(z0)
			x = int(1 + (L - 1)*x0)		! X coordinate
			y = int(1 + (L - 1)*y0)		! Y coordinate
			z = int(1 + (L - 1)*z0)		! Z coordinate
			s = spin(x, y, z)	! Spin
			! Periodic boundary conditions
			if (x == 1) then
				left = spin(L, y, z)
				right = spin(2, y, z)
			else if (x == L) then
				left = spin(L - 1, y, z)
				right = spin(1, y, z)
			else
				left = spin(x - 1, y, z)
				right = spin(x + 1, y, z)
			end if
			if (y == 1) then
				up = spin(x, 2, z)
				down = spin(x, L, z)
			else if (y == L) then
				up = spin(x, 1, z)
				down = spin(x, L - 1, z)
			else
				up = spin(x, y + 1, z)
				down = spin(x, y - 1, z)
			end if
			if (z == 1) then
				above = spin(x, y, 2)
				below = spin(x, y, L)
			else if (z == L) then
				above = spin(x, y, 1)
				below = spin(x, y, L - 1)
			else
				above = spin(x, y, z + 1)
				below = spin(x, y, z - 1)
			end if 		! End of periodic boundary conditions
			R = up + down + right + left + above + below
			dE = 2*s*R		! Energy difference
			call random_number(rand)	!Throws the die
			if (dE < 0.) then
				s = -s		! Flips the spin
			elseif (rand < exp(-dE*beta)) then
				s = -s		! Flips the spin
			endif
			spin(x, y, z) = s	! Returns the spin without flipping
		enddo
	enddo
enddo
end subroutine METROPOLIS
!-----------------------------------
real function CALCULATE_ENERGY(spin)
use shared_data
implicit none
integer, intent(in) :: spin(L, L, L)
real En, R
En = 0. ! Energy initialization
do x = 1, L
	do y = 1, L 
		do z = 1, L
			! Periodic boundary conditions
			if (x == 1) then
				left = spin(l, y, z)
				right = spin(2, y, z)
			else if (x == L) then
				left = spin(l - 1, y, z)
				right = spin(1, y, z)
			else
				left = spin(x - 1, y, z)
				right = spin(x + 1, y, z)
			end if
			if (y == 1) then
				up = spin(x, 2, z)
				down = spin(x, L, z)
			else if (y == L) then
				up = spin(x, 1, z)
				down = spin(x, L - 1, z)
			else
				up = spin(x, y + 1, z)
				down = spin(x, y - 1, z)
			end if
			if (z == 1) then
				above = spin(x, y, 2)
				below = spin(x, y, L)
			else if (z == L) then
				above = spin(x, y, 1)
				below = spin(x, y, L - 1)
			else
				above = spin(x, y, z + 1)
				below = spin(x, y, z - 1)
			end if
			R = up + down + right + left + above + below
			En = En - spin(x, y, z)*R	! Energy accumulation
		enddo
	enddo
enddo
CALCULATE_ENERGY = En	! Energy
end function
!------------------------------------------
real function CALCULATE_MAGNETIZATION(spin)
use shared_data
implicit none
integer, intent(in) :: spin(L, L, L)	! Spin
CALCULATE_MAGNETIZATION = sum(spin)	! Magnetization
end function
!------------------------------
program ising	! Main program
use shared_data
implicit none
integer sweep	! Loop counter
real T, B, B2	! Temperature, beta = 1/k_{B}T, beta^2
real summ0, summ1, summ2, summ3, summ4	! Accumulators
real E, M	! Energy, magnetization
real CALCULATE_ENERGY, CALCULATE_MAGNETIZATION	! external functions
real mean_E, mean_M, C_v, chi	! Mean E, mean M, specific heat, susceptibility
integer :: spin(L, L, L) = 1 ! Initial configuration: all spins up
integer c	! Counter
open(1, file = 'output.txt')	! Output file
!--------------------------------------------
T = T_0		! Temperature initialization
c = 0
10 format (6f12.6)
write(1, *)"#   T^{-1}         T           E          M          C_v         Chi"
write(1, *)"#  ========    ========    ========    ========    ========    ========"
do while (T <= T_f)		! Temperature loop (main loop)
	T = T + dT
	c = c + 1
	! Initializing accumulators ------------------
	summ0 = 0.; summ1 = 0.; summ2 = 0.; summ3 = 0.
	!---------------------------------
	B = 1./T; B2 = B*B	! Beta, beta^2
	do sweep = 1, eqstp	! Equilibration loop
		call METROPOLIS(spin, B)
	enddo
	!---------------------------------
	do sweep = 1, mcstp ! Main MC loop 
		call METROPOLIS(spin, B)
		E = CALCULATE_ENERGY(spin)
		M = CALCULATE_MAGNETIZATION(spin)
		summ0 = summ0 + E      ! E accumulator
		summ1 = summ1 + E*E    ! E^2 accumulator
		summ2 = summ2 + M      ! M accumulator
		summ3 = summ3 + M*M    ! M^2 accumulator
	enddo
	!--------------------------------------------------
	! Thermodynamic averages at each temperature value
	mean_E = summ0*norm		! Mean energy
	mean_M = summ2*norm		! Mean magnetization 
	C_v = (norm*summ1 - nr2*summ0*summ0)*B2		! Specific heat
	chi = (norm*summ3 - nr2*summ2*summ2)*B		! Susceptibility
	write(1, 10)B, T, mean_E, mean_M, C_v, chi
	write(1, 10)B, T, mean_E, mean_M, C_v, chi
	write(*, *)" >>> Temperature step", c
enddo
end program
