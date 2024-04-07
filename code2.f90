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
! PURPOSE: Ising model in 2 dimensions.
! 
module SHARED_DATA
implicit none
save
integer i, j	! Loop counters
integer x, y	! Position coordinates
integer up, down, right, left
integer, parameter :: L = 5         ! Size of the quare lattice
integer, parameter :: D = 2         ! Lattice dimension
integer, parameter :: N = L**D      ! Number of spins
integer, parameter :: mcstp = 100   ! Monte Carlo sweeps
integer, parameter :: eqstp = 100   ! Monte Carlo sweeps for equilibration
integer, parameter :: NN = 4        ! Number of nearest-neighbors
real, parameter :: T_0 = 0.5        ! Initial temperature
real, parameter :: T_f = 10.0       ! Final temperature
real, parameter :: norm = 1./(mcstp*N*NN)
real, parameter :: nr2 = norm/mcstp
real, parameter :: dT = 0.1         ! Temperature step
end module
!-------------------------------
subroutine METROPOLIS(sp, beta)
! The Metropolis algorithm"""
use SHARED_DATA
implicit none
integer, intent(inout) :: sp(L, L)	! Spin
real, intent(in) :: beta	! = 1/k_{B}T
integer s	! Box for spin
real R
integer dE	! Energy difference
real rand, x0, y0	! Random numbers
do i = 1, L 
	do j = 1, L 
		call random_number(x0)
		call random_number(y0)
		x = int(1 + (L - 1)*x0)
		y = int(1 + (L - 1)*y0)
		s = sp(x, y)
		! Periodic boundary conditions
		if (x == 1) then
			left = sp(L, y)
			right = sp(2, y)
		else if (x == L) then
			left = sp(L - 1, y)
			right = sp(1, y)
		else
			left = sp(x - 1, y)
			right = sp(x + 1, y)
		end if
		if (y == 1) then
			up = sp(x, 2)
			down = sp(x, L)
		else if (y == L) then
			up = sp(x, 1)
			down = sp(x, L - 1)
		else
			up = sp(x, y + 1)
			down = sp(x, y - 1)
		end if ! End of periodic boundary conditions
		R = up + down + right + left
		dE = 2*s*R
		call random_number(rand)
		if (dE < 0.) then
			s = -s  ! Flips the spin
		else if (rand < exp(-dE*beta)) then	! Throws the die
			s = -s  ! Flips the spin
		end if
		sp(x, y) = s  ! Returns the spin without flipping
	enddo
enddo
end subroutine METROPOLIS
!-----------------------------------
real function CALCULATE_ENERGY(sp1)
use SHARED_DATA
implicit none
integer, intent(in) :: sp1(L, L)	! Spin
real En, R
En = 0.	! Energy initialization
do x = 1, L
	do y = 1, L 
		if (x == 1) then
			left = sp1(L, y)
			right = sp1(2, y)
		else if (x == L) then
			left = sp1(L - 1, y)
			right = sp1(1, y)
		else
			left = sp1(x - 1, y)
			right = sp1(x + 1, y)
		end if
		if (y == 1) then
			up = sp1(x, 2)
			down = sp1(x, L)
		else if (y == L) then
			up = sp1(x, 1)
			down = sp1(x, L - 1)
		else
			up = sp1(x, y + 1)
			down = sp1(x, y - 1)
		end if
		R = up + down + right + left
		En = En - sp1(x, y)*R
	enddo
enddo
CALCULATE_ENERGY = En	! Energy
end function
!------------------------------------------
real function CALCULATE_MAGNETIZATION(sp2)
use SHARED_DATA
implicit none
integer, intent(in) :: sp2(L, L)	! Spin
CALCULATE_MAGNETIZATION = sum(sp2)	! Magnetization
end function
!------------------------------
program ising	! Main program
use SHARED_DATA
implicit none
integer sweep	! Loop counter
real T, B, B2	! Temperature, beta = 1/k_{B}T, beta^2
real summ0, summ1, summ2, summ3, summ4	! Accumulators
real E, M ! energy, magnetization
real CALCULATE_ENERGY, CALCULATE_MAGNETIZATION ! External functions
real mean_E, mean_M, C_v, chi ! Mean E, mean M, specific heat, susceptibility
integer :: spin(L, L) = 1 ! Initial configuration: all spins up
integer c	! Counter
open(1, file = 'output.txt') ! Output file: beta, T, E, M, Cv, chi
open(2, file = 'eqltst.txt') ! Output file for equilibration test
T = T_0		! Temperature initialization
c = 0
10 format (6f12.6)
write(1, *)"#   T^{-1}         T           E          M          C_v         Chi"
write(1, *)"#  ========    ========    ========    ========    ========    ========"
do while (T <= T_f)		! Temperature loop (main loop)
	T = T + dT
	c = c + 1
	summ0 = 0.; summ1 = 0.; summ2 = 0.; summ3 = 0. ! Initializing accumulators
	B = 1./T; B2 = B*B ! Beta, beta^2
	do sweep = 1, eqstp ! Equilibration loop
		call METROPOLIS(spin, B)
		! Activate the following 4 comment lines
		! (from the line 142 to 145) to check whether 
		! equilibration is achieved for a given 
		! temperature (here, T = 2). 
		! To this end, plot eqltst.txt.
		!E = CALCULATE_ENERGY(spin)
		!if (int(T) == 2) then
		!	write(2,*) step, E
		!end if
	enddo
	!-------------------------------------
	do sweep = 1, mcstp     ! Main MC loop 
		call METROPOLIS(spin, B)
		E = CALCULATE_ENERGY(spin)
		M = CALCULATE_MAGNETIZATION(spin)
		summ0 = summ0 + E       ! E accumulator
		summ1 = summ1 + E*E     ! E^2 accumulator
		summ2 = summ2 + M       ! M accumulator
		summ3 = summ3 + M*M     ! M^2 accumulator
	end do
	!--------------------------------------------------
	! Thermodynamic averages at each temperature value
	mean_E = summ0*norm       ! Mean energy
	mean_M = summ2*norm       ! Mean magnetization 
	C_v = (norm*summ1 - nr2*summ0*summ0)*B2    ! Specific heat
	chi = (norm*summ3 - nr2*summ2*summ2)*B     ! Susceptibility
	write(1, 10)B, T, mean_E, mean_M, C_v, chi
	write(*, *)" >>> Temperature step", c
end do
end program
