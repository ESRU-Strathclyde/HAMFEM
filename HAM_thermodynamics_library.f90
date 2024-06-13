! ***************************************************************************************
! ***************************************************************************************
!
								MODULE HAM_thermodynamics_library	  ! updated: mar 2005
!
! ***************************************************************************************
! ***************************************************************************************
	implicit none
!
		double precision,public::rho_l,L_v,R_v,cap_l,cap_v,cap_a,acc_g
!
!	rho_l:		density of liquid water
!	L_v:		vaporisation heat of water
!	R_v:		gas constant of water vapour
!	cap_l:		thermal capacity of liquid water
!	cap_v:		thermal capacity of water vapour
!	cap_a:		thermal capacity of air
!	acc_g:		acceleration of gravity
!
CONTAINS
!
!   _____________________________________________________________________________________
!
	subroutine thedyn_initialise(orie)
!
!
!   _____________________________________________________________________________________
!	initialisation of standard physical properties of water, vapour and air		
		double precision,intent(in)::orie(:)
!		---------------------------------------------------------------------------------
		rho_l=1.0d3; cap_l=4.18d3; cap_v=2.05d3; cap_a=1.024d3; L_v=2.5d6; R_v=4.6189d2
		acc_g=9.81d0; if (orie(3).LT.0.0d0.AND.orie(4).LT.0.0d0) acc_g=0.0d0
!
	end subroutine thedyn_initialise
!   _____________________________________________________________________________________
!
	pure function thedyn_psat(T)
!
!
!   _____________________________________________________________________________________
!	Calculates the saturation pressure [Pa] for a temperature T in [K]
		double precision,intent(in)::T
		double precision::thedyn_psat
!		---------------------------------------------------------------------------------
		thedyn_psat=dexp(6.58094d1-7.06627d3/T-5.976d0*dlog(T))
!
	end function thedyn_psat
!   _____________________________________________________________________________________
!
	pure function thedyn_dpsatdt(T)
!
!
!   _____________________________________________________________________________________
!	Calculates the derivative of saturation pressure [Pa] for a temperature T in [K]
		double precision,intent(in)::T
		double precision::thedyn_dpsatdt
!		---------------------------------------------------------------------------------
		thedyn_dpsatdt=(7.06627d3/(T*T)-5.976d0/T)*thedyn_psat(T)
!
	end function thedyn_dpsatdt
!   _____________________________________________________________________________________
!
	pure function thedyn_phi(T,pc)
!
!
!   _____________________________________________________________________________________
!	Calculates the relative humidity [-] for a capillary pressure in [Pa] and a tempera-
!	ture T in [K]
		double precision,intent(in)::T,pc
		double precision::thedyn_phi
!		---------------------------------------------------------------------------------
		thedyn_phi=dexp(pc/(rho_l*R_v*T))
!
	end function thedyn_phi
!   _____________________________________________________________________________________
!
	pure function thedyn_gam(T)
!
!
!   _____________________________________________________________________________________
!	Calculates the normalised thermal derivative of surface tension [1/K] for a tempera-
!	ture in [K]
		double precision,intent(in)::T
		double precision::thedyn_gam
!		---------------------------------------------------------------------------------
		thedyn_gam=0.0d0 !-0.1682d-3/(1.2157d-1-0.1682d-3*T)
!
	end function thedyn_gam
!   _____________________________________________________________________________________
!
	pure function thedyn_airden(T)
!
!
!   _____________________________________________________________________________________
!	Returns the air density [kg/m³]
		double precision,intent(in)::T
		double precision::thedyn_airden
!		---------------------------------------------------------------------------------
		thedyn_airden=1.01325d5/(2.87d2*T)
!       
	end function thedyn_airden
!   _____________________________________________________________________________________
!
	pure function thedyn_grav(orie)
!
!
!   _____________________________________________________________________________________
!	Returns the gravity vector [m/s²] depending on the input angles in grav(:)
		double precision,intent(in)::orie(:)
		double precision::thedyn_grav(3)
!		---------------------------------------------------------------------------------
		thedyn_grav(1)=-acc_g*dcos(orie(3))*dcos(orie(4))
		thedyn_grav(2)=-acc_g*dsin(orie(3))*dcos(orie(4))
		thedyn_grav(3)=-acc_g*dsin(orie(4))
!       
	end function thedyn_grav
!   _____________________________________________________________________________________


END MODULE HAM_thermodynamics_library
