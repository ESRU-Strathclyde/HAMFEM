
! This file is part of HAMFEM.
! HAMFEM is free software. You can redistribute it and/or modify it under the terms of the 
! GNU General Public License as published by the Free Software Foundation version 3 or later.

! HAMFEM is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.

! Copyright Hans Janssen, Katholic University Leuven, and 
! the Energy Systems Research Unit - University of Strathclyde.

! ***************************************************************************************
! ***************************************************************************************
!
								MODULE HAM_material_library			  
!
! ***************************************************************************************
! ***************************************************************************************

!
!
	use HAM_thermodynamics_library
	implicit none
!
!	*** private variable declarations ***
!
        integer,private::Kl_offsets(99),nret(100)
		double precision,private::Kl_data(9999,2),matdata(100,10)
		double precision,private::wi(100),Ci(100),vdc(100,7),keq(100),ka(100),kb(100)
		double precision,private::reta(100,3),retm(100,3),retn(100,3),retw(100,3)
		double precision,private::kreta(100,3),kretm(100,3),kretn(100,3),kretw(100,3)
		COMMON/data2/matdata
!
!	Kl_data is used to store linear conductivity splines (data points) of all materials
!	Kl_data(i,1) contains the log pc or w values, Kl_data(i,2) contains the log Kl valu-
!	es.  Make sure, that the array is large enough for all the spline data.
!	Kl_offsets stores the starting indices for the splines of the different materials.
!
!
CONTAINS
!
!   _____________________________________________________________________________________
!
	subroutine matlib_initialise(main_dir,xts,xtslenght)
!
!
!   _____________________________________________________________________________________
!	Initializes the library by loading data points from files if needed. The character 
!	string is hereby the name of the directory (including the trailing backslash) whe-
!	re the material data files can be found.
		character*(*),intent(in)::main_dir
		character(100),intent(in)::xts(:)
		integer,intent(in)::xtslenght(:)
		character(100)::dum1,dum2   
		character(xtslenght(1))::matdb
		character(xtslenght(4))::permdb
		integer::n,i,j,nmat
!		---------------------------------------------------------------------------------
	    Kl_offsets=0; Kl_offsets(1)=1

		permdb=TRIM(xts(4))
		open (unit=1001,file=TRIM(main_dir)//permdb//'.pdb'&
		&,action='read')

		do j=2,size(Kl_offsets) ! number of materials
			read (1001,*) n ! read number of data points
			Kl_offsets(j)=Kl_offsets(j-1)+n
			do i=Kl_offsets(j-1),Kl_offsets(j)-1
				read(1001,*) Kl_data(i,:)
			enddo
		enddo

		close(1001)

		matdb=TRIM(xts(1))
		open (unit=1009,file=TRIM(main_dir)//matdb//'.mtr',action='read')
		read (1009,*) nmat ! read number of material in the database
		read (1009,*);read (1009,*) ! read (skip) header line
		do j=1,nmat ! read each line
		read (1009,*) matdata(j,1),matdata(j,2),dum1,dum2,matdata(j,5),matdata(j,6),&
		&matdata(j,7),matdata(j,8),matdata(j,9),matdata(j,10),vdc(j,1),vdc(j,2),vdc(j,3),&
		&vdc(j,4),vdc(j,5),vdc(j,6),vdc(j,7),nret(j),reta(j,1),retn(j,1),retm(j,1),&
		&retw(j,1),Ci(j),wi(j),reta(j,2),retn(j,2),retm(j,2),retw(j,2),reta(j,3),retn(j,3),&
		&retm(j,3),retw(j,3),keq(j),ka(j),kb(j),kreta(j,1),kretn(j,1),kretm(j,1),kretw(j,1),&
		&kreta(j,2),kretn(j,2),kretm(j,2),kretw(j,2),kreta(j,3),kretn(j,3),kretm(j,3),&
		&kretw(j,3)

!		matdata(j,1) ID number
!		matdata(j,2) case number in old HAMFEM versions. only for documentation purpose.
!		dum1         material name
!		dum2         data source
!		matdata(j,5) 'bulk density'	
!		matdata(j,6) 'thermal capacity'	
!		matdata(j,7) 'thermal conductivity - dry' 
!		matdata(j,8) 'thermal conductivity - moisture dependant term'
!		matdata(j,9) 'capillary moisture content'
!		matdata(j,10) 'absorption for short wave radiation [-]'
!		vdc(1)       'vapour diffusion coefficient - equation type'
!		vdc(2 to 7)  'coefficients to the vapour diffusion coefficient calculation'
!       nret to retw 'coefficients to the moisture storage parameters calculation'
!		keq,ka,kb,kreta to kretw	,coefficients to the calculation ofthe logarithm
!									 to the base 10 of the liquid permeability/conductivity,
		enddo

		close(1009)

!
	end subroutine matlib_initialise
!   _____________________________________________________________________________________
!
	subroutine materialproperties	(t0,t,p,m0,m,mat,tpm,tpmloc,orie,Cxx,Cxy,Cxz,&
									&Kxx,Kxy,Kxz,Kxxgr,Kxygr,Kxzgr,Quan0,Quan)
!
!
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::tpm,tpmloc,mat
		double precision,intent(in)::t0,t,p,m0,m,orie(:)
		double precision,intent(out)::Cxx,Cxy,Cxz,Kxx(3),Kxy(3),Kxz(3),Kxxgr(3),Kxygr(3),Kxzgr(3),&
			&Quan0,Quan
		double precision::rho_a,Kpp(3),w,w0,dwdm,phi,psat,dpsatdt,delta(3),grav(3),&
			&Kmm_l(3),Kmm_v(3),rho_m,cap_m
			integer::dum23=0

!		---------------------------------------------------------------------------------
!	calculation of material properties for t,p, m
		Cxx=0.0d0; Cxy=0.0d0; Cxz=0.0d0; Kxx=0.0d0; Kxy=0.0d0; Kxz=0.0d0 
		Kxxgr=0.0d0; Kxygr=0.0d0; Kxzgr=0.0d0; Quan0=0.0d0; Quan=0.0d0 
!
		Quantity:select case(tpmloc)
			case(010) ! air transfer
				call matlib_moisture_storage(mat,t,m,w,dwdm)				
					rho_a=thedyn_airden(t)
					grav=thedyn_grav(orie)
					Kpp=matlib_airperm(mat,t,m,w)
				Kxx=rho_a*Kpp
				Kxxgr=rho_a*rho_a*grav*Kpp

			case(001) ! moisture transfer
				call matlib_moisture_storage(mat,t0,m0,w0,dwdm)
				call matlib_moisture_storage(mat,t,m,w,dwdm)
				Quan=w
				Quan0=w0
				Cxx=dwdm
					phi=thedyn_phi(t,m)
					psat=thedyn_psat(t)
					dpsatdt=thedyn_dpsatdt(t)
					delta=matlib_delta(mat,t,m,w)
					grav=thedyn_grav(orie)
					Kmm_l=matlib_Kl(mat,t,m,dlog10(-m),w)
					Kmm_v=delta*phi*psat/(rho_l*R_v*t)
				Kxx=Kmm_l+Kmm_v
				Kxxgr=Kmm_l*grav*rho_l
				if (tpm/100.EQ.1) then
					Kxy=delta*phi*(dpsatdt-psat*m/(R_v*rho_l*t*t))
				endif
				if (mod(tpm,100)/10.EQ.1) then
						rho_a=thedyn_airden(t)
						Kpp=matlib_airperm(mat,t,m,w)
					Kxz=phi*psat/(R_v*t)*Kpp
					Kxzgr=phi*psat/(R_v*t)*Kpp*rho_a*grav
				endif
			case(100) ! heat transfer
				call matlib_moisture_storage(mat,t0,m0,w0,dwdm)
				call matlib_moisture_storage(mat,t,m,w,dwdm)
					rho_m=matlib_rho(mat) !density red from the material library
					cap_m=matlib_cet(mat) !thremal capacity red from the material library
				Cxx=rho_m*cap_m+cap_l*w
				Quan=rho_m*cap_m*t+cap_l*w*t
				Quan0=rho_m*cap_m*t0+cap_l*w0*t0
				Kxx=matlib_lambda(mat,t,w)
				if (mod(tpm,100)/10.EQ.1) then
						rho_a=thedyn_airden(t)
						grav=thedyn_grav(orie)
						Kpp=matlib_airperm(mat,t,m,w)
					Kxy=cap_a*t*rho_a*Kpp
					Kxygr=cap_a*t*rho_a*rho_a*grav*Kpp
					if (mod(tpm,10).EQ.1) then
							phi=thedyn_phi(t,m)
							psat=thedyn_psat(t)
							grav=thedyn_grav(orie)
						Kxy=Kxy+(L_v+cap_v*t)*phi*psat/(R_v*t)*Kpp
						Kxygr=Kxygr+(L_v+cap_v*t)*phi*psat/(R_v*t)*Kpp*rho_a*grav
					endif
				endif
				if (mod(tpm,10).EQ.1) then
					Cxz=cap_l*t*dwdm
						phi=thedyn_phi(t,m)
						psat=thedyn_psat(t)
						dpsatdt=thedyn_dpsatdt(t)
						delta=matlib_delta(mat,t,m,w)
						grav=thedyn_grav(orie)
						Kmm_l=matlib_Kl(mat,t,m,dlog10(-m),w)
						Kmm_v=delta*phi*psat/(rho_l*R_v*t)
					Kxx=Kxx+(L_v+cap_v*t)*delta*phi*(dpsatdt-psat*m/(R_v*rho_l*t*t))
					Kxz=cap_l*t*Kmm_l+(L_v+cap_v*t)*Kmm_v
					Kxzgr=cap_l*t*Kmm_l*grav*rho_l
				endif
		end select quantity 
!
	end subroutine materialproperties
!	_____________________________________________________________________________________
!
	pure function matlib_rho(mat_nr)
!   _____________________________________________________________________________________
!	Returns the bulk density of the material [kg/m3]
		integer,intent(in)::mat_nr
		double precision::matlib_rho
!		---------------------------------------------------------------------------------
		matlib_rho=matdata(mat_nr,5)
!
	end function matlib_rho
!   _____________________________________________________________________________________
!
!
!
	pure function matlib_cet(mat_nr)
!   _____________________________________________________________________________________
!	Returns the thermal capacity [J/kgK]
		integer,intent(in)::mat_nr
		double precision::matlib_cet
!		---------------------------------------------------------------------------------
		matlib_cet=matdata(mat_nr,6)
!
	end function matlib_cet
!   _____________________________________________________________________________________
!
!   _____________________________________________________________________________________



!
	pure function matlib_lambda(mat_nr,T,w)
!   _____________________________________________________________________________________
!	Returns the thermal conductivity [W/mK] for a given temperature T in [K] and moisture 
!	content w in [kg/m3]
		integer,intent(in)::mat_nr
		double precision,intent(in)::T,w
		double precision::matlib_lambda(3)
!		---------------------------------------------------------------------------------
		matlib_lambda=matdata(mat_nr,7)+matdata(mat_nr,8)*w 
!
	end function matlib_lambda
!   _____________________________________________________________________________________
!
!
!
!

	pure function matlib_wsat(mat_nr)
!
!
!   _____________________________________________________________________________________
!	Returns the capillary moisture content [kg/m3] 
		integer,intent(in)::mat_nr
		double precision::matlib_wsat
!		---------------------------------------------------------------------------------
		matlib_wsat=matdata(mat_nr,9)
!
	end function matlib_wsat
!   _____________________________________________________________________________________
!
	pure subroutine matlib_moisture_storage(mat_nr,T,pc,w,C)
!   _____________________________________________________________________________________



!	Calculates the moisture storage parameters, moisture content w [kg/m3] and moisture 
!	capacity C [kg/m3Pa], in dependence of the capillary pressure pc in [Pa] and tempe-
!	rature in [K]
		integer,intent(in)::mat_nr
		double precision,intent(in)::t,pc
		double precision,intent(out)::w,C
		integer::i
		double precision::tmp,tmp2
!		---------------------------------------------------------------------------------

                C=0d0
				w=0d0
				do i=1,nret(mat_nr)
					tmp = (reta(mat_nr,i)*pc)**retn(mat_nr,i)
					w = w + retw(mat_nr,i) / ( (1d0 + tmp)**retm(mat_nr,i))
					tmp2 = (1d0 + tmp)**retm(mat_nr,i)
					C = C - retw(mat_nr,i)/tmp2 * retm(mat_nr,i)*retn(mat_nr,i)&
					&*tmp/((1d0 + tmp)*pc)
				enddo
				C = C*Ci(mat_nr)
				w = w*wi(mat_nr)
!
	end subroutine matlib_moisture_storage
!   _____________________________________________________________________________________
!
	pure function matlib_delta(mat_nr,T,pc,w)
!
!
!   _____________________________________________________________________________________
!	Returns the vapour diffusion coefficient in [s] for a given temperature T in [K], 
!	capillary pressure pc in [Pa] and moisture mass density w in [kg/m3(REV)]
	    integer,intent(in)::mat_nr
        integer::dum1
		double precision,intent(in)::T, pc, w
		double precision::tmp, mew, matlib_delta
!		---------------------------------------------------------------------------------
		tmp = 1.0d0 - w/matlib_wsat(mat_nr)
		dum1=vdc(mat_nr,1)
		material:select case(dum1)
			case(1)
				matlib_delta = vdc(mat_nr,2) * tmp/(R_v*T*vdc(mat_nr,3)*(vdc(mat_nr,4)*&
				&tmp*tmp + vdc(mat_nr,5)))
			case (2)
				matlib_delta = vdc(mat_nr,2)*(t/vdc(mat_nr,3))**vdc(mat_nr,4)/(R_v*t)/vdc(mat_nr,5)
			case (3)
				matlib_delta = vdc(mat_nr,2)*(t/vdc(mat_nr,3))**vdc(mat_nr,4)/(R_v*t)*&
				&(vdc(mat_nr,5)+vdc(mat_nr,6)*dexp(vdc(mat_nr,7)*thedyn_phi(t,pc)))
			case (4)
				matlib_delta = vdc(mat_nr,2)
			case (5)
                                tmp=pc
                                if(pc.LT.vdc(mat_nr,5))tmp=vdc(mat_nr,5)
                                if(pc.GT.vdc(mat_nr,6))tmp=vdc(mat_nr,6)
				matlib_delta = vdc(mat_nr,2)/tmp+vdc(mat_nr,3)*dexp(vdc(mat_nr,3)/tmp)
			case (6)
				matlib_delta = vdc(mat_nr,2)*(1-thedyn_phi(t,pc))+vdc(mat_nr,3)
		end select material
!                
	end function matlib_delta
!   _____________________________________________________________________________________


	pure function matlib_logKl(mat_nr,T,pc,lpc,w)
!
!
!   _____________________________________________________________________________________
!	Returns the logarithm to the base 10 of the liquid permeability/conductivity, depending 
!	on the capillary pressures pc in [Pa], logarithm to the base 10 of the capillary 
!	pressures lpc [unitless] and the moisture mass density [kg/m3(REV)]. You need to 
!	specify all three moisture potential parameters  so that if a different potential then 
!	pc is used, the function does not need to call the appropriate material functions 
!	functions itself and thus calculation time is saved.
		integer,intent(in)::mat_nr
		double precision,intent(in)::T,pc,lpc,w
		integer::i,n,left,right,middle,dum23
		double precision::matlib_logKl(3),Ks,tau,dum1,dum2,dum3,dum4,dum5
		double precision::tmp,tmp2
!		---------------------------------------------------------------------------------
	dum23=keq(mat_nr)
		material:select case(dum23)
			case(1)
				matlib_logKl=ka(mat_nr)  ! ka(mat_nr)=-50 indicates Kl = 0
			case(2)
				tmp=w+ka(mat_nr)
				tmp=kb(mat_nr)+0.0704d0*tmp-1.742d-4*tmp*tmp-2.7953d-6*tmp*tmp*tmp-&
					&1.1566d-7*tmp*tmp*tmp*tmp+2.5969d-9*tmp*tmp*tmp*tmp*tmp
				matlib_logKl=tmp*0.4342944819d0
			case(3)
				Ks=ka(mat_nr); tau=kb(mat_nr)
				dum2=0.0d0; dum3=0.0d0; dum4=0.0d0
				do i=1,3
					dum1=(-kreta(mat_nr,i)*pc)**kretn(mat_nr,i)
					dum2=dum2+kretw(mat_nr,i)*(1.0d0+dum1)**(-kretm(mat_nr,i))
					dum3=dum3+kretw(mat_nr,i)*kreta(mat_nr,i)*(1.0d0-(dum1/(1.0d0+dum1))&
					&**kretm(mat_nr,i))
					dum4=dum4+kretw(mat_nr,i)*kreta(mat_nr,i)
				enddo
				matlib_logKl=dlog10(Ks*(dum2**tau)*((dum3/dum4)**2.0d0))
	
			case(4)
				! calculate left, right and middle offset
				left=Kl_offsets(mat_nr); right=Kl_offsets(mat_nr+1)-1 
				middle=(left+right)/2
				! if log pc bigger then maximum log pc value, return lowest conductivity value
				if (lpc.GT.Kl_data(right,1)) then
					matlib_logKl=Kl_data(right,2)
				else
					! perform a binary lookup to find the index of the log pc value that 
					! is just a bit smaller then lpc
					do while (middle.NE.left)
						if (lpc.LT.Kl_data(middle,1)) then
							right=middle
						else                            
							left=middle
						endif
						middle=(left+right)/2
					end do
					! return interpolated value
					matlib_logKl=Kl_data(middle,2)+(Kl_data(middle+1,2)-Kl_data(middle,2))*&
						(lpc-Kl_data(middle,1))/(Kl_data(middle+1,1)-Kl_data(middle,1))
				endif
		end select material

	end function matlib_logKl
!   _____________________________________________________________________________________
!
	pure function matlib_Kl(mat_nr,T,pc,lpc,w)
!
!
!   _____________________________________________________________________________________
!	Returns the liquid permeability/conductivity in [s], depending on the capillary pressures 
!	pc in [Pa], logarithm to the base 10 of the capillary pressures lpc [unitless] and the
!	moisture mass density [kg/m3(REV)]. You need to specify all three moisture potential parameters  
!	so that if a different potential then pc is used, the function does not need to call the 
!	appropriate material functions functions itself and thus calculation time is saved.
		integer,intent(in)::mat_nr
		double precision,intent(in)::T,pc,lpc,w
		double precision::matlib_Kl(3)
!		---------------------------------------------------------------------------------
		matlib_Kl=10d0**matlib_logKl(mat_nr,T,pc,lpc,w)
!       
	end function matlib_Kl
!   _____________________________________________________________________________________
!
!
!
!
	pure function matlib_airperm(mat_nr,T,pc,w)
!   _____________________________________________________________________________________
!	Returns the air permeability [s]
		integer,intent(in)::mat_nr
		double precision,intent(in)::T,pc,w
		double precision::matlib_airperm(3),dum
!		---------------------------------------------------------------------------------
		material:select case(mod(mat_nr,100))
			case default
				matlib_airperm(1)=1.0d-8
				matlib_airperm(2)=matlib_airperm(1)
				matlib_airperm(3)=matlib_airperm(1)
				dum=1.0d0-w/matlib_wsat(mat_nr)
				matlib_airperm=matlib_airperm/(5.0d-8*t+3.3425d-6)*dum/(0.503d0*(dum**&
					&2.0d0)+0.497d0)
		end select material
!       
	end function matlib_airperm
!   _____________________________________________________________________________________
!
!

!
END MODULE HAM_material_library
!
!	More information on the material properties can be found in:
!
!	CEN/TC89 WG10: Hygrothermal performance of building components and building elements 
!	- Assessment of moisture transfer by numerical simulation (prenormative text) 2003.
!
!	Hagentoft C-E, Kalagasidis AS, Adl-Zarrabi B, Roels S, Carmeliet C, Hens H, Grunewald 
!	J, Funk M, Becker R, Shamir D, Adan O, Brocken H, Kumaran K, Djebbar R. Assessment 
!	method of numerical prediction models for combined heat, air and moisture transfer in 
!	building components: benchmarks for one-dimensional cases. Journal of Thermal Envelope 
!	and Building Science 2004, 27: 327-352.
!
!	Peuhkuri R.  Moisture dynamics in building envelopes. PhD thesis, Technical University
!	of Denmark, 2003.
!
!	Kumaran K.  Heat, air and moisture transfer in insulated envelope parts, Task 3: mate-
!	rial properties, IEA Annex 24 final report, 1996.

