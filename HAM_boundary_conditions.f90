
! This file is part of HAMFEM.
! HAMFEM is free software. You can redistribute it and/or modify it under the terms of the 
! GNU General Public License as published by the Free Software Foundation version 3 or later.

! HAMFEM is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.

! Copyright Hans Janssen, Katholic University Leuven, and 
! the Energy Systems Research Unit - University of Strathclyde.!   
!
!
!****************************************************************************************
!****************************************************************************************
!
							MODULE HAM_boundary_conditions            
!
!****************************************************************************************
!****************************************************************************************
! Documentation can be found in references at bottom.
!
	use HAM_thermodynamics_library
	implicit none
!
	double precision,private::bdst,locat(2)
	double precision,private,allocatable::clidat(:,:),catch(:,:,:,:)
!
!	bdst:		interval between atmospheric boundary value transitions
!	clidat:		contains climate data (exterior and interior)
!	catch:		contains numerical global catch factors for specified positions on facade
!
CONTAINS
!   _____________________________________________________________________________________
!
	subroutine bndcon_initialise(init_dir,atmos,xts,xtslenght,facloc,bdstep,cliquant,cliunit)
!
!
!   _____________________________________________________________________________________
!	initialisation of climatefile and catch factors		
		double precision,intent(in)::facloc(:)
		logical,intent(in)::atmos
		character(100),intent(in)::xts(:)
		character(1),intent(in)::cliunit
		integer,intent(in)::cliquant,xtslenght(:)
		character*(*),intent(in)::init_dir     
		double precision,intent(inout)::bdstep
		integer::i,j,k,l,m,n,o,p,numb
		character(xtslenght(2))::clim
		character(xtslenght(3))::wdrcr
		character(1)::step
!		character(2)::numb
!		---------------------------------------------------------------------------------
		if (atmos) then
			clim=TRIM(xts(2)); step=cliunit; numb=cliquant ! initialisation of climate data
			stepsize: select case(step)
				case default; bdstep=0.0d0
				case('s'); bdstep=1.0d0
				case('m'); bdstep=6.0d1
				case('h'); bdstep=3.6d3
				case('d'); bdstep=8.64d4
				case('y'); bdstep=3.1536d7
			end select stepsize
			number: select case(numb)
				case (01); bdstep=bdstep
				case (05); bdstep=bdstep*5.0d0
				case (10); bdstep=bdstep*1.0d1
				case (20); bdstep=bdstep*2.0d1
				case (30); bdstep=bdstep*3.0d1
			end select number
			bdst=bdstep
			allocate(clidat(idint(3.1536d7/bdst),19))
			open(2001,file=init_dir//clim//'.cli',action='read')
			read(2001,*); read(2001,*)
			do i=1,size(clidat,1)
				read(2001,*) clidat(i,:)
			enddo
			close(2001);
!			
			locat=facloc ! initialisation of wind-driven-rain catch factors
			wdrcr=TRIM(xts(3))
			open(2002,file=init_dir//wdrcr//'.wdr',action='read')
			read(2002,*); read(2002,*) i,j,k,l
			allocate(catch(i+1,j+1,k+1,l+1)); catch=0.0d0
			read(2002,*); read(2002,*) catch(2:i+1,1,1,1); read(2002,*) catch(1,2:j+1,1,1)
			read(2002,*); read(2002,*) catch(1,1,2:k+1,1); read(2002,*) catch(1,1,1,2:l+1)
			read(2002,*)
			do n=2,j+1
				do m=2,i+1
					read(2002,*); read(2002,*) catch(m,n,1,1)
					catch(m,n,2:k+1,1)=catch(1,1,2:k+1,1)
					catch(m,n,1,2:l+1)=catch(1,1,1,2:l+1)
					do o=2,k+1
						read(2002,*) catch(m,n,o,2:l+1)
					enddo
				enddo
			enddo					
			close(2002)
		endif
!
	end subroutine bndcon_initialise
!   _____________________________________________________________________________________
!
	subroutine boundaryconditions	(t,p,m,c,orie,bdloc,mat,tpm,tpmloc,atmos,bdkind,Bxx,&
									&dBxxdx,ttot,dt)
!
!
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::bdloc,mat,tpm,tpmloc
		double precision,intent(in)::t,p,m,c(:),orie(:),ttot,dt
		logical,intent(in)::atmos
		integer,intent(out)::bdkind
		double precision,intent(out)::Bxx,dBxxdx
		double precision::time,Tae,Dir,Dif,RHe,Rav,Wsp,Wdi,Clo,Tai,RHi,clima(18),Hsurt,&
			&Hsurm,lati
		COMMON/data1/lati
		
!	        real extheatc/1.0/,extheati/1.0/,intheatc/1.0/,intheati/1.0/
!		real extmoistc,extmoisti,intmoistc,intmoisti,tso,tsi,testcalc
!                COMMON/data2/extheatc,extheati,intheatc,intheati,extmoistc,extmoisti,&
!			&intmoistc,intmoisti

! only for debugging:
        double precision::Bxx1,dBxxdx1,Bxx2,dBxxdx2,Bxx3,dBxxdx3,wdrm

	real extesprc(3),extespri(3),intesprc(3),intespri(3),&
			&tso,tsi,testcalc
        COMMON/data3/extesprc,extespri,intesprc,intespri


!		---------------------------------------------------------------------------------
!	calculation of boundary condition values
		bdkind=0; Bxx=0.0d0; dBxxdx=0.0d0
		if (.NOT.atmos) then  ! non atmospherical (defined by user hence)
			Tae=293.15d0; RHe=0.4d0; Hsurt=2.0d1; Hsurm=1.0d-8
			!boundnonatmos:select case(mod(bdloc,100))
			boundnonatmos:select case(bdloc)
				case(1) ! boundary location 1 ! prescribed fluxes
					quantity1:select case(tpmloc)
						case(100) ! temperature
							bdkind=2
							Bxx=Hsurt*(Tae-t)+(L_v+cap_v*t)*Hsurm*(RHe*thedyn_psat(Tae)-&
								&thedyn_phi(t,m)*thedyn_psat(t))
							dBxxdx=-Hsurt-(L_v+cap_v*t)*Hsurm*thedyn_dpsatdt(t)
						case(010) ! air
							bdkind=1
							Bxx=1.0d0
							dBxxdx=0.0d0
						case(001) ! moisture
							bdkind=2
							Bxx=Hsurm*(RHe*thedyn_psat(Tae)-thedyn_phi(t,m)*thedyn_psat(t))
							dBxxdx=-Hsurm*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
					end select quantity1
				case(2) ! boundary location 2
					Tae=285.15d0; RHe=0.9d0
					quantity2:select case(tpmloc)
						case(100)
							bdkind=2
							Bxx=1.5d1*(Tae-t)+(L_v+cap_v*t)*Hsurm*(RHe*thedyn_psat(Tae)-&
								&thedyn_phi(t,m)*thedyn_psat(t))
							dBxxdx=-2.5d1
						case(010)
							bdkind=1
							Bxx=1.0d0
							dBxxdx=0.0d0
						case(001)
							bdkind=2
							Bxx=Hsurm*(RHe*thedyn_psat(Tae)-thedyn_phi(t,m)*thedyn_psat(t))
							dBxxdx=-Hsurm*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
					end select quantity2
				case(3) ! boundary location 3 ! prescribed values
				Tae=285.15d0; RHe=0.2d0
					quantity102:select case(tpmloc)
						case(100)
							bdkind=1
							Bxx=2.73d2
							dBxxdx=0d0
						case(010)
							bdkind=1
							Bxx=1.0d0
							dBxxdx=0.0d0
						case(001)
							bdkind=1
							Bxx=Hsurm*RHe*thedyn_psat(Tae)
							dBxxdx=0.0d0
					end select quantity102
				case(4) ! boundary location 4
				Tae=305.15d0; RHe=0.9d0
					quantity112:select case(tpmloc)
						case(100)
							bdkind=1
							Bxx=3.13d2
							dBxxdx=0d0
						case(010)
							bdkind=1
							Bxx=1.0d0
							dBxxdx=0.0d0
						case(001)
							bdkind=1
							Bxx=Hsurm*RHe*thedyn_psat(Tae)
							dBxxdx=0
					end select quantity112
			end select boundnonatmos
		else  ! atmospherical
			time=dmod(ttot,3.1536d7); clima=clidat(idint(time/bdst)+1,2:19)
			Tae=clima(1)+2.7315d2; Dir=clima(2); Dif=clima(3) 
			RHe=clima(4); Rav=clima(5)/3.6d3; Wsp=clima(6); Wdi=clima(7)
			Clo=clima(8);Tai=clima(9)+2.7315d2; RHi=clima(10)
!			extesprc(1)=clima(11);extespri(1)=clima(12)
!			intesprc(1)=clima(13);intespri(1)=clima(14)
!                        extesprc(3)=clima(15);extespri(3)=clima(16)
!			intesprc(3)=clima(17);intespri(3)=clima(18)

!			boundatmos:select case(mod(bdloc,100))	
			!boundatmos:
			select case(bdloc)	
				case(1)	 ! outside
					Hsurt=-1.0d0; Hsurm=-1.0d0 !1.42d-7*0.5d0
					Hsurm=Hsurm/7.7d-9
					quantity3: select case(tpmloc)
						case(100)
							bdkind=2  ! convection + solar & long wave radiation
							Bxx=bndcon_conv_t(t,Tae,Wsp,Wdi,orie,Hsurt)+bndcon_solar(&
								&time,dir,dif,orie,mat,lati)+bndcon_Longw(t,Tae,orie,Clo,mat)
							if (mod(tpm,10).EQ.1) then ! + latent & sensible heat moisture
								Bxx=Bxx+bndcon_evap_t(t,m,Wsp,Wdi,Tae,RHe,orie,Hsurm)+&
									&bndcon_rain_t(Tae,Rav,Wsp,Wdi,orie,c)  
							endif
							dBxxdx=bndcon_dconv_tdt(Wsp,Wdi,orie,Hsurt)+bndcon_dlongwdt(t,&
								&mat)
							if (mod(tpm,10).EQ.1) then
								dBxxdx=dBxxdx+bndcon_devap_tdt(t,m,Wsp,Wdi,Tae,RHe,orie,&
									&Hsurm)
							endif
						case(010)
							bdkind=0 
						case(001)
							bdkind=2
							Bxx=bndcon_evap_m(t,m,Wsp,Wdi,RHe,Tae,orie,Hsurm)+&
								&bndcon_rain_m(Rav,Wsp,Wdi,orie,c)  ! evaporation+rain
							dBxxdx=bndcon_devap_mdm(t,m,Wsp,Wdi,orie,Hsurm)
					end select quantity3
				case(2)  ! inside
					Hsurt=0.70d0; Hsurm=1.0d-10 !(insulated wall)
					Hsurm=Hsurm/7.7d-9
					quantity488: select case(tpmloc)
						case(100)
							bdkind=2
							Bxx=bndcon_conv_t(t,Tai,Wsp,Wdi,orie,Hsurt)
							if (mod(tpm,10).EQ.1) then
								Bxx=Bxx+bndcon_evap_t(t,m,Wsp,Wdi,Tai,RHi,orie,Hsurm)
							endif
							dBxxdx=bndcon_dconv_tdt(Wsp,Wdi,orie,Hsurt)
							if (mod(tpm,10).EQ.1) then
								dBxxdx=dBxxdx+bndcon_devap_tdt(t,m,Wsp,Wdi,Tai,RHi,orie,Hsurm)
							endif
						case(010)
							bdkind=0 
						case(001)
							bdkind=2
							Bxx=bndcon_evap_m(t,m,Wsp,Wdi,RHi,Tai,orie,Hsurm)
							dBxxdx=bndcon_devap_mdm(t,m,Wsp,Wdi,orie,Hsurm)
					end select quantity488


				! case(3)	 ! outside
					! Hsurt=-1.0d0; Hsurm=-1.0d0 !1.42d-7*0.5d0
					! Hsurm=Hsurm/7.7d-9
					! select case(tpmloc)
						! case(100)
							! bdkind=2  ! convection + solar & long wave radiation
							! Bxx=bndcon_conv_t(t,Tae,Wsp,Wdi,orie,Hsurt)+bndcon_solar(&
								! &time,dir,dif,orie,mat,lati)+bndcon_Longw(t,Tae,orie,Clo,mat)
							! if (mod(tpm,10).EQ.1) then ! + latent & sensible heat moisture
								! Bxx=Bxx+bndcon_evap_t(t,m,Wsp,Wdi,Tae,RHe,orie,Hsurm)+&
									! &bndcon_rain_t(Tae,Rav,Wsp,Wdi,orie,c)  
							! endif
							! dBxxdx=bndcon_dconv_tdt(Wsp,Wdi,orie,Hsurt)+bndcon_dlongwdt(t,&
								! &mat)
							! if (mod(tpm,10).EQ.1) then
								! dBxxdx=dBxxdx+bndcon_devap_tdt(t,m,Wsp,Wdi,Tae,RHe,orie,&
									! &Hsurm)
							! endif
						! case(010)
							! bdkind=0 
						! case(001)
							! bdkind=2
							! Bxx=bndcon_evap_m(t,m,Wsp,Wdi,RHe,Tae,orie,Hsurm)+&
								! &bndcon_rain_m(Rav,Wsp,Wdi,orie,c)  ! evaporation+rain
							! dBxxdx=bndcon_devap_mdm(t,m,Wsp,Wdi,orie,Hsurm)
					! end select 
					

				! There is some issue with using certain values for the case. 4 and 457, for some reason, lead to memory erros. Befor implementing new cases, make sure the arbitrary value assigned does not lead to errors.
				
				
				! case(5)  ! inside
					! Hsurt=0.70d0; Hsurm=1.0d-10 !(insulated wall)
					! Hsurm=Hsurm/7.7d-9
					! !quantity456: 
					! select case(tpmloc)
						! case(100)
							! bdkind=2
							! Bxx=bndcon_conv_t(t,Tai,Wsp,Wdi,orie,Hsurt)
							! if (mod(tpm,10).EQ.1) then
								! Bxx=Bxx+bndcon_evap_t(t,m,Wsp,Wdi,Tai,RHi,orie,Hsurm)
							! endif
							! dBxxdx=bndcon_dconv_tdt(Wsp,Wdi,orie,Hsurt)
							! if (mod(tpm,10).EQ.1) then
								! dBxxdx=dBxxdx+bndcon_devap_tdt(t,m,Wsp,Wdi,Tai,RHi,orie,Hsurm)
							! endif
						! case(010)
							! bdkind=0 
						! case(001)
							! bdkind=2
							! Bxx=bndcon_evap_m(t,m,Wsp,Wdi,RHi,Tai,orie,Hsurm)
							! dBxxdx=bndcon_devap_mdm(t,m,Wsp,Wdi,orie,Hsurm)
					! end select 
					! !quantity456


		! HAMSTAD exercise 1
			case(3)	 ! outside
					Hsurt=2.5d1; Hsurm=0
				select case(tpmloc)
					case(100)
					bdkind=2
					Bxx=Hsurt*(Tae-t)+(L_v+cap_v*t)*Hsurm*(RHe*thedyn_psat(Tae)-&
						&thedyn_phi(t,m)*thedyn_psat(t))
					dBxxdx=-Hsurt-(L_v+cap_v*t)*Hsurm*thedyn_dpsatdt(t)
					case(010)
					bdkind=0 
					case(001)
					bdkind=2
					Bxx=Hsurm*(RHe*thedyn_psat(Tae)-thedyn_phi(t,m)*thedyn_psat(t))
					dBxxdx=-Hsurm*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)					
				end select
			case(5)  ! inside
					Hsurt=7.0d0; Hsurm=2.0d-8
				select case(tpmloc)
					case(100)
					bdkind=2
					Bxx=Hsurt*(Tai-t)+(L_v+cap_v*t)*Hsurm*(RHi*thedyn_psat(Tai)-&
						&thedyn_phi(t,m)*thedyn_psat(t))
					dBxxdx=-Hsurt-(L_v+cap_v*t)*Hsurm*thedyn_dpsatdt(t)
					case(010)
					bdkind=0 
					case(001)
					bdkind=2
					Bxx=Hsurm*(RHi*thedyn_psat(Tai)-thedyn_phi(t,m)*thedyn_psat(t))
					dBxxdx=-Hsurm*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
				end select 

		! HAMSTAD 4 - original folder exercises
				case(500)	 ! outside
				Hsurt=2.5d1; Hsurm=2.0d-7
					quantity7: select case(tpmloc)
						case(100)
							bdkind=2
							Bxx=Hsurt*(Tae-t)+(L_v+cap_v*t)*Hsurm*(RHe-&
							&thedyn_phi(t,m)*thedyn_psat(t))+Rav*cap_l*Tae
							dBxxdx=-Hsurt-(L_v+cap_v*t)*Hsurm*thedyn_dpsatdt(t)
						case(010)
							bdkind=0 
						case(001)
							bdkind=2
							Bxx=Hsurm*(RHe-thedyn_phi(t,m)*thedyn_psat(t))+Rav
						dBxxdx=-Hsurm*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
					end select quantity7
				case(600) ! inside
				! WRITE(*,*) bdloc
					Hsurt=8.0d0; Hsurm=3.0d-8
					quantity8: select case(tpmloc)
						case(100)
							bdkind=2
							Bxx=Hsurt*(Tai-t)+(L_v+cap_v*t)*Hsurm*(RHi-&
							&thedyn_phi(t,m)*thedyn_psat(t))
							dBxxdx=-Hsurt-(L_v+cap_v*t)*Hsurm*thedyn_dpsatdt(t)
						case(010)
							bdkind=0 
						case(001)
							bdkind=2
							Bxx=Hsurm*(RHi-thedyn_phi(t,m)*thedyn_psat(t))
					dBxxdx=-Hsurm*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
					end select quantity8
	
				
		! HAMSTAD exercise 4
			! case(5)	 ! outside
				! Hsurt=2.5d1; Hsurm=2.0d-7
			! quantity7: select case(tpmloc)
					! case(100)
					! bdkind=2
					! Bxx=Hsurt*(Tae-t)+(L_v+cap_v*t)*Hsurm*(RHe-&
					! &thedyn_phi(t,m)*thedyn_psat(t))+Rav*cap_l*Tae
					! dBxxdx=-Hsurt-(L_v+cap_v*t)*Hsurm*thedyn_dpsatdt(t)
					! case(010)
						! bdkind=0 
					! case(001)
					! bdkind=2
					! Bxx=Hsurm*(RHe-thedyn_phi(t,m)*thedyn_psat(t))+Rav
					! dBxxdx=-Hsurm*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
			! end select quantity7
			! case(6)  ! inside
			! Hsurt=8.0d0; Hsurm=3.0d-8
			! quantity8: select case(tpmloc)
					! case(100)
					! bdkind=2
					! Bxx=Hsurt*(Tai-t)+(L_v+cap_v*t)*Hsurm*(RHi-&
						! &thedyn_phi(t,m)*thedyn_psat(t))
					! dBxxdx=-Hsurt-(L_v+cap_v*t)*Hsurm*thedyn_dpsatdt(t)
					! case(010)
						! bdkind=0 
					! case(001)
					! bdkind=2
					! Bxx=Hsurm*(RHi-thedyn_phi(t,m)*thedyn_psat(t))
					! dBxxdx=-Hsurm*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
				! end select quantity8						
!

		! ONE-WAY coupling (only HEAT) - values exported from ESP-r
			case(11)  ! outside
				Hsurt=2.5d1; Hsurm=2.0d-7
			quantity9: select case(tpmloc)
					case(100)
					bdkind=2
					Bxx=extesprc(1)*(t-2.7315d2)+extespri(1)
					if (mod(tpm,10).EQ.1) then ! + latent & sensible heat moisture
						Bxx=Bxx+bndcon_evap_t(t,m,Wsp,Wdi,Tae,RHe,orie,Hsurm)+&
							&bndcon_rain_t(Tae,Rav,Wsp,Wdi,orie,c)  
					endif
					dBxxdx=extesprc(1)
					if (mod(tpm,10).EQ.1) then
						dBxxdx=dBxxdx+bndcon_devap_tdt(t,m,Wsp,Wdi,Tae,RHe,orie,&
						&Hsurm)
					endif					
					case(010)
					bdkind=0 
					case(001)
					bdkind=2
!					Bxx=extesprc(3)*dexp(m/(rho_l*R_v*t))&
!			&*thedyn_psat(t)+extespri(3)
!					dBxxdx=extesprc(3)*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
					Bxx=bndcon_evap_m(t,m,Wsp,Wdi,RHe,Tae,orie,Hsurm)+&
						&bndcon_rain_m(Rav,Wsp,Wdi,orie,c)  ! evaporation+rain
					dBxxdx=bndcon_devap_mdm(t,m,Wsp,Wdi,orie,Hsurm)
				end select quantity9
			case(12)  ! inside
					Hsurt=7.0d0; Hsurm=2.0d-8
				quantity10: select case(tpmloc)
					case(100)
					bdkind=2
					Bxx=intesprc(1)*(t-2.7315d2)+intespri(1)
						!intmoistc,intmoisti	
					!Bxx=Hsurt*(Tai-t)+(L_v+cap_v*t)*Hsurm*(RHi*thedyn_psat(Tae)-&
					!	&thedyn_phi(t,m)*thedyn_psat(t))
					!dBxxdx=-Hsurt-(L_v+cap_v*t)*Hsurm*thedyn_dpsatdt(t)
					dBxxdx=intesprc(1)
					case(010)
					bdkind=0 
					case(001)
					bdkind=2
					Bxx=Hsurm*(RHi*thedyn_psat(Tai)-thedyn_phi(t,m)*thedyn_psat(t))
					dBxxdx=-Hsurm*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
				end select quantity10
! TWO-WAY coupling (only HEAT) - values exported from ESP-r				
			case(44)  ! outside
					Hsurt=-1.0d0; Hsurm=-1.0d0 !1.42d-7*0.5d0
					Hsurm=Hsurm/7.7d-9
			quantity11: select case(tpmloc)
					case(100)
					bdkind=2
					Bxx=extesprc(1)*(t-2.7315d2)+extespri(1)
                                        dBxxdx=extesprc(1)
!					if (mod(tpm,10).EQ.7) then
					if (mod(tpm,10).EQ.1) then ! + latent & sensible heat moisture

!						Bxx=Bxx+bndcon_rain_t(Tae,Rav,Wsp,Wdi,orie,c)

    			if(extesprc(3).NE.0)then
             call bndcon_rain_m_alt(wdrm,Rav,Wsp,Wdi,orie,c)
				Bxx=Bxx+(extesprc(3)*dexp(m/(rho_l*R_v*t))&
				&*thedyn_psat(t)+extespri(3))*(L_v+cap_v*t)&
                                &+wdrm*cap_l*Tae!&
!				&+bndcon_rain_t(Tae,Rav,Wsp,Wdi,orie,c)

				dBxxdx=dBxxdx-extesprc(3)*&
                       		&((extespri(3)/(-extesprc(3))-thedyn_phi(t,m)*thedyn_psat(t))&
                      		&*cap_v-thedyn_phi(t,m)*&
		     		 &(thedyn_dpsatdt(t)-m/(rho_l*R_v*t*t)*thedyn_psat(t))*(L_v+cap_v*t))

!      write(35,'(E17.10,A,E10.3,A,E10.3,A,E10.3,A,E10.3)')&
!      &ttot,'  ',Bxx2,' ',dBxxdx2,' ',Bxx3,' ',dBxxdx3


!				Bxx=Bxx1
!             			dBxxdx=dBxxdx1

!      			else

!				Bxx1=0
!             			dBxxdx1=0
!      write(35,'(E17.10,A,E10.3,A,E10.3,A,E10.3,A,E10.3)')&
!      &ttot,'  ',Bxx1,' ',dBxxdx1,' ',Bxx,' ',dBxxdx

   			endif

!!						Bxx=Bxx+bndcon_evap_t(t,m,Wsp,Wdi,Tae,RHe,orie,Hsurm)+&
!!							&bndcon_rain_t(Tae,Rav,Wsp,Wdi,orie,c)  

!!						dBxxdx=dBxxdx+bndcon_devap_tdt(t,m,Wsp,Wdi,Tae,RHe,orie,&
!!						&Hsurm)

					endif	!end latent heat				
					case(010)
					bdkind=0 
					case(001)
					bdkind=2
					Bxx=extesprc(3)*dexp(m/(rho_l*R_v*t))&
			&*thedyn_psat(t)+extespri(3)
!         Bxx1=0
!  		WRITE(*,'(A,E17.10)')"antes",ttot
!        write(*,'(E17.10,A,E10.3,A,E10.3,A,E10.3,A,E10.3)')&
!         &Rav,'  ',Wsp,'  ',Wdi,'  ',orie(1),'  ',c(1)
             call bndcon_rain_m_alt(wdrm,Rav,Wsp,Wdi,orie,c)
         Bxx=Bxx+wdrm
!        write(*,'(E17.10,A,E10.3,A,E10.3)')&
!         &ttot,'  ',Bxx,' ',Bxx1
!        write(*,'(A)')"done"
!            Bxx1=0
!            Bxx1=size(catch,3)
					dBxxdx=extesprc(3)*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)
!      write(*,'(E17.10,A,E10.3,A,E10.3)')&
!      &ttot,'  ',Bxx,' ',Bxx1
!      write(36,'(E17.10,A,E10.3,A,E10.3)')&
!      &ttot,'  ',Bxx,' ',Bxx1

!					Bxx=bndcon_evap_m(t,m,Wsp,Wdi,RHe,Tae,orie,Hsurm)+&
!						&bndcon_rain_m(Rav,Wsp,Wdi,orie,c)  ! evaporation+rain
!					dBxxdx=bndcon_devap_mdm(t,m,Wsp,Wdi,orie,Hsurm)
				end select quantity11
			case(45)  ! inside
					Hsurt=5.0d0; Hsurm=7.7d-9*Hsurt
				quantity12: select case(tpmloc)
					case(100)
					bdkind=2
					Bxx=intesprc(1)*(t-2.7315d2)+intespri(1)
					dBxxdx=intesprc(1)

!					if (mod(tpm,10).EQ.7) then 
					if (mod(tpm,10).EQ.1) then ! + latent & sensible heat moisture
!			Bxx=Bxx1+7.7d-9*Hsurt*(RHi*thedyn_psat(Tai)-&
!			&dexp(m/(rho_l*R_v*t))*thedyn_psat(t))*(L_v+cap_v*t)
 
!    			dBxxdx=dBxxdx1+7.7d-9*Hsurt*((RHi*thedyn_psat(Tai)-&
!			&dexp(m/(rho_l*R_v*t))*thedyn_psat(t))*cap_v-dexp(m/(rho_l*R_v*t))*&
!			&(thedyn_dpsatdt(t)-m/(rho_l*R_v*t*t)*thedyn_psat(t))*(L_v+cap_v*t))

			! if Hsurm from coupling is equal to 0, calculate latente heat
    			if(intesprc(3).NE.0)then
				Bxx=Bxx+(intesprc(3)*dexp(m/(rho_l*R_v*t))&
				&*thedyn_psat(t)+intespri(3))*(L_v+cap_v*t)

				dBxxdx=dBxxdx-intesprc(3)*&
                       		&((intespri(3)/(-intesprc(3))-thedyn_phi(t,m)*thedyn_psat(t))&
                      		&*cap_v-thedyn_phi(t,m)*&
		     		 &(thedyn_dpsatdt(t)-m/(rho_l*R_v*t*t)*thedyn_psat(t))*(L_v+cap_v*t))

!      write(36,'(E17.10,A,E10.3,A,E10.3,A,E10.3,A,E10.3)')&
!      &ttot,'  ',Bxx1,' ',dBxxdx1,' ',Bxx,' ',dBxxdx
!				Bxx=Bxx1
!             			dBxxdx=dBxxdx1
!      			else
!				Bxx1=0
!             			dBxxdx1=0
!      write(36,'(E17.10,A,E10.3,A,E10.3,A,E10.3,A,E10.3)')&
!      &ttot,'  ',Bxx1,' ',dBxxdx1,' ',Bxx,' ',dBxxdx

   			endif

					endif	!end latent heat	

					case(010)
					bdkind=0 
					case(001)
					bdkind=2
					Bxx=intesprc(3)*thedyn_phi(t,m)&
						&*thedyn_psat(t)+intespri(3)
					dBxxdx=intesprc(3)*thedyn_phi(t,m)*thedyn_psat(t)/(rho_l*R_v*t)

				end select quantity12	

			case(49)  ! outside
				Hsurt=2.5d1; Hsurm=2.0d-7
			quantity13: select case(tpmloc)
					case(100)
					bdkind=2
					Bxx=extesprc(1)*(t-2.7315d2)+extespri(1)
					dBxxdx=extesprc(1)
					case(010)
					bdkind=0 
					case(001)
					bdkind=2
					Bxx=0
					dBxxdx=0
!					Bxx=extesprc(3)*dexp(m/(rho_l*R_v*t))&
!			&*thedyn_psat(t)+extespri(3)
!					dBxxdx=extesprc(3)*thedyn_psat(t)*thedyn_phi(t,m)/(rho_l*R_v*t)

!					Bxx=bndcon_evap_m(t,m,Wsp,Wdi,RHe,Tae,orie,Hsurm)+&
!						&bndcon_rain_m(Rav,Wsp,Wdi,orie,c)  ! evaporation+rain
!					dBxxdx=bndcon_devap_mdm(t,m,Wsp,Wdi,orie,Hsurm)
				end select quantity13

			case(51)  ! BESTEST Floor. fixed temp, humidity

			quantity14: select case(tpmloc)
					case(100)
					bdkind=1
					Bxx=283
					dBxxdx=0
					case(010)
					bdkind=0 
					case(001)
					bdkind=1
					Bxx=-1.44d04
					dBxxdx=0

				end select quantity14

			case(52)  ! fixed 30Ctemp, 30%humidity

			quantity15: select case(tpmloc)
					case(100)
					bdkind=1
					Bxx=303.15
					dBxxdx=0
					case(010)
					bdkind=0 
					case(001)
					bdkind=1
					Bxx=-1.6858d+08
					dBxxdx=0

				end select quantity15

			case(53)  ! fixed 50Ctemp, 30%humidity

			quantity16: select case(tpmloc)
					case(100)
					bdkind=1
					Bxx=323.15
					dBxxdx=0
					case(010)
					bdkind=0 
					case(001)
					bdkind=1
					Bxx=-1.797d+08
					dBxxdx=0

				end select quantity16

			
			end select 
			!boundatmos
		endif				
!
	end subroutine BOUNDARYCONDITIONS
!	____
!
!	pure function bndcon_rain_mm(Ra,Ws,Wd,ori,c)
	subroutine bndcon_rain_m_alt(mm,Ra,Ws,Wd,ori,c)
!
!
!  THIS ROUTINE IS A COPY OF THE FUNCTION bndcon_rain_mm
!it is used because there is a runtime bug when the arrays catc(:,:),catc5(:,:),catc6(:,:)
! are defined as allocatable _____________________________________________________________________________________
!	computation of arriving rain (horizontal+wind-driven) [kg/m2s] by interpolation in
!	the catch factor matrix given in the input file [Blocken, 2004]
		double precision,intent(out)::mm
		double precision,intent(in)::Ra,Ws,Wd,ori(:),c(:)
		integer::inm,jnm,i,j
		double precision::bndcon_rain_mm,verti,horiz,Vloc,cat1,cat2,cat3,cat4,cat5,cat6,&
			&place(2)
!		double precision,allocatable::catc(:,:),catc5(:,:),catc6(:,:)
		double precision::catc(10,17),catc5(10,17),catc6(10,17)
!		---------------------------------------------------------------------------------
!
		horiz=Ra; Vloc=Ws*dcos(Wd-ori(1))
!        write(*,'(E17.10,A,E17.3)')Ra,'  ',Vloc

!			allocate(catc(size(catch,3),size(catch,4)),catc5(size(catch,3),size(catch,4)),&
!				&catc6(size(catch,3),size(catch,4)))
!                        catc=0.0d0;catc5=0.0d0;catc6=0.0d0;

		if (Ra*Vloc.GT.0.0d0) then
			place=locat ! location on facade
			if (size(c).EQ.1) place(1)=place(1)+c(1)
			if (size(c).EQ.2) place(1:2)=place(1:2)+c(1:2)
			i=size(catch,1)-1 ! interpolation based on vertical and horizontal location
			do inm=2,size(catch,1)-1
				if (place(1).GE.catch(inm,1,1,1).AND.place(1).LT.catch(inm+1,1,1,1)) i=inm
			enddo
			j=size(catch,2)-1
			do jnm=2,size(catch,2)-1
				if (place(2).GE.catch(1,jnm,1,1).AND.place(2).LT.catch(1,jnm+1,1,1)) j=jnm
			enddo
!			allocate(catc(size(catch,3),size(catch,4)),catc5(size(catch,3),size(catch,4)),&
!				&catc6(size(catch,3),size(catch,4)))

                        catc=0.0d0;catc5=0.0d0;catc6=0.0d0;

!        write(*,'(E17.10,A,E17.3)')&
!         &size(catch,3),'  ',size(catch,3)!,'  ',Wdi,'  ',orie(1),'  ',c(1)

!        write(*,'(I5,A,I5)')&
!         &size(catch,3),'  ',size(catch,4)

			catc5=catch(i,j,:,:)+(catch(i,j+1,:,:)-catch(i,j,:,:))*(place(2)-&
				&catch(1,j,1,1))/(catch(1,j+1,1,1)-catch(1,j,1,1))
			catc6=catch(i+1,j,:,:)+(catch(i+1,j+1,:,:)-catch(i+1,j,:,:))*(place(2)-&
				&catch(1,j,1,1))/(catch(1,j+1,1,1)-catch(1,j,1,1))
			catc=catc5+(catc6-catc5)*(place(1)-catch(i,1,1,1))/(catch(i+1,1,1,1)-&
				&catch(i,1,1,1))
			i=size(catc,1)-1 ! interpolation based on rain intensity and wind speed

			do inm=2,size(catc,1)-1
				if (Vloc.GE.catc(inm,1).AND.Vloc.LT.catc(inm+1,1)) i=inm
			enddo
			j=size(catc,2)-1
			do jnm=2,size(catc,2)-1
				if (Ra*3.6d3.GE.catc(1,jnm).AND.Ra*3.6d3.LT.catc(1,jnm+1)) j=jnm
			enddo

			cat1=catc(i,j); cat2=catc(i,j+1) ! bounds for interpolation
			cat3=catc(i+1,j); cat4=catc(i+1,j+1) ! bounds for interpolation
			cat5=cat1+(cat2-cat1)*(Ra*3.6d3-catc(1,j))/(catc(1,j+1)-catc(1,j))
			cat6=cat3+(cat4-cat3)*(Ra*3.6d3-catc(1,j))/(catc(1,j+1)-catc(1,j))
			verti=Ra*(cat5+(cat6-cat5)*(Vloc-catc(i,1))/(catc(i+1,1)-catc(i,1)))
!			deallocate(catc,catc5,catc6)		
		else
			verti=0.0d0
		endif
		mm=horiz*dcos(ori(2))+verti*dsin(ori(2))


!			deallocate(catc,catc5,catc6)
!       mm=0.1E-7
!
!
	end subroutine bndcon_rain_m_alt


!
!____________________________________________________________________________
!
	pure function bndcon_rain_m(Ra,Ws,Wd,ori,c)
!
!
!   ____________________________________________________________________________________
!	computation of arriving rain (horizontal+wind-driven) [kg/m/2s] by interpolation in
!	the catch factor matrix given in the input file [Blocken, 2004]
		double precision,intent(in)::Ra,Ws,Wd,ori(:),c(:)
		integer::inm,jnm,i,j
		double precision::bndcon_rain_m,verti,horiz,Vloc,cat1,cat2,cat3,cat4,cat5,cat6,&
			&place(2)
		double precision,allocatable::catc(:,:),catc5(:,:),catc6(:,:)
!		---------------------------------------------------------------------------------
!
		horiz=Ra; Vloc=Ws*dcos(Wd-ori(1))
		if (Ra*Vloc.GT.0.0d0) then
			place=locat ! location on facade
			if (size(c).EQ.1) place(1)=place(1)+c(1)
			if (size(c).EQ.2) place(1:2)=place(1:2)+c(1:2)
			i=size(catch,1)-1 ! interpolation based on vertical and horizontal location
			do inm=2,size(catch,1)-1
				if (place(1).GE.catch(inm,1,1,1).AND.place(1).LT.catch(inm+1,1,1,1)) i=inm
			enddo
			j=size(catch,2)-1
			do jnm=2,size(catch,2)-1
				if (place(2).GE.catch(1,jnm,1,1).AND.place(2).LT.catch(1,jnm+1,1,1)) j=jnm
			enddo
			allocate(catc(size(catch,3),size(catch,4)),catc5(size(catch,3),size(catch,4)),&
				&catc6(size(catch,3),size(catch,4)))
			catc5=catch(i,j,:,:)+(catch(i,j+1,:,:)-catch(i,j,:,:))*(place(2)-&
				&catch(1,j,1,1))/(catch(1,j+1,1,1)-catch(1,j,1,1))
			catc6=catch(i+1,j,:,:)+(catch(i+1,j+1,:,:)-catch(i+1,j,:,:))*(place(2)-&
				&catch(1,j,1,1))/(catch(1,j+1,1,1)-catch(1,j,1,1))
			catc=catc5+(catc6-catc5)*(place(1)-catch(i,1,1,1))/(catch(i+1,1,1,1)-&
				&catch(i,1,1,1))
			i=size(catc,1)-1 ! interpolation based on rain intensity and wind speed
			do inm=2,size(catc,1)-1
				if (Vloc.GE.catc(inm,1).AND.Vloc.LT.catc(inm+1,1)) i=inm
			enddo
			j=size(catc,2)-1
			do jnm=2,size(catc,2)-1
				if (Ra*3.6d3.GE.catc(1,jnm).AND.Ra*3.6d3.LT.catc(1,jnm+1)) j=jnm
			enddo
			cat1=catc(i,j); cat2=catc(i,j+1) ! bounds for interpolation
			cat3=catc(i+1,j); cat4=catc(i+1,j+1) ! bounds for interpolation
			cat5=cat1+(cat2-cat1)*(Ra*3.6d3-catc(1,j))/(catc(1,j+1)-catc(1,j))
			cat6=cat3+(cat4-cat3)*(Ra*3.6d3-catc(1,j))/(catc(1,j+1)-catc(1,j))
			verti=Ra*(cat5+(cat6-cat5)*(Vloc-catc(i,1))/(catc(i+1,1)-catc(i,1)))
			deallocate(catc,catc5,catc6)		
		else
			verti=0.0d0
		endif
		bndcon_rain_m=horiz*dcos(ori(2))+verti*dsin(ori(2))
!
!
!
	end function bndcon_rain_m
!   _____________________________________________________________________________________
!
	pure function bndcon_hsurf_t(Ws,Wd,ori,Hs)
!
!
!   _____________________________________________________________________________________
!	computation of surface heat transfer coefficient [W/m2K] [Sharples, 1984].  The mois-
!	ture heat transfer coefficient is derived from this with the Lewis analogy.
		double precision,intent(in)::Ws,Wd,ori(:),Hs
		double precision::bndcon_hsurf_t,Vloc
!		---------------------------------------------------------------------------------
		if (Hs.LT.0.0d0) then 
			! no given value for hsurf
			if (Wd.GE.0.0d0) then
				if (dcos(Wd-ori(1)).GE.0.0d0) then
					Vloc=1.8d0*Ws+0.2d0
				else
					Vloc=0.4d0*Ws+1.7d0
				endif
			else
				Vloc=1.1d0*Ws+0.95d0
			endif
			bndcon_hsurf_t=1.7d0*Vloc+5.1d0
		else 
			! given value for hsurf (from subroutine "boundaryconditions" in mainprog)
			bndcon_hsurf_t=Hs
		endif
!
	end function bndcon_hsurf_t
!   _____________________________________________________________________________________
!
	pure function bndcon_evap_m(t,m,Ws,Wd,RH,Ta,ori,Hs)
!
!
!   _____________________________________________________________________________________
!	computation of vapour exchange [kg/m2s]
		double precision,intent(in)::t,m,Ws,Wd,ori(:),RH,Ta,Hs
		double precision::bndcon_evap_m
!		---------------------------------------------------------------------------------
		bndcon_evap_m=7.7d-9*bndcon_hsurf_t(Ws,Wd,ori,Hs)*(RH*thedyn_psat(Ta)-dexp(m/&
			&(rho_l*R_v*t))*thedyn_psat(t))
!
	end function bndcon_evap_m
!   _____________________________________________________________________________________
!
	pure function bndcon_devap_mdm(t,m,Ws,Wd,ori,Hs)
!
!
!   _____________________________________________________________________________________
!	computation of derivative of vapour exchange [kg/m2sPa]
		double precision,intent(in)::t,m,Ws,Wd,ori(:),Hs
		double precision::bndcon_devap_mdm
!		---------------------------------------------------------------------------------
		bndcon_devap_mdm=-7.7d-9*bndcon_hsurf_t(Ws,Wd,ori,Hs)*dexp(m/(rho_l*R_v*t))/&
			&(rho_l*R_V*t)*thedyn_psat(t)
!
	end function bndcon_devap_mdm
!   _____________________________________________________________________________________
!
	pure function bndcon_conv_t(t,Ta,Ws,Wd,ori,Hs)
!
!
!   _____________________________________________________________________________________
!	computation of convective heat exchange [W/m2]
		double precision,intent(in)::t,Ta,Ws,Wd,ori(:),Hs
		double precision::bndcon_conv_t
!		---------------------------------------------------------------------------------
		bndcon_conv_t=bndcon_hsurf_t(Ws,Wd,ori,Hs)*(Ta-t)
!
	end function bndcon_conv_t
!   _____________________________________________________________________________________
!
	pure function bndcon_dconv_tdt(Ws,Wd,ori,Hs)
!
!
!   _____________________________________________________________________________________
!	computation of derivative of convective heat exchange [W/m2K]
		double precision,intent(in)::Ws,Wd,ori(:),Hs
		double precision::bndcon_dconv_tdt
!		---------------------------------------------------------------------------------
		bndcon_dconv_tdt=-bndcon_hsurf_t(Ws,Wd,ori,Hs)
!
	end function bndcon_dconv_tdt
!   _____________________________________________________________________________________
!
	pure function bndcon_longw(t,Ta,ori,Cl,mat_nr)
!
!
!   _____________________________________________________________________________________
!	computation of long wave radiation exchange with sky and soil surface [W/m2]		
		integer,intent(in)::mat_nr
		double precision,intent(in)::t,Ta,ori(:),Cl
		double precision::bndcon_longw,Fwg,Fws
!		---------------------------------------------------------------------------------
		Fwg=ori(2)/1.8d2; Fws=1.0d0-Fwg
		bndcon_longw=5.67d-8*bndcon_epsil(mat_nr)*((Ta-(2.38d1-0.2025d0*(Ta-2.7315d2))*&
			&(1.0d0-0.87d0*Cl))**4.0d0*Fws+Ta**4.0d0*Fwg-t**4.0d0)
!
	end function bndcon_longw
!   _____________________________________________________________________________________
!
	pure function bndcon_dlongwdt(t,mat_nr)
!
!
!   _____________________________________________________________________________________
!	computation of derivative of long wave radiation exchange [W/m2K]
		integer,intent(in)::mat_nr
		double precision,intent(in)::t
		double precision::bndcon_dlongwdt
!		---------------------------------------------------------------------------------
		bndcon_dlongwdt=-5.67d-8*bndcon_epsil(mat_nr)*4.0d0*t**3.0d0
!
	end function bndcon_dlongwdt
!   _____________________________________________________________________________________
!
	pure function bndcon_epsil(mat_nr)
!
!
!   _____________________________________________________________________________________
!	emissivity for long wave radiation [-]
		integer,intent(in)::mat_nr
		double precision::bndcon_epsil
!		---------------------------------------------------------------------------------
		mat_nrerial:select case(mat_nr)
			case default
				bndcon_epsil=0.9d0
		end select mat_nrerial
!
	end function bndcon_epsil
!   _____________________________________________________________________________________
!
	pure function bndcon_evap_t(t,m,Ws,Wd,Ta,RH,ori,Hs)
!
!
!   _____________________________________________________________________________________
!	computation of latent and sensible heat transfer by vapour exchange	[W/m2]
		double precision,intent(in)::t,m,Ws,Wd,ori(:),Hs,RH,Ta
		double precision::bndcon_evap_t
!		---------------------------------------------------------------------------------
		bndcon_evap_t=7.7d-9*bndcon_hsurf_t(Ws,Wd,ori,Hs)*(RH*thedyn_psat(Ta)-&
			&dexp(m/(rho_l*R_v*t))*thedyn_psat(t))*(L_v+cap_v*t)
!
	end function bndcon_evap_t
!   _____________________________________________________________________________________
!
	pure function bndcon_devap_tdt(t,m,Ws,Wd,Ta,RH,ori,Hs)
!
!
!   _____________________________________________________________________________________
!	computation of derivative of heat transfer by vapour exchange [W/2K]
		double precision,intent(in)::t,m,Ws,Wd,ori(:),Hs,RH,Ta
		double precision::bndcon_devap_tdt
!		---------------------------------------------------------------------------------
		bndcon_devap_tdt=7.7d-9*bndcon_hsurf_t(Ws,Wd,ori,Hs)*((RH*thedyn_psat(Ta)-&
			&dexp(m/(rho_l*R_v*t))*thedyn_psat(t))*cap_v-dexp(m/(rho_l*R_v*t))*&
			&(thedyn_dpsatdt(t)-m/(rho_l*R_v*t*t)*thedyn_psat(t))*(L_v+cap_v*t))
!
	end function bndcon_devap_tdt
!   _____________________________________________________________________________________
!
	pure function bndcon_rain_t(Ta,Ra,Ws,Wd,ori,c)
!
!
!   _____________________________________________________________________________________
!	computation of sensible heat transfer by rain [W/m2]	
		double precision,intent(in)::Ta,Ra,Ws,Wd,ori(:),c(:)
		double precision::bndcon_rain_t
!		---------------------------------------------------------------------------------
		bndcon_rain_t=bndcon_rain_m(Ra,Ws,Wd,ori,c)*cap_l*Ta
!
	end function bndcon_rain_t
!   _____________________________________________________________________________________
!
	pure function bndcon_solar(tim,Dr,Df,ori,mat_nr,lati)
!
!
!   _____________________________________________________________________________________
!	computation of arriving short wave radiation [W/m2] from direct and diffuse radiation.
!	Direct radiation is considered the radiation on a plane orthogonal to the sun beam di-
!	rection. [Kittler, 1981]
		integer,intent(in)::mat_nr
		double precision,intent(in)::tim,Dr,Df,ori(:),lati
		double precision::bndcon_solar,decl,heig,azim,matdata(100,10)
		COMMON/data2/matdata
!		---------------------------------------------------------------------------------
		bndcon_solar=Df*(1.0d0+dcos(ori(2)))/2.0d0
		if (Dr.GT.0.0d0) then
			decl=23.45d0*dsin(3.6d2/3.1536d7*tim-109.0d0) ! declination
			heig=dasin(dsin(lati)*dsin(decl)-dcos(lati)*dcos(decl)*dcos(tim/2.4d2))
			if (heig.GT.0.0d0) then
				azim=dacos(dmin1(dmax1((dcos(lati)*dsin(decl)+dsin(lati)*dcos(decl)*&
					&dcos(tim/2.4d2))/dcos(heig),-1.0d0),1.0d0))
				if (dmod(tim,8.64d4)/8.64d4.GT.0.5d0) azim=360.0d0-azim
				if (dcos(ori(1)-azim).GT.0.0d0) bndcon_solar=bndcon_solar+Dr*dmax1(dcos(&
					&ori(2))*dsin(heig)+dsin(ori(2))*dcos(azim-ori(1))*dcos(heig),&
					&0.0d0)
			endif
		endif

		bndcon_solar=matdata(mat_nr,10)*bndcon_solar ! multiply by the material absorption for shortwave
!
	end function bndcon_solar
!    _____________________________________________________________________________________
!
      SUBROUTINE  EXTCOUPLING_INITIALISE(fnhostport,servername)      
! Common block to hold socket IDs
      integer hsocket(8)
      common/skt/hsocket
      character servername*21
      integer fnhostport

!	real extesprc(3)/0,0,0.1E+02/,extespri(3)/3*0/&
!         &,intesprc(3)/0,0,-0.1E+02/,intespri(3)/3*0/
!        COMMON/data2/extesprc,extespri,intesprc,intespri

! Open the sockets and get socket ID
      call socketclientopen(fnhostport,servername)
      write (*,*) " "
      write (*,*) "Socket opened. Socket id."
      write (*,*) hsocket(8)
!
	end subroutine EXTCOUPLING_INITIALISE
!	
!       ________________________________________________________________________________
! 
        SUBROUTINE  EXTCOUPLING_SENDnRECEIVE(temp0,pres0,mois0,nnm)
!      
!       Common block to hold socket IDs
        integer hsocket(8),hamj
        common/skt/hsocket
!	implicit none
	integer,intent(in)::nnm
	real testers/0/,testerss/0/      
	double precision,intent(in)::temp0(:),pres0(:),mois0(:)
	real extesprc(3),extespri(3),intesprc(3),intespri(3),&
			&tso,tsi,testcalc
        real Gextheatc,Gextheati,Gintheatc,Gintheati,Gextmoistc,Gextmoisti,&
			&Gintmoistc,Gintmoisti,mso,msi
        COMMON/data2/extesprc,extespri,intesprc,intespri

        real RHi,RHo,PVsati,PVsato
			
!	---------------------------------------------------------------------------------

! Send information to ESP-r
!      write (33,*) "sending"
      tso=temp0(1)-273.15
      tsi=temp0(nnm)-273.15

! Partial vapor pressure = RH * saturation vapor pressure
      mso=thedyn_phi(temp0(1),mois0(1))*thedyn_psat(temp0(1))
      msi=thedyn_phi(temp0(nnm),mois0(nnm))*thedyn_psat(temp0(nnm))

!      write (33,*) tso
!      write (*,*) "Client sending information to the Server." 
      call socketread(testers,hsocket(8))
!		WRITE(*,*)"read"
      call socketwrite(tsi,hsocket(8))
!		WRITE(*,*)"write"
      call socketread(testers,hsocket(8))
      call socketwrite(tso,hsocket(8))
      call socketread(testers,hsocket(8))
      call socketwrite(msi,hsocket(8))
      call socketread(testers,hsocket(8))
      call socketwrite(mso,hsocket(8))
      RHi=thedyn_phi(temp0(nnm),mois0(nnm))
      RHo=thedyn_phi(temp0(1),mois0(1))
      PVsati=thedyn_psat(temp0(nnm))
      PVsato=thedyn_psat(temp0(1))

! LATEST VERSION
!      write(33,'(E14.7,A,E14.7,A,E14.7,A,E14.7&
!     &,A,E14.7,A,E14.7,A,E14.7,A,E14.7)')&
!     & tsi,' ',tso,' ',msi,' ',mso&
!     & ,' ',RHi,' ',RHo,' ',PVsati,' ',PVsato

      do hamj=1,3
! Receive information from ESP-r
!      write (*,*) "receive"
      call socketread(extesprc(hamj),hsocket(8))
      call socketwrite(testerss,hsocket(8))
      call socketread(extespri(hamj),hsocket(8))
      call socketwrite(testerss,hsocket(8))
      call socketread(intesprc(hamj),hsocket(8))
      call socketwrite(testerss,hsocket(8))
      call socketread(intespri(hamj),hsocket(8))
      call socketwrite(testerss,hsocket(8))
!      write (*,*) "Data received from the server." 
      enddo

! obsolete...?
      testcalc=Gextheatc+Gextheati+Gintheatc+Gintheati      
!      write(34,'(F10.2,A,F10.2,A,F10.2,A,F10.2,A,F10.2)')&
!     & Gextheatc,' ',Gextheati,' ',Gintheatc,' ',Gintheati,' ',testcalc

!write values received from esp-r
!      write(34,'(F10.2,A,F10.2,A,F10.2,A,F10.2,A,&
!     &F10.2,A,F10.2,A,F10.2,A,F10.2,A,&
!     &F10.2,A,F10.2,A,F10.2,A,F10.2)')&

!    LATEST VERSION
!     write(34,'(E,A,E,A,E,A,E,A,&
!     &E,A,E,A,E,A,E,A,&
!     &E,A,E,A,E,A,E)')&
!     &extesprc(1),' ',extespri(1),' ',intesprc(1),' ',intespri(1),' ',&
!     &extesprc(2),' ',extespri(2),' ',intesprc(2),' ',intespri(2),' ',&
!     &extesprc(3),' ',extespri(3),' ',intesprc(3),' ',intespri(3)

      end subroutine    EXTCOUPLING_SENDnRECEIVE
!      
!	
!       ________________________________________________________________________________
! 
	

END MODULE HAM_boundary_conditions
!
!	More information on the boundary conditions can be found in:
!	
!	Janssen H, Blocken B, Carmeliet C. Simulation of moisture and heat transfer in building 
!	components, part 1: conservative modelling. International Journal of Numerical Methods 
!	in Engineering 2005; submitted for publication.
!
!	Blocken B. Wind-driven rain on buildings: measurements, numerical modelling and appli-
!	cations.  Doctoral dissertation 2004, Katholieke Universiteit Leuven, Belgium.
!
!	Sharples S.  Full scale measurements of convective energy losses from exterior building 
!	surfaces.  Building and Environment 1984; 19:31-39.
!
!	Kittler R. A universal calculation method for simple predetermination of natural radia-
!	tion on building surfaces and solar collectors. Building and Environment 1981; 16:177-
!	
!	Costola, D., Blocken, B., & Hensen, J. 2009. "External coupling between BES and HAM programs 
!	for whole-building simulation", Proceedings of the 11th IBPSA Building Simulation Conference, 
!	27-30 July, International Building Performance Simulation Association, Glasgow, pp. 316-323.
