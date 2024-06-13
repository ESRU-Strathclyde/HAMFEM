!
!		    		1D/2D/3D-FINITE ELEMENT MODEL FOR THE SIMULATION OF 
!			    	 HEAT, AIR AND MOISTURE TRANSFER IN BUILDING PARTS
!
!																 written by: Hans Janssen
!****************************************************************************************
!****************************************************************************************
!
									PROGRAM HAMFEM					
!
!****************************************************************************************
!****************************************************************************************
!																		updated: mar 2005
!
!
!		use msflib			! WINDOWS
!		use msimsl			! WINDOWS
!		use HAM_graphical_output	! WINDOWS
		use HAM_material_library
		use HAM_thermodynamics_library
		use HAM_boundary_conditions
		use FEM_shape_functions
		implicit none
!   
!   VARIABLE DECLARATION
!	_____________________________________________________________________________________
		integer::tpm,dim,nogpel(10),nem,nnm,nhbw,orig,bdnr,bdne,bdnn,itime,currentpathlength,&
		&inputfilelenght,cliquant,couplingmark
		integer::err,perc
		integer,allocatable::nod(:,:),nodinv(:),mater(:),bdele(:,:),bdnod(:,:),matout(:),&
			&bndout(:)
		double precision::orie(4),ttot,dt,dtold,dtavg(51),tmax,period,bdstep,facloc(2),lati
		double precision,allocatable::coord(:,:),wgp(:,:),Ni(:,:,:),dNi(:,:,:,:),&
			&Nimidq(:),temp0(:),pres0(:),mois0(:),thet0(:),valold(:,:),bdcoor(:,:),&
			&retval(:),noddra(:,:),nodsup(:,:),eledra(:,:,:),elesup(:,:,:),cumdra(:,:),&
			&cumsup(:,:,:)
		logical::lin,atmos,satur,graph
		character(100)::xts(4)
		integer::xtslenght(4)
		character(1) cliunit
		character(100)::inputfile,currentpath
		character(200)::inputfilenpath,promptoptions
!		character(10)::installpath 		! WINDOWS size is hardcoded
!		character(29)::installpath		! UNIX SOLARIS size is hardcoded
		character(12)::installpath		! LINUX size is hardcoded
		logical, allocatable::ret(:)
		COMMON/data1/lati


! Common block to hold socket IDs
      integer hsocket(8)/8*-1/
      common/skt/hsocket
! Port number. Hardcoded.      
      integer fnhostport
      character sfnhostport*5
! Servername
!      character servername*19 /"rietveld.bwk.tue.nl"/
      character servername*21 /"bps01-net0.bwk.tue.nl"/
! Data to exchange      
      real extheatc,extheati,intheatc,intheati,extmoistc,extmoisti,&
			&intmoistc,intmoisti
      COMMON/data3/extheatc,extheati,intheatc,intheati,extmoistc,extmoisti,&
			&intmoistc,intmoisti      
       


      
      	
!	atmos:		ATMOSpheric(T) or non-(F) boundary conditions to be used
!	bdcoor:		coordinates and location of boundary conditions
!	bdele:		ELEments with BounDary conditions, contains element number, location of 
!					boundary condition, direction of boundary face and the nodes involved
!	bdne:		Number of Elements with specified BounDary conditions
!	bdnn:		Number of Nodes with specified BounDary conditions
!	bdnod:		NODes with BounDary conditions, contains the node number and the location 
!					of the boundary condition
!	bdnr:		NumbeR of BounDary conditions to be read from input file
!	bdstep:		STEP interval for atmospheric climate BounDary value transition [s]
!	bndout:		OUTput of flows per BouNDary, indicates in which column the boundary flow 
!					for the boundary condition index have to be printed
!	cliunit		unit of the climate file: s(econds),h(ours), d(ays), m(onths) and y(ears) (usually 'y')
!	cliquant	amount of data per unit in the climate file ‘01’, ‘05’, ’10’, ‘20’ and ‘30’. (usually '01')
!	coord:		COORDinates for the node index [m]: x(1), y(2), z(3)
!	cumdra:		time-cumulated value of 'noddra' for output
!	cumsup:		time-cumulated value of 'elesup' for output
!	dim:		number of DIMensions
!	dNi:		values of Derivatives of shape functions N (idem as wgp)
!	dt:			Time step [s]
!	dtavg:		last 50 time steps and average of those [s]
!	dtold:		previous time step [s]
!	eledra:		DRAined quantities into component during time step (stored per element and 
!					per node per element, with structure similar to 'nod', and for the dif-
!					ferent quantities: heat = 1, air = 2, moisture = 3)
!	elesup:		SUPplied quantities to component during time step (stored per boundary ele-
!					ment and per node per element, with structure similar to 'bdele', and 
!					for the different quantities: heat = 1, air = 2, moisture = 3)
!	facloc:		location of origin on facade of building (vertical and horizontal distance
!					from lower left corner)
!	graph:		graphical output full(T) or time only(F)
!	itime	:	number of TIME step
!	lin:		LINear(T)/quadratic(F) elements
!	mater:		MATERial property for element index
!	matout:		OUTput of cumulated values per MATerial, indicates in which column the 
!					cumulated values for the material index have to be printed
!	mois0:		capillary pressures [Pa]
!	nem:		Number of Elements in Mesh
!	Ni:			values of shape functions N (idem as wgp)
!	Nimidq:		values of shape functions for center of quadrangular element (for graph. 
!					output)
!	nnm:		Number of Nodes in Mesh
!	nod:		connectivity matrix (which nodes belong to element index)
!	nodinv:		INVerse connectivity (node index belongs to which element)
!	noddra:		DRAined quantities into the component, stored per node and per quantity
!	nodsup:		SUPplied quantities to the component, stored per node and per quantity
!	nogpel:		NOdes and Gauss Points per ELement (see FEM_shape_functions.f90)
!	orie:		ORIEntation information (1) orientation (2) slope (3,4) gravity direction
!	orig:		node number of ORIGin (for graphical output)
!	period:		interval for numerical output [s]
!	pres0:		air pressures [Pa]
!	satur:		check and (de)activate(T) or not(F) for surface SATURation
!	sup:		SUPplied moisture by boundary conditions during time step (stored per boun-
!					dary element
!	ret:		surface saturated(T) or not(F) per node index
!	retval:		RETained surface moisture per node index [kg]
!	temp0:		temperatures at ttot [K]
!	thet0:		moisture contents [kg/m³]
!	tmax:		end time of calculation [s]
!	tpm:		quantities to calculate: Temp, Pres, Mois (1 for yes, 0 for no)
!	ttot:		TOTal Time [s]
!	valold:		previous values of temp (1), pres(2), mois(3) for 'restart'
!	wgp:		Weight of Gauss Points (see FEM_shape_functions.f90) for trian/qaud 
!					elements for the standard and refined scheme
!	xts:		extensions for material file(1), climate file(2), WDR catch ratio file(3), permeability(4)
!
!   MAIN PROGRAM
!	_____________________________________________________________________________________
!	_____________________________________________________________________________________
!
! PREPROCESSOR




!MODE: COMMAND LINE or FROM CVF INTERFACE
		!Get input file name from command line, current directory,define instaltion folder.
		CALL GET_COMMAND_ARGUMENT(1,inputfile)
		CALL GET_COMMAND_ARGUMENT(2,promptoptions)
		CALL GET_COMMAND_ARGUMENT(3,sfnhostport)
!               CALL GETARG(1,inputfile)		! Solaris
!		CALL GETARG(2,promptoptions)
!		CALL GETARG(3,sfnhostport)
		CALL GETCWD(currentpath)
		
		READ(sfnhostport,'(I5)') fnhostport
		WRITE(*,*) fnhostport
		
		currentpathlength=LEN_TRIM(currentpath)
!		installpath='D:\Hamfem\'			! WINDOWS hardcoded
!		installpath='$HOME/bwdcosto/hamfem/Hamfem/'	! UNIX SOLARIS hardcoded
		installpath='/opt/hamfem/'	! LINUX hardcoded
		inputfilelenght=LEN_TRIM(inputfile)-4


		CALL INPUT (ttot,tmax,dt,dtold,period,bdstep,dtavg,tpm,dim,nem,nnm,orie,bdnr,&
			&nogpel,lin,atmos,satur,xts,xtslenght,graph,facloc,inputfile, currentpath,&
			&currentpathlength,inputfilenpath,inputfilelenght,cliquant,cliunit)
		
		ALLOCATE (coord(nnm,dim),nod(nem,nogpel(2)),nodinv(nnm),mater(nem),bdcoor(bdnr,&
			&2*dim+1),temp0(nnm),pres0(nnm),mois0(nnm),thet0(nnm),valold(nnm,6),ret(nnm),&
			&retval(nnm),noddra(nnm,3),nodsup(nnm,3),eledra(nem,nogpel(2),3),elesup(nem,&
			&nogpel(4),3),cumdra(nnm,3),cumsup(nem,nogpel(4),3),wgp(maxval(nogpel(5:10)),&
			&6),Ni(maxval(nogpel(1:4)),maxval(nogpel(5:10)),6),dNi(maxval(nogpel(1:4)),&
			&maxval(nogpel(5:10)),dim,6),Nimidq(nogpel(2)))
		CALL PREPROS1 (dim,nogpel,nem,nnm,nod,nodinv,mater,nhbw,coord,orig,temp0,pres0,&
			&mois0,thet0,valold,ret,retval,bdnr,bdcoor,bdne,bdnn)
		ALLOCATE(bdele(bdne,3+nogpel(4)),bdnod(bdnn,2),matout(maxval(mater)),bndout(idint(&
			&maxval(bdcoor(:,2*dim+1)))))
 		CALL PREPROS2 (dim,nogpel,lin,nem,nnm,nod,coord,mater,bdnr,bdcoor,bdele,bdnod,&
			&matout,bndout)
		CALL GAUSS (nogpel,lin,dim,wgp,Ni,dNi,Nimidq)
		CALL THEDYN_INITIALISE(orie)
		CALL MATLIB_INITIALISE(installpath//'databases/',xts,xtslenght)
		CALL BNDCON_INITIALISE(installpath//'databases/',atmos,xts,xtslenght,facloc,bdstep,&
		&cliquant, cliunit)
!HAM to be generalized
		if (bdele(1,2).EQ.44.OR.bdele(1,2).EQ.45.OR.bdele(2,2).EQ.45) &
                  &CALL EXTCOUPLING_INITIALISE(fnhostport,servername)
		if (graph) then
!			CALL GRAFIC_INITIALISE()						! WINDOWS
!			CALL INIDRAW (dim,coord,nod,nodinv,nem,nogpel,nnm,mater,temp0,thet0)	! WINDOWS
		else
			perc=5 !percentage of time to report progress in text mode
		endif
!
! MAIN PROGRAM		
		itime=0 ! number of the time step
		CALL OUTPUT (ttot,period,nem,nogpel,wgp,Ni,dNi,coord,mater,dim,nod,temp0,pres0,&
			&mois0,thet0,cumdra,cumsup,bdnn,bdnod,bdne,bdele,matout,bndout,itime)
		DO WHILE (ttot.LT.tmax)
			itime=itime+1
!		WRITE(*,*)"bb"
			CALL COFEM (ttot,dt,dtold,dtavg,period,bdstep,nogpel,nem,nnm,nhbw,dim,orie,&
				&nod,nodinv,mater,coord,tpm,temp0,pres0,mois0,thet0,valold,ret,retval,&
				&atmos,satur,noddra,nodsup,eledra,elesup,cumdra,cumsup,wgp,Ni,dNi,bdne,&
				&bdele,bdnn,bdnod)
!		WRITE(*,*)"cc"
			CALL OUTPUT (ttot,period,nem,nogpel,wgp,Ni,dNi,coord,mater,dim,nod,temp0,&
				&pres0,mois0,thet0,cumdra,cumsup,bdnn,bdnod,bdne,bdele,matout,bndout,&
				&itime)
!		WRITE(*,*)"dd"
			if (graph) then
!				IF (mod(itime,1).eq.0) CALL TIMEDRAW (ttot,dim,coord,nod,mater,nodinv,nem,&
!				&nogpel,nnm,temp0,mois0,tpm,Nimidq,graph,lin)
			else
				IF (itime.EQ.1) write(*,*) "simulation started. elapsed time (%):"
				IF (ttot.GT.(tmax*perc/100)) then
						write(*,*) perc
						perc=perc+5
				endif
			endif
!		WRITE(*,*)"aa"
			if ((dmod(ttot,period).EQ.0.0d0)) then
!		WRITE(*,*)"bb"
!                if ((bdele(1,2).EQ.44).OR.(bdele(2,2).EQ.45))WRITE(*,*)"cc"
			if ((bdele(1,2).EQ.44).OR.(bdele(2,2).EQ.45)) CALL EXTCOUPLING_SENDnRECEIVE&
				&(temp0,pres0,mois0,nnm)

				endif
		ENDDO
		write(*,*) "simulation finished"

! POST-PROCESSOR
		DEALLOCATE (nod,nodinv,mater,coord,ret,retval,temp0,pres0,mois0,thet0,valold,&
			&noddra,nodsup,eledra,elesup,cumdra,cumsup,wgp,Ni,dNi,bdele,bdnod,bdcoor,&
			&Nimidq,matout,bndout)
		close(20); close(21); close(22); close(23); close(30); close(31); close(32) 
		close(50) 
		close(33);close(34);close(35);close(36)
!
CONTAINS
!
		include	'./Mainprog.f90'
		
END PROGRAM HAMFEM	


