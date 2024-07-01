
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
    SUBROUTINE INPUT	(ttot,tmax,dt,dtold,period,bdstep,dtavg,tpm,dim,nem,nnm,orie,&
						&bdnr,nogpel,lin,atmos,satur,xts,xtslenght,graph,facloc,inputfile,&
						&currentpath,currentpathlength,inputfilenpath,inputfilelenght,&
						&cliquant, cliunit)

!______________________________________________
		implicit none
		integer,intent(out)::tpm,dim,nem,nnm,bdnr,nogpel(:),cliquant,xtslenght(:)
		integer,intent(in)::currentpathlength,inputfilelenght
		double precision,intent(out)::ttot,tmax,dt,period,bdstep,dtavg(:),dtold,orie(:),&
			&facloc(:)
		character(100),intent(out)::xts(:)
		character(1),intent(out)::cliunit
		logical,intent(out)::lin,atmos,satur,graph
		integer::dum1,dum2,dum3
		character(currentpathlength)::currentpathtrim
		character(inputfilelenght)::inputfiletrim
		character(100),intent(in)::inputfile,currentpath
		character(200),intent(out)::inputfilenpath
		double precision::lati
		COMMON/data1/lati

!		---------------------------------------------------------------------------------
!	reading in of main input variables
		currentpathtrim=TRIM(currentpath)
		inputfilenpath=currentpathtrim//'/'//inputfile
		inputfiletrim=inputfile(1:inputfilelenght)

		write (*,*) "input file:"
		write (*,*) inputfilenpath
		write (*,*)
		write (*,*) "current folder:"
		write (*,*) currentpath
		open(10,file=inputfilenpath,action='read')

!		read(10,*); read(10,*) xts(1) ! extension for output
!		xts(1)=inputfiletrim
		read(10,*); read(10,*) ttot,tmax,dt,period; dtold=dt; dtavg=period; bdstep=period
		read(10,*); read(10,*) tpm ! quantities to calculate
		read(10,*); read(10,*) dim ! number of dimensions
		read(10,*); read(10,*) nem,nnm  ! number of elements and nodes in mesh
		read(10,*); read(10,*) dum1; 
		lin=.FALSE.; if (dum1.EQ.1) lin=.TRUE. ! linear/quanoddratic
		read(10,*); read(10,*) bdnr ! number of boundary conditions to read in
		read(10,*); read(10,*) dum1,dum2; 
		atmos=.FALSE.; if (dum1.EQ.1) atmos=.TRUE. ! atmospherical boundary conditions
		satur=.FALSE.; if (dum2.EQ.1) satur=.TRUE. ! surface saturation possible
		read(10,*); read(10,*) xts(1), xts(4) ! material databese file, permeability database file
		read(10,*); read(10,*) xts(2),lati, cliunit, cliquant  !climate file name, latitude, primary unit of data, quantity of data per unit
		read(10,*); read(10,*) xts(3), facloc(1:2) ! Wind-driven rain file, and location of the point in the facade
		read(10,*); read(10,*) orie(1:2) ! orientation of building part
		read(10,*); read(10,*) orie(3:4) ! orientation of building part (gravity)
		read(10,*); read(10,*) dum1; 
		graph=.FALSE.; if (dum1.EQ.1) graph=.TRUE. ! graphical output full or time only
		read(10,*); read(10,*) 
!	determination of number of nodes/gauss points per (boundary) element
		CALL NODESANDGAUSS (dim,lin,nogpel)
!   opening the files for output
		open(20,file='./'//inputfiletrim//'_nodetemperat.out',action='write')
		open(21,file='./'//inputfiletrim//'_nodeairpress.out',action='write')
		open(22,file='./'//inputfiletrim//'_nodecappress.out',action='write')
		open(23,file='./'//inputfiletrim//'_nodemoisture.out',action='write')
!
		open(30,file='./'//inputfiletrim//'_fluxheat.out',action='write')
		open(31,file='./'//inputfiletrim//'_fluxair.out',action='write')
		open(32,file='./'//inputfiletrim//'_fluxmois.out',action='write')
!
		open(50,file='./'//inputfiletrim//'_time.out',action='write')

! files to track external coupling variables.
		open(33,file='./'//inputfiletrim//'_send.txt',action='write')
		open(34,file='./'//inputfiletrim//'_rec.txt',action='write')
                open(35,file='./'//inputfiletrim//'_bd-out.txt',action='write')
                open(36,file='./'//inputfiletrim//'_bd-in.txt',action='write')
                write(33,*) "tso    tsi"
                write(34,*) "extheatc    extheati    intheatc    intheati"
                write(35,*) "ttot     Bxx1     dBxxdx1     Bxx     dBxxdx"
                write(36,*) "ttot     Bxx1     dBxxdx1     Bxx     dBxxdx"


!calculates the size of the file names
do dum3=1,4
		xtslenght(dum3)=LEN_TRIM(xts(dum3))
enddo
!
	end subroutine INPUT
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE OUTPUT	(ttot,period,nem,nogpel,wgp,Ni,dNi,coord,mater,dim,nod,temp0,&
						&pres0,mois0,thet0,cumdra,cumsup,bdnn,bdnod,bdne,bdele,matout,&
						&bndout,itime)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::nem,nogpel(:),dim,nod(:,:),bdnn,bdnod(:,:),bdne,bdele(:,:),&
			&mater(:),matout(:),bndout(:),itime
		double precision,intent(in)::ttot,period,wgp(:,:),Ni(:,:,:),dNi(:,:,:,:),&
			&coord(:,:),temp0(:),pres0(:),mois0(:),thet0(:),cumdra(:,:),cumsup(:,:,:)
		integer::iem,nne,ngp,ine,igp,nn,bdie,mat,imat,ibnd,bdin
		double precision::detjc,m,clima,time,tgau,pgau,mgau,wgau,dum
		double precision,allocatable::wg(:),N(:,:),dN(:,:,:),el(:,:),jc(:,:),jcinv(:,:),&
			&flowto(:,:),flowin(:,:),cumul(:,:),tval(:),pval(:),mval(:),vol(:)
		logical::trian
!		---------------------------------------------------------------------------------
		if (itime.EQ.0) then ! initialisation of output files
			write(20,'(a25,25000e25.15)') "coordinate x (m)",coord(:,1)
			if (dim.GE.2) write(20,'(a25,50000e25.15)') "coordinate y (m)",coord(:,2)
			if (dim.GE.3) write(20,'(a25,50000e25.15)') "coordinate z (m)",coord(:,3)
			write(20,'(2a25)') "time (s)","temperature (K)"
			write(21,'(a25,25000e25.15)') "coordinate x (m)",coord(:,1)
			if (dim.GE.2) write(21,'(a25,50000e25.15)') "coordinate y (m)",coord(:,2)
			if (dim.GE.3) write(21,'(a25,50000e25.15)') "coordinate z (m)",coord(:,3)
			write(21,'(2a25)') "time (s)","air pressure (Pa)"
			write(22,'(a25,25000e25.15)') "coordinate x (m)",coord(:,1)
			if (dim.GE.2) write(22,'(a25,50000e25.15)') "coordinate y (m)",coord(:,2)
			if (dim.GE.3) write(22,'(a25,50000e25.15)') "coordinate z (m)",coord(:,3)
			write(22,'(2a25)') "time (s)","capil. pressure (Pa)"
			write(23,'(a25,25000e25.15)') "coordinate x (m)",coord(:,1)
			if (dim.GE.2) write(23,'(a25,50000e25.15)') "coordinate y (m)",coord(:,2)
			if (dim.GE.3) write(23,'(a25,50000e25.15)') "coordinate z (m)",coord(:,3)
			write(23,'(2a25)') "time (s)","moi. content (kg/m³)"
			write(30,'(a35)') "cumulated temperatures + heat flows"
			write(31,'(a35)') "cumulated air pressures + air flows"
			write(32,'(a35)') "cumulated moisture contents + flows"
			write(30,'(a18,a7)',advance='no') "time (s)",""
			write(31,'(a18,a7)',advance='no') "time (s)",""
			write(32,'(a18,a7)',advance='no') "time (s)",""
			do iem=1,size(matout)
				if (matout(iem).NE.0) then
					write(30,'(a15,i5,a5)',advance='no') "material ",iem,""
					write(31,'(a15,i5,a5)',advance='no') "material ",iem,""
					write(32,'(a15,i5,a5)',advance='no') "material ",iem,""
				endif
			enddo
			do iem=1,size(bndout)
				if (bndout(iem).NE.0) then
					write(30,'(a19,i5,a1)',advance='no') "flow to boundary ",iem,""
					write(31,'(a19,i5,a1)',advance='no') "flow to boundary ",iem,""
					write(32,'(a19,i5,a1)',advance='no') "flow to boundary ",iem,""
				endif
			enddo
			do iem=1,size(bndout)
				if (bndout(iem).NE.0) then
					write(30,'(a20,i5)',advance='no') "flow from boundary ",iem
					write(31,'(a20,i5)',advance='no') "flow from boundary ",iem
					write(32,'(a20,i5)',advance='no') "flow from boundary ",iem
				endif
			enddo
			write(30,*); write(31,*); write(32,*)		
			write(50,'(4a25)') "time step start   ","time step length  ",&
				&"avg time step length","required iterations" 
		endif			
		if (dmod(ttot,period).EQ.0.0d0) then
!	output independent variables
			if (tpm/100.EQ.1) write(20,'(50000e25.15)') ttot,temp0
			if (mod(tpm,100)/10.EQ.1) write(21,'(50000e25.15)') ttot,pres0
			if (mod(tpm,10).EQ.1) then
				write(22,'(50000e25.15)') ttot,mois0
				write(23,'(50000e25.15)') ttot,thet0
			endif
!	output derived quantities (cumulative + flows)
			allocate(cumul(maxval(matout),3),vol(maxval(matout)),flowto(maxval(bndout),3),&
				&flowin(maxval(bndout),3))
			cumul=0.0d0; flowto=0.0d0; flowin=0.0d0; vol=0.0d0
			do iem=1,nem ! calculation cumulative temp, pres, mois, thet
				nne=nogpel(2); ngp=nogpel(8); trian=.FALSE.
				if (nod(iem,nne).EQ.0) then
					nne=nogpel(1); ngp=nogpel(7); trian=.TRUE.
				endif
				allocate (wg(ngp),N(nne,ngp),dN(nne,ngp,dim),el(nne,dim),jc(dim,dim),&
					&jcinv(dim,dim),tval(nne),pval(nne),mval(nne))
				if (trian) then
					wg=wgp(:,3); N=Ni(:,:,3); dN=dNi(:,:,:,3)
				else
					wg=wgp(:,4); N=Ni(:,:,4); dN=dNi(:,:,:,4)
				endif
				do ine=1,nne
					nn=nod(iem,ine); el(ine,:)=coord(nn,:)
					tval(ine)=temp0(nn); pval(ine)=pres0(nn); mval(ine)=mois0(nn)
				enddo
				mat=mater(iem); imat=matout(mat)
				do igp=1,ngp
					call jacobian(dim,igp,el,dN,jc,jcinv,detjc,.FALSE.) 
					detjc=detjc*wg(igp)
					vol(imat)=vol(imat)+detjc
					tgau=dmax1(dmin1(dot_product(N(:,igp),tval),4.0d2),2.0d2)
					pgau=dmax1(dmin1(dot_product(N(:,igp),pval),1.0d3),-1.0d3)
					mgau=dmax1(dmin1(dot_product(N(:,igp),mval),-1.0d-5),-1.0d15)
					cumul(imat,1)=cumul(imat,1)+tgau*detjc
					cumul(imat,2)=cumul(imat,2)+pgau*detjc
					call matlib_moisture_storage(mat,tgau,mgau,wgau,dum)
					cumul(imat,3)=cumul(imat,3)+wgau*detjc
				enddo
				deallocate(wg,N,dN,el,jc,jcinv,tval,pval,mval)
			enddo
			do bdie=1,bdne ! calculation cumulated flows
				ibnd=bndout(bdele(bdie,2))
				nne=nogpel(4)
				if (bdele(bdie,3+nne).EQ.0) nne=nogpel(3)			
				do ine=1,nne
					flowto(ibnd,:)=flowto(ibnd,:)+cumsup(bdie,ine,:)
				enddo
			enddo
			do bdin=1,bdnn
				flowin(bdnod(bdin,2),:)=flowin(bdnod(bdin,2),:)+cumdra(bdnod(bdin,1),:)
			enddo
			if (tpm/100.EQ.1) write(30,'(2500e25.15)') ttot,cumul(:,1),flowto(:,1),&
				&flowin(:,1)
			if (mod(tpm,100)/10.EQ.1) write(31,'(2500e25.15)') ttot,cumul(:,2),&
				&flowto(:,2),flowin(:,2)
			if (mod(tpm,10).EQ.1) write(32,'(2500e25.15)') ttot,cumul(:,3),flowto(:,3),&
				&flowin(:,3)
			deallocate(cumul,flowto,flowin,vol)
		endif
!
	end subroutine OUTPUT
!	_____________________________________________________________________________________
!
!
!
    SUBROUTINE PREPROS1	(dim,nogpel,nem,nnm,nod,nodinv,mater,nhbw,coord,orig,temp0,pres0,&
						&mois0,thet0,valold,ret,retval,bdnr,bdcoor,bdne,bdnn)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::dim,nogpel(:),nem,nnm,bdnr
		integer,intent(out)::nod(:,:),nodinv(:),mater(:),nhbw,bdne,bdnn,orig
		double precision,intent(out)::coord(:,:),temp0(:),pres0(:),mois0(:),thet0(:),&
			&valold(:,:),retval(:),bdcoor(:,:)
		logical,intent(out)::ret(:)
		integer::iem,ine,nne,nw,inm,idim,i,nn,nnebdloc,nnebd,bddir
		double precision::dum
		logical::constant,origin,bdelem
!		---------------------------------------------------------------------------------
!	reading global coordinates of nodes and initialisation of calculation variables
		coord=0.0d0; constant=.FALSE.
		do inm=1,nnm
			if (.NOT.constant) then
				read(10,*) coord(inm,1:dim),temp0(inm),pres0(inm),mois0(inm)
			else
				read(10,*) coord(inm,1:dim)
			endif
			if (temp0(inm).EQ.-1.0d0) constant=.TRUE.
			if (constant) then
				temp0(inm)=temp0(inm-1); pres0(inm)=pres0(inm-1); mois0(inm)=mois0(inm-1)
			endif
			origin=.TRUE.
			do idim=1,dim
				if (coord(inm,idim).NE.0.0d0) origin=.FALSE.
			enddo
			if (origin) orig=inm
		end do
		valold(:,1)=temp0; valold(:,2)=pres0; valold(:,3)=mois0
		ret=.FALSE.; retval=-1.0d0; dtavg=period		
		read(10,*)
!	reading connectivity matrix & calculation of half-of-band-width
		nhbw=0; nodinv=0
		do iem=1,nem
			nne=nogpel(2)
			read(10,*) nod(iem,1:nne),mater(iem)
			if (nod(iem,nne).EQ.0) nne=nogpel(1)
			do ine=1,nne
				nn=nod(iem,ine)
				if (nodinv(nn).EQ.0) nodinv(nn)=iem
			enddo
			nw=maxval(nod(iem,1:nne))-minval(nod(iem,1:nne))+1
			if (nhbw.LT.nw) nhbw=nw
		enddo
		read(10,*); read(10,*)
!	calculation of moisture contents
		do inm=1,nnm
			call matlib_moisture_storage(mater(nodinv(inm)),temp0(inm),mois0(inm),&
				&thet0(inm),dum)
		enddo
!	reading data on boundary conditions	
		do i=1,bdnr
			read(10,*) bdcoor(i,1:2*dim+1)			
		enddo
		close(10)
		bdne=0; bdnn=0 ! total number of elements/nodes at specified boundaries		 
		do i=1,bdnr ! determine bdne and bdnn
			do idim=1,dim
				if (bdcoor(i,idim).EQ.bdcoor(i,idim+dim)) bddir=idim 
			enddo ! bddir = 1/2/3 if boundary perpendicular to x/y/z-axis
			do iem=1,nem ! determines nnebdloc (number of nodes of the current el. that 
						 ! are situated at the current bound.loc.)
				nnebdloc=0
				nne=nogpel(2); nnebd=nogpel(4)
				if (nod(iem,nne).EQ.0) then
					nne=nogpel(1); nnebd=nogpel(3)
				endif
				do ine=1,nne
					nn=nod(iem,ine); bdelem=.TRUE.
					do idim=1,dim
						if (.NOT.(coord(nn,idim).GE.bdcoor(i,idim).AND.coord(&
							&nn,idim).LE.bdcoor(i,dim+idim))) bdelem=.FALSE.
					enddo
					if (bdelem) nnebdloc=nnebdloc+1 
				enddo
				if (nnebdloc.EQ.nnebd) bdne=bdne+1	
			enddo
			do inm=1,nnm
				bdelem=.TRUE.
				do idim=1,dim
					if (.NOT.(coord(inm,idim).GE.bdcoor(i,idim).AND.coord(&
						&inm,idim).LE.bdcoor(i,dim+idim))) bdelem=.FALSE.
				enddo
				if (bdelem) bdnn=bdnn+1
			enddo
		enddo
!
	end subroutine PREPROS1
!	_____________________________________________________________________________________
!
!
!
    SUBROUTINE PREPROS2	(dim,nogpel,lin,nem,nnm,nod,coord,mater,bdnr,bdcoor,bdele,bdnod,&
						&matout,bndout)
!	_____________________________________________________________________________________
!	Computation of the previously allocated matrices bdele and bdnod
		implicit none
		integer,intent(in)::dim,nogpel(:),nem,nnm,bdnr,nod(:,:),mater(:)
		double precision,intent(in)::coord(:,:),bdcoor(:,:)
		logical,intent(in)::lin
		integer,intent(out)::bdele(:,:),bdnod(:,:),matout(:),bndout(:)
		integer::iem,ine,jne,idim,jdim,bdneloc,bdnnloc,i,nn,nnebdloc,dum,nne,nnebd,bddir,&
			&inm,imat,prev,ibnd
		integer,allocatable::bdeleloc(:),bdelelocbis(:)
		double precision::atanmin
		double precision,allocatable::coordloc(:,:)
		logical::trian,bdelem,swap
		character(1)::quan
		character(4)::kind
		character(5)::clima
		character(7)::ngcr
!		---------------------------------------------------------------------------------
!	Constructs the matrices bdele (with data for each element at each of the specified boun-
!	daries) and bdnod (with data for the nodes of each bd element that are situated on the 
!	spec. bds).
		allocate(bdeleloc(nogpel(4))); bdele=0; bdnod=0; bdneloc=0; bdnnloc=0 
		do i=1,bdnr
			do iem=1,nem ! determines nnebdloc (number of nodes of the current el. that are 
					! at the current bound.loc.) and bdeleloc (node numbers of the nodes of
					! the current el. that are at the current bound.loc.)
				nnebdloc=0; bdeleloc=0
				nne=nogpel(2); nnebd=nogpel(4); trian=.FALSE.
				if (nod(iem,nne).EQ.0) then
					nne=nogpel(1); nnebd=nogpel(3); trian=.TRUE.
				endif
				do idim=1,dim
					if (bdcoor(i,idim).EQ.bdcoor(i,idim+dim)) bddir=idim
				enddo
				do ine=1,nne
					nn=nod(iem,ine); bdelem=.TRUE.
					do idim=1,dim
						if (.NOT.(coord(nn,idim).GE.bdcoor(i,idim).AND.coord(nn,idim)&
							&.LE.bdcoor(i,dim+idim))) bdelem=.FALSE.
					enddo
					if (bdelem) then
						nnebdloc=nnebdloc+1	! nnebdloc is number of nodes of the current el. 
								! that are situated at the current boundary location
						bdeleloc(nnebdloc)=ine
					endif
				enddo
				if (nnebdloc.EQ.nnebd) then	 
					if (dim.EQ.3) then ! to range nodes in correct order
						if (lin) then
							if (.NOT.trian) then
								if (bdeleloc(2).NE.2.AND.bdeleloc(2).NE.6) then
										! rearrange nodes 3-4
									 dum=bdeleloc(3)
									 bdeleloc(3)=bdeleloc(4)
									 bdeleloc(4)=dum
								endif
							endif
						else
							if (trian) then
								if (bdeleloc(3).NE.3) then ! rearrange nodes 5-6
									 dum=bdeleloc(5)
									 bdeleloc(5)=bdeleloc(6)
									 bdeleloc(6)=dum
								endif
							else
								if (bdeleloc(6).NE.10.AND.bdeleloc(6).NE.18) then 
										! rearrange nodes 3-4 & 6-7-8
									dum=bdeleloc(3)
									bdeleloc(3)=bdeleloc(4)
									bdeleloc(4)=dum
									dum=bdeleloc(6)
									bdeleloc(6)=bdeleloc(7)
									bdeleloc(7)=bdeleloc(8)
									bdeleloc(8)=dum
								endif
							endif
						endif
					endif						
					bdneloc=bdneloc+1; bdele(bdneloc,1)=iem
					bdele(bdneloc,2)=bdcoor(i,2*dim+1)
					bdele(bdneloc,3)=bddir 
					bdele(bdneloc,4:3+nnebd)=bdeleloc(:)
				endif
			enddo
			do inm=1,nnm
				bdelem=.TRUE.
				do idim=1,dim
					if (.NOT.(coord(inm,idim).GE.bdcoor(i,idim).AND.coord(&
						&inm,idim).LE.bdcoor(i,dim+idim))) bdelem=.FALSE.
				enddo
				if (bdelem) then
					bdnnloc=bdnnloc+1; bdnod(bdnnloc,1)=inm
					bdnod(bdnnloc,2)=bdcoor(i,2*dim+1)
				endif
			enddo
		enddo
		deallocate(bdeleloc)
!
		matout=0; imat=1 ! determination of column number for output cumulated values
		do iem=1,nem
			if (matout(mater(iem)).EQ.0) then
				matout(mater(iem))=imat; imat=imat+1
			endif
		enddo
		do ine=1,size(matout) ! ranging in correct order
			prev=1				
			do jne=1,size(matout)
				if (matout(jne).NE.0) then
					if (matout(prev).GT.matout(jne)) then
						iem=matout(prev); matout(prev)=matout(jne); matout(jne)=iem
					endif
					prev=jne
				endif
			enddo
		enddo
!
		bndout=0; ibnd=1 ! determination of column number for output boundary flows
		do iem=1,bdnr
			if (bndout(idint(bdcoor(iem,2*dim+1))).EQ.0) then
				bndout(idint(bdcoor(iem,2*dim+1)))=ibnd; ibnd=ibnd+1
			endif
		enddo
		do ine=1,size(bndout) ! ranging in correct order
			prev=1				
			do jne=1,size(bndout)
				if (bndout(jne).NE.0) then
					if (bndout(prev).GT.bndout(jne)) then
						iem=bndout(prev); bndout(prev)=bndout(jne); bndout(jne)=iem
					endif
					prev=jne
				endif
			enddo
		enddo
!
	end subroutine PREPROS2
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE COFEM	(ttot,dt,dtold,dtavg,period,bdstep,nogpel,nem,nnm,nhbw,dim,orie,&
						&nod,nodinv,mater,coord,tpm,temp0,pres0,mois0,thet0,valold,ret,&
						&retval,atmos,satur,noddra,nodsup,eledra,elesup,cumdra,cumsup,&
						&wgp,Ni,dNi,bdne,bdele,bdnn,bdnod)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::nogpel(:),nem,nnm,nhbw,dim,nod(:,:),nodinv(:),mater(:),tpm,&
			&bdne,bdele(:,:),bdnn,bdnod(:,:)
		double precision,intent(in)::period,bdstep,coord(:,:),wgp(:,:),Ni(:,:,:),&
			&dNi(:,:,:,:),orie(:)
		logical,intent(in)::atmos,satur
		double precision,intent(inout)::ttot,dt,dtold,dtavg(:),temp0(:),pres0(:),&
			&mois0(:),thet0(:),valold(:,:),retval(:),noddra(:,:),nodsup(:,:),&
			&eledra(:,:,:),elesup(:,:,:),cumdra(:,:),cumsup(:,:,:)
		logical,intent(inout)::ret(:)
		integer::iter,maxiter,moditer,tpmloc,inm,i,info
		double precision::unbal,tolunbal,devia,toldevia,dum
		double precision,allocatable::temp(:),pres(:),mois(:),thet(:),valnew(:),&
			&Kot(:,:),Kop(:,:),Kom(:,:),Kn(:,:),dF(:),Fe(:)
		logical::convt,convm,convp,diver,overl,retact,retdeact,Kassem
!		---------------------------------------------------------------------------------
!	iterative solution of the system of equations for each time step
 		allocate(temp(nnm),pres(nnm),mois(nnm),thet(nnm),Kot(nhbw,nnm),Kop(nhbw,nnm),&
			&Kom(nhbw,nnm),Kn(nhbw,nnm),dF(nnm),Fe(nnm),valnew(nnm))
		convt=.FALSE.; convp=.FALSE.; convm=.FALSE.; maxiter=30; moditer=3
		tolunbal=1.0d-6; toldevia=1.0d-5; diver=.FALSE.
!		
		do while ((.NOT.convt.OR..NOT.convp.OR..NOT.convm).AND.dt.GT.1.0d-18)
			iter=0; temp=temp0; pres=pres0; mois=mois0; thet=thet0
			convt=.FALSE.; if (tpm/100.NE.1) convt=.TRUE.
			convp=.FALSE.; if (mod(tpm,100)/10.NE.1) convp=.TRUE.
			convm=.FALSE.; if (mod(tpm,10).NE.1) convm=.TRUE.
			do while ((.NOT.convt.OR..NOT.convp.OR..NOT.convm).AND.iter.LE.maxiter)
				iter=iter+1; Kassem=.FALSE.; if (iter.LE.moditer) Kassem=.TRUE.
				if (.NOT.convp) then	
					tpmloc=010
					CALL FEM (ttot,dt,nogpel,nem,nnm,nhbw,dim,nod,mater,coord,orie,temp0,&
						&temp,pres,mois0,mois,tpm,tpmloc,ret,retval,atmos,satur,wgp,Ni,dNi,&
						&bdne,bdele,Kn,Fe,dF,noddra,nodsup,eledra,elesup,unbal,overl,Kassem)
					if (iter.GT.1.AND.overl) exit

!					Original version using IMSL
!					if (Kassem) CALL DLFTQS(nnm,Kn(1:nhbw,1:nnm),nhbw,nhbw-1,Kop(1:nhbw,&
!						&1:nnm),nhbw)
!					CALL DLFSQS(nnm,Kop,nhbw,nhbw-1,dF,valnew)
!					CALL UPDATE(valnew,pres,devia,nnm)

!					LAPACK version
					if (Kassem) THEN
					CALL DPBTRF('U',nnm,nhbw-1,Kn(1:nhbw,1:nnm),nhbw,info)
					Kop=Kn(1:nhbw,1:nnm)
					ENDIF
					CALL DPBTRS('U',nnm,nhbw-1,1,Kop,nhbw,df,nnm,info)
					CALL UPDATE(df,pres,devia,nnm)
!					End LAPACK version


					if (unbal.LE.tolunbal.OR.devia.LE.toldevia) convp=.TRUE.
				endif
				if (.NOT.convm) then	
					tpmloc=001
					CALL FEM (ttot,dt,nogpel,nem,nnm,nhbw,dim,nod,mater,coord,orie,temp0,&
						&temp,pres,mois0,mois,tpm,tpmloc,ret,retval,atmos,satur,wgp,Ni,dNi,&
						&bdne,bdele,Kn,Fe,dF,noddra,nodsup,eledra,elesup,unbal,overl,Kassem)
					if (iter.GT.1.AND.overl) exit

!					Original version using IMSL
!					if (Kassem) 
!					CALL DLFTQS(nnm,Kn(1:nhbw,1:nnm),nhbw,nhbw-1,Kom(1:nhbw,&
!						&1:nnm),nhbw)
!					CALL DLFSQS(nnm,Kom,nhbw,nhbw-1,dF,valnew)
!					CALL UPDATEM(valnew,mois,devia,nnm,thet,nodinv,mater,temp)

!					LAPACK version
					if (Kassem) then
					CALL DPBTRF('U',nnm,nhbw-1,Kn(1:nhbw,1:nnm),nhbw,info)
					Kom=Kn(1:nhbw,1:nnm)
					ENDIF
					CALL DPBTRS('U',nnm,nhbw-1,1,Kom,nhbw,df,nnm,info)
					CALL UPDATEM(df,mois,devia,nnm,thet,nodinv,mater,temp)
!					End LAPACK version

					if (unbal.LE.tolunbal.OR.devia.LE.toldevia) convm=.TRUE.
				endif
				if (.NOT.convt) then	
					tpmloc=100
					CALL FEM (ttot,dt,nogpel,nem,nnm,nhbw,dim,nod,mater,coord,orie,temp0,&
						&temp,pres,mois0,mois,tpm,tpmloc,ret,retval,atmos,satur,wgp,Ni,dNi,&
						&bdne,bdele,Kn,Fe,dF,noddra,nodsup,eledra,elesup,unbal,overl,Kassem)
					if (iter.GT.1.AND.overl) exit

!					Original version using IMSL
!					if (Kassem) CALL DLFTQS(nnm,Kn(1:nhbw,1:nnm),nhbw,nhbw-1,Kot(1:nhbw,&
!						&1:nnm),nhbw)
!					CALL DLFSQS(nnm,Kot,nhbw,nhbw-1,dF,valnew)
!					CALL UPDATE(valnew,temp,devia,nnm)

!					LAPACK version
					if (Kassem) THEN
					CALL DPBTRF('U',nnm,nhbw-1,Kn(1:nhbw,1:nnm),nhbw,info)
					Kot=Kn(1:nhbw,1:nnm)
					ENDIF
					CALL DPBTRS('U',nnm,nhbw-1,1,Kot,nhbw,df,nnm,info)
					CALL UPDATE(df,temp,devia,nnm)
!					End LAPACK version

					if (unbal.LE.tolunbal.OR.devia.LE.toldevia) convt=.TRUE.
				endif
			enddo
			retact=.FALSE.; retdeact=.FALSE.
			if (mod(tpm,10).EQ.1.AND.satur.AND.convm.AND.convt.AND.convp) then
				do i=1,bdnn
					inm=bdnod(i,1)
					if (.NOT.ret(inm)) then
						if (mois(inm).GT.1.0d3.AND.dt.GT.1.0d-7) retact=.TRUE. 
					else 
						if (retval(inm)+nodsup(inm,3)-noddra(inm,3).LT.-1.0d-8.AND.dt.GT.&
							&bdstep/1.0d3) then
							retdeact=.TRUE. 
						else
							retval(inm)=retval(inm)+nodsup(inm,3)-noddra(inm,3)
							if (retval(inm).GT.0.0d0) retval(inm)=0.0d0
						endif
					endif
				enddo
				if (retact.OR.retdeact) then
					convm=.FALSE.
				endif
			endif
			if (.NOT.convt.OR..NOT.convp.OR..NOT.convm) then ! if no convergence is reached
				if (.NOT.retact.AND..NOT.retdeact) then
					dt=dt/2.0d0; diver=.TRUE.
				else 
					dt=dt/4.0d0
				endif
				iter=0
			else
				CALL CONVERGED (ttot,dt,dtold,dtavg,period,bdstep,nnm,nodinv,mater,iter,&
					&maxiter,diver,satur,temp0,pres0,mois0,thet0,temp,pres,mois,thet,&
					&valold,ret,retval,eledra,elesup,cumdra,cumsup,bdnn,bdnod)
			endif
		enddo
		if (dt.LE.1.0d-18.OR.dtavg(51).LT.bdstep/1.0d16) then ! checking for too small steps
			if (dt.LE.1.0d-18) then
				print*,ttot,"MAYDAY dt too small",dt
				write(50,'(e25.15,a30)') ttot,"MAYDAY dt too small"
			else
				print*,ttot,"MAYDAY dtavg too small",dtavg(51)
				write(50,'(e25.15,a30)') ttot,"MAYDAY dtavg too small" 
			endif			! jump to next bdstep with the previous values tempold, etc.
			ttot=dint((ttot+bdstep)/bdstep)*bdstep; dt=bdstep/1.0d1; dtavg=bdstep
			temp0=valold(:,1); pres0=valold(:,2); mois0=valold(:,3); thet0=valold(:,4)
			ret=.FALSE.; retval=-1.0d0
			do inm=1,bdnn
				if (valold(bdnod(inm,1),5).EQ.1.0d0) then
					ret(bdnod(inm,1))=.TRUE.; retval(bdnod(inm,1))=valold(bdnod(inm,1),6)
				endif
			enddo
		endif
 		deallocate(temp,pres,mois,thet,Kot,Kop,Kom,Kn,dF,Fe,valnew)
!
	end subroutine cofem
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE UPDATE	(valnew,val,devia,nnm)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::nnm
		double precision,intent(inout)::valnew(:),val(:),devia
!		---------------------------------------------------------------------------------
		val=val+valnew; devia=dot_product(val,val)
		if (devia.LT.1.0d-50) then
			devia=1.0d0
		else
			devia=dsqrt(dot_product(valnew,valnew)/devia)
		endif
!
	end subroutine UPDATE
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE UPDATEM	(valnew,mois,devia,nnm,thet,nodinv,mater,temp)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::nnm,nodinv(:),mater(:)
		double precision,intent(inout)::valnew(:),mois(:),devia,thet(:),temp(:)
		integer::inm
		double precision::deviapc,deviath,t,m,w,dum
!		---------------------------------------------------------------------------------
		mois=mois+valnew; deviapc=dot_product(mois,mois)
		if (deviapc.LT.1.0d-50) then
			deviapc=1.0d0
		else
			deviapc=dsqrt(dot_product(valnew,valnew)/deviapc)
		endif
		valnew=thet
		do inm=1,nnm
			t=dmax1(dmin1(temp(inm),4.0d2),2.0d2)
			m=dmax1(dmin1(mois(inm),-1.0d-5),-1.0d15)
			call matlib_moisture_storage(mater(nodinv(inm)),t,m,thet(inm),dum)
		enddo
		deviath=dot_product(thet,thet)
		if (deviath.LT.1.0d-50) then
			deviath=1.0d0
		else
			valnew=thet-valnew
			deviath=dsqrt(dot_product(valnew,valnew)/deviath)
		endif
		devia=dmin1(deviapc,deviath)
!
	end subroutine UPDATEM
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE CONVERGED	(ttot,dt,dtold,dtavg,period,bdstep,nnm,nodinv,mater,iter,&
							&maxiter,diver,satur,temp0,pres0,mois0,thet0,temp,pres,&
							&mois,thet,valold,ret,retval,eledra,elesup,cumdra,cumsup,&
							&bdnn,bdnod)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::maxiter,iter,nnm,nodinv(:),mater(:),bdnn,bdnod(:,:)
		double precision,intent(in)::temp(:),pres(:),mois(:),thet(:),period,bdstep
		logical,intent(in)::diver,satur
		double precision,intent(inout)::ttot,dt,dtold,dtavg(:),temp0(:),pres0(:),&	
			&mois0(:),thet0(:),valold(:,:),retval(:),eledra(:,:,:),elesup(:,:,:),&
			&cumdra(:,:),cumsup(:,:,:)			
		logical,intent(inout)::ret(:)
		integer::i,inm
		double precision::dtfix,tfix,t,m,dum
		logical::retact,retdeact
!		---------------------------------------------------------------------------------
!	intermediate results 
		if (dmod(ttot+dt,bdstep).EQ.0.0d0) then
			valold(:,1)=temp0; valold(:,2)=pres0; valold(:,3)=mois0
			valold(:,4)=thet0; valold(:,5)=0.0d0; valold(:,6)=retval
			do inm=1,bdnn
				if (ret(bdnod(inm,1))) valold(inm,5)=1.0d0
			enddo
		endif
		temp0=temp; pres0=pres; mois0=mois; thet0=thet
!		
!	saturation conditions ?
		retact=.FALSE.; retdeact=.FALSE.
		if (satur) then
			do i=1,bdnn
				inm=bdnod(i,1)
				if (ret(inm)) then
					if (retval(inm).LT.-1.0d-10) then
						retdeact=.TRUE.; ret(inm)=.FALSE.; retval(inm)=-1.0d0
					endif
				else
					if (mois0(inm).GT.0.0d0) then
						retact=.TRUE.; ret(inm)=.TRUE.; retval(inm)=0.0d0
					endif
				endif
			enddo
		endif
!
!	cumulative heat and mass flows
		if (dmod(ttot,period).EQ.0.0d0) then
			cumdra=0.0d0; cumsup=0.0d0 
		endif
		cumdra=cumdra+noddra; cumsup=cumsup+elesup
!
!	new time step determination
		ttot=ttot+dt; dtold=dt; dtavg(51)=0.0d0
		dtavg(1:49)=dtavg(2:50); dtavg(50)=dt
		dtavg(51)=sum(dtavg(1:50))/5.0d1
		if (dmod(ttot,period).EQ.0.0d0) dtavg=bdstep
		dtfix=dmin1(period,bdstep)
		if (diver) then
			dt=dt
		else				
			dt=dmin1(dt*dmin1(0.5d0*real(maxiter)/real(iter),2.0d0),dtfix)
			if (retact.AND.dt.LT.1.0d3) dt=1.0d3
			if (retdeact.AND.dt.LT.1.0d2) dt=1.0d2
		endif
		tfix=dtfix*(idint(ttot/dtfix)+1)
		if (((tfix-(ttot+dt)).LT.dt/2).OR.((ttot+dt).GT.tfix)) dt=tfix-ttot
		write(50,'(3e25.15,i15)') ttot-dtold,dtold,dtavg(51),iter
!
	end subroutine CONVERGED
!	_____________________________________________________________________________________
!
!
!
!
	SUBROUTINE FEM	(ttot,dt,nogpel,nem,nnm,nhbw,dim,nod,mater,coord,orie,temp0,temp,&
					&pres,mois0,mois,tpm,tpmloc,ret,retval,atmos,satur,wgp,Ni,dNi,bdne,&
					&bdele,K,Fe,dF,noddra,nodsup,eledra,elesup,unbal,overl,Kassem)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::nogpel(:),nem,nnm,nhbw,dim,nod(:,:),mater(:),tpm,tpmloc,bdne,&
			&bdele(:,:)
		double precision,intent(in)::ttot,dt,coord(:,:),wgp(:,:),Ni(:,:,:),dNi(:,:,:,:),&
			&temp0(:),mois0(:),temp(:),pres(:),mois(:),retval(:),orie(:)
		logical,intent(in)::Kassem,atmos,satur
		double precision,intent(inout)::noddra(:,:),nodsup(:,:),eledra(:,:,:),elesup(:,:,:)
		logical,intent(inout)::ret(:)
		double precision,intent(inout)::K(:,:),Fe(:),dF(:)
		double precision,intent(out)::unbal
		logical,intent(out)::overl
		integer::iem,ine,idim,jne,j,nr,nrn,nc,mat,ibd,bdloc,bddir,inm,nn,nne,ngp,&
			&nrloc,bdkind
		double precision::dum,w,Kmax,Kmin,Kact(3)
		double precision,allocatable::tval0(:),tval(:),pval(:),mval0(:),mval(:),&
			&elstif(:,:),elfe(:),eldf(:),el(:,:),wg(:),N(:,:),dN(:,:,:)
		logical::extra,trian
!		---------------------------------------------------------------------------------
!	calculation of matrices for system of equations (transfer & boundary conditions)
		K=0.0d0; dF(:)=0.0d0; Fe=0.0d0; 
		quantity1: select case(tpmloc)
			case(100); noddra(:,1)=0.0d0; nodsup(:,1)=0.0d0
			case(010); noddra(:,2)=0.0d0; nodsup(:,2)=0.0d0
			case(001); noddra(:,3)=0.0d0; nodsup(:,3)=0.0d0
		end select quantity1
!	mass and stiffness matrices (collection, per element) 
		do iem=1,nem
			nne=nogpel(2); trian=.FALSE.
			if (nod(iem,nne).EQ.0) then
				nne=nogpel(1); trian=.TRUE.
			endif
			mat=mater(iem)
			extra=.FALSE.
			if (tpmloc.EQ.1) then
				Kmax=-2.5d1; Kmin=0.0d0
				do ine=1,nne
					nn=nod(iem,ine)
					dum=dmax1(dmin1(mois(nn),-1.0d-5),-1.0d10)
				Kact=matlib_logKl(mat,temp(nn),dum,dlog10(-dum),0.0d0)
					if (Kact(1).GT.Kmax) Kmax=Kact(1); if (Kact(1).LT.Kmin) Kmin=Kact(1)
				enddo
				if (Kmax-Kmin.GE.1.0d0) extra=.TRUE.
			endif
			if (extra) then	
				if (.NOT.trian) then
					ngp=nogpel(8)
				else
					ngp=nogpel(7)
				endif
			else					
				if (.NOT.trian) then
					ngp=nogpel(6)
				else
					ngp=nogpel(5)
				endif
			endif
			allocate (wg(ngp),N(nne,ngp),dN(nne,ngp,dim),el(nne,dim),tval0(nne),&
				&tval(nne),pval(nne),mval0(nne),mval(nne),elstif(nne,nne),elfe(nne),&
				&eldf(nne))
			if (.NOT.extra) then
				if (trian) then
					wg=wgp(:,1); N=Ni(:,:,1); dN=dNi(:,:,:,1)
				else
					wg=wgp(:,2); N=Ni(:,:,2); dN=dNi(:,:,:,2)
				endif
			else					
				if (trian) then
					wg=wgp(:,3); N=Ni(:,:,3); dN=dNi(:,:,:,3)
				else
					wg=wgp(:,4); N=Ni(:,:,4); dN=dNi(:,:,:,4)
				endif
			endif
			el=0.0d0
			do ine=1,nne
				nn=nod(iem,ine); el(ine,:)=coord(nn,:)				
				tval0(ine)=temp0(nn); tval(ine)=temp(nn)
				                      pval(ine)=pres(nn)
				mval0(ine)=mois0(nn); mval(ine)=mois(nn)
			enddo
			CALL MATRIX(nne,ngp,dim,tpm,tpmloc,dt,el,tval0,tval,pval,mval0,&
				&mval,mat,orie,wg,N,dN,elstif,elfe,eldf)
			do ine=1,nne
				nr=nod(iem,ine)
				Fe(nr)=Fe(nr)+elfe(ine)
				dF(nr)=dF(nr)+eldf(ine)
				do jne=1,nne
					nc=nod(iem,jne)
					if (nr.LE.nc) then 
						nrn=nr-nc+nhbw; K(nrn,nc)=K(nrn,nc)+elstif(ine,jne)
					endif
				enddo
				quantity2: select case(tpmloc)
					case(100)
						noddra(nr,1)=noddra(nr,1)-eldf(ine); eledra(iem,1:nne,1)=-eldf
					case(010)
						noddra(nr,2)=noddra(nr,2)-eldf(ine); eledra(iem,1:nne,2)=-eldf
					case(001)
						noddra(nr,3)=noddra(nr,3)-eldf(ine); eledra(iem,1:nne,3)=-eldf
				end select quantity2
			enddo
			deallocate(wg,N,dN,el,tval0,tval,pval,mval0,mval,elstif,elfe,eldf)
		enddo
!		
!	imposing boundary conditions
		do ibd=1,bdne
			nne=nogpel(4); ngp=nogpel(10); trian=.FALSE.
			if (bdele(ibd,3+nne).EQ.0) then
				nne=nogpel(3); ngp=nogpel(9); trian=.TRUE.
			endif
			allocate(wg(ngp),N(nne,ngp),dN(nne,ngp,dim-1),el(nne,dim-1),tval(nne),&
				&pval(nne),mval(nne),elstif(nne,nne),elfe(nne))
			if (trian) then
				wg=wgp(:,5); N=Ni(:,:,5); dN=dNi(:,:,:,5)
			else
				wg=wgp(:,6); N=Ni(:,:,6); dN=dNi(:,:,:,6)
			endif
			mat=mater(bdele(ibd,1)); bdloc=bdele(ibd,2); bddir=bdele(ibd,3); el=0.0d0
			do ine=1,nne
				nn=nod(bdele(ibd,1),bdele(ibd,3+ine)); j=0
				do idim=1,dim
					if (idim.NE.bddir) then
						j=j+1; el(ine,j)=coord(nn,idim)
					endif
				enddo
				tval(ine)=temp(nn); pval(ine)=pres(nn); mval(ine)=mois(nn)
			enddo
			CALL BOUND (ttot,dt,nne,ngp,dim,el,tval,pval,mval,mat,tpm,tpmloc,atmos,wg,N,&
				&dN,bdloc,orie,bdkind,elstif,elfe)
			do ine=1,nne
				nr=nod(bdele(ibd,1),bdele(ibd,3+ine))
				if (bdkind.EQ.2) then					
					Fe(nr)=Fe(nr)+elfe(ine)
					dF(nr)=dF(nr)+elfe(ine)
					do jne=1,nne
						nc=nod(bdele(ibd,1),bdele(ibd,3+jne))
						if (nr.LE.nc) then
							nrn=nr-nc+nhbw; K(nrn,nc)=K(nrn,nc)+elstif(ine,jne)
						endif
					enddo
				endif
				quantity3: select case(tpmloc)
					case(100)
						if (bdkind.EQ.2) then
							nodsup(nr,1)=nodsup(nr,1)+elfe(ine); elesup(ibd,1:nne,1)=elfe
						endif
						if (bdkind.EQ.1) CALL FIX(nnm,nhbw,K,dF,nr,elfe(ine)-tval(ine))
					case(010)
						if (bdkind.EQ.2) then
							nodsup(nr,2)=nodsup(nr,2)+elfe(ine); elesup(ibd,1:nne,2)=elfe
						endif
						if (bdkind.EQ.1) CALL FIX(nnm,nhbw,K,dF,nr,elfe(ine)-pval(ine))
					case(001)
						if (bdkind.EQ.2) then
							nodsup(nr,3)=nodsup(nr,3)+elfe(ine); elesup(ibd,1:nne,3)=elfe
						endif
						if (bdkind.EQ.1) CALL FIX(nnm,nhbw,K,dF,nr,elfe(ine)-mval(ine))
				end select quantity3
			enddo
			deallocate(wg,N,dN,el,tval,pval,mval,elstif,elfe)
		enddo
!	impose (internal) saturation where necessary + correct thermal boundary condition
		if (satur.AND.mod(tpmloc,10).EQ.1) then
			do inm=1,nnm 
				if (ret(inm)) then
					CALL FIX(nnm,nhbw,K,dF,inm,retval(inm)-mois(inm)) ! impose saturation
					if (tpmloc/100.EQ.1) dF(inm)=dF(inm)-(nodsup(inm,3)-noddra(inm,3))*cap_l*&
						&temp(inm) ! correction of thermal condition
				endif
			enddo
		endif
!	cumulating for convergence and divergence criteria
		overl=.FALSE.
		do inm=1,nnm
			if (dabs(dF(inm)).GT.1.0d12.OR.dabs(Fe(inm)).GT.1.0d12) overl=.TRUE. ! divergence checks
		enddo
		if (.NOT.overl) unbal=dsqrt(dot_product(dF,dF)/dot_product(Fe,Fe))
		if (unbal.GE.1.0d6) overl=.TRUE.
!
	end subroutine fem
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE MATRIX	(nne,ngp,dim,tpm,tpmloc,dt,el,tval0,tval,pval,mval0,mval,mat,&
						&orie,wg,N,dN,elstif,elfe,eldf)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::nne,ngp,dim,tpm,tpmloc,mat
		double precision,intent(in)::dt,el(:,:),tval0(:),tval(:),pval(:),mval0(:),&
			&mval(:),orie(:),wg(:),N(:,:),dN(:,:,:)
		double precision,intent(out)::elstif(:,:),elfe(:),eldf(:)
		integer::igp,ine,jne,idim
		double precision::tgau0,tgau,pgau,mgau0,mgau,Cxx,Cxy,Cxz,Kxx(3),Kxy(3),Kxz(3),&
			&Kxxgr(3),Kxygr(3),Kxzgr(3),Quan0,Quan,detjc,vol
		double precision,allocatable::Mxx(:,:),Mxy(:,:),Mxz(:,:),Sxx(:,:),Sxy(:,:),&
			&Sxz(:,:),Gxx(:),Gxy(:),Gxz(:),Q(:),U(:),V(:,:),kfac(:,:),jc(:,:),&
			&jcinv(:,:),xval(:),yval(:),zval(:)
!		---------------------------------------------------------------------------------
!	calculation of element matrices		
		allocate(U(nne),V(nne,nne),kfac(nne,dim),jc(dim,dim),jcinv(dim,dim),Mxx(nne,nne),&
			&Mxy(nne,nne),Mxz(nne,nne),Sxx(nne,nne),Sxy(nne,nne),Sxz(nne,nne),Gxx(nne),&
			&Gxy(nne),Gxz(nne),Q(nne),xval(nne),yval(nne),zval(nne))
!	initialize element matrices
		Mxx=0.0d0; Mxy=0.0d0; Mxz=0.0d0; Q=0.0d0
		Sxx=0.0d0; Sxy=0.0d0; Sxz=0.0d0; Gxx=0.0d0; Gxy=0.0d0; Gxz=0.0d0 
		quantity:select case(tpmloc)
			case(100); xval=tval; yval=pval; zval=mval
			case(010); xval=pval; yval=mval; zval=tval
			case(001); xval=mval; yval=tval; zval=pval 
		end select quantity
!	do loop on number of Gauss points
		vol=0.0d0
		do igp=1,ngp
			tgau0=dmax1(dmin1(dot_product(N(:,igp),tval0),4.0d2),2.0d2)
			mgau0=dmax1(dmin1(dot_product(N(:,igp),mval0),-1.0d-5),-1.0d10)
			tgau=dmax1(dmin1(dot_product(N(:,igp),tval),4.0d2),2.0d2)
			pgau=dmax1(dmin1(dot_product(N(:,igp),pval),1.0d3),-1.0d3)
			mgau=dmax1(dmin1(dot_product(N(:,igp),mval),-1.0d-5),-1.0d10)
			CALL MATERIALPROPERTIES(tgau0,tgau,pgau,mgau0,mgau,mat,tpm,tpmloc,&
				&orie,Cxx,Cxy,Cxz,Kxx,Kxy,Kxz,Kxxgr,Kxygr,Kxzgr,Quan0,Quan)
!
			call jacobian(dim,igp,el,dN,jc,jcinv,detjc,.FALSE.)
			detjc=detjc*wg(igp)
			vol=vol+detjc
			U=N(:,igp)*detjc
			do ine=1,nne	
				do jne=ine,nne
					V(ine,jne)=U(ine)*N(jne,igp)
				enddo
				do jne=1,ine-1
					V(ine,jne)=V(jne,ine)
				enddo
			enddo
			do ine=1,nne		 
				do idim=1,dim
					kfac(ine,idim)=dot_product(dN(ine,igp,1:dim),jcinv(1:dim,idim))
				enddo
			enddo	
!
			Q=Q+U*(Quan-Quan0); Mxx=Mxx+V*Cxx
			CALL STIFF(dim,nne,kfac,detjc,Kxx,Sxx) 
			if (sum(Kxxgr).NE.0.0d0) CALL GRAVI(dim,nne,kfac,detjc,Kxxgr,Gxx) 
			if (Cxy.NE.0.0d0) Mxy=Mxy+V*Cxy
			if (sum(Kxy).NE.0.0d0) CALL STIFF(dim,nne,kfac,detjc,Kxy,Sxy)
			if (sum(Kxygr).NE.0.0d0) CALL GRAVI(dim,nne,kfac,detjc,Kxygr,Gxy) 
			if (Cxz.NE.0.0d0) Mxz=Mxz+V*Cxz
			if (sum(Kxz).NE.0.0d0) CALL STIFF(dim,nne,kfac,detjc,Kxz,Sxz)
			if (sum(Kxzgr).NE.0.0d0) CALL GRAVI(dim,nne,kfac,detjc,Kxzgr,Gxz) 
		enddo
!	collection of different matrices
		elstif=Mxx+Sxx*dt
		elfe=matmul(Mxx,xval)-Q
		eldf=-Q-matmul(Sxx,xval)*dt
!					elfe=matmul(Mxx,mval0) !-Q
!					eldf=elfe-matmul(elstif,xval)
		if (sum(Kxxgr).NE.0.0d0) then
			elfe=elfe+Gxx*dt; eldf=eldf+Gxx*dt
		endif
		if (Cxy.NE.0.0d0) elfe=elfe+matmul(Mxy,yval)
		if (sum(Kxy).NE.0.0d0) eldf=eldf-matmul(Sxy,yval)*dt
		if (sum(Kxygr).NE.0.0d0) then
			elfe=elfe+Gxy*dt; eldf=eldf+Gxy*dt
		endif
		if (Cxz.NE.0.0d0) elfe=elfe+matmul(Mxz,zval)
		if (sum(Kxz).NE.0.0d0) eldf=eldf-matmul(Sxz,zval)*dt
		if (sum(Kxzgr).NE.0.0d0) then
			elfe=elfe+Gxz*dt; eldf=eldf+Gxz*dt
		endif
!
		deallocate(U,V,kfac,jc,jcinv,Mxx,Mxy,Mxz,Sxx,Sxy,Sxz,Gxx,Gxy,Gxz,xval,yval,zval)
!		
	end subroutine MATRIX
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE BOUND	(ttot,dt,nne,ngp,dim,el,tval,pval,mval,mat,tpm,tpmloc,atmos,wg,&
						&N,dN,bdloc,orie,bdkind,elstif,elfe)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::nne,ngp,dim,mat,tpm,tpmloc,bdloc
		double precision,intent(in)::ttot,dt,el(:,:),tval(:),pval(:),mval(:),wg(:),&
			&N(:,:),dN(:,:,:),orie(:)
		logical,intent(in)::atmos
		integer,intent(out)::bdkind
		double precision,intent(out)::elstif(:,:),elfe(:)
		integer::igp,ine,jne
		double precision::tgau,pgau,mgau,xgau,detjc,Bxx,dBxxdx
		double precision,allocatable::jc(:,:),jcinv(:,:),U(:),V(:,:),cgau(:)
!		---------------------------------------------------------------------------------
!	calculation of boundary conditions elements and forces
		allocate(jc(dim-1,dim-1),jcinv(dim-1,dim-1),U(nne),V(nne,nne),cgau(dim-1))
		bdkind=0; elfe=0.0d0; elstif=0.0d0
!		
		do igp=1,ngp
			tgau=dmax1(dmin1(dot_product(N(:,igp),tval),4.0d2),2.0d2)
			pgau=dmax1(dmin1(dot_product(N(:,igp),pval),1.0d3),-1.0d3)
			mgau=dmax1(dmin1(dot_product(N(:,igp),mval),-1.0d-5),-1.0d15)
			do ine=1,dim-1
				cgau(ine)=dot_product(N(:,igp),el(:,ine))
			enddo
			CALL BOUNDARYCONDITIONS(tgau,pgau,mgau,cgau,orie,bdloc,mat,tpm,tpmloc,atmos,&
				&bdkind,Bxx,dBxxdx,ttot,dt)
!	for flux/load conditions: distribution to nodal values
			if (bdkind.EQ.1) elfe=Bxx
			if (bdkind.EQ.2) then
				call jacobian(dim-1,igp,el,dN,jc,jcinv,detjc,.TRUE.); detjc=detjc*wg(igp)
				U=N(:,igp)*detjc
				do ine=1,nne
					do jne=ine,nne
						V(ine,jne)=U(ine)*N(jne,igp)
					enddo
					do jne=1,ine-1
						V(ine,jne)=V(jne,ine)
					enddo
				enddo
				elfe=elfe+U*Bxx*dt
				elstif=elstif-V*dBxxdx*dt
			endif
		enddo
!		
		deallocate(jc,jcinv,U,V,cgau)
!
	end subroutine BOUND
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE JACOBIAN		(dimen,igp,el,dN,jc,jcinv,detjc,bound)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::dimen,igp 
		double precision,intent(in)::el(:,:),dN(:,:,:)
		logical,intent(in)::bound
		double precision,intent(out)::jc(:,:),jcinv(:,:),detjc
		double precision::invdetjc
		integer::i,j,ia,ib,ja,jb
!		---------------------------------------------------------------------------------
!   calculation of jacobian, its determinant and its inverse
		do i=1,dimen
			do j=1,dimen
				jc(i,j)=dot_product(el(:,i),dN(:,igp,j))
			enddo
		enddo
!
		determinant: select case(dimen)
			case(0)
				detjc=1.0d0
			case(1)
				detjc=jc(1,1)
			case(2)
				detjc=jc(1,1)*jc(2,2)-jc(1,2)*jc(2,1)
			case(3)
				detjc=jc(1,1)*jc(2,2)*jc(3,3)+jc(1,2)*jc(2,3)*jc(3,1)+ &
					&jc(2,1)*jc(3,2)*jc(1,3)-jc(3,1)*jc(2,2)*jc(1,3)- &
					&jc(2,1)*jc(1,2)*jc(3,3)-jc(3,2)*jc(2,3)*jc(1,1)
		end select determinant
!	
		if (.NOT.bound) then
			invdetjc=1.0d0/detjc
			inverse_jacob: select case(dimen)
				case(1)
					jcinv(1,1)=invdetjc
				case(2)
					jcinv(1,1)=jc(2,2); jcinv(1,2)=-jc(1,2); jcinv(2,1)=-jc(2,1)
					jcinv(2,2)=jc(1,1); jcinv=jcinv*invdetjc
				case(3)
					do i=1,3
						ia=mod(i,3)+1
						ib=mod(i+1,3)+1
						do j=1,3
							ja=mod(j,3)+1
							jb=mod(j+1,3)+1
							jcinv(j,i)=(jc(ia,ja)*jc(ib,jb)-jc(ib,ja)*jc(ia,jb))
						enddo
					enddo
					jcinv=jcinv*invdetjc
			end select inverse_jacob
		endif
		detjc=dabs(detjc)
!
	end subroutine JACOBIAN
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE STIFF	(dim,nne,kfac,detjc,Kxx,Sxx)  

!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::dim,nne
		double precision,intent(in)::kfac(:,:),detjc,Kxx(3)
		double precision,intent(out)::Sxx(:,:)
		integer::ine,jne
!		---------------------------------------------------------------------------------
!	calculation of stiffness matrix for one element	and one gausspoint
		do ine=1,nne	
			do jne=1,nne
				Sxx(ine,jne)=Sxx(ine,jne)+dot_product(Kxx(1:dim)*kfac(ine,1:dim),&
					&kfac(jne,1:dim))*detjc
			enddo
		enddo
!	
	end subroutine STIFF
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE GRAVI	(dim,nne,kfac,detjc,Kxxgr,Gxx)

!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::dim,nne
		double precision,intent(in)::kfac(:,:),detjc,Kxxgr(3)
		double precision,intent(out)::Gxx(:)
		integer::ine
!		---------------------------------------------------------------------------------
!	calculation of the gravity matrix for one element and one gausspoint
		do ine=1,nne	
			Gxx(ine)=Gxx(ine)+dot_product(Kxxgr(1:dim),kfac(ine,1:dim))*detjc
		enddo
!	
	end subroutine GRAVI
!	_____________________________________________________________________________________
!
!
!
	SUBROUTINE FIX		(nnm,nhbw,s,sl,fixno,fixval)
!	_____________________________________________________________________________________
		implicit none
		integer,intent(in)::nnm,nhbw,fixno
		double precision,intent(in)::fixval
		double precision,intent(inout)::s(:,:),sl(:)
		integer::i,nc
!		---------------------------------------------------------------------------------
!   Including prescribed nodal values on banded matrix (s cfr. gstif, sl cfr. gf)
		do nc=1,nnm
			if (iabs(nc-fixno).LT.nhbw) then
				if (fixno.GT.nc) then
					i=nhbw+nc-fixno
					sl(nc)=sl(nc)-s(i,fixno)*fixval; s(i,fixno)=0.0d0
				endif
				if (nc.GT.fixno) then
					i=nhbw+fixno-nc
					sl(nc)=sl(nc)-s(i,nc)*fixval; s(i,nc)=0.0d0
				endif
			endif
		enddo
		s(nhbw,fixno)=1.0d0; sl(fixno)=fixval											
	end subroutine fix
!	_____________________________________________________________________________________
