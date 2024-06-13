!
!		    		1D/2D/3D-FINITE ELEMENT MODEL FOR THE SIMULATION OF 
!			    	 HEAT, AIR AND MOISTURE TRANSFER IN BUILDING PARTS
!
!															 written by: Daniel Cóstola
!****************************************************************************************
!****************************************************************************************
!
									PROGRAM HAMFEMTESTER					
!
!****************************************************************************************
!****************************************************************************************
!																		updated: jan 2009
!
!
		implicit none
!   
!   VARIABLE DECLARATION
!	_____________________________________________________________________________________
		integer::ntest,j,i,m,n,x,y,nlines(5),nnodes(5),tpm(5),rel,io
		character(14)::installpath
	    character(40)::commandline

!		ntest=1 define the size of those arrays
		character(len=20),dimension(5)::inputfile
!		character(len=20),dimension(5)::xts
		character(len=5),dimension(8)::restype
		character(len=25),dimension(8)::restypen

       double precision::mindif,maxdif,meandif,stdevdif
	   double precision::timedummynew,timenew,timeref,timedummyref
       double precision,allocatable::statesref(:,:),statesnew(:,:), statesdif(:,:)
!

!	_____________________________________________________________________________________



!   MAIN PROGRAM
!	_____________________________________________________________________________________
!	_____________________________________________________________________________________
!Inicialization
!		installpath='D:\00baseversion\hamfem'
		j=1
		ntest=4 ! change to 5 to run Hamstad 4. this simulation is very slow due to intense fluid flux.

		i=1

		restype(1)='outte'
		restype(2)='outpr'
		restype(3)='outmo'
		restype(4)='outth'
		restype(5)='flote'
		restype(6)='flopr'
		restype(7)='flomo'
		restype(8)='times'




restypen(1)='_nodetemperat.out'
restypen(2)='_nodeairpress.out'
restypen(3)='_nodecappress.out'
restypen(4)='_nodemoisture.out'
restypen(5)='_fluxheat.out'
restypen(6)='_fluxair.out'
restypen(7)='_fluxmois.out'
restypen(8)='_time.out'



       open(10,file='summarytest.txt')


! cases to test


		inputfile(1)='00basevers'
		nlines(1)=865
		nnodes(1)=51
!		xts(1)='00input'
		tpm(1)=001


		inputfile(2)='07AtmNOwdr'
		nlines(2)=2401
		nnodes(2)=51
!		xts(2)='07input'
		tpm(2)=101


		inputfile(3)='08atmsrain'
		nlines(3)=2401
		nnodes(3)=51
	!		xts(3)='08input'
		tpm(3)=101

		inputfile(4)='09Hamstad1'
		nlines(4)=8758
		nnodes(4)=30
!		xts(4)='09input'
		tpm(4)=101


		inputfile(5)='10Hamstad4'
		nlines(5)=121
		nnodes(5)=42
!		xts(5)='10input'
		tpm(5)=101

!start loop. run hamfem and compare with reference results
		do while (j.Le.ntest)

! Run HAMFEM		
       commandline='hamfem '//TRIM(inputfile(j))//'.ham'
 write(*,*) commandline
       call system (commandline)

! Compare results

       write(*,*) inputfile(j)
       write(10,*) inputfile(j)

       do while (i.LT.9)

! Check if the result file relevance, considering the quantities calculated (HAM).
		rel=0

		!heat calculation 
		if(tpm(j).GT.99.AND.i.eq.1.OR.tpm(j).GT.99.AND.i.eq.5) rel=1
		!mositure calculation 
		if(mod(tpm(j),10).eq.1.AND.i.eq.3.OR.mod(tpm(j),10).eq.1.AND.i.eq.4.OR.mod(tpm(j),10).eq.1.AND.i.eq.7) rel=1
		if(i.eq.5.or.i.eq.6.or.i.eq.7) rel=0 ! flows not included in the tester. statesref matrixes need to be redefined for that
		!Flow files, not suported by this tester.
		!if(i.eq.2.OR.i.eq.6.OR.i.eq.8) rel=0


! If results are relevant:
		if (rel.eq.1) then

! Open Referencefile


		allocate(statesref(nlines(j)+10,nnodes(j)+1))
		allocate(statesnew(nlines(j)+10,nnodes(j)+1))
		allocate(statesdif(nlines(j)+10,nnodes(j)+1))



       open(20,file='../reference-output/'//TRIM(inputfile(j))//TRIM(restypen(i)),action='read')
		read(20,*);read(20,*) !headers
			do m=1,nlines(j) !nnodes(j)+1 
			read(20,*) statesref(m,:)
			enddo
		close(20)		

       open(30,file='./'//TRIM(inputfile(j))//TRIM(restypen(i)),action='read')
!       open(30,file='./'//restype(i)//TRIM(xts(j))//'.out',action='read')
	   		read(30,*);read(30,*) !headers
			do m=1,nlines(j) !nnodes(j)+1
			read(30,*) statesnew(m,:)

			enddo
		close(30)

		statesdif=statesref-statesnew

		mindif=minval(statesdif)
		maxdif=maxval(statesdif)
		
		meandif=sum(statesdif)/((nlines(j)+1)*(nnodes(j)+1))


			do m=1,nlines(j)+1
				do n=1,nnodes(j)+1
					statesdif(m,n)=statesdif(m,n)**2 !the array statedif is used as a dummy here
				enddo
			enddo
		

      stdevdif=sqrt(sum(statesdif)/((nlines(j)+1)*(nnodes(j)+1)))


		deallocate(statesdif,statesref,statesnew)

       write(*,'(A,4E8.1)') restype(i),mindif,maxdif,meandif,stdevdif
       write(10,'(A,4E8.1)') restype(i),mindif,maxdif,meandif,stdevdif


		endif


		


		! Open timestep results and get run time		
		if(i.eq.8) then

				open(40,file='../reference-output/'//TRIM(inputfile(j))//TRIM(restypen(i)),action='read')
				timeref=0
				timedummyref=0
				io=0
					do while(io.ge.0)
					
					timeref=timedummyref

					READ(40, *, IOSTAT = IO ) timedummyref
					enddo
				close(40)



				open(50,file='./'//TRIM(inputfile(j))//TRIM(restypen(i)),action='read')
				timedummynew=0
				timenew=0
				io=0
					do while(io.ge.0)
					
					timenew=timedummynew

					READ(50, *, IOSTAT = IO ) timedummynew
						 
					enddo
				close(50)

				timenew=timenew-timeref

				write(*,*) timenew
				write(10,*) timenew
		
		endif


		i=i+1
		enddo

       j=j+1
	   i=1
       enddo

       close(10)


       END PROGRAM HAMFEMTESTER
