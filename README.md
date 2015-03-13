# cfd_assn2
elliptic heat conduction

	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k,choice,t1,t2,iter
	double precision l,h,delx,dely,beta,x,y,eva,om
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	open(10,file='elliptic_tpsor.dat')
	l = 1.0d0
	h = 2.0d0
	om = 1.d0
	iter=0
	delx = l/(M-1)
	dely = h/(N-1)
	beta = (delx**2)/(dely**2)
	eva = 1/(2*(1+beta**2))
	print *,beta
	
10 	do i = 1,M
		do j = 1,N
		TEMP(i,j) = 0.d0
		TEMP1(i,j) = TEMP(i,j)
		TEMP2(i,j) = TEMP(i,j)
		end do
	end do
11 	do i = 1,M
	TEMP(i,1) = 100.d0
	TEMP1(i,1)= 100.d0
	TEMP2(i,1)= 100.d0
	end do
	
12 	call lgsi
	val = TEMP1 - TEMP
	TEMP = TEMP1
13 	do i = 2,M-1
		do j = 2,N-1
14 			if (abs(val(i,j)) .le. 0.001) then
				exit
			else
				iter=iter+1
				goto 12
			end if
		end do
	end do
	TEMP = TEMP - val
	x = 0.d0
	y = 0.d0
	write(10,*) 'TITLE="TEMP DIST"'
	write(10,*) 'VARIABLES = "X","Y" "TEMP"'
	write(10,*) 'ZONE T="ONLY ZONE", I=21,J=41, F=POINT'
	do j = 1,N
		do i = 1,M
		write(10,*) x,y,TEMP1(i,j)
		x = x + delx
		end do
		x = 0.d0
		y = y + dely
	end do
	write(*,*) iter
	stop
	end
!-------------------------------------------------------------------------
	subroutine jacobi
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k
	double precision l,h,delx,dely,beta,x,y,eva,tend,om
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	do i = 2,M-1
		do j = 2,N-1
		tend = (beta**2)*((TEMP(i,j+1))+(TEMP(i,j-1)))
		TEMP1(i,j) = eva*(TEMP((i+1),j) + TEMP((i-1),j) + tend)
		end do
	end do
	return
	end subroutine
!-------------------------------------------------------------------------
	subroutine pgs
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k
	double precision l,h,delx,dely,beta,x,y,eva,tend,om
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	do i = 2,M-1
		do j = 2,N-1
		tend = (beta**2)*((TEMP(i,j+1))+(TEMP1(i,j-1)))
		TEMP1(i,j) = eva*(TEMP((i+1),j) + TEMP1((i-1),j) + tend)
		end do
	end do
	return
	end subroutine
!--------------------------------------------------------------------------
	subroutine psor
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k
	double precision l,h,delx,dely,beta,x,y,eva,tend,om,tat
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)

	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om
	do i = 2,M-1
		do j = 2,N-1
		tend = (beta**2)*((TEMP(i,j+1))+(TEMP1(i,j-1)))
		tat = eva*(TEMP((i+1),j) + TEMP1((i-1),j) + tend)
		TEMP1(i,j)=((1-om)*TEMP(i,j))+(om*tat)
		end do
	end do
	return
	end subroutine
!--------------------------------------------------------------------------
	subroutine lgsi
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k
	double precision l,h,delx,dely,beta,x,y,eva,tend,om,tat
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)
	double precision a(M),b(M),c(M),d(M),nan,cas
	
	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	do j=2,N-1

	
	d(1)=-(2*(1+beta**2))
	
	c(1)=-(beta**2)*(TEMP1(1,j-1)+TEMP(1,j+1))
	
	do i=2,M-1
	a(i-1)=1.d0
	b(i)=1.d0
	d(i)=-(2*(1+beta**2))
	d(i)=d(i)-(b(i)*a(i-1))/d(i-1)
	
	nan=-(beta**2)*(TEMP1(i,j-1)+TEMP(i,j+1))
	
	c(i)=nan-(b(i)*c(i-1))/d(i-1)
	
	
	end do
	
	TEMP1(21,j)=c(M)/d(M)  
	
	do i=M-1,2,-1

	TEMP1(i,j)=(c(i)-(a(i)*TEMP1(i+1,j)))/d(i)
	
	end do
	end do
	return
	end subroutine
!---------------------------------------------------------------------------
	subroutine lsor
	implicit none
	integer M,N
	parameter(M=21)
	parameter(N=41)
	integer i,j,k
	double precision l,h,delx,dely,beta,x,y,eva,tend,om,tat
	double precision TEMP(M,N),TEMP1(M,N),TEMP2(M,N),Val(M,N)
	double precision a(M),b(M),c(M),d(M),nan,cas,tak
	
	common /yval1/ i,j,beta,eva,TEMP,TEMP1,TEMP2,om

	do j=2,N-1

	
	d(1)=-(2*(1+beta**2))
	cas =om*(beta**2)*(TEMP1(1,j-1)+TEMP(1,j+1))
	c(1)=-(1-om)*((1+beta**2))*TEMP(1,j)-cas
	write(*,*) c(1)
	do i=2,M-1
	a(i-1)=om
	b(i)=om
	d(i)=-(2*(1+beta**2))
	d(i)=d(i)-(b(i)*a(i-1))/d(i-1)
	tak=om*(beta**2)*(TEMP1(i,j-1)+TEMP(i,j+1))
	nan=-(1-om)*((1+beta**2))*TEMP(i,j)-tak
	c(i)=nan-(b(i)*c(i-1))/d(i-1)
		
	end do
	
	TEMP1(21,j)=c(M)/d(M)  
	
	do i=M-1,2,-1

	TEMP1(i,j)=(c(i)-(a(i)*TEMP1(i+1,j)))/d(i)
	
	end do
	end do
	return
	end subroutine
	
