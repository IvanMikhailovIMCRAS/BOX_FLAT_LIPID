program begin
	
	real(4) dx, dy, dz, xt, yt, zt, BOX(3), l_bond
	real(4) x, y, z, x1, x2, y1, y2, z1, z2
	integer(4) NA, NB, NW, N, typ, d, M, ioer, tmp , i, j
	
	open(11,file='input.txt')
	
	read(11,*,iostat=ioer) BOX(1), BOX(2), BOX(3)  
	if (ioer.ne.0) then
		backspace(11)
		read(11,*,iostat=ioer) BOX(1)
		BOX(2) = BOX(1); BOX(3) = BOX(1)
	endif
	read(11,*) N
	read(11,*) NW
	read(11,*) l_bond
	close(11)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	open(22,file='COORD')
	write(22,*) 'num_atoms ', NW+N*11, ' BOX ', BOX
	
	open(33,file='BONDS')
	write(33,*) 'num_bonds ', N*10, ' num_atoms ', NW+N*11, ' BOX ', BOX
	
	open(55,file='ANGLS')
	write(55,*) 'num_angl ', N*6
	
	tmp = 0
	do i = 1, N
		typ = 1
		call Random_Number(xt); call Random_Number(yt); call Random_Number(zt)
		x = (0.5 - xt)*BOX(1)
		y = (0.5 - yt)*BOX(2)
		!do while (x**2 + y**2 > (BOX(1)/2.0-3.0)**2)
		!	call Random_Number(xt); call Random_Number(yt)
		!	x = (0.5 - xt)*(BOX(1)-4.)
		!	y = (0.5 - yt)*(BOX(2)-4.)
		!enddo
		z = 2.0 * (-1)**i
		!x = (0.5 - xt)*BOX(1)
		!y = (0.5 - yt)*BOX(2)
		!z = (0.5 - zt)*BOX(3)
		tmp = tmp + 1
		write(22,*) tmp, x, y, z, typ

		call rndm_vector(xt,yt,zt,i)
		x1 = x + xt*l_bond; y1 = y + yt*l_bond; z1 = z + zt*l_bond 
		call periodic_test(x1, BOX(1))
		call periodic_test(y1, BOX(2))
		call periodic_test(z1, BOX(3))
		tmp = tmp + 1
		write(22,*) tmp, x1, y1, z1, typ
		write(33,*) tmp-1, tmp
		
		call rndm_vector(xt,yt,zt,i)
		x2 = x1 + xt*l_bond; y2 = y1 + yt*l_bond; z2 = z1 + zt*l_bond 
		call periodic_test(x2, BOX(1))
		call periodic_test(y2, BOX(2))
		call periodic_test(z2, BOX(3))
		tmp = tmp + 1
		write(22,*) tmp, x2, y2, z2, typ
		write(33,*) tmp-1, tmp
		
		typ = 2
		
		call rndm_vector(xt,yt,zt,i)
		x = x1 + xt*l_bond; y = y1 + yt*l_bond; z = z1 + zt*l_bond 
		call periodic_test(x, BOX(1))
		call periodic_test(y, BOX(2))
		call periodic_test(z, BOX(3))
		tmp = tmp + 1
		write(22,*) tmp, x, y, z, typ
		write(33,*) tmp-2, tmp
		write(55,*) tmp-2, tmp, tmp+1
		
		do j = 1, 3
			call rndm_vector(xt,yt,zt,i)
			x = x + xt*l_bond; y = y + yt*l_bond; z = z + zt*l_bond 
			call periodic_test(x, BOX(1))
			call periodic_test(y, BOX(2))
			call periodic_test(z, BOX(3))
			tmp = tmp + 1
			write(22,*) tmp, x, y, z, typ
			write(33,*) tmp-1, tmp
			write(55,*) tmp-1, tmp, tmp+1	
		enddo
		backspace(55)
		
		call rndm_vector(xt,yt,zt,i)
		x = x2 + xt*l_bond; y = y2 + yt*l_bond; z = z2 + zt*l_bond 
		call periodic_test(x, BOX(1))
		call periodic_test(y, BOX(2))
		call periodic_test(z, BOX(3))
		tmp = tmp + 1
		write(22,*) tmp, x, y, z, typ
		write(33,*) tmp-5, tmp
		write(55,*) tmp-5, tmp, tmp+1
		
		do j = 1, 3
			call rndm_vector(xt,yt,zt,i)
			x = x + xt*l_bond; y = y + yt*l_bond; z = z + zt*l_bond 
			call periodic_test(x, BOX(1))
			call periodic_test(y, BOX(2))
			call periodic_test(z, BOX(3))
			tmp = tmp + 1
			write(22,*) tmp, x, y, z, typ
			write(33,*) tmp-1, tmp
			write(55,*) tmp-1, tmp, tmp+1		
		enddo		
		backspace(55)
	enddo
	
	write(55,*) ' '
	
	typ = 3
	do while (tmp.ne.N*11+NW)
		tmp = tmp + 1
		call Random_Number(xt); call Random_Number(yt); call Random_Number(zt) 
		xt = (0.5 - xt)*BOX(1)
		yt = (0.5 - yt)*BOX(2)
		zt = (0.5 - zt)*BOX(3)
		do while (zt > -1.5 .and. zt < 1.5 )
			call Random_Number(zt)
			zt = (0.5 - zt)*BOX(3)
		enddo
		write(22,*) tmp, xt, yt, zt, typ
	enddo
			
	close(22)
	close(33)
	close(55)

stop
end program begin
!#######################################################################################

!***************************************************************************************
subroutine rndm(x,y,z)
!***************************************************************************************
implicit none
	real(4), intent(out) :: x, y, z
	
	call Random_Number(x); call Random_Number(y); call Random_Number(z)
	x = 1.0 - 2.0*x; y = 1.0 - 2.0*y
	z = 1.0 - 2.0*z 

return
end subroutine rndm
!***************************************************************************************

!***************************************************************************************
subroutine rndm_vector(x,y,z,i)
!***************************************************************************************
implicit none
	real(4), intent(out) :: x, y, z
	integer(4), intent(in) :: i
	real(4) r 
	
	call rndm(x,y,z)
	z = z - 0.5*(-1.0)**i
	r = sqrt(x**2 + y**2 + z**2)
	x = x/r; y = y/r; z = z/r

return
end subroutine rndm_vector
!***************************************************************************************

!***************************************************************************************
subroutine periodic_test(x, Box) ! проверка периодических граничных условий
!***************************************************************************************
	real(4),intent(inout) :: x
	real(4),intent(in) :: Box

	if (abs(x).gt.0.5*Box) x = x - x/abs(x)*Box

return
end subroutine periodic_test
!***************************************************************************************
