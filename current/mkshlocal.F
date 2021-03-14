	subroutine mkshlocal
	implicit none
	include "global.h"
	integer iq,i,j
	real*8 sm1,sp1,tm1,tp1,nm1,np1,rm1,rp1


C...SET UP FOR Semi-Discrete and Spacetime-Jump terms
C...TETS
        if (elem_type.eq.elem_tetr_sd .or. elem_type.eq.elem_tetr_st) then
        	call quad3d4n(iquad, nquad3, xq3, wq3, nsdpad, nquadpad)
		do iq = 1,nquad3
		sq3(0,1,iq) = xq3(1,iq)
		sq3(0,2,iq) = xq3(2,iq)
		sq3(0,3,iq) = xq3(3,iq)
		sq3(0,4,iq) = 1 - xq3(1,iq) - xq3(2,iq) - xq3(3,iq)

c               local first derivatives
                sq3(1,1,iq) =  1.0
                sq3(1,2,iq) =  0.0
                sq3(1,3,iq) =  0.0
                sq3(1,4,iq) = -1.0

                sq3(2,1,iq) =  0.0
                sq3(2,2,iq) =  1.0
                sq3(2,3,iq) =  0.0
                sq3(2,4,iq) = -1.0

                sq3(3,1,iq) =  0.0
                sq3(3,2,iq) =  0.0
                sq3(3,3,iq) =  1.0
                sq3(3,4,iq) = -1.0

		enddo
C...BRICKS
        else if (elem_type.eq.elem_hexa_sd .or. elem_type.eq.elem_hexa_st) then
        call quad3d8n(2, nquad3, xq3, wq3, nsdpad, nquadpad)
		do iq = 1,nquad3
		sq3(0,1,iq) = (1 - xq3(1,iq))
     &			   * (1 - xq3(2,iq)) * (1 - xq3(3,iq)) / 8
		sq3(0,2,iq) = (1 + xq3(1,iq))
     &			   * (1 - xq3(2,iq)) * (1 - xq3(3,iq)) / 8
		sq3(0,3,iq) = (1 + xq3(1,iq))
     &	                   * (1 + xq3(2,iq)) * (1 - xq3(3,iq)) / 8
		sq3(0,4,iq) = (1 - xq3(1,iq))
     &	                   * (1 + xq3(2,iq)) * (1 - xq3(3,iq)) / 8
		sq3(0,5,iq) = (1 - xq3(1,iq))
     &	                   * (1 - xq3(2,iq)) * (1 + xq3(3,iq)) / 8
		sq3(0,6,iq) = (1 + xq3(1,iq))
     &	                   * (1 - xq3(2,iq)) * (1 + xq3(3,iq)) / 8
		sq3(0,7,iq) = (1 + xq3(1,iq))
     &	                   * (1 + xq3(2,iq)) * (1 + xq3(3,iq)) / 8
		sq3(0,8,iq) = (1 - xq3(1,iq))
     &	                   * (1 + xq3(2,iq)) * (1 + xq3(3,iq)) / 8

c		local first derivatives
		sq3(1,1,iq) = - (1 - xq3(2,iq)) * (1 - xq3(3,iq)) / 8
		sq3(1,2,iq) = + (1 - xq3(2,iq)) * (1 - xq3(3,iq)) / 8
		sq3(1,3,iq) = + (1 + xq3(2,iq)) * (1 - xq3(3,iq)) / 8
		sq3(1,4,iq) = - (1 + xq3(2,iq)) * (1 - xq3(3,iq)) / 8
		sq3(1,5,iq) = - (1 - xq3(2,iq)) * (1 + xq3(3,iq)) / 8
		sq3(1,6,iq) = + (1 - xq3(2,iq)) * (1 + xq3(3,iq)) / 8
		sq3(1,7,iq) = + (1 + xq3(2,iq)) * (1 + xq3(3,iq)) / 8
		sq3(1,8,iq) = - (1 + xq3(2,iq)) * (1 + xq3(3,iq)) / 8
	
		sq3(2,1,iq) = - (1 - xq3(1,iq)) * (1 - xq3(3,iq)) / 8
		sq3(2,2,iq) = - (1 + xq3(1,iq)) * (1 - xq3(3,iq)) / 8
		sq3(2,3,iq) = + (1 + xq3(1,iq)) * (1 - xq3(3,iq)) / 8
		sq3(2,4,iq) = + (1 - xq3(1,iq)) * (1 - xq3(3,iq)) / 8
		sq3(2,5,iq) = - (1 - xq3(1,iq)) * (1 + xq3(3,iq)) / 8
		sq3(2,6,iq) = - (1 + xq3(1,iq)) * (1 + xq3(3,iq)) / 8
		sq3(2,7,iq) = + (1 + xq3(1,iq)) * (1 + xq3(3,iq)) / 8
		sq3(2,8,iq) = + (1 - xq3(1,iq)) * (1 + xq3(3,iq)) / 8
	
		sq3(3,1,iq) = - (1 - xq3(1,iq)) * (1 - xq3(2,iq)) / 8
		sq3(3,2,iq) = - (1 + xq3(1,iq)) * (1 - xq3(2,iq)) / 8
		sq3(3,3,iq) = - (1 + xq3(1,iq)) * (1 + xq3(2,iq)) / 8
		sq3(3,4,iq) = - (1 - xq3(1,iq)) * (1 + xq3(2,iq)) / 8
		sq3(3,5,iq) = + (1 - xq3(1,iq)) * (1 - xq3(2,iq)) / 8
		sq3(3,6,iq) = + (1 + xq3(1,iq)) * (1 - xq3(2,iq)) / 8
		sq3(3,7,iq) = + (1 + xq3(1,iq)) * (1 + xq3(2,iq)) / 8
		sq3(3,8,iq) = + (1 - xq3(1,iq)) * (1 + xq3(2,iq)) / 8
	
		enddo
C...TRIANGLES
        else if (elem_type.eq.elem_tria_sd) then
       	call quad2d3n(iquad, nquad3, xq3, wq3, nsdpad, nquadpad)
		do iq = 1,nquad3
                        sq3(0,1,iq) = xq3(1,iq)
                        sq3(0,2,iq) = xq3(2,iq)
                        sq3(0,3,iq) = 1.0 -xq3(1,iq) -xq3(2,iq)
                enddo
C...QUADS
        else if (elem_type.eq.elem_quad_sd) then
       	call quad2d4n(iquad, nquad3, xq3, wq3, nsdpad, nquadpad)
		  do iq = 1,nquad3
                         sm1     = (1-xq3(1,iq))
                         sp1     = (1+xq3(1,iq))
                         tm1     = (1-xq3(2,iq))
                         tp1     = (1+xq3(2,iq))
C...SHAPE FUNCTIONS VALUES
                        sq3(0,1,iq) = 0.25*sm1*tm1
                        sq3(0,2,iq) = 0.25*sp1*tm1
                        sq3(0,3,iq) = 0.25*sp1*tp1
                        sq3(0,4,iq) = 0.25*sm1*tp1

C...LOCAL DERIVATIVES
                        sq3(1,1,iq) = - 0.25*tm1
                        sq3(1,2,iq) = + 0.25*tm1
                        sq3(1,3,iq) = + 0.25*tp1
                        sq3(1,4,iq) = - 0.25*tp1

                        sq3(2,1,iq) = - 0.25*sm1
                        sq3(2,2,iq) = - 0.25*sp1
                        sq3(2,3,iq) = + 0.25*sp1
                        sq3(2,4,iq) = + 0.25*sm1
		enddo
        end if


C...SETUP SPACETIME SHAPEFUNCTIONS
C...TETS
        if (elem_type.eq.elem_tetr_st) then
                call quad4d4n(iquad,nquad4,xq4,wq4,4,16)

                do iq = 1,nquad4
                        rm1     = (1-xq4(4,iq))
                        rp1     = (1+xq4(4,iq))
C...SHAPE FUNCTIONS VALUES
                        sq4(0,1,iq) = 0.5*xq4(1,iq)*rm1
                        sq4(0,2,iq) = 0.5*xq4(2,iq)*rm1
                        sq4(0,3,iq) = 0.5*xq4(3,iq)*rm1
                        sq4(0,4,iq) =
     .                     0.5*(1.0-xq4(1,iq)-xq4(2,iq)-xq4(3,iq))*rm1
                        sq4(0,5,iq) = 0.5*xq4(1,iq)*rp1
                        sq4(0,6,iq) = 0.5*xq4(2,iq)*rp1
                        sq4(0,7,iq) = 0.5*xq4(3,iq)*rp1
                        sq4(0,8,iq) =
     .                     0.5*(1.0-xq4(1,iq)-xq4(2,iq)-xq4(3,iq))*rp1

C...LOCAL DERIVATIVES
                        sq4(1,1,iq) =  0.5*rm1
                        sq4(1,2,iq) =  0.0
                        sq4(1,3,iq) =  0.0
                        sq4(1,4,iq) = -0.5*rm1
                        sq4(1,5,iq) =  0.5*rp1
                        sq4(1,6,iq) =  0.0
                        sq4(1,7,iq) =  0.0
                        sq4(1,8,iq) = -0.5*rp1

                        sq4(2,1,iq) =  0.0
                        sq4(2,2,iq) =  0.5*rm1
                        sq4(2,3,iq) =  0.0
                        sq4(2,4,iq) = -0.5*rm1
                        sq4(2,5,iq) =  0.0
                        sq4(2,6,iq) =  0.5*rp1
                        sq4(2,7,iq) =  0.0
                        sq4(2,8,iq) = -0.5*rp1

                        sq4(3,1,iq) =  0.0
                        sq4(3,2,iq) =  0.0
                        sq4(3,3,iq) =  0.5*rm1
                        sq4(3,4,iq) = -0.5*rm1
                        sq4(3,5,iq) =  0.0
                        sq4(3,6,iq) =  0.0
                        sq4(3,7,iq) =  0.5*rp1
                        sq4(3,8,iq) = -0.5*rp1

                        sq4(4,1,iq) = -0.5*xq4(1,iq)
                        sq4(4,2,iq) = -0.5*xq4(2,iq)
                        sq4(4,3,iq) = -0.5*xq4(3,iq)
                        sq4(4,4,iq) = -0.5*(1.0-xq4(1,iq)-xq4(2,iq)-xq4(3,iq))
                        sq4(4,5,iq) =  0.5*xq4(1,iq)
                        sq4(4,6,iq) =  0.5*xq4(2,iq)
                        sq4(4,7,iq) =  0.5*xq4(3,iq)
                        sq4(4,8,iq) =  0.5*(1.0-xq4(1,iq)-xq4(2,iq)-xq4(3,iq))

                enddo

C...BRICKS
        else if (elem_type.eq.elem_hexa_st) then
		call quad4d8n(iquad,nquad4,xq4,wq4,4,16)

		do iq = 1,nquad4
       			 sm1     = (1-xq4(1,iq))
       			 sp1     = (1+xq4(1,iq))
       			 tm1     = (1-xq4(2,iq))
       			 tp1     = (1+xq4(2,iq))
       			 nm1     = (1-xq4(3,iq))
       			 np1     = (1+xq4(3,iq))
       			 rm1     = (1-xq4(4,iq))
       			 rp1     = (1+xq4(4,iq))
C...SHAPE FUNCTIONS VALUES
        		sq4(0, 1,iq) = 0.0625*sm1*tm1*nm1*rm1
        		sq4(0, 2,iq) = 0.0625*sp1*tm1*nm1*rm1
        		sq4(0, 3,iq) = 0.0625*sp1*tp1*nm1*rm1
        		sq4(0, 4,iq) = 0.0625*sm1*tp1*nm1*rm1
        		sq4(0, 5,iq) = 0.0625*sm1*tm1*np1*rm1
        		sq4(0, 6,iq) = 0.0625*sp1*tm1*np1*rm1
        		sq4(0, 7,iq) = 0.0625*sp1*tp1*np1*rm1
        		sq4(0, 8,iq) = 0.0625*sm1*tp1*np1*rm1
        		sq4(0, 9,iq) = 0.0625*sm1*tm1*nm1*rp1
        		sq4(0,10,iq) = 0.0625*sp1*tm1*nm1*rp1
        		sq4(0,11,iq) = 0.0625*sp1*tp1*nm1*rp1
        		sq4(0,12,iq) = 0.0625*sm1*tp1*nm1*rp1
        		sq4(0,13,iq) = 0.0625*sm1*tm1*np1*rp1
        		sq4(0,14,iq) = 0.0625*sp1*tm1*np1*rp1
        		sq4(0,15,iq) = 0.0625*sp1*tp1*np1*rp1
        		sq4(0,16,iq) = 0.0625*sm1*tp1*np1*rp1

C...LOCAL DERIVATIVES
        		sq4(1, 1,iq) = - 0.0625*tm1*nm1*rm1
        		sq4(1, 2,iq) = + 0.0625*tm1*nm1*rm1
        		sq4(1, 3,iq) = + 0.0625*tp1*nm1*rm1
        		sq4(1, 4,iq) = - 0.0625*tp1*nm1*rm1
        		sq4(1, 5,iq) = - 0.0625*tm1*np1*rm1
        		sq4(1, 6,iq) = + 0.0625*tm1*np1*rm1
        		sq4(1, 7,iq) = + 0.0625*tp1*np1*rm1
        		sq4(1, 8,iq) = - 0.0625*tp1*np1*rm1
        		sq4(1, 9,iq) = - 0.0625*tm1*nm1*rp1
        		sq4(1,10,iq) = + 0.0625*tm1*nm1*rp1
        		sq4(1,11,iq) = + 0.0625*tp1*nm1*rp1
        		sq4(1,12,iq) = - 0.0625*tp1*nm1*rp1
        		sq4(1,13,iq) = - 0.0625*tm1*np1*rp1
        		sq4(1,14,iq) = + 0.0625*tm1*np1*rp1
        		sq4(1,15,iq) = + 0.0625*tp1*np1*rp1
        		sq4(1,16,iq) = - 0.0625*tp1*np1*rp1

        		sq4(2, 1,iq) = - 0.0625*sm1*nm1*rm1
        		sq4(2, 2,iq) = - 0.0625*sp1*nm1*rm1
        		sq4(2, 3,iq) = + 0.0625*sp1*nm1*rm1
        		sq4(2, 4,iq) = + 0.0625*sm1*nm1*rm1
        		sq4(2, 5,iq) = - 0.0625*sm1*np1*rm1
        		sq4(2, 6,iq) = - 0.0625*sp1*np1*rm1
        		sq4(2, 7,iq) = + 0.0625*sp1*np1*rm1
        		sq4(2, 8,iq) = + 0.0625*sm1*np1*rm1
        		sq4(2, 9,iq) = - 0.0625*sm1*nm1*rp1
        		sq4(2,10,iq) = - 0.0625*sp1*nm1*rp1
        		sq4(2,11,iq) = + 0.0625*sp1*nm1*rp1
        		sq4(2,12,iq) = + 0.0625*sm1*nm1*rp1
        		sq4(2,13,iq) = - 0.0625*sm1*np1*rp1
        		sq4(2,14,iq) = - 0.0625*sp1*np1*rp1
        		sq4(2,15,iq) = + 0.0625*sp1*np1*rp1
        		sq4(2,16,iq) = + 0.0625*sm1*np1*rp1

        		sq4(3, 1,iq) = - 0.0625*sm1*tm1*rm1
        		sq4(3, 2,iq) = - 0.0625*sp1*tm1*rm1
        		sq4(3, 3,iq) = - 0.0625*sp1*tp1*rm1
        		sq4(3, 4,iq) = - 0.0625*sm1*tp1*rm1
        		sq4(3, 5,iq) = + 0.0625*sm1*tm1*rm1
        		sq4(3, 6,iq) = + 0.0625*sp1*tm1*rm1
        		sq4(3, 7,iq) = + 0.0625*sp1*tp1*rm1
        		sq4(3, 8,iq) = + 0.0625*sm1*tp1*rm1
        		sq4(3, 9,iq) = - 0.0625*sm1*tm1*rp1
        		sq4(3,10,iq) = - 0.0625*sp1*tm1*rp1
        		sq4(3,11,iq) = - 0.0625*sp1*tp1*rp1
        		sq4(3,12,iq) = - 0.0625*sm1*tp1*rp1
        		sq4(3,13,iq) = + 0.0625*sm1*tm1*rp1
        		sq4(3,14,iq) = + 0.0625*sp1*tm1*rp1
        		sq4(3,15,iq) = + 0.0625*sp1*tp1*rp1
        		sq4(3,16,iq) = + 0.0625*sm1*tp1*rp1

        		sq4(4, 1,iq) = - 0.0625*sm1*tm1*nm1
        		sq4(4, 2,iq) = - 0.0625*sp1*tm1*nm1
        		sq4(4, 3,iq) = - 0.0625*sp1*tp1*nm1
        		sq4(4, 4,iq) = - 0.0625*sm1*tp1*nm1
        		sq4(4, 5,iq) = - 0.0625*sm1*tm1*np1
        		sq4(4, 6,iq) = - 0.0625*sp1*tm1*np1
        		sq4(4, 7,iq) = - 0.0625*sp1*tp1*np1
        		sq4(4, 8,iq) = - 0.0625*sm1*tp1*np1
        		sq4(4, 9,iq) = + 0.0625*sm1*tm1*nm1
        		sq4(4,10,iq) = + 0.0625*sp1*tm1*nm1
        		sq4(4,11,iq) = + 0.0625*sp1*tp1*nm1
        		sq4(4,12,iq) = + 0.0625*sm1*tp1*nm1
        		sq4(4,13,iq) = + 0.0625*sm1*tm1*np1
        		sq4(4,14,iq) = + 0.0625*sp1*tm1*np1
        		sq4(4,15,iq) = + 0.0625*sp1*tp1*np1
        		sq4(4,16,iq) = + 0.0625*sm1*tp1*np1

		enddo

	endif	



 100	format(3f12.8)
	return
	end
