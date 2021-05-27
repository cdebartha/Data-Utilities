        subroutine saxpy(n,a,x,incx,y,incy)
        implicit none
        integer n,incx,incy
        real*8 a
        real*8 x(1),y(1)
        integer i

        do i = 1,n
                y(i) = y(i) + a*x(i)
        enddo

        return

        end

        function sdot(n,x,incx,y,incy)
        integer n,incx,incy
        real*8 sdot
        real*8 x(1),y(1)

        sdot = 0.0

        do i = 1,n
                sdot = sdot + x(i)*y(i)
        enddo

        return
        end



        subroutine equal(f,d,m)
        implicit none
        integer i,m
        real* 8 f(m),d(m)

        do i=1,m
                d(i) = f(i)
        enddo

        return
        end


         subroutine equal1(f,d,m)
        implicit none
        integer i,m
        integer  f(m),d(m)
                                      
                                                                        
        do i=1,m
                d(i) = f(i)
                write(44,*)d(i),f(i),i
        enddo
                                                                                                                   
        return
        end

        subroutine fclear(f,m)
        implicit none
        integer i,m
        real* 8 f(m)

        do i=1,m
                f(i) = 0.0
        enddo

        return
        end

!**********************************************************************
        subroutine mkblk(ien,block1,selm)

        implicit none
        include "global.h"
        include "fmalloc.h"

        integer ien(nen,nec)
        integer iel,i,node,j,k,lt1,lt2
        real*8  selm(nee,nee,nec), block(ndf,ndf,nnl)
        real*8  block1(ndf,ndf,nnc)

        pointer (blockptr,block)

        blockptr = malloc(ndf*ndf*nnl*fsize)

        call fclear(block,ndf*ndf*nnl)
        call fclear(block1,ndf*ndf*nnc)
        do iel = 1,nec
            do i = 1,nen
                node = ien(i,iel)
                do j = 1,ndf
                    do k = 1,ndf
                        lt1 = (i-1)*ndf + k
                        lt2 = (i-1)*ndf + j
                        block(k,j,node)=block(k,j,node)+selm(lt1,lt2,iel)
!                       if (k.ne.j) then
!                          block(k,j,node) = 0.0e0
!                       endif
                    end do
                end do
            end do
        end do
        
        call scatter(block,block1,ndf*ndf,1)
        
        call free(blockptr)
        return
        end

!**********************************************************************
       subroutine facblk(id,block,nnc)
        IMPLICIT REAL*8 (A-H,O-Z)

        parameter (ndf=5)
        real*8 id(ndf,nnc),i,j,k
        real*8 block(ndf,ndf,nnc) 

        do i = 1,nnc
            do j = 1,ndf

                if(id(j,i).eq.0) then
                    do k = 1,ndf
                        block(j,k,i) = 0.0
                        block(k,j,i) = 0.0
                    end do
                    block(j,j,i) = 1.0
                end if
            end do
            do j = 1,ndf-1
                do k = j+1,ndf
                    block(k,j,i) = block(k,j,i)/block(j,j,i)
                        do l = j+1,ndf
                            block(k,l,i) = block(k,l,i) - 
     $                          block(k,j,i)*block(j,l,i)
                        end do
                end do
            end do

            do j = ndf-1,1,-1
                do k = j+1,ndf
                    block(j,k,i) = block(j,k,i)/block(j,j,i)
               end do
            end do
        end do
        
        return
       end

******************************************************
        subroutine sc1blk(b1,b2,block,id,ien)
        implicit real*8 (A-H,O-Z)

        include "global.h"
        include "fmalloc.h"
!        parameter (ndf=5)
        real*8 b1(nqc), b2(nqc)
        real*8 block(ndf,ndf,nnc)
        real*8 id(ndf,nnc),dl(ndf)
        integer i,j,k,ien(nen,nec),idf,inode

        pointer (dlptr,dl)

        dlptr = malloc(ndf*fsize)
        !call fclear(dl,ndf)
!        do k = 1,nen
!           node=ien(k,1)
!        if(its.eq.nts) then
!        do i = 1,nnc
!           do j = 1,ndf
!                   write(500,*)b1((i-1)*ndf+j)
!           enddo
!        enddo
!        endif
!        enddo


        do i = 1,nnc
           do j = 1,ndf
              if(id(j,i).eq.0) then
                dl(j) = 0.0
              else
                dl(j) = b1((i-1)*ndf+j)
              endif
           enddo
           do j = 2,ndf
              do k = 1,j-1
                 dl(j) = dl(j) - dl(k)*block(j,k,i)
              enddo
           enddo

           do j = 1,ndf
              dl(j) = dl(j)/block(j,j,i)
           enddo

           do j = ndf-1,1,-1
              do k = ndf,j+1,-1
                dl(j) = dl(j) - dl(k)*block(j,k,i)
              enddo
           enddo

           do j = 1,ndf
              if(id(j,i).ne.0) then
                b2((i-1)*ndf+j) = dl(j)
              endif
           enddo
        enddo

        call free(dlptr)
        return
        end

*******************************************************************                
