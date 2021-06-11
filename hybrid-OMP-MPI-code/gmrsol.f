        program gmrsol(lselm,idn,p,ndf,nen,munel,nee,neq,ikg,epsgm
     &    ,ngmr,block,sol,vv,w,wk1,wk2)

        include "global.h"

        dimension idn(ndf,numnp),ien(nen,numel)
        dimension lm(ndf,nen,numel)
        integer i,j,k


        do i = 1,ndf
           do j = 1, nen
              do k = 1, numel
                lm(i,j,k) = idn(i,ien(j,k))
              enddo
           enddo
        enddo
        call addrhs(b,p,lm(1,1,numel),nee)
        
        icode = 0
1       continue
        call cgmres(neq,ikg,rhs,sol,i,vv,w,wk1(1),wk2(1),epsgm,ngmr,iout
     &             ,icode)
                
        if(icode.eq.1) then
           wk2(0) = 0.0e0
           call sc1blk(wk(1),wk2(1),block,id,neq)
           goto 1
        else if(icode.eq.2) then
           do ii = 0,neq
               wk2(ii) = 0.0e0
           enddo
           call geres(numel,selm,wk1(1),nee,lm,wk2(1))
           goto 1
        else
           do ii=1,neq
              b(ii) = sol(ii)
           end do
        endif
        return
        end


*********************************************************        
        subroutine addrhs(b,p,lm,nee)
        IMPLICIT REAL*8(A-H,O-Z)
        dimension b(1),p(1),lm(1)
        
        do i = 1, nee
           if(lm(i).eq.0) b(lm(i) = b(lm(i)) + p(i)
        enddo
        return
        end        
