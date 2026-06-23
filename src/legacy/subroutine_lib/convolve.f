chelp+
!----------------------------------------------
! Code to convolve a data array with a template
! 1-D only now. 
!                        --wr, 03 Nov, 2010
!----------------------------------------------
chelp-


       subroutine convolve(A,nin,T,nt,OutArr)

       implicit none
       real*4      A(*), T(*), OutArr(*) 
       !real*4      weight
       integer*4   nin, nt 
       integer*4   i, k, k1

       do i = 1,nin+nt-1
          OutArr(i) = 0.0
       enddo
       
       do i = 1,nin
         do k = 1,nt
            k1 = i + nt - k
            OutArr(k1) = OutArr(k1) + A(i)*T(k)
         enddo
       enddo
       !weight = 0.0
       !do i = 1,nt
       !   weight = weight + T(i)
       !enddo
       !do i = 1,nin+nt-1
       !   OutArr(i) = OutArr(i)/weight
       !enddo

       return
       end
