 function max_(a,b)
        implicit none
        integer a,b,max_
        if(a.gt.b) then
            max_=a
        else
            max_=b
        end if
    end

    subroutine calc_lusq(ia,di,al,au)
        implicit none
        common/matrix/n,nn
        integer n,nn,ia(n+1),i,j,k,ind,ind1,ind2,num_elem_i,num_elem_j,max_
        real*4 al(nn),au(nn),di(n),sum_l,sum_u,sum_
        do i=1,n
            num_elem_i=ia(i+1)-ia(i) ! к-во эл-тов в i строке/столбце
            ind=ia(i) ! индекс в массиве al
            sum_=0d0 ! сумма для диагонали
            do j=i-num_elem_i,i-1
                sum_l=0d0
                sum_u=0d0
                num_elem_j=ia(j+1)-ia(j) ! к-во эл-тов в j строке/столбце
                k=max_(i-num_elem_i,j-num_elem_j)
                ind1=ia(i)-(i-num_elem_i)+k
                ind2=ia(j)-(j-num_elem_j)+k
                k=j-1-k+ind1
                do ind1=ind1,k
                    sum_l=sum_l+al(ind1)*au(ind2)
                    sum_u=sum_u+au(ind1)*al(ind2)
                    ind2=ind2+1
                end do
                al(ind)=(al(ind)-sum_l)/di(j)
                au(ind)=(au(ind)-sum_u)/di(j)
                sum_=sum_+al(ind)*au(ind)
                ind=ind+1
            end do
            if(di(i)-sum_.le.0d0) then
                print*, 'Matrix is NOT LU(sq) decomposable!'
                stop
            end if
            di(i)=sqrt(di(i)-sum_)
        end do
    end

    subroutine calc_y(ia,di,al,vect)
        implicit none
        common/matrix/n,nn
        integer n,nn,ia(n+1),i,j,ind,num_elem_i
        real*4 al(nn),vect(n),di(n),sum_
             do i=1,n
            sum_=0d0
            num_elem_i=ia(i+1)-ia(i)
            ind=ia(i)
            do j=i-num_elem_i,i-1
                sum_=sum_+al(ind)*vect(j)
                ind=ind+1
            end do
            vect(i)=(vect(i)-sum_)/di(i)
        end do
    end

    subroutine calc_x(ia,di,au,vect)
        implicit none
        common/matrix/n,nn
        integer n,nn,ia(n+1),i,j,ind,num_elem_i
        real*4 au(nn),vect(n),di(n),xi
        do i=n,1,-1
            xi=vect(i)
            xi=xi/di(i)
            vect(i)=xi
            ind=ia(i+1)-1
            num_elem_i=ia(i+1)-ia(i)
            do j=i-1,i-num_elem_i,-1
                vect(j)=vect(j)-au(ind)*xi
                ind=ind-1
            end do
        end do
    end

    subroutine read_(mem,ia,di,al,au,vect)
        implicit none
        common/matrix/n,nn
        integer n,nn,ia(4097),di,al,au,vect,i
        real*4 mem(16781312)
        read(10,*) n
        read(10,*) (ia(i), i=1,n+1)
        nn=ia(n+1)-1
        di=0
        read(10,*) (mem(di+i), i=1,n)
        al=di+n+1
        read(10,*) (mem(al+i), i=1,nn)
        au=al+nn
        read(10,*) (mem(au+i), i=1,nn)
        vect=au+nn
        read(20,*) (mem(vect+i), i=1,n)
    end

    subroutine write_(vect)
        implicit none
        common/matrix/n,nn
        integer n,nn,i
        real*4 vect(n)
        write(30,*) (vect(i), i=1,n)
    end

    program main
        implicit none
        common/matrix/n,nn
        integer n,nn,ia(4097),di,al,au,vect,i
        real*4 mem(16781312) !di(4096),al(8386560),au(8386560),vect(4096)
        open(10,file='../profile.txt',status='old',err=101)
        open(20,file='../vector.txt',status='old',err=101)
        call read_(mem,ia,di,al,au,vect)
        close(10)
        close(20)
        call calc_lusq(ia,mem(di+1),mem(al+1),mem(au+1))
        call calc_y(ia,mem(di+1),mem(al+1),mem(vect+1))
        call calc_x(ia,mem(di+1),mem(au+1),mem(vect+1))
        open(30,file='../tmp3.txt',status='unknown')
        call write_(mem(vect))
        close(30)
              do i=1,n
            print*,mem(vect+i)
        end do
        goto 300
101     print*,'ERROR: Can`t reading profile.txt'
300     continue
    end
