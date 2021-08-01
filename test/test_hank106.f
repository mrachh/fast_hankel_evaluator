        implicit real *8 (a-h,o-z)
        complex *16, allocatable :: z0(:),h0(:),h1(:),h00(:),h01(:)
        
        n = 1 000 000
        call prini(6,13)
        allocate(z0(n),h0(n),h1(n),h00(n),h01(n))
        
        do i=1,n
          z0(i) = (i-1.0d0)/(n+0.0d0)*20 + hkrand(0)/(n+0.0d0) 
        enddo
        z0(1) = 1.15d-6
        
        ifexpon = 1
        call cpu_time(t3)
        do i=1,n
          call hank103(z0(i),h00(i),h01(i),ifexpon)
        enddo
        call cpu_time(t4)
        
        print *, "time taken hank103=",t4-t3
        
        r = 0
        do i=1,n
          r = r + abs(h00(i))*abs(h01(i))
        enddo

        ifexpon = 1
        ier = 0
        call hank106datagen(z0(1),ier)
        
        call cpu_time(t1)
        do i=1,n
          call hank106(z0(i),h0(i),h1(i),ifexpon)
        enddo
        call cpu_time(t2)
        
        print *, "time taken hank106=",t2-t1
        
        r = 0
        erra = 0
        ra = 0
        do i=1,n
          r = r + abs(h0(i))*abs(h1(i))
          erra = erra + abs(h0(i)-h00(i))**2
          ra = ra + abs(h0(i))**2
        enddo
        erra = sqrt(erra/ra)
        print *, "erra=",erra

        stop
        end
