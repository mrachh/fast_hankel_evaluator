        implicit real *8 (a-h,o-z)
        complex *16, allocatable :: z0(:),h0(:),h1(:),h00(:),h01(:)
        real *8, allocatable :: rz0(:)
        real *8, allocatable :: w(:)

        call prini(6,13)
        
        n = 1000 000
        allocate(z0(n),h0(n),h1(n),h00(n),h01(n),rz0(n))
         
        do i=1,n
          z0(i) = (i-1.0d0)/(n+0.0d0)*20 + hkrand(0)/(n+0.0d0)
          z0(i) = 20.0d0*hkrand(0)
          rz0(i) = real(z0(i))
        enddo

        z0(1) = 1.15d-6
        rz0(1)= real(z0(1))

        
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
        lw = 50000
        allocate(w(lw))
        call hank106datagen_r(w,lw,ier)
   
        intnum = 1
        call cpu_time(t1)
        do i=1,n
          call hank106_r_h0(rz0(i),h0(i),ifexpon,intnum,w,lw)
        enddo
        call cpu_time(t2)

        
        print *, "time taken hank106_fast=",t2-t1
        
        r = 0
        erra = 0
        ra = 0
        do i=1,n
          erra = erra + abs(h0(i)-h00(i))**2
          ra = ra + abs(h0(i))**2
        enddo
        erra = sqrt(erra/ra)
        print *, "erra h0=",erra

        do i=1,n
          call hank106_r_h1(rz0(i),h1(i),ifexpon,intnum,w,lw)
        enddo

        r = 0
        erra = 0
        ra = 0
        do i=1,n
          erra = erra + abs(h1(i)-h01(i))**2
          ra = ra + abs(h1(i))**2
        enddo
        erra = sqrt(erra/ra)
        print *, "erra h1=",erra

        do i=1,n
          call hank106_r_h01(rz0(i),h0(i),h1(i),ifexpon,intnum,w,lw)
        enddo

        r = 0
        erra = 0
        ra = 0
        do i=1,n
          erra = erra + abs(h1(i)-h01(i))**2
          erra = erra + abs(h0(i)-h00(i))**2
          ra = ra + abs(h1(i))**2
          ra = ra + abs(h0(i))**2
        enddo
        erra = sqrt(erra/ra)
        print *, "erra h1 and h0=",erra

        stop
        end
