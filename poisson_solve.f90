
subroutine fftstp(mm,nfft,m,nn,n,zin,zout,trig,after,now,before,isign)
  implicit real(kind=8) (a-h,o-z)
  integer, parameter :: nfft_max=24000
! Arguments
  integer :: mm,nfft,m,nn,n,after,before,isign,now,ia,ias,ib,itrig,itt,j
  integer :: nin1,nin2,nin3,nin4,nin5,nin6,nin7,nin8,nout1,nout2,nout3,nout4,nout5
  integer :: nout6,nout7,nout8
  real(kind=8) :: trig(2,nfft_max),zin(2,mm,m),zout(2,nn,n)
! Local variables
  integer :: atn,atb
  atn=after*now
  atb=after*before

!         sqrt(.5d0)
  rt2i=0.7071067811865475d0

! First now == 2
  if (now.eq.2) then
     ia=1
     nin1=ia-after
     nout1=ia-atn
     do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nout1=nout1+atn
        nout2=nout1+after
        do j=1,nfft
           r1=zin(1,j,nin1)
           s1=zin(2,j,nin1)
           r2=zin(1,j,nin2)
           s2=zin(2,j,nin2)
           zout(1,j,nout1)= r2 + r1
           zout(2,j,nout1)= s2 + s1
           zout(1,j,nout2)= r1 - r2
           zout(2,j,nout2)= s1 - s2
        end do
     end do

   ! Big loop
     Big_Loop: do ia=2,after
        ias=ia-1
        if (2*ias.eq.after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r2=zin(2,j,nin2)
                    s2=zin(1,j,nin2)
                    zout(1,j,nout1)= r1 - r2
                    zout(2,j,nout1)= s2 + s1
                    zout(1,j,nout2)= r2 + r1
                    zout(2,j,nout2)= s1 - s2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r2=zin(2,j,nin2)
                    s2=zin(1,j,nin2)
                    zout(1,j,nout1)= r2 + r1
                    zout(2,j,nout1)= s1 - s2
                    zout(1,j,nout2)= r1 - r2
                    zout(2,j,nout2)= s2 + s1
                 end do
              end do
           endif
        else if (4*ias.eq.after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r - s)*rt2i
                    s2=(r + s)*rt2i
                    zout(1,j,nout1)= r2 + r1
                    zout(2,j,nout1)= s2 + s1
                    zout(1,j,nout2)= r1 - r2
                    zout(2,j,nout2)= s1 - s2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(s - r)*rt2i
                    zout(1,j,nout1)= r2 + r1
                    zout(2,j,nout1)= s2 + s1
                    zout(1,j,nout2)= r1 - r2
                    zout(2,j,nout2)= s1 - s2
                 end do
              end do
           endif
        else if (4*ias.eq.3*after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(r - s)*rt2i
                    zout(1,j,nout1)= r1 - r2
                    zout(2,j,nout1)= s2 + s1
                    zout(1,j,nout2)= r2 + r1
                    zout(2,j,nout2)= s1 - s2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(s - r)*rt2i
                    s2=(r + s)*rt2i
                    zout(1,j,nout1)= r2 + r1
                    zout(2,j,nout1)= s1 - s2
                    zout(1,j,nout2)= r1 - r2
                    zout(2,j,nout2)= s2 + s1
                 end do
              end do
           endif
        else
           itrig=ias*before+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nout1=nout1+atn
              nout2=nout1+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 zout(1,j,nout1)= r2 + r1
                 zout(2,j,nout1)= s2 + s1
                 zout(1,j,nout2)= r1 - r2
                 zout(2,j,nout2)= s1 - s2
              end do
           end do
        endif

     end do Big_Loop
!    End of Big_loop
! End of if (now.eq.2)

  else if (now.eq.4) then
     if (isign.eq.1) then
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
           nin1=nin1+after
           nin2=nin1+atb
           nin3=nin2+atb
           nin4=nin3+atb
           nout1=nout1+atn
           nout2=nout1+after
           nout3=nout2+after
           nout4=nout3+after
           do j=1,nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r2=zin(1,j,nin2)
              s2=zin(2,j,nin2)
              r3=zin(1,j,nin3)
              s3=zin(2,j,nin3)
              r4=zin(1,j,nin4)
              s4=zin(2,j,nin4)
              r=r1 + r3
              s=r2 + r4
              zout(1,j,nout1) = r + s
              zout(1,j,nout3) = r - s
              r=r1 - r3
              s=s2 - s4
              zout(1,j,nout2) = r - s
              zout(1,j,nout4) = r + s
              r=s1 + s3
              s=s2 + s4
              zout(2,j,nout1) = r + s
              zout(2,j,nout3) = r - s
              r=s1 - s3
              s=r2 - r4
              zout(2,j,nout2) = r + s
              zout(2,j,nout4) = r - s
           end do
        end do
        do ia=2,after
           ias=ia-1
           if (2*ias.eq.after) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r-s)*rt2i
                    s2=(r+s)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=(r + s)*rt2i
                    s4=(r - s)*rt2i
                    r=r1 - r3
                    s=r2 - r4
                    zout(1,j,nout1) = r + s
                    zout(1,j,nout3) = r - s
                    r=r1 + r3
                    s=s2 - s4
                    zout(1,j,nout2) = r - s
                    zout(1,j,nout4) = r + s
                    r=s1 + s3
                    s=s2 + s4
                    zout(2,j,nout1) = r + s
                    zout(2,j,nout3) = r - s
                    r=s1 - s3
                    s=r2 + r4
                    zout(2,j,nout2) = r + s
                    zout(2,j,nout4) = r - s
                 end do
              end do
           else
              itt=ias*before
              itrig=itt+1
              cr2=trig(1,itrig)
              ci2=trig(2,itrig)
              itrig=itrig+itt
              cr3=trig(1,itrig)
              ci3=trig(2,itrig)
              itrig=itrig+itt
              cr4=trig(1,itrig)
              ci4=trig(2,itrig)
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=r*cr2 - s*ci2
                    s2=r*ci2 + s*cr2
                    r=zin(1,j,nin3)
                    s=zin(2,j,nin3)
                    r3=r*cr3 - s*ci3
                    s3=r*ci3 + s*cr3
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=r*cr4 - s*ci4
                    s4=r*ci4 + s*cr4
                    r=r1 + r3
                    s=r2 + r4
                    zout(1,j,nout1) = r + s
                    zout(1,j,nout3) = r - s
                    r=r1 - r3
                    s=s2 - s4
                    zout(1,j,nout2) = r - s
                    zout(1,j,nout4) = r + s
                    r=s1 + s3
                    s=s2 + s4
                    zout(2,j,nout1) = r + s
                    zout(2,j,nout3) = r - s
                    r=s1 - s3
                    s=r2 - r4
                    zout(2,j,nout2) = r + s
                    zout(2,j,nout4) = r - s
                 end do
              end do
           endif
        end do
     else
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
           nin1=nin1+after
           nin2=nin1+atb
           nin3=nin2+atb
           nin4=nin3+atb
           nout1=nout1+atn
           nout2=nout1+after
           nout3=nout2+after
           nout4=nout3+after
           do j=1,nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r2=zin(1,j,nin2)
              s2=zin(2,j,nin2)
              r3=zin(1,j,nin3)
              s3=zin(2,j,nin3)
              r4=zin(1,j,nin4)
              s4=zin(2,j,nin4)
              r=r1 + r3
              s=r2 + r4
              zout(1,j,nout1) = r + s
              zout(1,j,nout3) = r - s
              r=r1 - r3
              s=s2 - s4
              zout(1,j,nout2) = r + s
              zout(1,j,nout4) = r - s
              r=s1 + s3
              s=s2 + s4
              zout(2,j,nout1) = r + s
              zout(2,j,nout3) = r - s
              r=s1 - s3
              s=r2 - r4
              zout(2,j,nout2) = r - s
              zout(2,j,nout4) = r + s
           end do
        end do
        do ia=2,after
           ias=ia-1
           if (2*ias.eq.after) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(s - r)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=(s - r)*rt2i
                    s4=(r + s)*rt2i
                    r=r1 + r3
                    s=r2 + r4
                    zout(1,j,nout1) = r + s
                    zout(1,j,nout3) = r - s
                    r=r1 - r3
                    s=s2 + s4
                    zout(1,j,nout2) = r + s
                    zout(1,j,nout4) = r - s
                    r=s1 - s3
                    s=s2 - s4
                    zout(2,j,nout1) = r + s
                    zout(2,j,nout3) = r - s
                    r=s1 + s3
                    s=r2 - r4
                    zout(2,j,nout2) = r - s
                    zout(2,j,nout4) = r + s
                 end do
              end do
           else
              itt=ias*before
              itrig=itt+1
              cr2=trig(1,itrig)
              ci2=trig(2,itrig)
              itrig=itrig+itt
              cr3=trig(1,itrig)
              ci3=trig(2,itrig)
              itrig=itrig+itt
              cr4=trig(1,itrig)
              ci4=trig(2,itrig)
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=r*cr2 - s*ci2
                    s2=r*ci2 + s*cr2
                    r=zin(1,j,nin3)
                    s=zin(2,j,nin3)
                    r3=r*cr3 - s*ci3
                    s3=r*ci3 + s*cr3
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=r*cr4 - s*ci4
                    s4=r*ci4 + s*cr4
                    r=r1 + r3
                    s=r2 + r4
                    zout(1,j,nout1) = r + s
                    zout(1,j,nout3) = r - s
                    r=r1 - r3
                    s=s2 - s4
                    zout(1,j,nout2) = r + s
                    zout(1,j,nout4) = r - s
                    r=s1 + s3
                    s=s2 + s4
                    zout(2,j,nout1) = r + s
                    zout(2,j,nout3) = r - s
                    r=s1 - s3
                    s=r2 - r4
                    zout(2,j,nout2) = r - s
                    zout(2,j,nout4) = r + s
                 end do
              end do
           endif
        end do
     endif
! End of else if (now.eq.4)
  else if (now.eq.8) then
     if (isign.eq.-1) then
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
           nin1=nin1+after
           nin2=nin1+atb
           nin3=nin2+atb
           nin4=nin3+atb
           nin5=nin4+atb
           nin6=nin5+atb
           nin7=nin6+atb
           nin8=nin7+atb
           nout1=nout1+atn
           nout2=nout1+after
           nout3=nout2+after
           nout4=nout3+after
           nout5=nout4+after
           nout6=nout5+after
           nout7=nout6+after
           nout8=nout7+after
           do j=1,nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r2=zin(1,j,nin2)
              s2=zin(2,j,nin2)
              r3=zin(1,j,nin3)
              s3=zin(2,j,nin3)
              r4=zin(1,j,nin4)
              s4=zin(2,j,nin4)
              r5=zin(1,j,nin5)
              s5=zin(2,j,nin5)
              r6=zin(1,j,nin6)
              s6=zin(2,j,nin6)
              r7=zin(1,j,nin7)
              s7=zin(2,j,nin7)
              r8=zin(1,j,nin8)
              s8=zin(2,j,nin8)
              r=r1 + r5
              s=r3 + r7
              ap=r + s
              am=r - s
              r=r2 + r6
              s=r4 + r8
              bp=r + s
              bm=r - s
              r=s1 + s5
              s=s3 + s7
              cp=r + s
              cm=r - s
              r=s2 + s6
              s=s4 + s8
              dp=r + s
              dm=r - s
              zout(1,j,nout1) = ap + bp
              zout(2,j,nout1) = cp + dp
              zout(1,j,nout5) = ap - bp
              zout(2,j,nout5) = cp - dp
              zout(1,j,nout3) = am + dm
              zout(2,j,nout3) = cm - bm
              zout(1,j,nout7) = am - dm
              zout(2,j,nout7) = cm + bm
              r=r1 - r5
              s=s3 - s7
              ap=r + s
              am=r - s
              r=s1 - s5
              s=r3 - r7
              bp=r + s
              bm=r - s
              r=s4 - s8
              s=r2 - r6
              cp=r + s
              cm=r - s
              r=s2 - s6
              s=r4 - r8
              dp=r + s
              dm=r - s
              r = ( cp + dm)*rt2i
              s = ( dm - cp)*rt2i
              cp= ( cm + dp)*rt2i
              dp = ( cm - dp)*rt2i
              zout(1,j,nout2) = ap + r
              zout(2,j,nout2) = bm + s
              zout(1,j,nout6) = ap - r
              zout(2,j,nout6) = bm - s
              zout(1,j,nout4) = am + cp
              zout(2,j,nout4) = bp + dp
              zout(1,j,nout8) = am - cp
              zout(2,j,nout8) = bp - dp
           end do
        end do
        do ia=2,after
           ias=ia-1
           itt=ias*before
           itrig=itt+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           itrig=itrig+itt
           cr3=trig(1,itrig)
           ci3=trig(2,itrig)
           itrig=itrig+itt
           cr4=trig(1,itrig)
           ci4=trig(2,itrig)
           itrig=itrig+itt
           cr5=trig(1,itrig)
           ci5=trig(2,itrig)
           itrig=itrig+itt
           cr6=trig(1,itrig)
           ci6=trig(2,itrig)
           itrig=itrig+itt
           cr7=trig(1,itrig)
           ci7=trig(2,itrig)
           itrig=itrig+itt
           cr8=trig(1,itrig)
           ci8=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nin5=nin4+atb
              nin6=nin5+atb
              nin7=nin6+atb
              nin8=nin7+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              nout5=nout4+after
              nout6=nout5+after
              nout7=nout6+after
              nout8=nout7+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 r=zin(1,j,nin3)
                 s=zin(2,j,nin3)
                 r3=r*cr3 - s*ci3
                 s3=r*ci3 + s*cr3
                 r=zin(1,j,nin4)
                 s=zin(2,j,nin4)
                 r4=r*cr4 - s*ci4
                 s4=r*ci4 + s*cr4
                 r=zin(1,j,nin5)
                 s=zin(2,j,nin5)
                 r5=r*cr5 - s*ci5
                 s5=r*ci5 + s*cr5
                 r=zin(1,j,nin6)
                 s=zin(2,j,nin6)
                 r6=r*cr6 - s*ci6
                 s6=r*ci6 + s*cr6
                 r=zin(1,j,nin7)
                 s=zin(2,j,nin7)
                 r7=r*cr7 - s*ci7
                 s7=r*ci7 + s*cr7
                 r=zin(1,j,nin8)
                 s=zin(2,j,nin8)
                 r8=r*cr8 - s*ci8
                 s8=r*ci8 + s*cr8
                 r=r1 + r5
                 s=r3 + r7
                 ap=r + s
                 am=r - s
                 r=r2 + r6
                 s=r4 + r8
                 bp=r + s
                 bm=r - s
                 r=s1 + s5
                 s=s3 + s7
                 cp=r + s
                 cm=r - s
                 r=s2 + s6
                 s=s4 + s8
                 dp=r + s
                 dm=r - s
                 zout(1,j,nout1) = ap + bp
                 zout(2,j,nout1) = cp + dp
                 zout(1,j,nout5) = ap - bp
                 zout(2,j,nout5) = cp - dp
                 zout(1,j,nout3) = am + dm
                 zout(2,j,nout3) = cm - bm
                 zout(1,j,nout7) = am - dm
                 zout(2,j,nout7) = cm + bm
                 r=r1 - r5
                 s=s3 - s7
                 ap=r + s
                 am=r - s
                 r=s1 - s5
                 s=r3 - r7
                 bp=r + s
                 bm=r - s
                 r=s4 - s8
                 s=r2 - r6
                 cp=r + s
                 cm=r - s
                 r=s2 - s6
                 s=r4 - r8
                 dp=r + s
                 dm=r - s
                 r = ( cp + dm)*rt2i
                 s = ( dm - cp)*rt2i
                 cp= ( cm + dp)*rt2i
                 dp = ( cm - dp)*rt2i
                 zout(1,j,nout2) = ap + r
                 zout(2,j,nout2) = bm + s
                 zout(1,j,nout6) = ap - r
                 zout(2,j,nout6) = bm - s
                 zout(1,j,nout4) = am + cp
                 zout(2,j,nout4) = bp + dp
                 zout(1,j,nout8) = am - cp
                 zout(2,j,nout8) = bp - dp
              end do
           end do
        end do
     else ! else for isign.eq.1
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
           nin1=nin1+after
           nin2=nin1+atb
           nin3=nin2+atb
           nin4=nin3+atb
           nin5=nin4+atb
           nin6=nin5+atb
           nin7=nin6+atb
           nin8=nin7+atb
           nout1=nout1+atn
           nout2=nout1+after
           nout3=nout2+after
           nout4=nout3+after
           nout5=nout4+after
           nout6=nout5+after
           nout7=nout6+after
           nout8=nout7+after
           do j=1,nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r2=zin(1,j,nin2)
              s2=zin(2,j,nin2)
              r3=zin(1,j,nin3)
              s3=zin(2,j,nin3)
              r4=zin(1,j,nin4)
              s4=zin(2,j,nin4)
              r5=zin(1,j,nin5)
              s5=zin(2,j,nin5)
              r6=zin(1,j,nin6)
              s6=zin(2,j,nin6)
              r7=zin(1,j,nin7)
              s7=zin(2,j,nin7)
              r8=zin(1,j,nin8)
              s8=zin(2,j,nin8)
              r=r1 + r5
              s=r3 + r7
              ap=r + s
              am=r - s
              r=r2 + r6
              s=r4 + r8
              bp=r + s
              bm=r - s
              r=s1 + s5
              s=s3 + s7
              cp=r + s
              cm=r - s
              r=s2 + s6
              s=s4 + s8
              dp=r + s
              dm=r - s
              zout(1,j,nout1) = ap + bp
              zout(2,j,nout1) = cp + dp
              zout(1,j,nout5) = ap - bp
              zout(2,j,nout5) = cp - dp
              zout(1,j,nout3) = am - dm
              zout(2,j,nout3) = cm + bm
              zout(1,j,nout7) = am + dm
              zout(2,j,nout7) = cm - bm
              r= r1 - r5
              s=-s3 + s7
              ap=r + s
              am=r - s
              r=s1 - s5
              s=r7 - r3
              bp=r + s
              bm=r - s
              r=-s4 + s8
              s= r2 - r6
              cp=r + s
              cm=r - s
              r=-s2 + s6
              s= r4 - r8
              dp=r + s
              dm=r - s
              r = ( cp + dm)*rt2i
              s = ( cp - dm)*rt2i
              cp= ( cm + dp)*rt2i
              dp= ( dp - cm)*rt2i
              zout(1,j,nout2) = ap + r
              zout(2,j,nout2) = bm + s
              zout(1,j,nout6) = ap - r
              zout(2,j,nout6) = bm - s
              zout(1,j,nout4) = am + cp
              zout(2,j,nout4) = bp + dp
              zout(1,j,nout8) = am - cp
              zout(2,j,nout8) = bp - dp
           end do
        end do
        do ia=2,after
           ias=ia-1
           itt=ias*before
           itrig=itt+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           itrig=itrig+itt
           cr3=trig(1,itrig)
           ci3=trig(2,itrig)
           itrig=itrig+itt
           cr4=trig(1,itrig)
           ci4=trig(2,itrig)
           itrig=itrig+itt
           cr5=trig(1,itrig)
           ci5=trig(2,itrig)
           itrig=itrig+itt
           cr6=trig(1,itrig)
           ci6=trig(2,itrig)
           itrig=itrig+itt
           cr7=trig(1,itrig)
           ci7=trig(2,itrig)
           itrig=itrig+itt
           cr8=trig(1,itrig)
           ci8=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nin5=nin4+atb
              nin6=nin5+atb
              nin7=nin6+atb
              nin8=nin7+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              nout5=nout4+after
              nout6=nout5+after
              nout7=nout6+after
              nout8=nout7+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 r=zin(1,j,nin3)
                 s=zin(2,j,nin3)
                 r3=r*cr3 - s*ci3
                 s3=r*ci3 + s*cr3
                 r=zin(1,j,nin4)
                 s=zin(2,j,nin4)
                 r4=r*cr4 - s*ci4
                 s4=r*ci4 + s*cr4
                 r=zin(1,j,nin5)
                 s=zin(2,j,nin5)
                 r5=r*cr5 - s*ci5
                 s5=r*ci5 + s*cr5
                 r=zin(1,j,nin6)
                 s=zin(2,j,nin6)
                 r6=r*cr6 - s*ci6
                 s6=r*ci6 + s*cr6
                 r=zin(1,j,nin7)
                 s=zin(2,j,nin7)
                 r7=r*cr7 - s*ci7
                 s7=r*ci7 + s*cr7
                 r=zin(1,j,nin8)
                 s=zin(2,j,nin8)
                 r8=r*cr8 - s*ci8
                 s8=r*ci8 + s*cr8
                 r=r1 + r5
                 s=r3 + r7
                 ap=r + s
                 am=r - s
                 r=r2 + r6
                 s=r4 + r8
                 bp=r + s
                 bm=r - s
                 r=s1 + s5
                 s=s3 + s7
                 cp=r + s
                 cm=r - s
                 r=s2 + s6
                 s=s4 + s8
                 dp=r + s
                 dm=r - s
                 zout(1,j,nout1) = ap + bp
                 zout(2,j,nout1) = cp + dp
                 zout(1,j,nout5) = ap - bp
                 zout(2,j,nout5) = cp - dp
                 zout(1,j,nout3) = am - dm
                 zout(2,j,nout3) = cm + bm
                 zout(1,j,nout7) = am + dm
                 zout(2,j,nout7) = cm - bm
                 r= r1 - r5
                 s=-s3 + s7
                 ap=r + s
                 am=r - s
                 r=s1 - s5
                 s=r7 - r3
                 bp=r + s
                 bm=r - s
                 r=-s4 + s8
                 s= r2 - r6
                 cp=r + s
                 cm=r - s
                 r=-s2 + s6
                 s= r4 - r8
                 dp=r + s
                 dm=r - s
                 r = ( cp + dm)*rt2i
                 s = ( cp - dm)*rt2i
                 cp= ( cm + dp)*rt2i
                 dp= ( dp - cm)*rt2i
                 zout(1,j,nout2) = ap + r
                 zout(2,j,nout2) = bm + s
                 zout(1,j,nout6) = ap - r
                 zout(2,j,nout6) = bm - s
                 zout(1,j,nout4) = am + cp
                 zout(2,j,nout4) = bp + dp
                 zout(1,j,nout8) = am - cp
                 zout(2,j,nout8) = bp - dp
              end do
           end do
        end do
     endif  !end if of isign

  else if (now.eq.3) then
!         .5d0*sqrt(3.d0)
     bb=real(isign,kind=8)*0.8660254037844387d0
     ia=1
     nin1=ia-after
     nout1=ia-atn
     do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do j=1,nfft
           r1=zin(1,j,nin1)
           s1=zin(2,j,nin1)
           r2=zin(1,j,nin2)
           s2=zin(2,j,nin2)
           r3=zin(1,j,nin3)
           s3=zin(2,j,nin3)
           r=r2 + r3
           s=s2 + s3
           zout(1,j,nout1) = r + r1
           zout(2,j,nout1) = s + s1
           r1=r1 - .5d0*r
           s1=s1 - .5d0*s
           r2=bb*(r2-r3)
           s2=bb*(s2-s3)
           zout(1,j,nout2) = r1 - s2
           zout(2,j,nout2) = s1 + r2
           zout(1,j,nout3) = r1 + s2
           zout(2,j,nout3) = s1 - r2
        end do
     end do
     loop_3000: do ia=2,after
        ias=ia-1
        if (4*ias.eq.3*after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r2=zin(2,j,nin2)
                    s2=zin(1,j,nin2)
                    r3=zin(1,j,nin3)
                    s3=zin(2,j,nin3)
                    r=r3 + r2
                    s=s2 - s3
                    zout(1,j,nout1) = r1 - r
                    zout(2,j,nout1) = s + s1
                    r1=r1 + .5d0*r
                    s1=s1 - .5d0*s
                    r2=bb*(r2-r3)
                    s2=bb*(s2+s3)
                    zout(1,j,nout2) = r1 - s2
                    zout(2,j,nout2) = s1 - r2
                    zout(1,j,nout3) = r1 + s2
                    zout(2,j,nout3) = s1 + r2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r2=zin(2,j,nin2)
                    s2=zin(1,j,nin2)
                    r3=zin(1,j,nin3)
                    s3=zin(2,j,nin3)
                    r=r2 - r3
                    s=s2 + s3
                    zout(1,j,nout1) = r + r1
                    zout(2,j,nout1) = s1 - s
                    r1=r1 - .5d0*r
                    s1=s1 + .5d0*s
                    r2=bb*(r2+r3)
                    s2=bb*(s2-s3)
                    zout(1,j,nout2) = r1 + s2
                    zout(2,j,nout2) = s1 + r2
                    zout(1,j,nout3) = r1 - s2
                    zout(2,j,nout3) = s1 - r2
                 end do
              end do
           endif
        else if (8*ias.eq.3*after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r - s)*rt2i
                    s2=(r + s)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=r2 - r3
                    s=s2 + s3
                    zout(1,j,nout1) = r + r1
                    zout(2,j,nout1) = s + s1
                    r1=r1 - .5d0*r
                    s1=s1 - .5d0*s
                    r2=bb*(r2+r3)
                    s2=bb*(s2-s3)
                    zout(1,j,nout2) = r1 - s2
                    zout(2,j,nout2) = s1 + r2
                    zout(1,j,nout3) = r1 + s2
                    zout(2,j,nout3) = s1 - r2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(s - r)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=r2 + r3
                    s=s2 - s3
                    zout(1,j,nout1) = r + r1
                    zout(2,j,nout1) = s + s1
                    r1=r1 - .5d0*r
                    s1=s1 - .5d0*s
                    r2=bb*(r2-r3)
                    s2=bb*(s2+s3)
                    zout(1,j,nout2) = r1 - s2
                    zout(2,j,nout2) = s1 + r2
                    zout(1,j,nout3) = r1 + s2
                    zout(2,j,nout3) = s1 - r2
                 end do
              end do
           endif
        else
           itt=ias*before
           itrig=itt+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           itrig=itrig+itt
           cr3=trig(1,itrig)
           ci3=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 r=zin(1,j,nin3)
                 s=zin(2,j,nin3)
                 r3=r*cr3 - s*ci3
                 s3=r*ci3 + s*cr3
                 r=r2 + r3
                 s=s2 + s3
                 zout(1,j,nout1) = r + r1
                 zout(2,j,nout1) = s + s1
                 r1=r1 - .5d0*r
                 s1=s1 - .5d0*s
                 r2=bb*(r2-r3)
                 s2=bb*(s2-s3)
                 zout(1,j,nout2) = r1 - s2
                 zout(2,j,nout2) = s1 + r2
                 zout(1,j,nout3) = r1 + s2
                 zout(2,j,nout3) = s1 - r2
              end do
           end do
        endif
     end do loop_3000
! End of if (now.eq.3)

  else if (now.eq.5) then
!     cos(2.d0*pi/5.d0)
     cos2=0.3090169943749474d0
!     cos(4.d0*pi/5.d0)
     cos4=-0.8090169943749474d0
!     sin(2.d0*pi/5.d0)
     sin2=real(isign,kind=8)*0.9510565162951536d0
!      sin(4.d0*pi/5.d0)
     sin4=real(isign,kind=8)*0.5877852522924731d0
     ia=1
     nin1=ia-after
     nout1=ia-atn
     do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        do j=1,nfft
           r1=zin(1,j,nin1)
           s1=zin(2,j,nin1)
           r2=zin(1,j,nin2)
           s2=zin(2,j,nin2)
           r3=zin(1,j,nin3)
           s3=zin(2,j,nin3)
           r4=zin(1,j,nin4)
           s4=zin(2,j,nin4)
           r5=zin(1,j,nin5)
           s5=zin(2,j,nin5)
           r25 = r2 + r5
           r34 = r3 + r4
           s25 = s2 - s5
           s34 = s3 - s4
           zout(1,j,nout1) = r1 + r25 + r34
           r = r1 + cos2*r25 + cos4*r34
           s = sin2*s25 + sin4*s34
           zout(1,j,nout2) = r - s
           zout(1,j,nout5) = r + s
           r = r1 + cos4*r25 + cos2*r34
           s = sin4*s25 - sin2*s34
           zout(1,j,nout3) = r - s
           zout(1,j,nout4) = r + s
           r25 = r2 - r5
           r34 = r3 - r4
           s25 = s2 + s5
           s34 = s3 + s4
           zout(2,j,nout1) = s1 + s25 + s34
           r = s1 + cos2*s25 + cos4*s34
           s = sin2*r25 + sin4*r34
           zout(2,j,nout2) = r + s
           zout(2,j,nout5) = r - s
           r = s1 + cos4*s25 + cos2*s34
           s = sin4*r25 - sin2*r34
           zout(2,j,nout3) = r + s
           zout(2,j,nout4) = r - s
        end do
     end do
     loop_5000: do ia=2,after
        ias=ia-1
        if (8*ias.eq.5*after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nin5=nin4+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 nout5=nout4+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r - s)*rt2i
                    s2=(r + s)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=(r + s)*rt2i
                    s4=(r - s)*rt2i
                    r5=zin(1,j,nin5)
                    s5=zin(2,j,nin5)
                    r25 = r2 - r5
                    r34 = r3 + r4
                    s25 = s2 + s5
                    s34 = s3 - s4
                    zout(1,j,nout1) = r1 + r25 - r34
                    r = r1 + cos2*r25 - cos4*r34
                    s = sin2*s25 + sin4*s34
                    zout(1,j,nout2) = r - s
                    zout(1,j,nout5) = r + s
                    r = r1 + cos4*r25 - cos2*r34
                    s = sin4*s25 - sin2*s34
                    zout(1,j,nout3) = r - s
                    zout(1,j,nout4) = r + s
                    r25 = r2 + r5
                    r34 = r4 - r3
                    s25 = s2 - s5
                    s34 = s3 + s4
                    zout(2,j,nout1) = s1 + s25 + s34
                    r = s1 + cos2*s25 + cos4*s34
                    s = sin2*r25 + sin4*r34
                    zout(2,j,nout2) = r + s
                    zout(2,j,nout5) = r - s
                    r = s1 + cos4*s25 + cos2*s34
                    s = sin4*r25 - sin2*r34
                    zout(2,j,nout3) = r + s
                    zout(2,j,nout4) = r - s
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nin5=nin4+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 nout5=nout4+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(s - r)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=(s - r)*rt2i
                    s4=(r + s)*rt2i
                    r5=zin(1,j,nin5)
                    s5=zin(2,j,nin5)
                    r25 = r2 - r5
                    r34 = r3 + r4
                    s25 = s2 + s5
                    s34 = s4 - s3
                    zout(1,j,nout1) = r1 + r25 + r34
                    r = r1 + cos2*r25 + cos4*r34
                    s = sin2*s25 + sin4*s34
                    zout(1,j,nout2) = r - s
                    zout(1,j,nout5) = r + s
                    r = r1 + cos4*r25 + cos2*r34
                    s = sin4*s25 - sin2*s34
                    zout(1,j,nout3) = r - s
                    zout(1,j,nout4) = r + s
                    r25 = r2 + r5
                    r34 = r3 - r4
                    s25 = s2 - s5
                    s34 = s3 + s4
                    zout(2,j,nout1) = s1 + s25 - s34
                    r = s1 + cos2*s25 - cos4*s34
                    s = sin2*r25 + sin4*r34
                    zout(2,j,nout2) = r + s
                    zout(2,j,nout5) = r - s
                    r = s1 + cos4*s25 - cos2*s34
                    s = sin4*r25 - sin2*r34
                    zout(2,j,nout3) = r + s
                    zout(2,j,nout4) = r - s
                 end do
              end do
           endif
        else !if of (ias...
           ias=ia-1
           itt=ias*before
           itrig=itt+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           itrig=itrig+itt
           cr3=trig(1,itrig)
           ci3=trig(2,itrig)
           itrig=itrig+itt
           cr4=trig(1,itrig)
           ci4=trig(2,itrig)
           itrig=itrig+itt
           cr5=trig(1,itrig)
           ci5=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nin5=nin4+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              nout5=nout4+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 r=zin(1,j,nin3)
                 s=zin(2,j,nin3)
                 r3=r*cr3 - s*ci3
                 s3=r*ci3 + s*cr3
                 r=zin(1,j,nin4)
                 s=zin(2,j,nin4)
                 r4=r*cr4 - s*ci4
                 s4=r*ci4 + s*cr4
                 r=zin(1,j,nin5)
                 s=zin(2,j,nin5)
                 r5=r*cr5 - s*ci5
                 s5=r*ci5 + s*cr5
                 r25 = r2 + r5
                 r34 = r3 + r4
                 s25 = s2 - s5
                 s34 = s3 - s4
                 zout(1,j,nout1) = r1 + r25 + r34
                 r = r1 + cos2*r25 + cos4*r34
                 s = sin2*s25 + sin4*s34
                 zout(1,j,nout2) = r - s
                 zout(1,j,nout5) = r + s
                 r = r1 + cos4*r25 + cos2*r34
                 s = sin4*s25 - sin2*s34
                 zout(1,j,nout3) = r - s
                 zout(1,j,nout4) = r + s
                 r25 = r2 - r5
                 r34 = r3 - r4
                 s25 = s2 + s5
                 s34 = s3 + s4
                 zout(2,j,nout1) = s1 + s25 + s34
                 r = s1 + cos2*s25 + cos4*s34
                 s = sin2*r25 + sin4*r34
                 zout(2,j,nout2) = r + s
                 zout(2,j,nout5) = r - s
                 r = s1 + cos4*s25 + cos2*s34
                 s = sin4*r25 - sin2*r34
                 zout(2,j,nout3) = r + s
                 zout(2,j,nout4) = r - s
              end do
           end do
        endif
     end do loop_5000
! end of if now.eq.5

  else if (now.eq.6) then
!     .5d0*sqrt(3.d0)
     bb=real(isign,kind=8)*0.8660254037844387d0
     ia=1
     nin1=ia-after
     nout1=ia-atn
     do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nin6=nin5+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        nout6=nout5+after
        do j=1,nfft
           r2=zin(1,j,nin3)
           s2=zin(2,j,nin3)
           r3=zin(1,j,nin5)
           s3=zin(2,j,nin5)
           r=r2 + r3
           s=s2 + s3
           r1=zin(1,j,nin1)
           s1=zin(2,j,nin1)
           ur1 = r + r1
           ui1 = s + s1
           r1=r1 - .5d0*r
           s1=s1 - .5d0*s
           r=r2-r3
           s=s2-s3
           ur2 = r1 - s*bb
           ui2 = s1 + r*bb
           ur3 = r1 + s*bb
           ui3 = s1 - r*bb

           r2=zin(1,j,nin6)
           s2=zin(2,j,nin6)
           r3=zin(1,j,nin2)
           s3=zin(2,j,nin2)
           r=r2 + r3
           s=s2 + s3
           r1=zin(1,j,nin4)
           s1=zin(2,j,nin4)
           vr1 = r + r1
           vi1 = s + s1
           r1=r1 - .5d0*r
           s1=s1 - .5d0*s
           r=r2-r3
           s=s2-s3
           vr2 = r1 - s*bb
           vi2 = s1 + r*bb
           vr3 = r1 + s*bb
           vi3 = s1 - r*bb

           zout(1,j,nout1)=ur1+vr1
           zout(2,j,nout1)=ui1+vi1
           zout(1,j,nout5)=ur2+vr2
           zout(2,j,nout5)=ui2+vi2
           zout(1,j,nout3)=ur3+vr3
           zout(2,j,nout3)=ui3+vi3
           zout(1,j,nout4)=ur1-vr1
           zout(2,j,nout4)=ui1-vi1
           zout(1,j,nout2)=ur2-vr2
           zout(2,j,nout2)=ui2-vi2
           zout(1,j,nout6)=ur3-vr3
           zout(2,j,nout6)=ui3-vi3
        end do
     end do

  else
     stop 'error fftstp'
  endif !end of now

end subroutine fftstp

subroutine fourier_dim(n,n_next)
  implicit none
  !Arguments
  integer, intent(in) :: n
  integer, intent(out) :: n_next

  !Local variables
  integer, parameter :: ndata = 180
  !Multiple of 2,3,5
  integer, dimension(ndata), parameter :: idata = (/   &
          3,     4,     5,     6,     8,     9,    12,    15,    16,    18, &
         20,    24,    25,    27,    30,    32,    36,    40,    45,    48, &
         54,    60,    64,    72,    75,    80,    81,    90,    96,   100, &
        108,   120,   125,   128,   135,   144,   150,   160,   162,   180, &
        192,   200,   216,   225,   240,   243,   256,   270,   288,   300, &
        320,   324,   360,   375,   384,   400,   405,   432,   450,   480, &
        486,   500,   512,   540,   576,   600,   625,   640,   648,   675, &
        720,   729,   750,   768,   800,   810,   864,   900,   960,   972, &
       1000,  1024,  1080,  1125,  1152,  1200,  1215,  1280,  1296,  1350, &
       1440,  1458,  1500,  1536,  1600,  1620,  1728,  1800,  1875,  1920, &
       1944,  2000,  2025,  2048,  2160,  2250,  2304,  2400,  2430,  2500, &
       2560,  2592,  2700,  2880,  3000,  3072,  3125,  3200,  3240,  3375, &
       3456,  3600,  3750,  3840,  3888,  4000,  4050,  4096,  4320,  4500, &
       4608,  4800,  5000,  5120,  5184,  5400,  5625,  5760,  6000,  6144, &
       6400,  6480,  6750,  6912,  7200,  7500,  7680,  8000,  8192,  8640, &
       9000,  9216,  9375,  9600, 10000, 10240, 10368, 10800, 11250, 11520, &
      12000, 12288, 12500, 12800, 13824, 14400, 15000, 15360, 15625, 16000, &
      16384, 17280, 18000, 18432, 18750, 19200, 20000, 20480, 23040, 24000   /)
  integer :: i

  loop_data: do i=1,ndata
     if (n <= idata(i)) then
        n_next = idata(i)
        return
     end if
  end do loop_data
  write(unit=*,fmt=*) "fourier_dim: ",n," is bigger than ",idata(ndata)
  stop
end subroutine fourier_dim

subroutine ctrig(n,trig,after,before,now,isign,ic)
  implicit none
! Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
! Arguments
  integer :: n,isign,ic,i,itt,j,nh
  integer :: after(7),before(7),now(7)
  real(kind=8) :: trig(2,nfft_max), angle,trigc,trigs,twopi
! Local variables
  integer, parameter :: ndata = 180
  integer, dimension(7,ndata) :: idata
! The factor 6 is only allowed in the first place!
        data ((idata(i,j),i=1,7),j=1,ndata) /                     &
            3,   3, 1, 1, 1, 1, 1,       4,   4, 1, 1, 1, 1, 1,   &
            5,   5, 1, 1, 1, 1, 1,       6,   6, 1, 1, 1, 1, 1,   &
            8,   8, 1, 1, 1, 1, 1,       9,   3, 3, 1, 1, 1, 1,   &
           12,   4, 3, 1, 1, 1, 1,      15,   5, 3, 1, 1, 1, 1,   &
           16,   4, 4, 1, 1, 1, 1,      18,   6, 3, 1, 1, 1, 1,   &
           20,   5, 4, 1, 1, 1, 1,      24,   8, 3, 1, 1, 1, 1,   &
           25,   5, 5, 1, 1, 1, 1,      27,   3, 3, 3, 1, 1, 1,   &
           30,   6, 5, 1, 1, 1, 1,      32,   8, 4, 1, 1, 1, 1,   &
           36,   4, 3, 3, 1, 1, 1,      40,   8, 5, 1, 1, 1, 1,   &
           45,   5, 3, 3, 1, 1, 1,      48,   4, 4, 3, 1, 1, 1,   &
           54,   6, 3, 3, 1, 1, 1,      60,   5, 4, 3, 1, 1, 1,   &
           64,   8, 8, 1, 1, 1, 1,      72,   8, 3, 3, 1, 1, 1,   &
           75,   5, 5, 3, 1, 1, 1,      80,   5, 4, 4, 1, 1, 1,   &
           81,   3, 3, 3, 3, 1, 1,      90,   6, 5, 3, 1, 1, 1,   &
           96,   8, 4, 3, 1, 1, 1,     100,   5, 5, 4, 1, 1, 1,   &
          108,   4, 3, 3, 3, 1, 1,     120,   8, 5, 3, 1, 1, 1,   &
          125,   5, 5, 5, 1, 1, 1,     128,   8, 4, 4, 1, 1, 1,   &
          135,   5, 3, 3, 3, 1, 1,     144,   6, 8, 3, 1, 1, 1,   &
          150,   6, 5, 5, 1, 1, 1,     160,   8, 5, 4, 1, 1, 1,   &
          162,   6, 3, 3, 3, 1, 1,     180,   5, 4, 3, 3, 1, 1,   &
          192,   6, 8, 4, 1, 1, 1,     200,   8, 5, 5, 1, 1, 1,   &
          216,   8, 3, 3, 3, 1, 1,     225,   5, 5, 3, 3, 1, 1,   &
          240,   6, 8, 5, 1, 1, 1,     243,   3, 3, 3, 3, 3, 1,   &
          256,   8, 8, 4, 1, 1, 1,     270,   6, 5, 3, 3, 1, 1,   &
          288,   8, 4, 3, 3, 1, 1,     300,   5, 5, 4, 3, 1, 1,   &
          320,   5, 4, 4, 4, 1, 1,     324,   4, 3, 3, 3, 3, 1,   &
          360,   8, 5, 3, 3, 1, 1,     375,   5, 5, 5, 3, 1, 1,   &
          384,   8, 4, 4, 3, 1, 1,     400,   5, 5, 4, 4, 1, 1,   &
          405,   5, 3, 3, 3, 3, 1,     432,   4, 4, 3, 3, 3, 1,   &
          450,   6, 5, 5, 3, 1, 1,     480,   8, 5, 4, 3, 1, 1,   &
          486,   6, 3, 3, 3, 3, 1,     500,   5, 5, 5, 4, 1, 1,   &
          512,   8, 8, 8, 1, 1, 1,     540,   5, 4, 3, 3, 3, 1,   &
          576,   4, 4, 4, 3, 3, 1,     600,   8, 5, 5, 3, 1, 1,   &
          625,   5, 5, 5, 5, 1, 1,     640,   8, 5, 4, 4, 1, 1,   &
          648,   8, 3, 3, 3, 3, 1,     675,   5, 5, 3, 3, 3, 1,   &
          720,   5, 4, 4, 3, 3, 1,     729,   3, 3, 3, 3, 3, 3,   &
          750,   6, 5, 5, 5, 1, 1,     768,   4, 4, 4, 4, 3, 1,   &
          800,   8, 5, 5, 4, 1, 1,     810,   6, 5, 3, 3, 3, 1,   &
          864,   8, 4, 3, 3, 3, 1,     900,   5, 5, 4, 3, 3, 1,   &
          960,   5, 4, 4, 4, 3, 1,     972,   4, 3, 3, 3, 3, 3,   &
         1000,   8, 5, 5, 5, 1, 1,    1024,   4, 4, 4, 4, 4, 1,   &
         1080,   6, 5, 4, 3, 3, 1,    1125,   5, 5, 5, 3, 3, 1,   &
         1152,   6, 4, 4, 4, 3, 1,    1200,   6, 8, 5, 5, 1, 1,   &
         1215,   5, 3, 3, 3, 3, 3,    1280,   8, 8, 5, 4, 1, 1,   &
         1296,   6, 8, 3, 3, 3, 1,    1350,   6, 5, 5, 3, 3, 1,   &
         1440,   6, 5, 4, 4, 3, 1,    1458,   6, 3, 3, 3, 3, 3,   &
         1500,   5, 5, 5, 4, 3, 1,    1536,   6, 8, 8, 4, 1, 1,   &
         1600,   8, 8, 5, 5, 1, 1,    1620,   5, 4, 3, 3, 3, 3,   &
         1728,   6, 8, 4, 3, 3, 1,    1800,   6, 5, 5, 4, 3, 1,   &
         1875,   5, 5, 5, 5, 3, 1,    1920,   6, 5, 4, 4, 4, 1,   &
         1944,   6, 4, 3, 3, 3, 3,    2000,   5, 5, 5, 4, 4, 1,   &
         2025,   5, 5, 3, 3, 3, 3,    2048,   8, 4, 4, 4, 4, 1,   &
         2160,   6, 8, 5, 3, 3, 1,    2250,   6, 5, 5, 5, 3, 1,   &
         2304,   6, 8, 4, 4, 3, 1,    2400,   6, 5, 5, 4, 4, 1,   &
         2430,   6, 5, 3, 3, 3, 3,    2500,   5, 5, 5, 5, 4, 1,   &
         2560,   8, 5, 4, 4, 4, 1,    2592,   6, 4, 4, 3, 3, 3,   &
         2700,   5, 5, 4, 3, 3, 3,    2880,   6, 8, 5, 4, 3, 1,   &
         3000,   6, 5, 5, 5, 4, 1,    3072,   6, 8, 4, 4, 4, 1,   &
         3125,   5, 5, 5, 5, 5, 1,    3200,   8, 5, 5, 4, 4, 1,   &
         3240,   6, 5, 4, 3, 3, 3,    3375,   5, 5, 5, 3, 3, 3,   &
         3456,   6, 4, 4, 4, 3, 3,    3600,   6, 8, 5, 5, 3, 1,   &
         3750,   6, 5, 5, 5, 5, 1,    3840,   6, 8, 5, 4, 4, 1,   &
         3888,   6, 8, 3, 3, 3, 3,    4000,   8, 5, 5, 5, 4, 1,   &
         4050,   6, 5, 5, 3, 3, 3,    4096,   8, 8, 4, 4, 4, 1,   &
         4320,   6, 5, 4, 4, 3, 3,    4500,   5, 5, 5, 4, 3, 3,   &
         4608,   6, 8, 8, 4, 3, 1,    4800,   6, 8, 5, 5, 4, 1,   &
         5000,   8, 5, 5, 5, 5, 1,    5120,   8, 8, 5, 4, 4, 1,   &
         5184,   6, 8, 4, 3, 3, 3,    5400,   6, 5, 5, 4, 3, 3,   &
         5625,   5, 5, 5, 5, 3, 3,    5760,   6, 8, 8, 5, 3, 1,   &
         6000,   6, 8, 5, 5, 5, 1,    6144,   6, 8, 8, 4, 4, 1,   &
         6400,   8, 8, 5, 5, 4, 1,    6480,   6, 8, 5, 3, 3, 3,   &
         6750,   6, 5, 5, 5, 3, 3,    6912,   6, 8, 4, 4, 3, 3,   &
         7200,   6, 5, 5, 4, 4, 3,    7500,   5, 5, 5, 5, 4, 3,   &
         7680,   6, 8, 8, 5, 4, 1,    8000,   8, 8, 5, 5, 5, 1,   &
         8192,   8, 8, 8, 4, 4, 1,    8640,   8, 8, 5, 3, 3, 3,   &
         9000,   8, 5, 5, 5, 3, 3,    9216,   6, 8, 8, 8, 3, 1,   &
         9375,   5, 5, 5, 5, 5, 3,    9600,   8, 5, 5, 4, 4, 3,   &
        10000,   5, 5, 5, 5, 4, 4,   10240,   8, 8, 8, 5, 4, 1,   &
        10368,   6, 8, 8, 3, 3, 3,   10800,   6, 8, 5, 5, 3, 3,   &
        11250,   6, 5, 5, 5, 5, 3,   11520,   8, 8, 5, 4, 3, 3,   &
        12000,   8, 5, 5, 5, 4, 3,   12288,   8, 8, 8, 8, 3, 1,   &
        12500,   5, 5, 5, 5, 5, 4,   12800,   8, 8, 8, 5, 5, 1,   &
        13824,   8, 8, 8, 3, 3, 3,   14400,   8, 8, 5, 5, 3, 3,   &
        15000,   8, 5, 5, 5, 5, 3,   15360,   6, 8, 8, 8, 5, 1,   &
        15625,   5, 5, 5, 5, 5, 5,   16000,   8, 5, 5, 5, 4, 4,   &
        16384,   8, 8, 8, 8, 4, 1,   17280,   6, 8, 8, 5, 3, 3,   &
        18000,   6, 8, 5, 5, 5, 3,   18432,   8, 8, 8, 4, 3, 3,   &
        18750,   6, 5, 5, 5, 5, 5,   19200,   8, 8, 5, 5, 4, 3,   &
        20000,   8, 5, 5, 5, 5, 4,   20480,   8, 8, 8, 8, 5, 1,   &
        23040,   8, 8, 8, 5, 3, 3,   24000,   8, 8, 5, 5, 5, 3    /

  do i=1,ndata
     if (n.eq.idata(1,i)) then
        ic=0
        do j=1,6
           itt=idata(1+j,i)
           if (itt.gt.1) then
              ic=ic+1
              now(j)=idata(1+j,i)
           else
              goto 1000
           endif
        end do
        goto 1000
     endif
  end do
  print *,'VALUE OF',n,'NOT ALLOWED FOR FFT, ALLOWED VALUES ARE:'
  write(6,"(15(i5))") (idata(1,i),i=1,ndata)
  stop
1000 continue

  after(1)=1
  before(ic)=1
  do i=2,ic
     after(i)=after(i-1)*now(i-1)
     before(ic-i+1)=before(ic-i+2)*now(ic-i+2)
  end do

!  write(6,"(6(i3))") (after(i),i=1,ic)
!  write(6,"(6(i3))") (now(i),i=1,ic)
!  write(6,"(6(i3))") (before(i),i=1,ic)

  twopi=6.283185307179586d0
  angle=real(isign,kind=8)*twopi/real(n,kind=8)
  if (mod(n,2).eq.0) then
     nh=n/2
     trig(1,1)=1.d0
     trig(2,1)=0.d0
     trig(1,nh+1)=-1.d0
     trig(2,nh+1)=0.d0
     do i=1,nh-1
        trigc=cos(real(i,kind=8)*angle)
        trigs=sin(real(i,kind=8)*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
     end do
  else
     nh=(n-1)/2
     trig(1,1)=1.d0
     trig(2,1)=0.d0
     do i=1,nh
        trigc=cos(real(i,kind=8)*angle)
        trigs=sin(real(i,kind=8)*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
      end do
  endif

end subroutine ctrig

subroutine test_functions(n01,n02,n03,acell,a_gauss,hx,hy,hz,&
     density,potential)
  implicit none
  integer, intent(in) :: n01,n02,n03
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz
  real(kind=8), dimension(n01,n02,n03), intent(out) :: density,potential

  !local variables
  integer :: i1,i2,i3,nu,ifx,ify,ifz
  real(kind=8) :: x,x1,x2,x3,y,length,denval,pi,a2,derf,hgrid,factor,r,r2
  real(kind=8) :: fx,fx2,fy,fy2,fz,fz2,a,ax,ay,az,bx,by,bz,tt,potion_fac



     !grid for the free BC case
     !hgrid=max(hx,hy,hz)

     pi = 4.d0*atan(1.d0)
     a2 = a_gauss**2

     !Normalization
     factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
     !gaussian function
     do i3=1,n03
        x3 = hz*real(i3-n03/2,kind=8)
        do i2=1,n02
           x2 = hy*real(i2-n02/2,kind=8)
           do i1=1,n01
              x1 = hx*real(i1-n01/2,kind=8)
              r2 = x1*x1+x2*x2+x3*x3
              density(i1,i2,i3) = factor*exp(-r2/a2)
              r = sqrt(r2)
              !Potential from a gaussian
              if (r == 0.d0) then
                 potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
              else
                 potential(i1,i2,i3) = derf(r/a_gauss)/r
              end if
           end do
        end do
     end do
end subroutine test_functions

subroutine gequad(nterms,p,w,urange,drange,acc)
!
  implicit none
  real(kind=8) :: urange,drange,acc
  integer :: nterms
  real(kind=8) :: p(*),w(*)
!
!       range [10^(-9),1] and accuracy ~10^(-8);
!
  p(1)=4.96142640560223544d19
  p(2)=1.37454269147978052d19
  p(3)=7.58610013441204679d18
  p(4)=4.42040691347806996d18
  p(5)=2.61986077948367892d18
  p(6)=1.56320138155496681d18
  p(7)=9.35645215863028402d17
  p(8)=5.60962910452691703d17
  p(9)=3.3666225119686761d17
  p(10)=2.0218253197947866d17
  p(11)=1.21477756091902017d17
  p(12)=7.3012982513608503d16
  p(13)=4.38951893556421099d16
  p(14)=2.63949482512262325d16
  p(15)=1.58742054072786174d16
  p(16)=9.54806587737665531d15
  p(17)=5.74353712364571709d15
  p(18)=3.455214877389445d15
  p(19)=2.07871658520326804d15
  p(20)=1.25064667315629928d15
  p(21)=7.52469429541933745d14
  p(22)=4.5274603337253175d14
  p(23)=2.72414006900059548d14
  p(24)=1.63912168349216752d14
  p(25)=9.86275802590865738d13
  p(26)=5.93457701624974985d13
  p(27)=3.5709554322296296d13
  p(28)=2.14872890367310454d13
  p(29)=1.29294719957726902d13
  p(30)=7.78003375426361016d12
  p(31)=4.68148199759876704d12
  p(32)=2.8169955024829868d12
  p(33)=1.69507790481958464d12
  p(34)=1.01998486064607581d12
  p(35)=6.13759486539856459d11
  p(36)=3.69320183828682544d11
  p(37)=2.22232783898905102d11
  p(38)=1.33725247623668682d11
  p(39)=8.0467192739036288d10
  p(40)=4.84199582415144143d10
  p(41)=2.91360091170559564d10
  p(42)=1.75321747475309216d10
  p(43)=1.0549735552210995d10
  p(44)=6.34815321079006586d9
  p(45)=3.81991113733594231d9
  p(46)=2.29857747533101109d9
  p(47)=1.38313653595483694d9
  p(48)=8.32282908580025358d8
  p(49)=5.00814519374587467d8
  p(50)=3.01358090773319025d8
  p(51)=1.81337994217503535d8
  p(52)=1.09117589961086823d8
  p(53)=6.56599771718640323d7
  p(54)=3.95099693638497164d7
  p(55)=2.37745694710665991d7
  p(56)=1.43060135285912813d7
  p(57)=8.60844290313506695d6
  p(58)=5.18000974075383424d6
  p(59)=3.116998193057466d6
  p(60)=1.87560993870024029d6
  p(61)=1.12862197183979562d6
  p(62)=679132.441326077231d0
  p(63)=408658.421279877969d0
  p(64)=245904.473450669789d0
  p(65)=147969.568088321005d0
  p(66)=89038.612357311147d0
  p(67)=53577.7362552358895d0
  p(68)=32239.6513926914668d0
  p(69)=19399.7580852362791d0
  p(70)=11673.5323603058634d0
  p(71)=7024.38438577707758d0
  p(72)=4226.82479307685999d0
  p(73)=2543.43254175354295d0
  p(74)=1530.47486269122675d0
  p(75)=920.941785160749482d0
  p(76)=554.163803906291646d0
  p(77)=333.46029740785694d0
  p(78)=200.6550575335041d0
  p(79)=120.741366914147284d0
  p(80)=72.6544243200329916d0
  p(81)=43.7187810415471025d0
  p(82)=26.3071631447061043d0
  p(83)=15.8299486353816329d0
  p(84)=9.52493152341244004d0
  p(85)=5.72200417067776041d0
  p(86)=3.36242234070940928d0
  p(87)=1.75371394604499472d0
  p(88)=0.64705932650658966d0
  p(89)=0.072765905943708247d0
!
  w(1)=47.67445484528304247d10
  w(2)=11.37485774750442175d9
  w(3)=78.64340976880190239d8
  w(4)=46.27335788759590498d8
  w(5)=24.7380464827152951d8
  w(6)=13.62904116438987719d8
  w(7)=92.79560029045882433d8
  w(8)=52.15931216254660251d8
  w(9)=31.67018011061666244d8
  w(10)=1.29291036801493046d8
  w(11)=1.00139319988015862d8
  w(12)=7.75892350510188341d7
  w(13)=6.01333567950731271d7
  w(14)=4.66141178654796875d7
  w(15)=3.61398903394911448d7
  w(16)=2.80225846672956389d7
  w(17)=2.1730509180930247d7
  w(18)=1.68524482625876965d7
  w(19)=1.30701489345870338d7
  w(20)=1.01371784832269282d7
  w(21)=7.86264116300379329d6
  w(22)=6.09861667912273717d6
  w(23)=4.73045784039455683d6
  w(24)=3.66928949951594161d6
  w(25)=2.8462050836230259d6
  w(26)=2.20777394798527011d6
  w(27)=1.71256191589205524d6
  w(28)=1.32843556197737076d6
  w(29)=1.0304731275955989d6
  w(30)=799345.206572271448d0
  w(31)=620059.354143595343d0
  w(32)=480986.704107449333d0
  w(33)=373107.167700228515d0
  w(34)=289424.08337412132d0
  w(35)=224510.248231581788d0
  w(36)=174155.825690028966d0
  w(37)=135095.256919654065d0
  w(38)=104795.442776800312d0
  w(39)=81291.4458222430418d0
  w(40)=63059.0493649328682d0
  w(41)=48915.9040455329689d0
  w(42)=37944.8484018048756d0
  w(43)=29434.4290473253969d0
  w(44)=22832.7622054490044d0
  w(45)=17711.743950151233d0
  w(46)=13739.287867104177d0
  w(47)=10657.7895710752585d0
  w(48)=8267.42141053961834d0
  w(49)=6413.17397520136448d0
  w(50)=4974.80402838654277d0
  w(51)=3859.03698188553047d0
  w(52)=2993.51824493299154d0
  w(53)=2322.1211966811754d0
  w(54)=1801.30750964719641d0
  w(55)=1397.30379659817038d0
  w(56)=1083.91149143250697d0
  w(57)=840.807939169209188d0
  w(58)=652.228524366749422d0
  w(59)=505.944376983506128d0
  w(60)=392.469362317941064d0
  w(61)=304.444930257324312d0
  w(62)=236.162932842453601d0
  w(63)=183.195466078603525d0
  w(64)=142.107732186551471d0
  w(65)=110.23530215723992d0
  w(66)=85.5113346705382257d0
  w(67)=66.3325469806696621d0
  w(68)=51.4552463353841373d0
  w(69)=39.9146798429449273d0
  w(70)=30.9624728409162095d0
  w(71)=24.018098812215013d0
  w(72)=18.6312338024296588d0
  w(73)=14.4525541233150501d0
  w(74)=11.2110836519105938d0
  w(75)=8.69662175848497178d0
  w(76)=6.74611236165731961d0
  w(77)=5.23307018057529994d0
  w(78)=4.05937850501539556d0
  w(79)=3.14892659076635714d0
  w(80)=2.44267408211071604d0
  w(81)=1.89482240522855261d0
  w(82)=1.46984505907050079d0
  w(83)=1.14019261330527007d0
  w(84)=0.884791217422925293d0
  w(85)=0.692686387080616483d0
  w(86)=0.585244576897023282d0
  w(87)=0.576182522545327589d0
  w(88)=0.596688817388997178d0
  w(89)=0.607879901151108771d0
!
  urange = 1.d0
  drange=1d-08
  acc   =1d-08
!
end subroutine gequad

subroutine scramble_unpack(i1,j2,lot,nfft,n1,n3,md2,nd3,zw,zmpi2,cosinarr)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nd3
  real(kind=8), dimension(2,lot,n3/2), intent(in) :: zw
  real(kind=8), dimension(2,n3/2), intent(in) :: cosinarr
  real(kind=8), dimension(2,n1,md2/1,nd3), intent(out) :: zmpi2
  !Local variables
  integer :: i3,i,ind1,ind2
  real(kind=8) ::  a,b,c,d,cp,sp,feR,feI,foR,foI,fR,fI

  !Body

  !case i3=1 and i3=n3/2+1
  do i=0,nfft-1
     a=zw(1,i+1,1)
     b=zw(2,i+1,1)
     zmpi2(1,i1+i,j2,1)=a+b
     zmpi2(2,i1+i,j2,1)=0.d0
     zmpi2(1,i1+i,j2,n3/2+1)=a-b
     zmpi2(2,i1+i,j2,n3/2+1)=0.d0
  end do
  !case 2<=i3<=n3/2
  do i3=2,n3/2
     ind1=i3
     ind2=n3/2-i3+2
     cp=cosinarr(1,i3)
     sp=cosinarr(2,i3)
     do i=0,nfft-1
        a=zw(1,i+1,ind1)
        b=zw(2,i+1,ind1)
        c=zw(1,i+1,ind2)
        d=zw(2,i+1,ind2)
        feR=.5d0*(a+c)
        feI=.5d0*(b-d)
        foR=.5d0*(a-c)
        foI=.5d0*(b+d)
        fR=feR+cp*foI-sp*foR
        fI=feI-cp*foR-sp*foI
        zmpi2(1,i1+i,j2,ind1)=fR
        zmpi2(2,i1+i,j2,ind1)=fI
     end do
  end do

end subroutine scramble_unpack

subroutine switch(nfft,n2,lot,n1,lzt,zt,zw)
  implicit none
  integer :: nfft, n2, lot, n1, lzt, i, j
  real(kind=8) :: zw(2,lot,n2),zt(2,lzt,n1)

  do j=1,nfft
     do i=1,n2
        zw(1,j,i)=zt(1,i,j)
        zw(2,j,i)=zt(2,i,j)
     end do
  end do

end subroutine switch

subroutine realcopy(lot,nfft,n2,nk1,nk2,zin,zout)
  implicit none
  integer, intent(in) :: nfft,lot,n2,nk1,nk2
  real(kind=8), dimension(2,lot,n2), intent(in) :: zin
  real(kind=8), dimension(nk1,nk2), intent(out) :: zout
  !local variables
  integer :: i,j

  do i=1,nk2
     do j=1,nfft
        zout(j,i)=zin(1,j,i)
     end do
  end do

end subroutine realcopy

subroutine mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,zmpi1,zw)
  implicit none
  integer      :: nfft,lot,n1,nd2,nd3, I1, J2, Jp2, mfft, j3, Jp2st,J2st
  real(kind=8) :: zmpi1(2,n1,nd2/1,nd3/1,1),zw(2,lot,n1)

  mfft=0
  do Jp2=Jp2st,1
     do J2=J2st,nd2/1
        mfft=mfft+1
        if (mfft.gt.nfft) then
           Jp2st=Jp2
           J2st=J2
           return
        endif
        do I1=1,n1
           zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
           zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
        end do
     end do
     J2st=1
  end do

end subroutine mpiswitch

subroutine inserthalf(n1,n3,lot,nfft,i1,zf,zw)
  implicit none
  integer, intent(in) :: n1,n3,lot,nfft,i1
  real(kind=8), dimension(n1/2+1,n3/2+1), intent(in) :: zf
  real(kind=8), dimension(2,lot,n3/2), intent(out) :: zw
  !local variables
  integer :: l1,l3,i01,i03r,i03i,i3

  i3=0
  do l3=1,n3,2
     i3=i3+1
     i03r=abs(l3-n3/2-1)+1
     i03i=abs(l3-n3/2)+1
     do l1=1,nfft
        i01=abs(l1-1+i1-n1/2-1)+1
        zw(1,l1,i3)=zf(i01,i03r)
        zw(2,l1,i3)=zf(i01,i03i)
     end do
  end do

end subroutine inserthalf

subroutine kernelfft(n1,n2,n3,nd1,nd2,nd3,nk1,nk2,nk3,zf,zr)
    implicit none
!value of ncache for the FFT routines
integer, parameter :: ncache_optimal=8*1024
!flag that states if we must perform the timings or not
!definitions of the timing variambles just in case
integer, parameter :: timing_flag=0
integer :: count_time1,count_time2,count_rate,count_max,number_time,index_time
real(kind=8) :: time_b0,time_b1,serial_time,parallel_time

  !Arguments
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,nk1,nk2,nk3
  real(kind=8), dimension(n1/2+1,n3/2+1,nd2/1), intent(in) :: zf
  real(kind=8), dimension(nk1,nk2,nk3/1), intent(inout) :: zr
  !Local variables
  !Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2st,J2st
  integer :: j2,j3,i1,i3,i,j,inzee,ierr,i_all,i_stat,k
  real(kind=8) :: twopion
  !work arrays for transpositions
  real(kind=8), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(kind=8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(kind=8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(kind=8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(kind=8), dimension(:,:), allocatable :: trig1,trig2,trig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, &
       after2,now2,before2,after3,now3,before3

  !Body

  !check input
  if (nd1.lt.n1) stop 'ERROR:nd1'
  if (nd2.lt.n2) stop 'ERROR:nd2'
  if (nd3.lt.n3/2+1) stop 'ERROR:nd3'
  if (mod(nd3,1).ne.0) stop 'ERROR:nd3'
  if (mod(nd2,1).ne.0) stop 'ERROR:nd2'

  !defining work arrays dimensions
  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4
  lzt=n2
  if (mod(n2,2).eq.0) lzt=lzt+1
  if (mod(n2,4).eq.0) lzt=lzt+1

  !Allocations
  allocate(trig1(2,nfft_max),stat=i_stat)
  allocate(after1(7),stat=i_stat)
  allocate(now1(7),stat=i_stat)
  allocate(before1(7),stat=i_stat)
  allocate(trig2(2,nfft_max),stat=i_stat)
  allocate(after2(7),stat=i_stat)
  allocate(now2(7),stat=i_stat)
  allocate(before2(7),stat=i_stat)
  allocate(trig3(2,nfft_max),stat=i_stat)
  allocate(after3(7),stat=i_stat)
  allocate(now3(7),stat=i_stat)
  allocate(before3(7),stat=i_stat)
  allocate(zw(2,ncache/4,2),stat=i_stat)
  allocate(zt(2,lzt,n1),stat=i_stat)
  allocate(zmpi2(2,n1,nd2/1,nd3),stat=i_stat)
  allocate(cosinarr(2,n3/2),stat=i_stat)

  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3/2,trig3,after3,before3,now3,1,ic3)
  call ctrig(n1,trig1,after1,before1,now1,1,ic1)
  call ctrig(n2,trig2,after2,before2,now2,1,ic2)

  !Calculating array of phases for HalFFT decoding
  twopion=8.d0*datan(1.d0)/real(n3,kind=8)
  do i3=1,n3/2
     cosinarr(1,i3)= dcos(twopion*real(i3-1,kind=8))
     cosinarr(2,i3)=-dsin(twopion*real(i3-1,kind=8))
  end do

  lot=ncache/(2*n3)
  if (lot.lt.1) stop 'kernelfft:enlarge ncache for z'

  do j2=1,nd2
     !this condition ensures that we manage only the interesting part for the FFT
     if (j2.le.n2) then
        do i1=1,n1
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           !input: I1,I3,J2,(Jp2)

           call inserthalf(n1,n3,lot,nfft,i1,zf(1,1,j2),zw(1,1,1))

           !performing FFT
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo
           !output: I1,i3,J2,(Jp2)

           !unpacking FFT in order to restore correct result,
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_unpack(i1,j2,lot,nfft,n1,n3,nd2,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        end do
     endif
  end do


  do j3=1,nd3
     !this condition ensures that we manage only the interesting part for the FFT
     if (j3.le.n3/2+1) then
        Jp2st=1
        J2st=1

        !transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) stop 'kernelfft:enlarge ncache for x'

        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !reverse ordering
           !input: I1,J2,j3,Jp2,(jp3)
           call mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,zmpi2,zw(1,1,1))
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo
           !storing the last step into zt
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), &
                trig1,after1(i),now1(i),before1(i),1)
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis, and taking only the first half
        lot=ncache/(4*n2)
        if (lot.lt.1) stop 'kernelfft:enlarge ncache for y'

        do j=1,nk1,lot
           ma=j
           mb=min(j+(lot-1),nk1)
           nfft=mb-ma+1

           !reverse ordering
           !input: I2,i1,j3,(jp3)
           call switch(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)

           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo

           call realcopy(lot,nfft,n2,nk1,nk2,zw(1,1,inzee),zr(j,1,j3))

        end do
        !output: i1,i2,j3,(jp3)
     endif
  end do

end subroutine kernelfft

subroutine F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3)
 implicit none
 integer, intent(in) :: n01,n02,n03
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space, inverted for convenience
 m1=n01
 m2=n03
 m3=n02
 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=2*m2
 l3=m3 !beware of the half dimension
 do
    call fourier_dim(l1,n1)
    if (modulo(n1,2) == 0) then
       exit
    end if
    l1=l1+1
 end do
 do
    call fourier_dim(l2,n2)
    if (modulo(n2,2) == 0) then
       exit
    end if
    l2=l2+1
 end do
 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2/2
 md3=n3/2
151 if (1*(md2/1).lt.n2/2) then
    md2=md2+1
    goto 151
 endif

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with 1
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1

250 if (modulo(nd3,1) .ne. 0) then
    nd3=nd3+1
    goto 250
 endif

end subroutine F_FFT_dimensions

subroutine back_trans_14(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

    integer, parameter :: m=16
    real(kind=8), dimension(-m:m) :: ch = (/ &
         0.d0,0.d0,0.d0,0.0000275373458862304687D0,0.D0,-0.000423073768615722656D0,0.D0,&
         0.00310254096984863281D0,0.D0,-0.0146262645721435547D0,0.D0,&
         0.0511919260025024414D0,0.D0,-0.153575778007507324D0,0.D0,0.614303112030029297D0,&
         1.D0,0.614303112030029297D0,0.D0,-0.153575778007507324D0,0.D0,&
         0.0511919260025024414D0,0.D0,-0.0146262645721435547D0,0.D0,&
         0.00310254096984863281D0,0.D0,-0.000423073768615722656D0,0.D0,&
         0.0000275373458862304687D0,0.d0,0.d0,0.d0&
         /)
    real(kind=8), dimension(-m:m) ::  cg,cht,cgt

    !******** coefficients for wavelet transform *********************
    do i=-m,m
       cht(i)=0.d0
       cg(i)=0.d0
       cgt(i)=0.d0
    enddo

    ! the normalization is chosen such that a constant function remains the same constant
    ! on each level of the transform

    cht( 0)=1.D0

    ! g coefficients from h coefficients
    do i=-m,m-1
       cg(i+1)=cht(-i)*(-1.d0)**(i+1)
       cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
    enddo

  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0

     do j=-m/2,m/2-1

        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)

     end do
  end do

end subroutine back_trans_14


subroutine scaling_function(itype,nd,nrange,a,x)
    implicit none

    !Type of interpolating functions
    integer, intent(in) :: itype
    !Number of points: must be 2**nex
    integer, intent(in) :: nd
    integer, intent(out) :: nrange
    real(kind=8), dimension(0:nd), intent(out) :: a,x

    !Local variables
    real(kind=8), dimension(0:nd) :: y
    integer :: i,nt,ni

    ni=2*itype
    nrange = ni
    x(0:nd)=0
    y(0:nd)=0
    nt=ni
    x(nt/2-1)=1.d0

    loop1: do
        nt=2*nt
        call back_trans_14(nd,nt,x,y)
!        call dcopy(nt,y,1,x,1)
        x(1:nt)=y(1:nt)
        if (nt.eq.nd) then
            exit loop1
        end if
    end do loop1

    do i=0,nd
        a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
    end do

end subroutine scaling_function


subroutine scf_recursion(itype,n_iter,n_range,kernel_scf,kern_1_scf)
    implicit none
    !Arguments
    integer, intent(in) :: itype,n_iter,n_range
    real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
    real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
    !Local variables
    real(kind=8) :: kern,kern_tot
    integer :: i_iter,i,j,ind

    integer, parameter :: m=16
    real(kind=8), dimension(-m:m) :: ch = (/ &
         0.d0,0.d0,0.d0,0.0000275373458862304687D0,0.D0,-0.000423073768615722656D0,0.D0,&
         0.00310254096984863281D0,0.D0,-0.0146262645721435547D0,0.D0,&
         0.0511919260025024414D0,0.D0,-0.153575778007507324D0,0.D0,0.614303112030029297D0,&
         1.D0,0.614303112030029297D0,0.D0,-0.153575778007507324D0,0.D0,&
         0.0511919260025024414D0,0.D0,-0.0146262645721435547D0,0.D0,&
         0.00310254096984863281D0,0.D0,-0.000423073768615722656D0,0.D0,&
         0.0000275373458862304687D0,0.d0,0.d0,0.d0&
         /)
    real(kind=8), dimension(-m:m) ::  cg,cht,cgt

    !******** coefficients for wavelet transform *********************
    do i=-m,m
       cht(i)=0.d0
       cg(i)=0.d0
       cgt(i)=0.d0
    enddo

    ! the normalization is chosen such that a constant function remains the same constant
    ! on each level of the transform

    cht( 0)=1.D0

    ! g coefficients from h coefficients
    do i=-m,m-1
       cg(i+1)=cht(-i)*(-1.d0)**(i+1)
       cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
    enddo

    !Start the iteration to go from p0gauss to pgauss
    loop_iter_scf: do i_iter=1,n_iter
        kern_1_scf(:) = kernel_scf(:)
        kernel_scf(:) = 0.d0
        loop_iter_i: do i=0,n_range
            kern_tot = 0.d0
            do j=-m,m
               ind = 2*i-j
               if (abs(ind) > n_range) then
                  kern = 0.d0
               else
                  kern = kern_1_scf(ind)
               end if
               kern_tot = kern_tot + ch(j)*kern
            end do
            if (kern_tot == 0.d0) then
               !zero after (be sure because strictly == 0.d0)
               exit loop_iter_i
            else
               kernel_scf( i) = 0.5d0*kern_tot
               kernel_scf(-i) = kernel_scf(i)
            end if
        end do loop_iter_i
    end do loop_iter_scf

end subroutine scf_recursion

subroutine Free_Kernel(n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,&
     hx,hy,hz,itype_scf,karray)

 implicit none

 !Arguments
 integer, intent(in) :: n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,itype_scf
 real(kind=8), intent(in) :: hx,hy,hz
 real(kind=8), dimension(n1k,n2k,n3k), intent(out) :: karray

 !Local variables
 !Do not touch !!!!
 integer, parameter :: n_gauss = 89
 !Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
 integer, parameter :: n_points = 2**6

 !Better p_gauss for calculation
 !(the support of the exponential should be inside [-n_range/2,n_range/2])
 real(kind=8), parameter :: p0_ref = 1.d0
 real(kind=8), dimension(n_gauss) :: p_gauss,w_gauss

 real(kind=8), dimension(:), allocatable :: kern_1_scf,x_scf ,y_scf
 real(kind=8), dimension(:,:), allocatable :: kernel_scf
 real(kind=8), dimension(:,:,:), allocatable :: kp


 real(kind=8) :: ur_gauss,dr_gauss,acc_gauss,pgauss,kern,a_range,kern_tot
 real(kind=8) :: pi,factor,factor2,urange,dx,absci,p0gauss,weight,p0_cell,u1,u2,u3
 real(kind=8) :: a1,a2,a3,amax,ratio,hgrid,pref1,pref2,pref3,p01,p02,p03,kern1,kern2,kern3
 integer :: n_scf,nker1,nker2,nker3
 integer :: i_gauss,n_range,n_cell,istart,iend,istart1
 integer :: i,j,n_iter,i_iter,ind,i1,i2,i3,i_kern,i_stat,i_all
 integer :: i01,i02,i03,n1h,n2h,n3h,nit1,nit2,nit3

 !grid spacing
 hgrid=max(hx,hy,hz)
 !Number of integration points : 2*itype_scf*n_points
 n_scf=2*itype_scf*n_points
 !Set karray

 !here we must set the dimensions for the fft part, starting from the nfft
 !remember that actually nfft2 is associated to n03 and viceversa

 !dimensions that define the center of symmetry
 n1h=nfft1/2
 n2h=nfft2/2
 n3h=nfft3/2

 !Auxiliary dimensions only for building the FFT part
 nker1=nfft1
 nker2=nfft2
 nker3=nfft3/2+1

 !this will be the array of the kernel in the real space
 allocate(kp(n1h+1,n3h+1,nker2/1),stat=i_stat)

 !defining proper extremes for the calculation of the
 !local part of the kernel

 istart=1
 iend=min((1)*nker2/1,n2h+n03)

 istart1=n2h-n03+2

 !Allocations
 allocate(x_scf(0:n_scf),stat=i_stat)
 allocate(y_scf(0:n_scf),stat=i_stat)

 !Build the scaling function
 call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)

 !Step grid for the integration
 dx = real(n_range,kind=8)/real(n_scf,kind=8)
 !Extend the range (no more calculations because fill in by 0.d0)
 n_cell = max(n01,n02,n03)
 n_range = max(n_cell,n_range)

 !Lengthes of the box (use FFT dimension)
 a1 = hx * real(n01,kind=8)
 a2 = hy * real(n02,kind=8)
 a3 = hz * real(n03,kind=8)


 !Initialization of the gaussian (Beylkin)
 call gequad(n_gauss,p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
 !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3)
 !(biggest length in the cube)
 !We divide the p_gauss by a_range**2 and a_gauss by a_range
 a_range = sqrt(a1*a1+a2*a2+a3*a3)
 factor = 1.d0/a_range
 !factor2 = factor*factor
 factor2 = 1.d0/(a1*a1+a2*a2+a3*a3)
 do i_gauss=1,n_gauss
    p_gauss(i_gauss) = factor2*p_gauss(i_gauss)
 end do
 do i_gauss=1,n_gauss
    w_gauss(i_gauss) = factor*w_gauss(i_gauss)
 end do

 kp(:,:,:)=0.d0

 !Allocations
 allocate(kern_1_scf(-n_range:n_range),stat=i_stat)
 !add the treatment for inhomogeneous hgrids
    allocate(kernel_scf(-n_range:n_range,1),stat=i_stat)

    hgrid=hx

    !To have a correct integration
    p0_cell = p0_ref/(hgrid*hgrid)

    !Use in this order (better for accuracy).
    loop_gauss1: do i_gauss=n_gauss,1,-1
       !Gaussian
       pgauss = p_gauss(i_gauss)
       !We calculate the number of iterations to go from pgauss to p0_ref
       n_iter = nint((log(pgauss) - log(p0_cell))/log(4.d0))
       if (n_iter <= 0)then
          n_iter = 0
          p0gauss = pgauss
       else
          p0gauss = pgauss/4.d0**n_iter
       end if

       !Stupid integration
       !Do the integration with the exponential centered in i_kern
       kernel_scf(:,1) = 0.d0
       do i_kern=0,n_range
          kern = 0.d0
          do i=0,n_scf
             absci = x_scf(i) - real(i_kern,kind=8)
             absci = absci*absci*hgrid**2
             kern = kern + y_scf(i)*dexp(-p0gauss*absci)
          end do
          kernel_scf(i_kern,1) = kern*dx
          kernel_scf(-i_kern,1) = kern*dx
          if (abs(kern) < 1.d-18) then
             !Too small not useful to calculate
             exit
          end if
       end do

       !Start the iteration to go from p0gauss to pgauss
       call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,kern_1_scf)

       !Add to the kernel (only the local part)

       do i3=istart1,iend
          i03 = i3 - n2h -1
          do i2=1,n02
             i02 = i2-1
             do i1=1,n01
                i01 = i1-1
                kp(i1,i2,i3-istart+1) = kp(i1,i2,i3-istart+1) + w_gauss(i_gauss)* &
                     kernel_scf(i01,1)*kernel_scf(i02,1)*kernel_scf(i03,1)
             end do
          end do
       end do

    end do loop_gauss1

 !De-allocations
 i_all=-product(shape(kernel_scf))*kind(kernel_scf)
 deallocate(kernel_scf,stat=i_stat)
 i_all=-product(shape(kern_1_scf))*kind(kern_1_scf)
 deallocate(kern_1_scf,stat=i_stat)
 i_all=-product(shape(x_scf))*kind(x_scf)
 deallocate(x_scf,stat=i_stat)
 i_all=-product(shape(y_scf))*kind(y_scf)
 deallocate(y_scf,stat=i_stat)

!!!!END KERNEL CONSTRUCTION

!!$ if(0 .eq. 0) print *,"Do a 3D PHalFFT for the kernel"

 call kernelfft(nfft1,nfft2,nfft3,nker1,nker2,nker3,n1k,n2k,n3k,&
      kp,karray)

 !De-allocations
 i_all=-product(shape(kp))*kind(kp)
 deallocate(kp,stat=i_stat)

end subroutine Free_Kernel

subroutine PS_dim4allocation(n01,n02,n03,n3d,n3p,n3pi,i3xcsh,i3s)
  implicit none
  integer, intent(in) :: n01,n02,n03
  integer, intent(out) :: n3d,n3p,n3pi,i3xcsh,i3s
  !local variables
  integer, parameter :: nordgr=4
  integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
  integer :: istart,iend,nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr


  !calculate the dimensions wrt the geocode
  call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3)

  !formal start and end of the slice
  istart=0
  iend=min((1)*md2,m2)

     n3d=n03
     n3p=n03
     i3xcsh=0
     i3s=min(istart,m2-1)+1
     n3pi=max(iend-istart,0)

!!$  print *,'P4,0',0,'nxc,ncxl,ncxr,nwbl,nwbr',nxc,nxcl,nxcr,nwbl,nwbr,&
!!$       '0,n3d,n3p,i3xcsh,i3s',0,n3d,n3p,i3xcsh,i3s

end subroutine PS_dim4allocation

subroutine halfill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
  implicit none
  integer :: md1,md3,lot,nfft,n3,i1,i3
  real(kind=8) :: zw(2,lot,n3/2),zf(md1,md3)
! WARNING: Assuming that high frequencies are in the corners
!          and that n3 is multiple of 4
!in principle we can relax this condition

  do i3=1,n3/4
     do i1=1,nfft
        zw(1,i1,i3)=0.d0
        zw(2,i1,i3)=0.d0
     end do
  end do
  do i3=n3/4+1,n3/2
     do i1=1,nfft
        zw(1,i1,i3)=zf(i1,2*i3-1-n3/2)
        zw(2,i1,i3)=zf(i1,2*i3-n3/2)
     end do
  end do

end subroutine halfill_upcorn

subroutine mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,zmpi1,zw)
  implicit none
  integer :: j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3
  real(kind=8) :: zmpi1(2,n1/2,md2/1,nd3/1,1),zw(2,lot,n1)
  integer :: mfft, Jp2, J2, I1

! WARNING: Assuming that high frequencies are in the corners
!          and that n1 is multiple of 2

  mfft=0
  do Jp2=Jp2stb,1
     do J2=J2stb,md2/1
        mfft=mfft+1
        if (mfft.gt.nfft) then
        Jp2stb=Jp2
        J2stb=J2
        return
        endif
        do I1=1,n1/2
           zw(1,mfft,I1)=0.d0
           zw(2,mfft,I1)=0.d0
        end do
        do I1=n1/2+1,n1
           zw(1,mfft,I1)=zmpi1(1,I1-n1/2,J2,j3,Jp2)
           zw(2,mfft,I1)=zmpi1(2,I1-n1/2,J2,j3,Jp2)
        end do
     end do
     J2stb=1
  end do

end subroutine mpiswitch_upcorn

subroutine multkernel(nd1,nd2,n1,n2,lot,nfft,jS,pot,zw)
  implicit none
  !Argments
  integer, intent(in) :: nd1,nd2,n1,n2,lot,nfft,jS
  real(kind=8), dimension(nd1,nd2), intent(in) :: pot
  real(kind=8), dimension(2,lot,n2), intent(inout) :: zw
  !Local variables
  integer :: j,j1,i2,j2,isign

  !Body

  !case i2=1
  do j=1,nfft
     j1=j+jS-1
     !isign=(j1/(n1/2+2))
     !j1=(1-2*isign)*j1+isign*(n1+2) !n1/2+1-abs(n1/2+2-jS-i1)
     j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!$     j1=n1/2+1-abs(n1/2+2-jS-j)!this stands for j1=min(jS-1+j,n1+3-jS-j)
     zw(1,j,1)=zw(1,j,1)*pot(j1,1)
     zw(2,j,1)=zw(2,j,1)*pot(j1,1)
  end do

  !generic case
  do i2=2,n2/2
     do j=1,nfft
        j1=j+jS-1
        j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!$        j1=n1/2+1-abs(n1/2+2-jS-j)
        j2=n2+2-i2
        zw(1,j,i2)=zw(1,j,i2)*pot(j1,i2)
        zw(2,j,i2)=zw(2,j,i2)*pot(j1,i2)
        zw(1,j,j2)=zw(1,j,j2)*pot(j1,i2)
        zw(2,j,j2)=zw(2,j,j2)*pot(j1,i2)
     end do
  end do

  !case i2=n2/2+1
  do j=1,nfft
     j1=j+jS-1
     j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!$     j1=n1/2+1-abs(n1/2+2-jS-j)
     j2=n2/2+1
     zw(1,j,j2)=zw(1,j,j2)*pot(j1,j2)
     zw(2,j,j2)=zw(2,j,j2)*pot(j1,j2)
  end do

end subroutine multkernel

subroutine switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
  implicit none
  integer :: nfft,n2,lot,n1,lzt,i,j
  real(kind=8) :: zw(2,lot,n2),zt(2,lzt,n1)
! WARNING: Assuming that high frequencies are in the corners
!          and that n2 is multiple of 2

! Low frequencies
  do j=1,nfft
     do i=n2/2+1,n2
        zw(1,j,i)=zt(1,i-n2/2,j)
        zw(2,j,i)=zt(2,i-n2/2,j)
     end do
  end do

! High frequencies
  do i=1,n2/2
     do j=1,nfft
        zw(1,j,i)=0.d0
        zw(2,j,i)=0.d0
     end do
  end do

end subroutine switch_upcorn

subroutine unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf,scal)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(2,lot,n3/2), intent(in) :: zw
  real(kind=8), dimension(md1,md3),intent(inout) :: zf
  real(kind=8), intent(in) :: scal
  !real(kind=8), intent(out) :: ehartreetmp
  !Local variables
  integer :: i3,i1

  !Body
  !ehartreetmp=0.d0
  do i3=1,n3/4
     do i1=1,nfft
      zf(i1,2*i3-1)= scal*zw(1,i1,i3)
      zf(i1,2*i3)= scal*zw(2,i1,i3)
     enddo
  end do

end subroutine unfill_downcorn

subroutine unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,zw,zmpi1)
  implicit none
  integer j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,I1,J2,Jp2,mfft
  real(kind=8) :: zmpi1(2,n1/2,md2,nd3,1),zw(2,lot,n1)
! WARNING: Assuming that high frequencies are in the corners
!          and that n1 is multiple of 2

  mfft=0
  do Jp2=Jp2stf,1
     do J2=J2stf,md2/1
        mfft=mfft+1
        if (mfft.gt.nfft) then
           Jp2stf=Jp2
           J2stf=J2
           return
        endif
        do I1=1,n1/2
           zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
           zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
        end do
     end do
     J2stf=1
  end do

end subroutine unmpiswitch_downcorn

subroutine unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
  implicit none
  integer :: nfft,n2,lot,n1,lzt,j,i
  real(kind=8) :: zw(2,lot,n2),zt(2,lzt,n1)
! WARNING: Assuming that high frequencies are in the corners
!          and that n2 is multiple of 2

! Low frequencies
  do j=1,nfft
     do i=1,n2/2
        zt(1,i,j)=zw(1,j,i)
        zt(2,i,j)=zw(2,j,i)
     end do
  end do

end subroutine unswitch_downcorn

subroutine unscramble_pack(i1,j2,lot,nfft,n1,n3,md2,nd3,zmpi2,zw,cosinarr)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nd3
  real(kind=8), dimension(2,lot,n3/2), intent(out) :: zw
  real(kind=8), dimension(2,n3/2), intent(in) :: cosinarr
  real(kind=8), dimension(2,n1,md2/1,nd3), intent(in) :: zmpi2
  !Local variables
  integer :: i3,i,indA,indB
  real(kind=8) ::  a,b,c,d,cp,sp,re,ie,ro,io,rh,ih

  !Body

  do i3=1,n3/2
     indA=i3
     indB=n3/2+2-i3
     cp=cosinarr(1,i3)
     sp=cosinarr(2,i3)
     do i=0,nfft-1
        a=zmpi2(1,i1+i,j2,indA)
        b=zmpi2(2,i1+i,j2,indA)
        c= zmpi2(1,i1+i,j2,indB)
        d=-zmpi2(2,i1+i,j2,indB)
        re=(a+c)
        ie=(b+d)
        ro=(a-c)*cp-(b-d)*sp
        io=(a-c)*sp+(b-d)*cp
        rh=re-io
        ih=ie+ro
        zw(1,i+1,indA)=rh
        zw(2,i+1,indA)=ih
     end do
  end do

end subroutine unscramble_pack

subroutine F_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,pot,zf&
             ,scal)!,hgrid)!,ehartree)
  implicit none
!value of ncache for the FFT routines
integer, parameter :: ncache_optimal=8*1024
!flag that states if we must perform the timings or not
!definitions of the timing variambles just in case
integer, parameter :: timing_flag=0
integer :: count_time1,count_time2,count_rate,count_max,number_time,index_time
real(kind=8) :: time_b0,time_b1,serial_time,parallel_time

  !Arguments
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3
  real(kind=8), intent(in) :: scal!,hgrid
  !real(kind=8), intent(out) :: ehartree
  real(kind=8), dimension(nd1,nd2,nd3), intent(in) :: pot
  real(kind=8), dimension(md1,md3,md2), intent(inout) :: zf
  !Local variables
  !Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
  integer :: j1,j2,j3,i1,i2,i3,i,j,inzee,ierr,i_all,i_stat
  real(kind=8) :: twopion!,ehartreetmp
  !work arrays for transpositions
  real(kind=8), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(kind=8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(kind=8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(kind=8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(kind=8), dimension(:,:), allocatable :: btrig1,btrig2,btrig3,&
       ftrig1,ftrig2,ftrig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, &
       after2,now2,before2,after3,now3,before3


  !Body
  ! check input
  if (mod(n1,2).ne.0) stop 'Parallel convolution:ERROR:n1'
  if (mod(n2,2).ne.0) stop 'Parallel convolution:ERROR:n2'
  if (mod(n3,2).ne.0) stop 'Parallel convolution:ERROR:n3'
  if (nd1.lt.n1/2+1) stop 'Parallel convolution:ERROR:nd1'
  if (nd2.lt.n2/2+1) stop 'Parallel convolution:ERROR:nd2'
  if (nd3.lt.n3/2+1) stop 'Parallel convolution:ERROR:nd3'
  if (md1.lt.n1/2) stop 'Parallel convolution:ERROR:md1'
  if (md2.lt.n2/2) stop 'Parallel convolution:ERROR:md2'
  if (md3.lt.n3/2) stop 'Parallel convolution:ERROR:md3'
  if (mod(nd3,1).ne.0) stop 'Parallel convolution:ERROR:nd3'
  if (mod(md2,1).ne.0) stop 'Parallel convolution:ERROR:md2'

  !defining work arrays dimensions

  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4
  lzt=n2/2
  if (mod(n2/2,2).eq.0) lzt=lzt+1
  if (mod(n2/2,4).eq.0) lzt=lzt+1

  !Allocations
  allocate(btrig1(2,nfft_max),stat=i_stat)
  allocate(ftrig1(2,nfft_max),stat=i_stat)
  allocate(after1(7),stat=i_stat)
  allocate(now1(7),stat=i_stat)
  allocate(before1(7),stat=i_stat)
  allocate(btrig2(2,nfft_max),stat=i_stat)
  allocate(ftrig2(2,nfft_max),stat=i_stat)
  allocate(after2(7),stat=i_stat)
  allocate(now2(7),stat=i_stat)
  allocate(before2(7),stat=i_stat)
  allocate(btrig3(2,nfft_max),stat=i_stat)
  allocate(ftrig3(2,nfft_max),stat=i_stat)
  allocate(after3(7),stat=i_stat)
  allocate(now3(7),stat=i_stat)
  allocate(before3(7),stat=i_stat)
  allocate(zw(2,ncache/4,2),stat=i_stat)
  allocate(zt(2,lzt,n1),stat=i_stat)
  allocate(zmpi2(2,n1,md2/1,nd3),stat=i_stat)
  allocate(cosinarr(2,n3/2),stat=i_stat)

  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3/2,btrig3,after3,before3,now3,1,ic3)
  call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
  call ctrig(n2,btrig2,after2,before2,now2,1,ic2)
  do  j=1,n1
     ftrig1(1,j)= btrig1(1,j)
     ftrig1(2,j)=-btrig1(2,j)
  enddo
  do  j=1,n2
     ftrig2(1,j)= btrig2(1,j)
     ftrig2(2,j)=-btrig2(2,j)
  enddo
  do  j=1,n3
     ftrig3(1,j)= btrig3(1,j)
     ftrig3(2,j)=-btrig3(2,j)
  enddo

  !Calculating array of phases for HalFFT decoding
  twopion=8.d0*datan(1.d0)/real(n3,kind=8)
  do i3=1,n3/2
     cosinarr(1,i3)= dcos(twopion*real(i3-1,kind=8))
     cosinarr(2,i3)=-dsin(twopion*real(i3-1,kind=8))
  end do

  !initializing integrals
  !ehartree=0.d0

  ! transform along z axis
  lot=ncache/(2*n3)
  if (lot.lt.1) then
     write(6,*) &
          'convolxc_on:ncache has to be enlarged to be able to hold at' // &
          'least one 1-d FFT of this size even though this will' // &
          'reduce the performance for shorter transform lengths'
     stop
  endif

print *,"md2",md2
print *,"lot",lot
print *,"(n1/2)",(n1/2)
print *,"zf",size(zf)

  do j2=1,md2
     !this condition ensures that we manage only the interesting part for the FFT
     if (j2.le.n2/2) then
        do i1=1,(n1/2),lot
           ma=i1
           mb=min(i1+(lot-1),(n1/2))
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           call halfill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

           !performing FFT
           !input: I1,I3,J2,(Jp2)
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo
           !output: I1,i3,J2,(Jp2)

           !unpacking FFT in order to restore correct result,
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_unpack(i1,j2,lot,nfft,n1/2,n3,md2,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        end do
     endif
  end do

  !now each process perform complete convolution of its planes
  do j3=1,nd3/1
     !this condition ensures that we manage only the interesting part for the FFT
     if (0*(nd3/1)+j3.le.n3/2+1) then
      Jp2stb=1
      J2stb=1
      Jp2stf=1
      J2stf=1

        ! transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) then
           write(6,*) &
                'convolxc_on:ncache has to be enlarged to be able to hold at' // &
                'least one 1-d FFT of this size even though this will' // &
                'reduce the performance for shorter transform lengths'
           stop
        endif

        do j=1,n2/2,lot
           ma=j
           mb=min(j+(lot-1),n2/2)
           nfft=mb-ma+1

           !reverse index ordering, leaving the planes to be transformed at the end
           !input: I1,J2,j3,Jp2,(jp3)
            call mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,zmpi2,zw(1,1,1))
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo

           !storing the last step into zt array
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), &
                btrig1,after1(i),now1(i),before1(i),1)
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis
        lot=ncache/(4*n2)
        if (lot.lt.1) then
           write(6,*) &
                'convolxc_on:ncache has to be enlarged to be able to hold at' // &
                'least one 1-d FFT of this size even though this will' // &
                'reduce the performance for shorter transform lengths'
           stop
        endif

        do j=1,n1,lot
           ma=j
           mb=min(j+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering
           !input: I2,i1,j3,(jp3)
           call switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)

           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo
           !output: i1,i2,j3,(jp3)

           !Multiply with kernel in fourier space
           call multkernel(nd1,nd2,n1,n2,lot,nfft,j,pot(1,1,j3),zw(1,1,inzee))

           !TRANSFORM BACK IN REAL SPACE

           !transform along y axis
           !input: i1,i2,j3,(jp3)
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig2,after2(i),now2(i),before2(i),-1)
              inzee=3-inzee
           enddo

           !reverse ordering
           !input: i1,I2,j3,(jp3)
           call unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
           !output: I2,i1,j3,(jp3)
        end do

        !transform along x axis
        !input: I2,i1,j3,(jp3)
        lot=ncache/(4*n1)
        do j=1,n2/2,lot
           ma=j
           mb=min(j+(lot-1),n2/2)
           nfft=mb-ma+1

           !performing FFT
           i=1
           call fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                ftrig1,after1(i),now1(i),before1(i),-1)

           inzee=1
           do i=2,ic1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig1,after1(i),now1(i),before1(i),-1)
              inzee=3-inzee
           enddo
           !output: I2,I1,j3,(jp3)

           !reverse ordering
           !input: J2,Jp2,I1,j3,(jp3)
              call unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,zw(1,1,inzee),zmpi2)
           ! output: I1,J2,j3,Jp2,(jp3)
        end do
     endif
  end do

print *,"md2",md2
print *,"(n1/2)",(n1/2)
print *,"zf",size(zf)

  !transform along z axis
  !input: I1,J2,i3,(Jp2)
  lot=ncache/(2*n3)
print *,"lot",lot
  do j2=1,md2
     !this condition ensures that we manage only the interesting part for the FFT
     if (j2.le.n2/2) then
        do i1=1,(n1/2),lot
           ma=i1
           mb=min(i1+(lot-1),(n1/2))
           nfft=mb-ma+1

           !reverse ordering and repack the FFT data in order to be backward HalFFT transformed
           !input: I1,J2,i3,(Jp2)
           call unscramble_pack(i1,j2,lot,nfft,n1/2,n3,md2,nd3,zmpi2,zw(1,1,1),cosinarr)
           !output: I1,i3,J2,(Jp2)

           !performing FFT
           !input: I1,i3,J2,(Jp2)
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig3,after3(i),now3(i),before3(i),-1)
              inzee=3-inzee
           enddo
           !output: I1,I3,J2,(Jp2)

           !calculates the exchange correlation terms locally and rebuild the output array
           call unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2)&
                ,scal)!,ehartreetmp)

           !integrate local pieces together
           !ehartree=ehartree+0.5d0*ehartreetmp*hgrid**3
        end do
     endif
  end do
end subroutine F_PoissonSolver

subroutine xc_energy(m1,m2,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,&
     nxcl,nxcr,hx,hy,hz,rhopot,sumpion,zf,zfionxc)

  implicit none

  !Arguments----------------------
  logical, intent(in) :: sumpion
  integer, intent(in) :: m1,m2,m3,nxc,nwb,nxcl,nxcr,nxt,md1,md2,md3
  integer, intent(in) :: nwbl,nwbr
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), dimension(m1,m3,nxt,1), intent(in) :: rhopot
  real(kind=8), dimension(md1,md3,md2), intent(out) :: zf
  real(kind=8), dimension(md1,md3,md2,1), intent(out) :: zfionxc

  !Local variables----------------
  real(kind=8), dimension(:,:,:), allocatable :: exci,d2vxci
  real(kind=8), dimension(:,:,:,:), allocatable :: vxci,dvxci,dvxcdgr
  real(kind=8), dimension(:,:,:,:,:), allocatable :: gradient
  real(kind=8) :: elocal,vlocal,rho,pot,potion,hgrid,facpotion,sfactor
  integer :: npts,i_all,order,offset,i_stat,ispden
  integer :: i1,i2,i3,j1,j2,j3,jp2,jpp2,jppp2
  integer :: ndvxc,nvxcdgr,ngr2


  !Body

  !check for the dimensions
  if (  nwb/=nxcl+nxc+nxcr-2 .or. nxt/=nwbr+nwb+nwbl) then
     print *,'the XC dimensions are not correct'
     print *,'nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr',nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr
     stop
  end if

  !these are always the same
!  nspden=1
  order=1

  !useful for the freeBC case
  hgrid=max(hx,hy,hz)

  !starting point of the density array for the GGA cases in parallel
  offset=nwbl+1

     !case without XC terms
     !distributing the density in the zf array
     print *,"m1",m1
     print *,"m3",m3
     print *,"nxc",nxc
     print *,"zf",size(zf)

     do jp2=1,nxc
        j2=offset+jp2+nxcl-2
        jpp2=jp2
        do j3=1,m3
           do j1=1,m1
              zf(j1,j3,jp2)=rhopot(j1,j3,j2,1)
           end do
           do j1=m1+1,md1
              zf(j1,j3,jp2)=0.d0
           end do
        end do
        do j3=m3+1,md3
           do j1=1,md1
              zf(j1,j3,jp2)=0.d0
           end do
        end do
     end do
     do jp2=nxc+1,md2
        do j3=1,md3
           do j1=1,md1
              zf(j1,j3,jp2)=0.d0
           end do
        end do
     end do

end subroutine xc_energy

subroutine PSolver(n01,n02,n03,hx,hy,hz,&
     rhopot,karray,offset,n1k,n2k,n3k)
  implicit none
  integer, intent(in) :: n01,n02,n03,n1k,n2k,n3k
  real(kind=8), intent(in) :: hx,hy,hz,offset
!  real(kind=8), dimension(*), intent(in) :: karray
  real(kind=8), dimension(n1k,n2k,n3k), intent(in) :: karray
  real(kind=8), dimension(*), intent(inout) :: rhopot
  !local variables
  integer, parameter :: nordgr=4 !the order of the finite-difference gradient (fixed)
  integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
  integer :: i_all,i_stat,ierr,ind,ind2,ind3,ind4
  integer :: i1,i2,i3,j2,istart,iend,i3start,jend,jproc,i3xcsh,is_step,i4,ind2nd
  integer :: nxc,nwbl,nwbr,nxt,nwb,nxcl,nxcr,nlim,j1,j3,jp2,jpp2,offset2
  real(kind=8) :: scal,newoffset,factor,hgrid
  real(kind=8), dimension(:,:,:), allocatable :: zf
  real(kind=8), dimension(:,:,:,:), allocatable :: zfionxc
  integer, dimension(:,:), allocatable :: gather_arr
  real(kind=8), dimension(:), allocatable :: energies_mpi,rhopot_G

  !calculate the dimensions wrt the geocode
          write(*,'(1x,a,3(i5),a,i5,a,i3,a)',advance='no')&
          'PSolver, free  BC, dimensions: ',n01,n02,n03,'   proc',1,'   0:',0,' ...'
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3)

  !array allocations
  i_all=0
  allocate(zf(md1,md3,md2),stat=i_stat)
  allocate(zfionxc(md1,md3,md2,1),stat=i_stat)

  istart=0*(md2/1)
  iend=min((0+1)*md2/1,m2)
  nxc=iend-istart

    nwbl=0
    nwbr=0
    nxcl=1
    nxcr=1

  nwb=nxcl+nxc+nxcr-2
  nxt=nwbr+nwb+nwbl

     !starting address of rhopot in the case of global i/o
     i3start=istart+2-nxcl-nwbl

  !calculate the actual limit of the array for the zero padded FFT
     nlim=n2/2

!!$  print *,'density must go from',min(istart+1,m2),'to',iend,'with n2/2=',n2/2
!!$  print *,'        it goes from',i3start+nwbl+nxcl-1,'to',i3start+nxc-1

      call xc_energy(m1,m2,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,nxcl,nxcr,&
           hx,hy,hz,rhopot(1+n01*n02*(i3start-1)),.true.,zf,zfionxc)

  !this routine builds the values for each process of the potential (zf), multiplying by scal
    !hgrid=max(hx,hy,hz)
    scal=hx*hy*hz/real(n1*n2*n3,kind=8)
    call F_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,karray,zf(1,1,1),scal)
    factor=0.5d0*hx*hy*hz!hgrid**3


  !the value of the shift depends on the distributed i/o or not
    i3xcsh=istart
    is_step=n01*n02*n03

  !recollect the final data
     do j2=1,nxc
        i2=j2+i3xcsh !in this case the shift is always zero for a parallel run
        ind3=(i2-1)*n01*n02
        do i3=1,m3
           ind2=(i3-1)*n01+ind3
           do i1=1,m1
              ind=i1+ind2
              rhopot(ind)=zf(i1,i3,j2)
           end do
        end do
     end do

  write(*,'(a)')'done.'

end subroutine PSolver

subroutine solve( n, rho, V )
    integer, intent(in)                                 :: n
    real ( kind = 8 ), dimension(n,n,n), intent(inout)      :: rho
    real ( kind = 8 ), dimension(n,n,n), intent(out)        :: V

  real(kind=8), dimension(:,:,:), allocatable :: density,potential,karray
  real(kind=8) :: hx,hy,hz,max_diff,length,eh,exc,vxc,hgrid,diff_parser,offset
  real(kind=8) :: ehartree,eexcu,vexcu,diff_par,diff_ser
  integer :: m1,m2,m3,md1,md2,md3,nd1,nd2,nd3,n1,n2,n3,itype_scf,i_all,i_stat
  integer :: i1,i2,i3,j1,j2,j3,i1_max,i2_max,i3_max,ierr
  integer :: i_allocated,l1,nsp1,nsp2,nsp3,n3d,n3p,n3pi,i3xcsh,i3s

    V(:,:,:)=0

    hx=10.0/real(n,kind=8)
    hy=10.0/real(n,kind=8)
    hz=10.0/real(n,kind=8)

    hgrid=max(hx,hy,hz)

    call F_FFT_dimensions(n,n,n,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3)

    allocate(karray(nd1,nd2,nd3))

    itype_scf=14

    call Free_Kernel(n,n,n,n1,n2,n3,nd1,nd2,nd3,hx,hy,hz,itype_scf,karray)

     !Allocations
     !Density
     allocate(density(n,n,n),stat=i_stat)
     !Density then potential
     allocate(potential(n,n,n),stat=i_stat)

     call test_functions(n,n,n,10.0,1.0,hx,hy,hz,density,potential)

     print *,"max_diff",maxval(abs(potential(:,:,:)-density(:,:,:)))

     offset=potential(1,1,1)

     !dimension needed for allocations
     call PS_dim4allocation(n,n,n,n3d,n3p,n3pi,i3xcsh,i3s)

     !apply the Poisson Solver
     call PSolver(n,n,n,hx,hy,hz,&
            density,karray,offset,nd1,nd2,nd3)

     print *,"max_diff",maxval(abs(potential(:,:,:)-density(:,:,:)))

end subroutine solve

