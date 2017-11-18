subroutine fast_GammaG_RPA(GammaGNew, GammaG, G, W0, Beta, rIndex, SpinIndex, Spin2Index, Vol, MaxTauBin)
!OrderSign=-1, FermiLoopSign=-1, therefore TotalSign=1
  implicit none
  integer :: Vol, MaxTauBin, UP, DOWN, UPUP, DOWNDOWN, DOWNUP, UPDOWN
!f2py intent(in) Vol, MaxTauBin
  double precision :: Beta
!f2py intent(in) Beta
  integer :: rIndex(0:Vol-1, 0:Vol-1)
!f2py intent(in) rIndex
  integer :: SpinIndex(0:2-1)
!f2py intent(in) SpinIndex
  integer :: Spin2Index(0:4-1)
!f2py intent(in) Spin2Index
  Complex*16 :: GammaG(0:2-1, 0:Vol-1, 0:MaxTauBin-1, 0:MaxTauBin-1)
!f2py intent(in) GammaG
  Complex*16 :: W0(0:4-1, 0:4-1, 0:Vol-1)
  Complex*16 :: G(0:2-1, 0:2-1, 0:MaxTauBin-1)
!f2py intent(in) G, W0
  Complex*16 :: GammaGNew(0:2-1, 0:Vol-1, 0:MaxTauBin-1, 0:MaxTauBin-1)
!f2py intent(in, out, copy) GammaGNew
  integer :: t3, tin, dtin, tout, dtout, r1, r2, dr, sign
  Complex*16 :: GG, G1, G2

  UP=SpinIndex(0)
  DOWN=SpinIndex(1)
  UPUP=Spin2Index(0)
  DOWNDOWN=Spin2Index(1)
  UPDOWN=Spin2Index(2)
  DOWNUP=Spin2Index(3)

  GammaGNew=(0.0,0.0)

  do t3=0, MaxTauBin-1
    do tin=0, MaxTauBin-1
      dtin=t3-tin
      sign=1
      if(dtin<0) then
        dtin=dtin+MaxTauBin
        sign=-sign
      endif
      G1=0.5*sign*G(UP,UP,dtin)
      sign=1
      dtin=t3-tin-1
      if(dtin<0) then
        dtin=dtin+MaxTauBin
        sign=-sign
      endif
      G1=G1+0.5*sign*G(UP,UP,dtin)

      do tout=0, MaxTauBin-1
        dtout=tout-t3-1
        sign=1
        if(dtout<0) then
            dtout=dtout+MaxTauBin
            sign=-sign
        endif
        G2=0.5*sign*G(UP,UP,dtout)
        dtout=tout-t3
        sign=1
        if(dtout<0) then
            dtout=dtout+MaxTauBin
            sign=-sign
        endif
        G2=G2+0.5*sign*G(UP,UP,dtout)
        GG=G1*G2

        do r1=0, Vol-1
          do r2=0, Vol-1
            dr=rIndex(r1,r2)
            GammaGNew(UP,r1,tout,tin)=GammaGNew(UP, r1, tout,tin)+GG*W0(UPUP,UPUP,dr)*GammaG(UP,r2,t3,t3)
            GammaGNew(UP,r1,tout,tin)=GammaGNew(UP, r1, tout,tin)+GG*W0(UPUP,DOWNDOWN,dr)*GammaG(DOWN,r2,t3,t3)
            GammaGNew(DOWN,r1,tout,tin)=GammaGNew(DOWN,r1,tout,tin)+GG*W0(DOWNDOWN,UPUP,dr)*GammaG(UP,r2,t3,t3)
            GammaGNew(DOWN,r1,tout,tin)=GammaGNew(DOWN,r1,tout,tin)+GG*W0(DOWNDOWN,DOWNDOWN,dr)*GammaG(DOWN,r2,t3,t3)
          enddo
        enddo
      enddo
    enddo
  enddo
  GammaGNew=GammaGNew*Beta/MaxTauBin
end subroutine

subroutine fast_WWGammaW(GammaW, W0, W, Beta, rIndex, SpinIndex, Spin2Index, Vol, MaxTauBin)
! receive a GammaW object, multiple by two W, then return the same GammaW object
  implicit none
  integer :: Vol, MaxTauBin, UP, DOWN, UPUP, DOWNDOWN, DOWNUP, UPDOWN
!f2py intent(in) Vol, MaxTauBin
  double precision :: Beta
!f2py intent(in) Beta
  integer :: rIndex(0:Vol-1, 0:Vol-1)
!f2py intent(in) rIndex
  Complex*16 :: GammaW(0:6-1, 0:Vol-1, 0:Vol-1, 0:MaxTauBin-1, 0:MaxTauBin-1)
!f2py intent(in, out, copy) GammaW
  integer :: SpinIndex(0:2-1)
!f2py intent(in) SpinIndex
  integer :: Spin2Index(0:4-1)
!f2py intent(in) Spin2Index
  Complex*16 :: W0(0:4-1, 0:4-1, 0:Vol-1)
  Complex*16 :: W(0:4-1, 0:4-1, 0:Vol-1, 0:MaxTauBin-1)
!f2py intent(in) W0, W
  Complex :: WGammaW(0:6-1, 0:Vol-1, 0:Vol-1, 0:MaxTauBin-1, 0:MaxTauBin-1)
  integer :: r, rout, t, tout, dt_out, rin, tin, dt_in, dr_out, dr_in
  !complex :: Wout(0:4-1, 0:4-1)
  !complex :: Win(0:4-1, 0:4-1)
  complex*16 :: Wuuuu, Wuudd, Wuddu
  double precision :: deltaT

  UP=SpinIndex(0)
  DOWN=SpinIndex(1)
  UPUP=Spin2Index(0)
  DOWNDOWN=Spin2Index(1)
  UPDOWN=Spin2Index(2)
  DOWNUP=Spin2Index(3)
  WGammaW=(0.0, 0.0)
  deltaT=Beta/MaxTauBin

  print *, "calculating WGammaW with f2py..."
  do r=0, Vol-1
    do rout=0, Vol-1
      dr_out = rIndex(r, rout)
      do t=0, MaxTauBin-1
        do tout=0, MaxTauBin-1
            !if(real(GammaW(0, rout, rin, tout, tin)<1e-10) .and. aimag(GammaW(0, rout, rin, tout, tin)<1e-10)) then
              !WGammaW(:, 

            dt_out = t - tout -1
            if(dt_out<0) dt_out=dt_out+MaxTauBin
            Wuuuu=0.5*W(UPUP, UPUP, dr_out, dt_out)*deltaT

            dt_out = t - tout
            if(dt_out<0) dt_out=dt_out+MaxTauBin
            Wuuuu=Wuuuu+0.5*W(UPUP,UPUP,dr_out,dt_out)*deltaT

            if(t == tout) Wuuuu = Wuuuu+ W0(UPUP,UPUP,dr_out)

            Wuudd=-Wuuuu
            Wuddu=2.0*Wuuuu

            do rin=0, Vol-1
              do tin=0, MaxTauBin-1
                ! UPUP UPUP
                WGammaW(0, r, rin, t, tin)=WGammaW(0, r, rin, t, tin)+Wuuuu*GammaW(0, rout, rin, tout, tin)

                ! DOWNDOWN DOWNDOWN
                WGammaW(1, r, rin, t, tin)=WGammaW(1, r, rin, t, tin)+Wuuuu*GammaW(1, rout, rin, tout, tin)

                ! out:UPUP in:DOWNDOWN 
                WGammaW(2, r, rin, t, tin)=WGammaW(2, r, rin, t, tin)+Wuudd*GammaW(1, rout, rin, tout, tin)

                ! out:DOWNDOWN in:UPUP
                WGammaW(3, r, rin, t, tin)=WGammaW(3, r, rin, t, tin)+Wuudd*GammaW(0, rout, rin, tout, tin)

                ! out:UPDOWN in:DOWNUP
                WGammaW(4, r, rin, t, tin)=WGammaW(4, r, rin, t, tin)+Wuddu*GammaW(5, rout, rin, tout, tin)

                ! out:DOWNUP in:UPDOWN
                WGammaW(5, r, rin, t, tin)=WGammaW(5, r, rin, t, tin)+Wuddu*GammaW(4, rout, rin, tout, tin)
              enddo
            enddo
        enddo
      enddo
    enddo
  enddo

  GammaW=(0.0, 0.0)

  print *,"calculating WWGammaW with f2py..."
  do r=0, Vol-1
    do rin=0, Vol-1
      dr_in = rIndex(rin, r)
      do t=0, MaxTauBin-1
        do tin=0, MaxTauBin-1
          dt_in=tin-t
          if(dt_in < 0) dt_in =dt_in+MaxTauBin
          !Win = 0.5*W(:,:,dr_in,dt_in)*(Beta/MaxTauBin)
          Wuuuu = 0.5*W(UPUP,UPUP,dr_in,dt_in)*deltaT
          dt_in = tin - t-1
          if(dt_in < 0) dt_in =dt_in+MaxTauBin
          Wuuuu = Wuuuu + 0.5*W(UPUP,UPUP,dr_in,dt_in)*deltaT
          !Win=Win+0.5*W(:,:,dr_in,dt_in)*(Beta/MaxTauBin)

          if(t==tin) Wuuuu=Wuuuu+W0(UPUP,UPUP,dr_in)

          Wuudd=-Wuuuu
          Wuddu=2.0*Wuuuu

          do rout=0, Vol-1
            do tout=0, MaxTauBin-1

              ! out:UPUP in:UPUP
              GammaW(0, rout, r, tout, t) =GammaW(0, rout, r, tout, t)+Wuuuu*WGammaW(0, rout, rin, tout, tin)
              GammaW(0, rout, r, tout, t) =GammaW(0, rout, r, tout, t)+Wuudd*WGammaW(2, rout, rin, tout, tin)

              ! out:DOWNDOWN in:DOWNDOWN
              GammaW(1, rout, r, tout, t)= GammaW(1, rout, r, tout, t)+Wuudd*WGammaW(3, rout, rin, tout, tin)
              GammaW(1, rout, r, tout, t)= GammaW(1, rout, r, tout, t)+Wuuuu*WGammaW(1, rout, rin, tout, tin)

              ! out:UPUP in:DOWNDOWN
              GammaW(2, rout, r, tout, t)=GammaW(2, rout, r, tout, t)+Wuudd*WGammaW(0, rout, rin, tout, tin)
              GammaW(2, rout, r, tout, t)=GammaW(2, rout, r, tout, t)+Wuuuu*WGammaW(2, rout, rin, tout, tin)

              ! out:DOWNDOWN in:UPUP
              GammaW(3, rout, r, tout, t)=GammaW(3, rout, r, tout, t)+Wuuuu*WGammaW(3, rout, rin, tout, tin)
              GammaW(3, rout, r, tout, t)=GammaW(3, rout, r, tout, t)+Wuudd*WGammaW(1, rout, rin, tout, tin)

              ! out:UPDOWN in:DOWNUP
              GammaW(4, rout, r, tout, t)=GammaW(4, rout, r, tout, t)+Wuddu*WGammaW(4, rout, rin, tout, tin)

              ! out:DOWNUP in:UPDOWN
              GammaW(5, rout, r, tout, t)=GammaW(5, rout, r, tout, t)+Wuddu*WGammaW(5, rout, rin, tout, tin)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  GammaW=-1.0*GammaW
end subroutine

subroutine fast_fouier_WWGammaW(GammaW, W, Beta, rIndex, SpinIndex, Spin2Index, Vol, MaxTauBin)
! receive a GammaW object, multiple by two W, then return the same GammaW object
! all objects are in Matsubara frequencies
  implicit none
  integer :: Vol, MaxTauBin, UP, DOWN, UPUP, DOWNDOWN, DOWNUP, UPDOWN
!f2py intent(in) Vol, MaxTauBin
  double precision :: Beta
!f2py intent(in) Beta
  integer :: rIndex(0:Vol-1, 0:Vol-1)
!f2py intent(in) rIndex
  Complex*16 :: GammaW(0:6-1, 0:Vol-1, 0:Vol-1, 0:MaxTauBin-1, 0:MaxTauBin-1)
!f2py intent(in, out, copy) GammaW
  integer :: SpinIndex(0:2-1)
!f2py intent(in) SpinIndex
  integer :: Spin2Index(0:4-1)
!f2py intent(in) Spin2Index
  Complex*16 :: W(0:4-1, 0:4-1, 0:Vol-1, 0:MaxTauBin-1)
!f2py intent(in) W0, W
  Complex*16 :: WGammaW(0:6-1, 0:Vol-1, 0:Vol-1, 0:MaxTauBin-1, 0:MaxTauBin-1)
  double precision :: deltaT
  integer :: kout, wout, kin, win

  UP=SpinIndex(0)
  DOWN=SpinIndex(1)
  UPUP=Spin2Index(0)
  DOWNDOWN=Spin2Index(1)
  UPDOWN=Spin2Index(2)
  DOWNUP=Spin2Index(3)
  WGammaW=(0.0, 0.0)
  deltaT=Beta/MaxTauBin

  print *, "calculating WGammaW with fouier and f2py..."
  do kout=0, Vol-1
    do wout=0, Vol-1
      do kin=0, Vol-1
        do win=0, Vol-1
          ! UPUP UPUP
          WGammaW(0, kout, kin, wout, win)  = W(UPUP, UPUP, kout, wout) * GammaW(0, kout, kin, wout, win)

          ! DOWNDOWN DOWNDOWN
          WGammaW(1, kout, kin, wout, win)  = W(DOWNDOWN, DOWNDOWN, kout, wout) * GammaW(1, kout, kin, wout, win)

          ! out:UPUP in:DOWNDOWN 
          WGammaW(2, kout, kin, wout, win)  = W(UPUP, DOWNDOWN, kout, wout) * GammaW(1, kout, kin, wout, win)

          ! out:DOWNDOWN in:UPUP
          WGammaW(3, kout, kin, wout, win)  = W(DOWNDOWN, UPUP, kout, wout) * GammaW(0, kout, kin, wout, win)

          ! out:UPDOWN in:DOWNUP
          WGammaW(4, kout, kin, wout, win)  = W(UPDOWN, DOWNUP, kout, wout) * GammaW(5, kout, kin, wout, win)

          ! out:DOWNUP in:UPDOWN
          WGammaW(5, kout, kin, wout, win)  = W(DOWNUP, UPDOWN, kout, wout) * GammaW(4, kout, kin, wout, win)
        enddo
      enddo
    enddo
  enddo
  
  print *, "calculating WWGammaW with fouier and f2py..."
  do kout=0, Vol-1
    do wout=0, Vol-1
      do kin=0, Vol-1
        do win=0, Vol-1
          ! out:UPUP in:UPUP
          GammaW(0, kout, kin, wout, win)=W(UPUP, UPUP, kin, win)*WGammaW(0,kout,kin,wout,win) &
            +W(UPUP,DOWNDOWN,kin, win)*WGammaW(2, kout, kin, wout, win)

          ! out:DOWNDOWN in:DOWNDOWN
          GammaW(1, kout, kin, wout, win)=W(DOWNDOWN, UPUP,kin, win)*WGammaW(3, kout, kin, wout, win) &
            +W(DOWNDOWN, DOWNDOWN,kin,win)*WGammaW(1, kout, kin, wout, win)

          ! out:UPUP in:DOWNDOWN
          GammaW(2, kout, kin, wout, win) = W(DOWNDOWN, UPUP,kin, win) * WGammaW(0, kout, kin, wout, win) &
            + W(DOWNDOWN,DOWNDOWN,kin, win)*WGammaW(2,kout, kin,wout, win)

          ! out:DOWNDOWN in:UPUP
          GammaW(3, kout, kin, wout, win) = W(UPUP, UPUP,kin, win) * WGammaW(3, kout, kin, wout, win) &
            + W(UPUP, DOWNDOWN, kin, win)*WGammaW(1,kout, kin,wout, win)

          ! out:UPDOWN in:DOWNUP
          GammaW(4, kout, kin, wout, win) = W(UPDOWN, DOWNUP, kin, win) * WGammaW(4, kout, kin, wout, win)

          ! out:DOWNUP in:UPDOWN
          GammaW(5, kout, kin, wout, win) = W(DOWNUP, UPDOWN, kin, win) * WGammaW(5, kout, kin, wout, win)
        enddo
      enddo
    enddo
  enddo

  GammaW=-1.0*GammaW
end subroutine
