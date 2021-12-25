program PhysProc
implicit none
real*8::pi=3.141592653589793d0
real*8::T10,B12,Q23,Q23_1,Q23_10,Q23_100,Q23_1000,ro,x,n30,n30_,log_TkeV,log_n_e
real*8::n_e30,n_p30,n_e_ini30,TkeV,n30app,T_res,T_bil12,T_bil13,T_bil14,beta,mdot10,Ei,res,res_,res_app,theta,b,P23
real*8::n_e30_,n_p30_

real*8::ComptPressure00,ComptPressure01,ComptonEf,ComptPressureM,mu,ComptPressure_Max,ComptPressure_mag00
real*8::laguer_general_new,laguer_general

real*8::bessI0,bessK0,bessI1,bessK1,bessKn,dSigmadOmega_pl2   !==function==!
integer::det,n_max,n_max_e,n_max_p,n_sum_max
  200 format (4(es11.4,"   "),1(I6,"   "),6(es11.4,"   "))
  !call test_Ai()
  !call Test_Compt_Ave()
  !call test_Fnn()
  call Test_NeutrinoEmAnnih()
  !call test_ee_pairs()
  !call Test_ComptPressure_Max()
  !call Test_ComptPressure01()
  goto 122
  !call test_IntAnn(); goto 122
  !call test_Ai(); goto 122
  !goto 333
  !goto 123

  !goto 122

  !call TestSigmaTot()
  !call TestAmp()
  !n_e30=dSigmadOmega_pl2(1,2,1.1d0,1.d0,pi/3,pi/4,0.d0,1.d0)
  !goto 122

  !call T_max_PairInfluence(1.d0/4.412d0,1.d-1,100); goto 122

  !n_e_ini30=1.d-2 !-2
  !b=.5d0!/4.412d0
  !TkeV=10.d0 !973.09753417968750 !5.d0
  !call Test_EE_distribution(n_e_ini30,b,TkeV)
  !goto 122

  ! internalTemp(T0_keV,Tmax_keV,n30)
  !write(*,*)
  !call internalTemp(5.d0,1.d4,1.d4,9.d-6)
  !goto 122

  !n_e_ini30=1.d-5 !-2
  !b=1.d0/4.412d0
  !TkeV=10.d0 !973.09753417968750 !5.d0
  !do while(TkeV.le.1000.d0)
  !  !call NeutrinoEmAnnih(TkeV*1.1606d-3,b*44.12d0,Q23)
  !  call EE_pairs(b,n_e_ini30,TkeV,n_e30,n_p30,mu,n_max)
  !  !call EE_pairs_PressureZ(b,n_e_ini30,TkeV,mu,n_max,P23)
  !  write(*,200)b,Tkev,n_e_ini30,mu,n_max,n_e30,n_p30,P23,1.602d-2*TkeV*n_e_ini30,4.6d-10*TkeV**4
  !  !write(*,*)b,TkeV,Q23
  !  TkeV=TkeV*1.1d0
  !end do
  !!internalTemp(T0_keV,Tmax_keV,n30)
  !goto 122

  !goto 333  !==pairs code==!

  !goto 123  !==Compton pressure==!

  201 format (4(es11.4,"   "),2(I6,"   "))
  202 format (4(es11.4,"   "),1(I6,"   "))
  !open(unit = 23, file = './res/res_Ann', status = 'old')
  !write(*,*) "# b,n_e30,TkeV,res_Q23,n_max_e,n_max_p"; write(*,*)
  !write(23,*)"# b,n_e30,TkeV,res_Q23,n_max_e,n_max_p"; write(23,*)

  write(*,*) "# b,n_e30,TkeV,res_Q23,n_sum_max"; write(*,*)
  write(23,*)"# b,n_e30,TkeV,res_Q23,n_sum_max"; write(23,*)

  !close(23)
  
  TkeV=500.d0 !973.09753417968750 !5.d0

  b=0.6d0
  log_n_e=25.d0
  do while(log_n_e.le.29.1d0) 
    n_e30= 10**(log_n_e-30) !1.d-5 

    log_TkeV=1.d0
    do while(log_TkeV.le.2.01d0)
    !do while(TkeV.le.500.d0)
      TkeV=10.d0**log_TkeV
      !call NeutrinoSyn(b,TkeV,n_e30,res,n_max_e,n_max_p)
      !open(unit = 23, file = './res/res_Syn00', status = 'old',form='formatted',position="append")
      !write(*,201) b,n_e30,TkeV,res,n_max_e,n_max_p  !,res_,res_app!/1.d3
      !write(23,201)b,n_e30,TkeV,res,n_max_e,n_max_p  !,res_,res_app!/1.d3
      !close(23)

!!!call EE_pairs(b,n_e30,TkeV,n_e30_,n_p30_,mu,n_max,n_max_e,n_max_p)
      write(*,*)b,n_e30,TkeV,n_e30_,n_p30_,mu,n_max
      !read(*,*)
      !call NeutrinoAnn(b,TkeV,n_e30,res_,n_sum_max)
      !open(unit = 23, file = './res/res_Ann', status = 'old',form='formatted',position="append")
      !write(*,202) b,n_e30,TkeV,res_,n_sum_max !,n_max_e,n_max_p  !,res_,res_app!/1.d3
      !write(23,202)b,n_e30,TkeV,res_,n_sum_max !,n_max_e,n_max_p  !,res_,res_app!/1.d3
      !close(23)

      !call NeutrinoEmAnnih(TkeV*1.1606d-3,b*44.13d0,res_app)

      !write(*,*)b,n_e30,TkeV,res_app/1.d3
      log_TkeV=log_TkeV+0.05d0
      !TkeV=TkeV*2.d0
    end do
    write(*,*)
    !open(unit = 23, file = './res/res_Ann', status = 'old',form='formatted',position="append"); write(23,*)
    log_n_e=log_n_e+0.1d0
  end do

  goto 122


  123 Ei=.05d0
  theta=1.57d0
  beta=-0.7d0
  write(*,*)
  do while(beta.le.0.7d0)!0.1d0)
    !res=ComptPressure01(Ei,beta,theta)
    !res=ComptPressure_mag00(Ei,beta,0.d0,1)
    !res=ComptPressure_mag00(Ei,beta,pi/4.d0)
   !write(*,*)Ei,beta,theta,ComptPressure_mag00(Ei,beta,theta,1),ComptPressure_mag00(Ei,beta,theta,2),ComptPressure01(Ei,beta,theta)!,ComptPressure_mag00(Ei,beta,pi/4.d0),&
              ! ComptPressure_mag00(Ei,beta,pi/2.d0)
    !write(*,*)Ei,beta,theta,ComptPressure01(Ei,beta,0.d0),ComptPressure01(Ei,beta,pi/4.d0),&
    !          ComptPressure01(Ei,beta,pi/2.d0)
    !write(*,*)Ei,beta,theta,ComptPressure01(Ei,beta,theta),&
    !          ComptPressure_Max(Ei,theta,beta,10.d0),ComptPressure_Max(Ei,theta,beta,100.d0)
    !write(*,*)Ei,beta,theta,ComptPressure_mag00(Ei,beta,theta,1),&
    !          ComptPressure_Max(Ei,theta,beta,100.d0,1)!,ComptPressure_Max(Ei,theta,beta,100.d0,1   )

    !read(*,*)
    beta=beta+0.05d0
  end do
  goto 122

  beta=1.e-4
  mdot10=1.d0
  do while(beta.le.1.e-3)
    call calcT10(mdot10,1.4d0,beta,T_res,T_bil12,T_bil13,T_bil14)
    !write(*,*)mdot10,beta,T_res,T_bil12,T_bil13,T_bil14
    !read(*,*)
    beta=beta*1.2d0
  end do
  goto 122

  333 write(*,*)
  b=0.d0 !1.d0/4.413d0 !0.1d0
  TkeV=1.d1
  do while(TkeV.le.2.d2)
    n_e_ini30=1.d-3 !1.d0
    !call EE_pairs_app1(TkeV/511.d0,n30app)
!!!call EE_pairs(b,n_e_ini30,TkeV,n_e30,n_p30,mu,n_max,n_max_e,n_max_p)
    call EE_pairs_app1(TkeV/511.d0,res)
    call EE_pairs_app2(TkeV/511.d0,res_)
    !write(*,*)TkeV,n_e_ini30,n_e30,n_p30,n_p30/n_e_ini30
    !,n30app/n_e_ini30!, (3.5d-8*exp(-1022.d0/TkeV))/(n_e30*n_p30)
    write(*,*)b,TkeV,mu,n_max,n_e_ini30,n_e30,res,res_  !,1.d0/(1.d0-exp(-2.d0*mu/TkeV))
    TkeV=TkeV*1.1d0
  end do
  goto 122

  write(*,*)
  B12=10.d0
  ro=10.d0
  T10=0.1d0
  do while(T10.le.10.d0)
    call NeutrinoEmAnnih(T10,B12,Q23)
    !call NeutrinoGammaE(T10,B12,ro,Q23)
    write(*,*)T10,Q23/100.d0
    T10=T10*1.1d0
  end do
  !write(*,*)"dssds",pi
122 return
end program PhysProc


subroutine internalTemp(T0_keV,Tmax_keV,x_max,n30)
implicit none
real*8,intent(in)::T0_keV,Tmax_keV,x_max,n30
real*8::t,t_0,t_max,x,dt
integer::n,i
  n=20000
  t_0=T0_keV/511.d0
  t_max=Tmax_keV/511.d0
  dt=(t_max-t_0)/(n-1)
  t=t_0
  x=0
  i=1
  do while((i.le.n).and.(x.le.x_max))
    write(*,*)i,t,x
    i=i+1
    t=t+dt
    x=x+f(n30,t,T0_keV)*dt  !*x_max/(x_max-x)
  end do
return
contains
  real*8 function f(n30,t,T0_keV)
  implicit none
  real*8,intent(in)::n30,t,T0_keV
    f=1.d6/1.82d0/T0_keV**4/(n30+5.35d0*exp(-1.d0/t)*t**1.5d0)*t**3
    !   f=1.d6/1.82d0/T0_keV**4/( n30 )*t**3
  return
  end function f
end subroutine internalTemp



subroutine calcT10(dotM10,m,beta,res_out,res_bil12,res_bil13,res_bil14)
implicit none
real*8,intent(in)::dotM10,m,beta
real*8::T10,res,n_e_ini30,n_e30,n_p30,TkeV,mu,R6,pairs,P24,P24crit
real*8::res_out,res_bil12,res_bil13,res_bil14,B12
integer::det,n_max,n_max_e,n_max_p
  n_e_ini30=2.d-7*dotM10/beta
  R6=1.d0
  T10=1.d-2
  res_out=T10
  res_bil12=T10
  res_bil13=T10
  res_bil14=T10
  B12=1.d0
  det=13
  do while(T10.le.3.2d0)
    TkeV=T10*1.d7/11600.d0
!!!call EE_pairs(0.d0,n_e_ini30,TkeV,n_e30,n_p30,n_max,n_max_e,n_max_p)
    mu=0.5d0-n_p30/(2.d0*(n_e_ini30+n_e30))
    pairs=2.d0*n_p30/(n_e_ini30+n_e30-n_p30)
    !res=2.21d0*m/R6-7.5d0*beta**2-5.053d6*T10**4*beta/dotM10-3.45d-2*T10/mu-8.2d-3*pairs
    res=2.21d0*m/R6-7.5d0*beta**2-5.053d6*T10**4*beta/dotM10-3.45d-2*T10/mu-8.2d-3*pairs
    P24=18.9d0*T10**4+2.75d-7*dotM10*T10/beta/mu
    if((res.lt.0.d0).and.(det.eq.13))then; res_out=T10; det=11; end if
    if((0.32d0*1.d0**2-P24).gt.0.d0)then; res_bil12=T10; end if
    if((0.32d0*1.d1**2-P24).gt.0.d0)then; res_bil13=T10; end if
    if((0.32d0*5.d1**2-P24).gt.0.d0)then; res_bil14=T10; end if
    !write(*,*)dotM10,beta,T10,res,18.9d0*T10**4/(2.75d-7*dotM10*T10/beta/mu)
    write(*,*)dotM10,beta,T10,18.9d0*T10**4/(2.75d-7*dotM10*T10/beta/mu)
    T10=T10+0.002d0!*1.05d0
    !T10=T10+0.02d0!*1.05d0
  end do
  write(*,*)
return
end subroutine



