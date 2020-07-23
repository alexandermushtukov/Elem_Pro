!=======================================================================================
! The subroutine calculates approximate(!) neutrino emissivity due to pair annihilation.
! The calculations are based on the approximation given in Kaminker+ 1992,Ph.Rev.D,46,10
! See eq.(26) there. The temperature should be well above Fermi temperature.
!=======================================================================================
subroutine NeutrinoEmAnnih(T10,B12,Q23)
implicit none
real*8,intent(in)::T10,B12
real*8::Q23,Qc23,C2V,C2A,t,b
real*8::C_Ve,C_Ae,C_Vo,C_Ao
real*8::pi=3.141592653589793d0

  t=T10/0.593d0  !==see Kaminker 1992==!
  b=B12/44.12d0
  Qc23=1.015d0
  C_Ve=2*0.223d0+0.5d0
  C_Ae=0.5d0
  C_Vo=2*0.223d0-0.5d0
  C_Ao=-0.5d0
  C2V= C_Ve**2 +2*C_Vo**2 !(2*0.223d0+0.5d0)**2 !0.93d0
  C2A= C_Ae**2 +2*C_Ao**2 !0.75d0

  Q23=Qc23/pi**4*( (t**3*C2V*(1.d0+3.75d0*t)+t**4*(C2V+C2A)*P(t))*F(t,b) &
      + t*b**2/6.d0/(1.d0+b)*(C2V+C2A*b/(1.d0+b))*S(t) )*exp(-2.d0/t)
return
contains
  real*8 function S(t)
  implicit none
  real*8,intent(in)::t
  real*8::b1,b2,b3
    b1=1.058d0; b2=0.6701d0; b3=0.9143d0
    S=1.d0+b1*t+b2*t**2+b3*t**3+0.472*t**4
  return
  end function S

  real*8 function P(t)
  implicit none
  real*8,intent(in)::t
  real*8::a1,a2,a3,a4
    a1=3.581d0; a2=39.64d0; a3=24.43d0; a4=36.49d0
    P=1.d0+a1*t+a2*t**2+a3*t**3+a4*t**4+18.75d0*t**5
  return
  end function P

  real*8 function F(t,b)
  implicit none
  real*8,intent(in)::t,b
  real*8::R,c
  dimension R(3),c(3)
    c(1)=3.106d-6
    c(2)=1.491d-3
    c(3)=4.839d-6
    R(1:3)=1.d0+c(1:3)*b/t**2*exp(sqrt(2*b)/3/t)
    F=1.d0/R(1)/R(2)/R(3)
  return
  end function F
end subroutine NeutrinoEmAnnih


!==Subroutine testing Annihilation==!
subroutine Test_NeutrinoEmAnnih()
implicit none
real*8::T_kev,b,T10,lg_TkeV,B12,Q23_app,n_e_ini30,lg_n_e_ini30,Q23,mu_keV
integer::n_sum_max,n_max_e,n_max_p
character(len=100):: file_name

  200 format (8(es11.4,"   "),1(I6,"   "))
  201 format (7(es11.4,"   "),3(I6,"   "))
  202 format (1(A12),(es8.2))
  b=0.2d0
  write(file_name,202)"./res/res_b_",b
  write(*,*)"# file_name=",file_name; write(*,*)
  !!open (unit = 20, file = file_name)
  !!write(20,*)"# b=",b
  !!write(20,*)"# Format: b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,Q23_app,res/Q23_app,n_sum_max"
  !!write(20,*)
  !!close(20)

  lg_TkeV=1.78d0
  do while(lg_TkeV.le.2.d0)
    T_keV=10.d0**lg_TkeV
    T10=T_keV*1.1602d-3 !*1.1602d-3
    B12=b*44.12d0
    if(lg_TkeV.eq.1.78d0)then
      lg_n_e_ini30=2.d0
    else
      lg_n_e_ini30=-5.d0
    end if
    do while(lg_n_e_ini30.le.2.d0)
      n_e_ini30=10.d0**lg_n_e_ini30
      !call EE_pairs(b,n_e_ini30,TkeV,n_e30,n_p30,mu_keV,n_max)
      call NeutrinoEmAnnih(T10,B12,Q23_app)
      call NeutrinoAnn(b,T_keV,n_e_ini30,Q23,mu_keV,n_sum_max,n_max_e,n_max_p)
      open(unit = 20, file = file_name, status = 'old',form='formatted',position="append")
      !write(*,200)b,lg_TkeV,T_keV,lg_n_e_ini30,n_e_ini30,res,Q23_app,res/Q23_app,n_sum_max
      !write(*,201)b,T_keV,n_e_ini30,mu_keV,res,Q23_app,res/Q23_app,n_sum_max
      write(*,201)b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,Q23_app,Q23/Q23_app,n_sum_max,n_max_e,n_max_p
      write(20,201)b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,Q23_app,Q23/Q23_app,n_sum_max,n_max_e,n_max_p
      close(20)
      lg_n_e_ini30=lg_n_e_ini30+0.1d0
    end do
    lg_TkeV=lg_TkeV+2.d-2
    write(*,*)
  end do


return
end subroutine Test_NeutrinoEmAnnih
!==================================================================================================




!==================================================================================================
! The subroutine calculates accurately Annihilation emission of neutrinos in a strong magnetic field.
! The calculations are based on the paper written by Kaminker+ 1992, PhRevD.
! b - dimensionless magnetic filed strength
! TkeV - temperature in [keV]
! res is represented in units of [1.e23 erg/s/cm^3]
! n_sum_max - the maximal n_e+n_p taken into accpount in calculations.
!==================================================================================================
subroutine NeutrinoAnn(b,TkeV,n_e_ini30,res,mu_keV,n_sum_max,n_max_e,n_max_p)
implicit none
real*8,intent(in)::b,TkeV,n_e_ini30
integer::n_e,n_p,n_max,n_max_e,n_max_p,det,n_sum,det_n,i,n_sum_max
real*8::res,res_add,Z1,Z2,T,mu,mu_keV,f1,f2,eps,f,n_e30,n_p30,mas
dimension mas(1000)
real*8::pi=3.141592653589793d0
integer::ifEven    !==function==!
real*8::sum_array  !==function==!

  eps=1.d-2
  !==getting the actual chemical potential==!
  call EE_pairs(b,n_e_ini30,TkeV,n_e30,n_p30,mu_keV,n_max,n_max_e,n_max_p)   
  T=TkeV/511
  mu=mu_keV/511
  !write(*,*)b,TkeV,mu,n_max
  !read(*,*)

  !n_max=1000   !==the actual n_max is found later authomatically==!

  res=0.d0

  !==limits of electron/pozitron momentum, where we integrate==!
  Z1= -sqrt(T**2+2*T)*4 !-2.d0
  Z2=  sqrt(T**2+2*T)*4 !+2.d0

  mas(1:1000)=0.d0; det_n=11
  n_sum=0
  do while(det_n.eq.11)
  !do while(n_sum.le.n_max)
    if(n_sum.le.8)then
      call sum_over_diagonal_par(n_sum,n_max_e,n_max_p,b,T,Z1,Z2,mu,res_add)
      !write(*,*)res_add
    else
      !call sum_over_diagonal(n_sum,b,T,Z1,Z2,mu,res_add)
      call sum_over_diagonal_approx_par(n_sum,n_max_e,n_max_p,b,T,Z1,Z2,mu,res_add)
    end if
    res=res+res_add
    mas(n_sum+1)=res_add
    !write(*,*)"# ",n_sum,res,res_add
    !write(*,*)n_sum,res,res_add

    if(((n_sum.gt.2).and.(ifEven(n_sum+1).eq.0)).or.(res_add.lt.(1.d-5*res)))then
      !write(*,*)"qq",n_sum,(n_sum+1)/2+1,n_sum+1
      if((sum_array(mas,1000,1,n_sum+1).eq.0.d0).or.&
         ((sum_array(mas,1000,(n_sum+1)/2+1,n_sum+1)/sum_array(mas,1000,1,n_sum+1)).lt.1.d-2))then
        det_n=13
        n_sum_max=n_sum
      end if
    end if

    n_sum=n_sum+1
  end do
  res=res*b/3/(2*pi)**5   !==10^(23)==!
return
contains


  subroutine sum_over_diagonal(n_sum,b,T,Z1,Z2,mu,res)
  implicit none
  integer,intent(in)::n_sum
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer::n_e,n_p
  real*8::res
    res=0.d0
    n_e=0
    do while(n_e.le.n_sum)
      n_p=n_sum-n_e
      !==integration over electron momentum==!

      f=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)
      res=res+f
      !write(*,*)" ",n_e,n_p,res,f,f/res
      n_e=n_e+1
    end do
  return
  end subroutine sum_over_diagonal



  subroutine sum_over_diagonal_par(n_sum,n_max_e,n_max_p,b,T,Z1,Z2,mu,res)
  use omp_lib 
  implicit none
  integer,intent(in)::n_sum,n_max_e,n_max_p
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer::n_e,n_p,n_e_task,n_p_task
  real*8::res
    res=0.d0
    n_e=0

    f=0.d0
    call omp_set_num_threads(n_sum+1)
    !$omp parallel private(n_e_task,n_p_task) firstprivate(n_sum) reduction(+: f)
    n_e_task=omp_get_thread_num()
    n_p_task=n_sum-n_e_task
    if((n_e_task.le.n_max_e).and.(n_p_task.le.n_max_p))then
      f=NeutrinoSyn_intZe(b,T,n_e_task,n_p_task,Z1,Z2,mu)
    else
      f=0.d0
    end if
    !$omp end parallel
    res=f
  return
  end subroutine sum_over_diagonal_par



  !==test==!
  subroutine sum_over_diagonal_approx(n_sum,b,T,Z1,Z2,mu,res)
  implicit none
  integer,intent(in)::n_sum
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer::n_e,n_p,M1
  integer::num,num_new,add_num,added_num
  real*8::fun,fun_new
  dimension num(n_sum+1),fun(n_sum+1),num_new(n_sum+1),fun_new(n_sum+1),added_num(n_sum+1)
  integer::i,N1,step_num
  real*8::res,sum_array,sum_array2,sum_array2_new,eps
    !write(*,*)"#0 :",n_sum

    eps=1.d-2
    N1=3
    num(1)=1; n_e=num(1); n_p=n_sum-n_e      
    fun(1)=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)

    num(2)=(n_sum+1)/2; n_e=num(2); n_p=n_sum-n_e
    fun(2)=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)

    num(3)=n_sum-1; n_e=num(3); n_p=n_sum-n_e 
    fun(3)=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)

    sum_array2=13.d0
    sum_array2_new=2*sum_array2 
    step_num=0
    M1=11
    do while((M1.gt.0).and.( (abs(sum_array2_new-sum_array2)/sum_array2_new).gt.eps ) )
      !write(*,*)n_sum+1,N1,(abs(sum_array2_new-sum_array2)/sum_array2_new)
      step_num=step_num+1
      sum_array2=sum_array2_new
      call get_new_points_for_sum(n_sum+1,N1,num,fun,num_new,fun_new,added_num,M1)
      !write(*,*)"#3 M1=",M1
      i=1
      do while(i.le.M1)
        n_e=num_new(added_num(i)) 
        n_p=n_sum-n_e
        fun_new(added_num(i))=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)
        i=i+1
      end do
      num(1:n_sum+1)=num_new(1:n_sum+1)
      fun(1:n_sum+1)=fun_new(1:n_sum+1)
      N1=N1+M1
      call sum_over_incomplete_array2(n_sum+1,N1,num,fun,sum_array2_new)
      !write(*,*)"#4 ",N1,M1,sum_array,sum_array2_new
    end do
    res=sum_array +NeutrinoSyn_intZe(b,T,0,n_sum,Z1,Z2,mu) +NeutrinoSyn_intZe(b,T,n_sum,0,Z1,Z2,mu)
  return
  end subroutine sum_over_diagonal_approx




  subroutine sum_over_diagonal_approx_par(n_sum,n_max_e,n_max_p,b,T,Z1,Z2,mu,res)
  use omp_lib
  implicit none
  integer,intent(in)::n_sum,n_max_e,n_max_p
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer::n_e,n_p,M1,n_e_task,n_p_task
  integer::num,num_new,add_num,added_num
  real*8::fun,fun_new
  dimension num(n_sum+1),fun(n_sum+1),num_new(n_sum+1),fun_new(n_sum+1),added_num(n_sum+1)
  integer::i,N1,step_num
  real*8::res,sum_array,sum_array2,sum_array2_new,eps,fun_task
    !write(*,*)"#0 :",n_sum,n_max_e,n_max_p

    eps=1.d-2
    N1=3
    num(1)=1; n_e=num(1); n_p=n_sum-n_e 
    if((n_e.le.n_max_e).and.(n_p.le.n_max_p))then; fun(1)=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu); else; fun(1)=0.d0; end if

    num(2)=(n_sum+1)/2; n_e=num(2); n_p=n_sum-n_e
    if((n_e.le.n_max_e).and.(n_p.le.n_max_p))then; fun(2)=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu); else; fun(2)=0.d0; end if

    num(3)=n_sum-1; n_e=num(3); n_p=n_sum-n_e 
    if((n_e.le.n_max_e).and.(n_p.le.n_max_p))then; fun(3)=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu); else; fun(3)=0.d0; end if

    sum_array2=13.d0
    sum_array2_new=2*sum_array2 
    step_num=0
    M1=11
    do while((M1.gt.0).and.( (abs(sum_array2_new-sum_array2)/sum_array2_new).gt.eps ) )
      !write(*,*)n_sum+1,N1,(abs(sum_array2_new-sum_array2)/sum_array2_new)
      step_num=step_num+1
      sum_array2=sum_array2_new
      call get_new_points_for_sum(n_sum+1,N1,num,fun,num_new,fun_new,added_num,M1)
      !write(*,*)"#3 M1=",M1

      if(M1.gt.0)then
        call omp_set_num_threads(M1)
        !$omp parallel private(i,n_e_task,n_p_task,fun_task) firstprivate(n_sum,b,T,Z1,Z2,mu)
        i=omp_get_thread_num()+1
        n_e_task=num_new(added_num(i)) 
        n_p_task=n_sum-n_e_task
        if((n_e_task.le.n_max_e).and.(n_p_task.le.n_max_p))then
          fun_task=NeutrinoSyn_intZe(b,T,n_e_task,n_p_task,Z1,Z2,mu)
        else 
          fun_task=0.d0
        end if
        !write(*,*)omp_get_thread_num(),fun_task,M1
        fun_new(added_num(i))=fun_task
        !$omp end parallel
      end if
      num(1:n_sum+1)=num_new(1:n_sum+1)
      fun(1:n_sum+1)=fun_new(1:n_sum+1)
      N1=N1+M1
      call sum_over_incomplete_array2(n_sum+1,N1,num,fun,sum_array2_new)
    end do
    res=sum_array 
    if(n_sum.le.n_max_e)then; res=res+NeutrinoSyn_intZe(b,T,n_sum,0,Z1,Z2,mu); end if 
    if(n_sum.le.n_max_p)then; res=res+NeutrinoSyn_intZe(b,T,0,n_sum,Z1,Z2,mu); end if 
  return
  end subroutine sum_over_diagonal_approx_par



  !==the function integrates over the Z-momentum of the electron/positron==!
  !==todo: it would be better to use res_min as a parameter of a function==!
  recursive real*8 function NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu) result(res)
  implicit none
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer,intent(in)::n_e,n_p
  integer::nn,nn_max
  real*8::res1,res2,res_min
    res_min=1.d-90
    eps=1.d-2
    nn_max=800
    nn=100; det=13
    res1=NeutrinoSyn_intZe_n(b,T,n_e,n_p,Z1,Z2,mu,nn)
    do while(det.gt.11)
      nn=nn*2+1
      res2=NeutrinoSyn_intZe_n(b,T,n_e,n_p,Z1,Z2,mu,nn)
      if(abs((res2-res1)/res2).lt.eps)then
        det=det-1
      end if
      if(nn.gt.n_max)then
        det=11
      end if
      res1=res2
    end do

    !==If we reach the maximum amount of the intervals,
    !==we devide the region into 2 parts and repeate the procedure.
    if((nn.lt.nn_max).or.(res1.lt.res_min))then
      res=res1
    else
      res=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,(Z1+Z2)/2,mu)+NeutrinoSyn_intZe(b,T,n_e,n_p,(Z1+Z2)/2,Z2,mu)
    end if
  return
  end function NeutrinoSyn_intZe



  !====================================================================================
  ! Here we integrate over the Z-momentum of electrons using given amount of intervals.
  ! Simpson method is in use.
  !====================================================================================
  real*8 function NeutrinoSyn_intZe_n(b,T,n_e,n_p,Z1,Z2,mu,nn)
  implicit none
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer,intent(in)::n_e,n_p,nn
  real*8::res,coeff,dZ,Z,f_i,f_i_,Zp1,Zp2
  integer::i,det
  real*8::SimpsonCoeff  !==function==!
    res=0.d0
    dZ=(Z2-Z1)/(nn-1)
    res=0.d0
    i=1
    do while(i.le.nn)
      Z=Z1+dZ*(i-1)
      Zp1= -sqrt(T**2+2*T)
      Zp2=  sqrt(T**2+2*T)
      f_i=NeutrinoSyn_intZp(Zp1,Zp2,Z,b,n_e,n_p,mu,T)
      det=11
      do while(det.eq.11)
        Zp1=2*Zp1
        Zp2=2*Zp2
        f_i_=NeutrinoSyn_intZp(Zp1,Zp2,Z,b,n_e,n_p,mu,T)
        if(((abs(f_i_)+abs(f_i)).eq.0.d0).or.((abs(f_i-f_i_)/max(abs(f_i),abs(f_i_))).le.1.d-2))then
          det=12
        else
          f_i=f_i_
        end if
      end do
      coeff=SimpsonCoeff(i,nn)
      res=res+coeff*f_i*dZ/3
      i=i+1
    end do
    NeutrinoSyn_intZe_n=res
    !write(*,*)"#NeutrinoSyn_intZe_n: ",nn,res
  return
  end function NeutrinoSyn_intZe_n


  recursive real*8 function NeutrinoSyn_intZp(Zp1,Zp2,Z_i,b,n_e,n_p,mu,T) result(res)
  implicit none
  real*8,intent(in)::Zp1,Zp2,Z_i,b,mu,T
  integer,intent(in)::n_e,n_p
  integer::nn,nn_max
  real*8::res1,res2
    eps=1.d-2
    nn_max=200
    nn=20; det=13
    res1=NeutrinoSyn_intZp_n(Zp1,Zp2,Z_i,b,n_e,n_p,mu,T,nn) !int_simpson_n(f,a,b,n)
    do while(det.gt.11)
      nn=nn*2+1
      res2=NeutrinoSyn_intZp_n(Zp1,Zp2,Z_i,b,n_e,n_p,mu,T,nn)
      if( (abs((res2-res1)/res2).lt.eps).or.((res2-res1).eq.0.d0) )then
        det=det-1
      end if
      if(nn.gt.nn_max)then
        det=11
      end if
      res1=res2
    end do

    if(nn.lt.nn_max)then
      res=res1
    else
      res=NeutrinoSyn_intZp(Zp1,(Zp1+Zp2)/2,Z_i,b,n_e,n_p,mu,T)+NeutrinoSyn_intZp((Zp1+Zp2)/2,Zp2,Z_i,b,n_e,n_p,mu,T)
    end if
  return
  end function NeutrinoSyn_intZp



  !=================================================================================================
  ! The function calculates approximately the integral over q_{z} using given number of intervals
  ! and Simpson method.
  !=================================================================================================
  real*8 function NeutrinoSyn_intZp_n(Zp1,Zp2,Z_i,b,n_e,n_p,mu,T,nn)
  implicit none
  real*8,intent(in)::Zp1,Zp2,Z_i,b,mu,T
  integer,intent(in)::n_e,n_p,nn
  real*8::dZp,Zp,coeff,res,f_i,q_p_min,q_p_max,f1,f2,eps,omega,Z_f,q_z
  integer::i,nn_,nn_max,det
  real*8::SimpsonCoeff  !==functions==!
    eps=1.d-2
    nn_max=2000
    dZp=(Zp2-Zp1)/(nn-1)
    res=0.d0
    i=1
    do while(i.le.nn)
      Zp=Zp1+(i-1)*dZp
      q_z=Z_i+Zp
      omega=sqrt(1.d0+2*n_e*b+Z_i**2)+sqrt(1.d0+2*b*n_p+Zp**2)
      coeff=SimpsonCoeff(i,nn)

      !==integration over q_p==!
      q_p_min=0.d0
      q_p_max=sqrt(omega**2-q_z**2)
      nn_=20
      f2=NeutrinoSyn_intQp_n(q_p_min,q_p_max,q_z,Z_i,Zp,omega,b,n_e,n_p,mu,T,nn_)
      det=13
      do while(det.eq.13)
        f1=f2
        nn_=nn_*2
        f2=NeutrinoSyn_intQp_n(q_p_min,q_p_max,q_z,Z_i,Zp,omega,b,n_e,n_p,mu,T,nn_)
        if((abs(f2)).gt.1.d-90)then
          if(((abs((f2-f1)/f2)).lt.eps).or.(nn_.gt.nn_max))then
            det=1
            f_i=f2
          end if
        else
          det=1
          f_i=0.d0
        end if
      end do
      !==end of the integration over q_p==!

      res=res+coeff*f_i*dZp/3
      i=i+1
    end do
    NeutrinoSyn_intZp_n=res
  return
  end function NeutrinoSyn_intZp_n


  !=================================================================================================
  ! The function calculates approximately the integral over q_{\perp} using given number of intervals
  ! and Simpson method.
  !=================================================================================================
  real*8 function NeutrinoSyn_intQp_n(q_p1,q_p2,q_z,Ze,Zp,omega,b,n_e,n_p,mu,T,nn)
  implicit none
  real*8,intent(in)::q_p1,q_p2,q_z,Ze,Zp,omega,b,mu,T
  integer,intent(in)::n_e,n_p,nn
  real*8::dq_p,q_p,f_i,res,coeff
  integer::i
  real*8::SimpsonCoeff  !==function==!
  real*8::IntNeutrinoAnn  !==function==!
    dq_p=(q_p2-q_p1)/(nn-1)
    res=0.d0
    i=1
    do while(i.le.nn)
      q_p=q_p1+(i-1)*dq_p
      coeff=SimpsonCoeff(i,nn)
      f_i=IntNeutrinoAnn(n_e,n_p,Ze,Zp,q_z,q_p,b,mu,T)
      res=res+coeff*f_i*dq_p/3
      i=i+1
    end do
    NeutrinoSyn_intQp_n=res
  return
  end function NeutrinoSyn_intQp_n

end subroutine NeutrinoAnn



Subroutine test_IntAnn()
implicit none
real*8::q_p1,q_p2,q_z,Ze,Zp,omega,b,mu,T,dq_p,res,q_p,coeff,f_i
integer::nn,n_e,n_p,i
real*8::SimpsonCoeff,IntNeutrinoAnn
  q_p1=3.d0 !0.d0
  q_p2=6.d0 !12.727695838029952
  q_z= -3.4142980019670275d0
  Ze=-1.7071490009835137d0
  Zp=-1.7071490009835137d0
  omega=13.177696012264093d0
  b=0.5d0
  n_e=40
  n_p=39
  mu=2.0962818003913899d-4
  T=0.97847358121330719d0
  nn=2000
  dq_p=(q_p2-q_p1)/(nn-1)
  res=0.d0
  i=1
  do while(i.le.nn)
    q_p=q_p1+(i-1)*dq_p
    coeff=SimpsonCoeff(i,nn)
    f_i=IntNeutrinoAnn(n_e,n_p,Ze,Zp,q_z,q_p,b,mu,T)
    !if(nn.gt.4.e2)then
    !  write(*,*)i,q_p,f_i
    !end if
    res=res+coeff*f_i*dq_p/3
    i=i+1
  end do

return
end subroutine test_IntAnn


!=====================================================================================================================
! The function which is used in integration to get the intencity if neutrino emission due to the annihilation process.
!=====================================================================================================================
real*8 function IntNeutrinoAnn(n_e,n_p,Z_e,Z_p,q_z,q_p,b,mu,T)
implicit none
integer,intent(in)::n_e,n_p
real*8,intent(in)::Z_e,Z_p,q_z,q_p,b,mu,T
real*8::E_e,E_p,omega
real*8::res
  if(q_p.ne.0.d0)then
    E_e=sqrt(1.d0+2*b*n_e+Z_e**2)
    E_p=sqrt(1.d0+2*b*n_p+Z_p**2)
    omega=E_e+E_p
    res=q_p*A(b,n_e,n_p,Z_e,Z_p,E_e,E_p,q_z,q_p,omega)*omega*fe(E_e,mu,T)*fp(E_p,mu,T)
    !res=A(b,n_e,n_p,Z_e,Z_p,E_e,E_p,q_z,q_p,omega)
  else
    res=0.d0
  end if
  IntNeutrinoAnn=res
  !write(*,*)res
return
contains
  !==Electron distribution function for a given chemical potential mu, see eq (7) in Kaminker 1992==!
  real*8 function fe(E,mu,T)
  implicit none
  real*8,intent(in)::E,mu,T
    fe=1.d0/(exp((E-mu)/T)+1.d0)
  return
  end function fe

  !==Positron distribution function for a given chemical potential mu, see eq (7) in Kaminker 1992==!
  real*8 function fp(E,mu,T)
  implicit none
  real*8,intent(in)::E,mu,T
    fp=1.d0/(exp((E+mu)/T)+1.d0)
  return
  end function fp


  !==function A given by (24) in Kaminker 1992==!
  real*8 function A(b,n_i,n_f,Z_i,Z_f,E_i,E_f,q_z,q_p,omega)
  implicit none
  real*8,intent(in)::b,Z_i,Z_f,E_i,E_f,q_z,q_p,omega
  integer,intent(in)::n_i,n_f
  real*8::funPhi,funPsi   !==functions==!
  real*8::C_V_e,C_A_e,C_V_o,C_A_o,u,Psi_min,Psi_plus,Phi_min,Phi_plus,p_p_i2,p_p_f2
  real*8::A1,A2,A3
    C_V_e=2*0.223d0+0.5d0
    C_A_e= 0.5d0
    C_V_o=2*0.223d0-0.5d0
    C_A_o=-0.5d0

    u=q_p**2/2/b
    p_p_i2=2*b*n_i
    p_p_f2=2*b*n_f
    Phi_plus=funPhi(n_f,n_i,u,1)
    Phi_min =funPhi(n_f,n_i,u,-1)
    Psi_plus=funPsi(n_f,n_i,u,1)
    Psi_min =funPsi(n_f,n_i,u,-1)

    A1=Psi_plus*( 2*(E_i*E_f-Z_i*Z_f)**2+(E_i*E_f-Z_i*Z_f)*(2.d0+p_p_i2+p_p_f2-2*q_p**2)+q_p**2/2*(q_p**2-p_p_i2-p_p_f2-3.d0) )
    A1=A1+Phi_plus*( (E_i*E_f-Z_i*Z_f)*(1.d0+p_p_i2+p_p_f2) +1.d0+ 0.5d0*(p_p_i2+p_p_f2)*(3.d0-q_p**2+p_p_i2+p_p_f2) )
    A1=A1*(C_V_e**2+C_A_e**2)/E_i/E_f +2*A1*(C_V_o**2+C_A_o**2)/E_i/E_f

    A2=( omega**2-q_z**2-q_p**2/2 )*Psi_plus + ( omega**2/2-q_z**2/2-q_p**2 )*Phi_plus
    A2=A2*(C_V_e**2-C_A_e**2)/E_i/E_f +2*A2*(C_V_o**2-C_A_o**2)/E_i/E_f

    A3=( 2*omega**2-2*q_z**2-3*q_p**2 )*Psi_min + (p_p_i2-p_p_f2)*Phi_min
    A3=A3*(E_i*Z_f-E_f*Z_i)*C_V_e*C_A_e/E_i/E_f +2*A3*(E_i*Z_f-E_f*Z_i)*C_V_o*C_A_o/E_i/E_f

    A=A1+A2-A3
    !write(*,*)A
  return
  end function A

end function IntNeutrinoAnn




!==================================================================================
! Kaminker 1992 (16)
!==================================================================================
real*8 function funPhi(n1,n2,u,sign)
implicit none
integer,intent(in)::n1,n2,sign
real*8,intent(in)::u
real*8::Fnn   !==function==!
  funPhi=(Fnn(n1-1,n2-1,u))**2+sign*(Fnn(n1,n2,u))**2
return
end function funPhi



!==================================================================================
! Kaminker 1992 (16)
!==================================================================================
real*8 function funPsi(n1,n2,u,sign)
implicit none
integer,intent(in)::n1,n2,sign
real*8,intent(in)::u
real*8::Fnn   !==function==!
  funPsi=(Fnn(n1-1,n2,u))**2+sign*(Fnn(n1,n2-1,u))**2
return
end function funPsi




!================================================================================================================
! N  - total number of points for summation
! N1 - number of points which is taken to calculate approximate summ
! num(N) - array containing numbers of taken pointa
! fun(N) - array containing function in a given points
! M1 - a new number of points to be considered to calculate approximate sum (output)
! num_new(N) - array containing numbers of new points, actually we use only first (M1-N1) positions of this array
!================================================================================================================
subroutine get_new_points_for_sum(N,N1,num,fun,num_new,fun_new,added_num,M1)
implicit none
integer,intent(in)::N,N1,num
real*8,intent(in)::fun
dimension num(N),fun(N),num_new(N),fun_new(N),added_num(N)
integer::M1,num_new,added_num,i,j,k,new_1,new_2,add_1,add_2,help1,help2
real*8::eps,sum_1,sum_2,sum_tot,fun_new
  
  eps=1.d-2
  M1=0
  i=1; j=0; k=1
  do while(i.lt.(N1-1))
    !write(*,*)"a1"
    call mf_sum_fn_approx_(num(i+1)-num(i)+1  ,fun(i)  ,fun(i+1),sum_1,new_1,add_1)
    new_1=(new_1-1)+num(i)
    call mf_sum_fn_approx_(num(i+2)-num(i+1)+1,fun(i+1),fun(i+2),sum_2,new_2,add_2)
    new_2=(new_2-1)+num(i+1)   
    call mf_sum_fn_approx_(num(i+2)-num(i)+1  ,fun(i)  ,fun(i+2),sum_tot,help1,help2) 
    if(( abs((sum_1+sum_2-sum_tot)/(sum_1+sum_2)) ).lt.eps)then
      !==we will not add points==!
      !write(*,*)"a1.1",abs((sum_1+sum_2-sum_tot)/(sum_1+sum_2))
      M1=M1+0
      num_new(k)  =num(i)  ; fun_new(k)  =fun(i)
      num_new(k+1)=num(i+1); fun_new(k+1)=fun(i+1)
      num_new(k+2)=num(i+2); fun_new(k+2)=fun(i+2)
      k=k+2
    else
      !==we have to add points==!
      M1=M1+add_1+add_2
      num_new(k)=num(i)  ; fun_new(k)=fun(i)
      if(add_1.eq.1)then; 
        j=j+1; added_num(j)=k+add_1
        num_new(k+add_1)=new_1
      end if
      num_new(k+add_1+1)=num(i+1)  ; fun_new(k+add_1+1)=fun(i+1)
      if(add_2.eq.1)then; 
        j=j+1; added_num(j)=k+add_1+1+add_2
        num_new(k+add_1+1+add_2)=new_2
      end if  
      num_new(k+add_1+add_2+2)=num(i+2); fun_new(k+add_1+add_2+2)=fun(i+2)
      k=k+2+add_1+add_2
      !==we have added a new point to the massinve==!
    end if  
    i=i+2
  end do
  !==if there is only one interval left==!
  if(i.eq.(N1-1))then
    !write(*,*)"a2"
    call mf_sum_fn_approx_(num(i)-num(i-1)+1,fun(i-1),fun(i)  ,sum_1,new_1,add_1) 
    call mf_sum_fn_approx_(num(i+1)-num(i)+1,fun(i)  ,fun(i+1),sum_2,new_2,add_2)   
    call mf_sum_fn_approx_(num(i+1)-num(i-1)+1,fun(i-1),fun(i+1),sum_tot,help1,help2) 
    if(( (sum_1+sum_2-sum_tot)/(sum_1+sum_2) ).lt.eps)then
      !we will not add points
      M1=M1+0
      num_new(k+1)=num(i+1); fun_new(k+1)=fun(i+1)
    else
      ! we have to add points
      M1=M1+add_2
      if(add_2.eq.1)then
        j=j+1; added_num(j)=k+add_2
        num_new(k+add_2)=new_2
      end if  
      num_new(k+add_2+1)=num(i+1); fun_new(k+add_2+1)=fun(i+1)
      !==we have added a new point to the massinve==!
    end if  
  end if
  !==now we have all numbers to add==!
  !write(*,*)"a_res: ",M1
return
end subroutine get_new_points_for_sum



!================================================================================================
! N - a number of elements in a set
! f_1=f(1)
! f_N=f(N)
! number_add - number of points to add
!================================================================================================
subroutine mf_sum_fn_approx_(N,f_1,f_N,sum_app,new_N,number_add)
implicit none
integer,intent(in)::N
real*8,intent(in)::f_1,f_N
real*8::sum_app,df,f_i
integer::new_N,i,number_add
  if(N.gt.2)then
    df=(f_N-f_1)/(N-1)
    i=1; sum_app=0.d0
    do while(i.le.N) 
      f_i=f_1+(i-1)*df
      sum_app=sum_app+f_i
      i=i+1
    end do
    new_N=1+N/2
    number_add=1
  else
    if(N.eq.2)then
      sum_app=f_1+f_N
    end if
    if(N.eq.1)then
      sum_app=f_1
    end if
    new_N=1
    number_add=0
  end if
return
end subroutine mf_sum_fn_approx_



!===================================================================================================
! The subroutine calculates the summ over the elements of array fun(N), which total length is N, but where
! only N1<N elements are known. The numbers of these elements and fun(i) there are reprennted in first N1
! elements of arrays num() and fun().
!===================================================================================================
subroutine sum_over_incomplete_array(N,N1,num,fun,sum_array)
implicit none
integer,intent(in)::N,N1,num
real*8,intent(in)::fun
dimension num(N),fun(N)
integer::i
real*8::sum_array,sum_add
integer::new_N,number_add  !==for help==!
  sum_array=0.d0
  i=1
  do while(i.le.(N1-1))
    call mf_sum_fn_approx_(num(i+1)-num(i)+1,fun(i),fun(i+1),sum_add,new_N,number_add)
    sum_add=sum_add-fun(i+1)
    sum_array=sum_array+sum_add
    i=i+1
  end do
  sum_array=sum_array+fun(N1)
return
end subroutine sum_over_incomplete_array




!========================================================================================================
! The subroutine calculates the summ over the elements of array fun(N), which total length is N, but where
! only N1<N elements are known. The numbers of these elements and fun(i) there are reprennted in first N1
! elements of arrays num() and fun().
! Advansed version of subroutine sum_over_incomplete_array.
!========================================================================================================
subroutine sum_over_incomplete_array2(N,N1,num,fun,sum_array)
implicit none
integer,intent(in)::N,N1,num
real*8,intent(in)::fun
dimension num(N),fun(N)
integer::i,j
real*8::sum_array,sum_add,a,b,c
  sum_array=0.d0
  i=1
  do while(i.le.(N1-1))
    if((num(i+1)-num(i)).le.0)then
      sum_add=0.d0
    else
      if((num(i+1)-num(i)).eq.1)then
        sum_add=fun(i)
      else
        !==there are some points b/w num(i) and num(i+1)==!
        sum_add=0.d0
        if(i.eq.1)then
          !==we are looking at the first point in a row==!
          if(i.le.(N1-2))then
            !==there are two points after i==!
            call get_parabola_coeff(1.d0*num(i),fun(i),1.d0*num(i+1),fun(i+1),1.d0*num(i+2),fun(i+2),a,b,c)
            j=num(i)
            do while(j.lt.num(i+1))
              sum_add=sum_add+( a*j**2 + b*j + c )
              j=j+1
            end do
          else
            !==there are NO two points after i==!
            j=num(i)
            do while(j.lt.num(i+1))
              sum_add=sum_add+( fun(i) +(fun(i+1)-fun(i))/(num(i+1)-num(i))*(j-num(i)) )
              j=j+1
            end do
          end if
        else
          !==the point is not the first in a row==!
          call get_parabola_coeff(1.d0*num(i-1),fun(i-1),1.d0*num(i),fun(i),1.d0*num(i+1),fun(i+1),a,b,c)
          j=num(i)
          do while(j.lt.num(i+1))
            sum_add=sum_add+( a*j**2 + b*j + c )
            j=j+1
          end do
        end if
        !write(*,*)"qqq"
      end if
    end if
    sum_array=sum_array+sum_add
    i=i+1
  end do
  sum_array=sum_array+fun(N1)

return
contains
  !===========================================================================================
  ! The subroutine calculates coeffisients of parabola using coordinates of 3 points. 
  !===========================================================================================
  subroutine get_parabola_coeff(x1,y1,x2,y2,x3,y3,a,b,c)
  implicit none
  real*8,intent(in)::x1,y1,x2,y2,x3,y3
  real*8::a,b,c
    a=( y3 - (x3*(y2-y1)+x2*y1-x1*y2)/(x2-x1) )/( x3*(x3-x1-x2) + x1*x2 )
    b=(y2-y1)/(x2-x1) - a*(x1+x2)
    c=(x2*y1-x1*y2)/(x2-x1) + a*x1*x2
  return
  end subroutine get_parabola_coeff
end subroutine sum_over_incomplete_array2





