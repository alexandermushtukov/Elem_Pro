!====================================================================================================
! Subroutine constructs the Data Base of neutrino energy losses due to the annihilation of e-e-pairs.
!====================================================================================================
subroutine Test_NeutrinoEmSyn()
implicit none
real*8::T_kev,b,T10,lg_TkeV,B12,Q23_app,n_e_ini30,lg_n_e_ini30,Q23,mu_keV
integer::n_sum_max,n_max_e,n_max_p,det_res,num_stream,i,det_omp
character(len=100):: file_name
real*8::Z_e_max(5000),Z_p_max(5000),res_omp(10,5)
real*8::IntNeutrinoAnn,x,Fnn  !==function==!
real*8::n_e_ini30_task,Q23_task,mu_keV_task
integer::n_max_e_task,n_max_p_task

  200 format (8(es11.4,"   "),1(I6,"   "))
  201 format (5(es11.4,"   "),4(I6,"   "))
  202 format (1(A22),(es8.2))
  b=2.d0
  write(file_name,202)"./res/res2_nu_EmSyn_bb",b
  write(*,*)"# file_name=",file_name; write(*,*)
  open (unit = 20, file = file_name)
  write(20,*)"# b=",b
  write(20,*)"# Format: b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,n_max_e,n_max_p"
  write(20,*)
  close(20)

  !==set OpenMP parameters====!
  det_omp=1         !== determins if we use OpenMP: det_opm=0 - no , det_opm=1 - yes, we use ==!
  num_stream=12 !10
  !===========================!

  !lg_TkeV=1.d0
  lg_TkeV=2.d0
  !do while(lg_TkeV.lt.2.d0)
  do while(lg_TkeV.le.3.d0)
    T_keV=10.d0**lg_TkeV
    T10=T_keV*1.1602d-3 !*1.1602d-3
    B12=b*44.12d0

    !if(lg_TkeV.eq.2.18d0)then
    !  lg_n_e_ini30=0.2d0
    !else
      lg_n_e_ini30=1.1d0
      !lg_n_e_ini30=-5.d0
    !end if
    do while(lg_n_e_ini30 .le. 2.0013d0)
      n_e_ini30=10.d0**lg_n_e_ini30
      call NeutrinoSyn(b,T_keV,n_e_ini30,Q23,mu_keV,n_max_e,n_max_p,det_omp,num_stream)
      open(unit = 20, file = file_name, status = 'old',form='formatted',position="append")
      !write(*,200)b,lg_TkeV,T_keV,lg_n_e_ini30,n_e_ini30,res,Q23_app,res/Q23_app,n_sum_max
      !write(*,201)b,T_keV,n_e_ini30,mu_keV,res,Q23_app,res/Q23_app,n_sum_max
      write(*,201)b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,n_max_e,n_max_p!,det_res
      write(20,201)b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,n_max_e,n_max_p!,det_res
      close(20)
      lg_n_e_ini30=lg_n_e_ini30+0.1d0
    end do

    lg_TkeV=lg_TkeV+2.d-2
    !write(*,*)
  end do

return
end subroutine Test_NeutrinoEmSyn
!==================================================================================================



!==================================================================================================
! The subroutine calculates accurately Syncrotron emission of neutrinos in a strong magnetic field.
! The calculations are based on the paper written by Kaminker+ 1992, PhRevD.
! b - dimensionless magnetic filed strength
! TkeV - temperature in [keV]
! res is represented in units of [1.e23 erg/s/cm^3]
!==================================================================================================
subroutine NeutrinoSyn(b,TkeV,n_e_ini30,res,mu_keV,n_max_e,n_max_p,det_omp,num_stream)
use omp_lib
implicit none
real*8,intent(in)::b,TkeV,n_e_ini30
integer,intent(in)::det_omp,num_stream    !==related to OpneMP==!
integer::n_i,n_f,det,n_max,n_max_e,n_max_p,det_part,det_n,det_sim,n_steps_omp,omp_step
real*8::res,Z1,Z2,T,mu,mu_keV,f1,f2,eps,f,f_i,n_e30,n_p30,mas(1000),Z_e_max(5000),Z_p_max(5000)
real*8::pi=3.141592653589793d0
integer::ifEven    !==function==!
real*8::sum_array  !==function==!
integer::det_task,n_f_task
real*8::Z1_task,Z2_task,f1_task,f2_task

  det_sim=1    !== det_sim=1 - we do simplify calculations ==!
  eps=2.d-2
  
  !==getting the actual chemical potential==!
  call EE_pairs(b,n_e_ini30,TkeV,n_e30,n_p30,mu_keV,n_max,n_max_e,n_max_p,Z_e_max,Z_p_max)
  T=TkeV/511
  mu=mu_keV/511
  !write(*,*)"mu=",mu_keV,b,n_e_ini30,TkeV,n_max,n_max_e,n_max_p

  res=0.d0
  !==making a sum over the particle type: electrons and positrons==!
  det_part=1         ! det_part.eq.1 - electron, det_part.ne.1
  do while(det_part.le.2)
    n_i=0
    mas(1:1000)=0.d0       !==mas(i) contains input of transitions from (i-1)-level into the total energy loses==!
    det_n=11
    do while(det_n.eq.11)
      f_i=0.d0

      if(det_omp.eq.0)then
        !==calculation without OpenMP==!
        n_f=0
        do while(n_f.lt.n_i)
          !==integration over Z==!
          det=11
          Z1= -sqrt(T**2+2*T)
          Z2=  sqrt(T**2+2*T)
          f2=NeutrinoSyn_intZ(b,T,n_i,n_f,Z1,Z2,mu,det_part)
          do while(det.eq.11)
            f1=f2
            Z1=2*Z1
            Z2=2*Z2
            f2=NeutrinoSyn_intZ(b,T,n_i,n_f,Z1,Z2,mu,det_part)
            if(((abs((f2-f1)/f2)).lt.eps)&
              .or.( (Z2.gt.Z_e_max(n_i+1)).and.(det_part.eq.1) )&
              .or.( (Z2.gt.Z_p_max(n_i+1)).and.(det_part.ne.1) ).or.(f2.eq.0.d0))then
              det=1
            end if
          end do
          !==end integration over Z==!
          f=f2
          f_i=f_i+f
          !res=res+f
          n_f=n_f+1
        end do
      else
        !==we use OpenMP==!
        n_steps_omp=n_i/num_stream
        omp_step=0
        do while(omp_step.le.n_steps_omp)
           if(omp_step.le.n_steps_omp)then
             call omp_set_num_threads(num_stream)
           else
             call omp_set_num_threads(mod(n_i,num_stream))
           end if
           !$omp parallel private(det_task,n_f_task,Z1_task,Z2_task,f1_task,f2_task) reduction(+: f_i)
               n_f_task = omp_step*num_stream + omp_get_thread_num()
               det_task=11
               Z1_task= -sqrt(T**2+2*T)
               Z2_task=  sqrt(T**2+2*T)
               f2_task=NeutrinoSyn_intZ(b,T,n_i,n_f_task,Z1_task,Z2_task,mu,det_part)
               do while(det_task.eq.11)
                 f1_task=f2_task
                 Z1_task=2*Z1_task
                 Z2_task=2*Z2_task
                 f2_task=NeutrinoSyn_intZ(b,T,n_i,n_f_task,Z1_task,Z2_task,mu,det_part)
                 if(((abs((f2_task-f1_task)/f2_task)).lt.eps)&
                   .or.( (Z2_task.gt.Z_e_max(n_i+1)).and.(det_part.eq.1) )&
                   .or.( (Z2_task.gt.Z_p_max(n_i+1)).and.(det_part.ne.1) ).or.(f2_task.eq.0.d0))then
                   det_task=1
                 end if
               end do
               !==end integration over Z==!
               f_i=f2_task
               !write(*,*)n_f_task, omp_step , num_stream , omp_get_thread_num(),f2_task
           !$omp end parallel
           omp_step=omp_step+1
        end do
      end if
      res=res+f_i
      mas(n_i+1)=f_i     !==saving the contribution of transitions from i-th level==!
      if((n_i.gt.4).and.(ifEven(n_i+1).eq.0))then
        if((sum_array(mas,1000,1,n_i+1).eq.0.d0).or.&
               ((sum_array(mas,1000,(n_i+1)/2,n_i+1)/sum_array(mas,1000,1,n_i+1)).lt.5.d-2))then
          det_n=13
        end if
      end if
      if(( (det_part.eq.1).and.(n_i.ge.n_max_e) ).or.( (det_part.ne.1).and.(n_i.ge.n_max_p) ))then
        det_n=14
      end if

      n_i=n_i+1
    end do
    det_part=det_part+1
  end do

  res=res*b/3/(2*pi)**5
return
contains

  !==the function integrates over the Z-momentum of the electron/positron==!
  recursive real*8 function NeutrinoSyn_intZ(b,T,n_i,n_f,Z1,Z2,mu,det_part) result(res)
  implicit none
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer,intent(in)::n_i,n_f,det_part
  integer::nn,nn_max
  real*8::res1,res2
    eps=2.d-2
    nn_max=3200
    nn=20; det=13
    res1=NeutrinoSyn_intZn(b,T,n_i,n_f,Z1,Z2,mu,nn,det_part)
    do while(det.gt.11)
      nn=nn*2+1
      res2=NeutrinoSyn_intZn(b,T,n_i,n_f,Z1,Z2,mu,nn,det_part)
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
    if(nn.lt.nn_max)then
      res=res1
    else
      res=NeutrinoSyn_intZ(b,T,n_i,n_f,Z1,(Z1+Z2)/2,mu,det_part)+NeutrinoSyn_intZ(b,T,n_i,n_f,(Z1+Z2)/2,Z2,mu,det_part)
    end if

  return
  end function NeutrinoSyn_intZ

  !=============================================================================================
  ! Here we integrate over the Z-momentum of electrons using a given amount of intervals.
  !=============================================================================================
  real*8 function NeutrinoSyn_intZn(b,T,n_i,n_f,Z1,Z2,mu,nn,det_part)
  implicit none
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer,intent(in)::n_i,n_f,nn,det_part
  real*8::res,coeff,dZ,Z,f_i,q_z_min,q_z_max,eps,f1,f2
  integer::i,nn_,nn_max,det
  real*8::SimpsonCoeff  !==function==!
    eps=2.d-2
    nn_max=int(1.e6)
    res=0.d0
    dZ=(Z2-Z1)/(nn-1)
    res=0.d0
    i=1
    do while(i.le.nn)
      Z=Z1+dZ*(i-1)

      !== checking if omega>0 ==!
      if((Z**2+2*b*(n_i-n_f)).ge.0.d0)then
        q_z_min=Z-sqrt(Z**2+2*b*(n_i-n_f))
        q_z_max=Z+sqrt(Z**2+2*b*(n_i-n_f))
        det=13
        nn_=20
        f2=NeutrinoSyn_intQz_n(q_z_min,q_z_max,Z,b,n_i,n_f,mu,T,nn_,det_part)
        do while(det.eq.13)
          f1=f2
          nn_=nn_*2
          f2=NeutrinoSyn_intQz_n(q_z_min,q_z_max,Z,b,n_i,n_f,mu,T,nn_,det_part)
          if((abs(f2)).gt.1.d-40)then
            if(((abs((f2-f1)/f2)).lt.eps).or.(nn_.gt.nn_max))then
              det=1
              f_i=f2
            end if
          else
            det=1
            f_i=0.d0
          end if
        end do
      else
        f_i=0.d0
      end if
      coeff=SimpsonCoeff(i,nn)
      res=res+coeff*f_i*dZ/3
      i=i+1
    end do
    NeutrinoSyn_intZn=res
  return
  end function NeutrinoSyn_intZn


  !=================================================================================================
  ! The function calculates approximately the integral over q_{z} using given number of intervals
  ! and Simpson method.
  !=================================================================================================
  real*8 function NeutrinoSyn_intQz_n(q_z1,q_z2,Z_i,b,n_i,n_f,mu,T,nn,det_part)
  implicit none
  real*8,intent(in)::q_z1,q_z2,Z_i,b,mu,T
  integer,intent(in)::n_i,n_f,nn,det_part
  real*8::dq_z,coeff,q_z,res,f_i,q_p_min,q_p_max,f1,f2,eps,omega,Z_f
  integer::i,nn_,nn_max,det
  real*8::SimpsonCoeff  !==functions==!
    eps=2.d-2
    nn_max=int(1.e6)
    dq_z=(q_z2-q_z1)/(nn-1)
    res=0.d0
    i=1
    do while(i.le.nn)
      q_z=q_z1+(i-1)*dq_z
      Z_f=Z_i-q_z
      omega= sqrt(1.d0+Z_i**2+2*b*n_i) - sqrt(1.d0+Z_f**2+2*b*n_f)

      coeff=SimpsonCoeff(i,nn)

      if(((omega**2-q_z**2).ge.0.d0).and.(omega.gt.0.d0))then
        q_p_min=0.d0
        q_p_max=sqrt(omega**2-q_z**2)
        nn_=20
        f2=NeutrinoSyn_intQp_n(q_p_min,q_p_max,q_z,Z_i,omega,b,n_i,n_f,mu,T,nn_,det_part)
        det=13
        do while(det.eq.13)
          f1=f2
          nn_=nn_*2
          f2=NeutrinoSyn_intQp_n(q_p_min,q_p_max,q_z,Z_i,omega,b,n_i,n_f,mu,T,nn_,det_part)
          if(abs(f2).gt.1.d-40)then
            if(((abs((f2-f1)/f2)).lt.eps).or.(nn_.gt.nn_max))then
              det=1
              f_i=f2
            end if
          else
            det=1
            f_i=0.d0
          end if
        end do
      else
        f_i=0.d0
      end if
      res=res+coeff*f_i*dq_z/3
      i=i+1
    end do
    NeutrinoSyn_intQz_n=res
  return
  end function NeutrinoSyn_intQz_n


  !=================================================================================================
  ! The function calculates approximately the integral over q_{\perp} using given number of intervals
  ! and Simpson method.
  !=================================================================================================
  real*8 function NeutrinoSyn_intQp_n(q_p1,q_p2,q_z,Z_i,omega,b,n_i,n_f,mu,T,nn,det_part)
  implicit none
  real*8,intent(in)::q_p1,q_p2,q_z,Z_i,omega,b,mu,T
  integer,intent(in)::n_i,n_f,nn,det_part
  real*8::dq_p,q_p,f_i,res,coeff
  integer::i
  real*8::SimpsonCoeff  !==function==!
  real*8::IntNeutrinoCyc  !==function==!
    dq_p=(q_p2-q_p1)/(nn-1)
    res=0.d0
    i=1
    do while(i.le.nn)
      q_p=q_p1+(i-1)*dq_p
      coeff=SimpsonCoeff(i,nn)
      f_i=IntNeutrinoCyc(n_i,n_f,Z_i,q_z,q_p,b,mu,T,det_part)
      res=res+coeff*f_i*dq_p/3
      i=i+1
    end do
    NeutrinoSyn_intQp_n=res
  return
  end function NeutrinoSyn_intQp_n

end subroutine NeutrinoSyn



!=====================================================================================================================
! The function which is used in integration to get the intencity if neutrino emission due to the cyclotrone mechanism.
! det_part.eq.1 - electron, det_part.ne.1 - positron.
!=====================================================================================================================
real*8 function IntNeutrinoCyc(n_i,n_f,Z_i,q_z,q_p,b,mu,T,det_part)
implicit none
integer,intent(in)::n_i,n_f,det_part
real*8,intent(in)::Z_i,q_z,q_p,b,mu,T
real*8::E_i,E_f,Z_f,omega
real*8::res
  if(q_p.ne.0.d0)then
    E_i=sqrt(1.d0+Z_i**2+2*n_i*b)
    Z_f=Z_i-q_z
    E_f=sqrt(1.d0+Z_f**2+2*n_f*b)
    omega=E_i-E_f
    if(det_part.eq.1)then
      res=q_p*omega*A(b,n_i,n_f,Z_i,Z_f,E_i,E_f,q_z,q_p,omega)*fe(E_i,mu,T)*(1.d0-fe(E_f,mu,T))
    else
      res=q_p*omega*A(b,n_f,n_i,-Z_f,-Z_i,-E_f,-E_i,q_z,q_p,omega)*fp(E_i,mu,T)*(1.d0-fp(E_f,mu,T))
    end if
  else
    res=0.d0
  end if
  IntNeutrinoCyc=res
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

  !==function A given by (17) in Kaminker 1992==!
  real*8 function A(b,n_i,n_f,Z_i,Z_f,E_i,E_f,q_z,q_p,omega)
  implicit none
  real*8,intent(in)::b,Z_i,Z_f,E_i,E_f,q_z,q_p,omega
  integer,intent(in)::n_i,n_f
  real*8::funPhi,funPsi   !==functions==!
  real*8::u,C_V_e,C_A_e,C_V_o,C_A_o,A1,A2,A3,Psi_plus,Psi_min,Phi_plus,Phi_min,p_p_i2,p_p_f2
    C_V_e=2*0.223d0+0.5d0
    C_A_e=0.5d0
    C_V_o=2*0.223d0-0.5d0
    C_A_o=-0.5d0
    u=q_p**2/2/b
    Phi_plus=funPhi(n_f,n_i,u,1)
    Phi_min =funPhi(n_f,n_i,u,-1)
    Psi_plus=funPsi(n_f,n_i,u,1)
    Psi_min =funPsi(n_f,n_i,u,-1)
    p_p_i2=2*b*n_i
    p_p_f2=2*b*n_f

    A1=Psi_plus*( -2*(E_i*E_f-Z_i*Z_f)**2 +(E_i*E_f-Z_i*Z_f)*(2.d0+p_p_i2+p_p_f2-2*q_p**2)&
                 -q_p**2/2*(q_p**2-p_p_i2-p_p_f2-3.d0) )
    A1=A1+Phi_plus*( (E_i*E_f-Z_i*Z_f)*(1.d0+p_p_i2+p_p_f2) -1.d0 -0.5d0*(p_p_i2+p_p_f2)*(3.d0-q_p**2+p_p_i2+p_p_f2) )
    A1=A1*(C_V_e**2+C_A_e**2)/E_i/E_f +2*A1*(C_V_o**2+C_A_o**2)/E_i/E_f

    A2=( Psi_plus*(-omega**2+q_z**2+q_p**2/2) +Phi_plus*(-omega**2/2 + q_z**2/2 + q_p**2) )
    A2=A2*(C_V_e**2-C_A_e**2)/E_i/E_f +2*A2*(C_V_o**2-C_A_o**2)/E_i/E_f

    A3=( (2*omega**2-2*q_z**2-3*q_p**2)*Psi_min +(p_p_i2-p_p_f2)*Phi_min )
    A3=A3*C_V_e*C_A_e/E_i/E_f*(E_i*Z_f-E_f*Z_i) +2*A3*C_V_o*C_A_o/E_i/E_f*(E_i*Z_f-E_f*Z_i)

    A=A1+A2-A3
  return
  end function A
end function IntNeutrinoCyc



!====================================================================================================
! Summation of the elements from array mas(n) starting from element i_1 and up to the element i_2
!====================================================================================================
real*8 function sum_array(mas,n,i_1,i_2)
implicit none
real*8,intent(in)::mas(n)
integer,intent(in)::n,i_1,i_2
integer::i
  sum_array=0.d0
  if((i_1.ge.1).or.(i_2.ge.1).or.(i_2.gt.i_1))then
    i=i_1
    do while(i.le.i_2)
      sum_array=sum_array+mas(i)
      i=i+1
    end do
  end if
return
end function sum_array

