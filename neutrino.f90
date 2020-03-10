


!==================================================================================
! ro - mass density
!==================================================================================
subroutine NeutrinoGammaE(T10,B12,ro,Q23)
implicit none
real*8,intent(in)::T10,B12,ro
real*8::Q23
real*8::Bcr12,n_e23,n30,TkeV,p_Fc_keV
  Bcr12=44.13d0
  n_e23=6.d0*ro  !==for the case of hydrogen complitely ionised plasma==!
  call EE_pairs_app2(T10/0.593d0,n30)  !==due to pair crestion==!
  n_e23=max(n_e23,n30*1.e7)
  TkeV=T10*1.d10/1.14d4/1.d3   !==check==!
  p_Fc_keV=2.57d-3*n_e23/B12
  Q23=1.3d-4*(B12/Bcr12)*(TkeV/511.d0)**8*(TkeV/p_Fc_keV)
return
end subroutine NeutrinoGammaE








