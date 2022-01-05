objects =	./obj/c.o \
		./obj/neutrino.o \
		./obj/electron_in_BF.o \
		./obj/neutrino_Ann.o \
		./obj/neutrino_Syn.o #\
		#./obj/ee_pairs.o #\
		#./obj/pressure.o

objects_math_fnlib = ./obj/mf_fnlib_ai.o \
		./obj/mf_fnlib_ai.o \
		./obj/mf_fnlib_aie.o \
		./obj/mf_fnlib_csevl.o \
		./obj/mf_fnlib_inits.o \
		./obj/mf_fnlib_r9aimp.o 

objects_math =	./obj/mf_algebra_polinom.o \
		./obj/mf_bessel.o \
		./obj/mf_c.o \
		./obj/mf_coord_sys.o \
		./obj/mf_equations.o \
		./obj/mf_ermit_pol.o \
		./obj/mf_erm_link.o \
		./obj/mf_fact.o \
		./obj/mf_fi_function.o \
		./obj/mf_for_mass.o \
		./obj/mf_gamma-function.o \
		./obj/mf_gener_fun.o \
		./obj/mf_geom_2d.o \
		./obj/mf_geom_3d.o \
		./obj/mf_g.o \
		./obj/mf_integrals.o \
		./obj/mf_integral_2d.o \
		./obj/mf_int_gauss_100.o \
		./obj/mf_int_gauss.o \
		./obj/mf_intGauss.o \
		./obj/mf_intGauss_p.o \
		./obj/mf_int_G.o \
		./obj/mf_intSimps.o \
		./obj/mf_intSimp_p.o \
		./obj/mf_interpol.o \
		./obj/mf_ksi.o \
		./obj/mf_laguer_exp.o \
		./obj/mf_laguer.o \
		./obj/mf_laguer_gen.o \
		./obj/mf_laguer_n.o \
		./obj/mf_matrix_matrix.o \
		./obj/mf_numbers.o \
		./obj/mf_random.o \
		./obj/mf_vector_sk_mult.o \
		./obj/lib_specfun.o \
		./obj/mf_intExp.o


objects_phys =	./obj/PhFun_eMaxvell.o \
		./obj/PhFun_ComptonBF_cons_laws.o \
		./obj/PhFun_blackbody.o \
		./obj/PhFun_ee_pairs.o 


objects_astro =	./obj/af_ref_frame.o



c : $(objects) $(objects_math) $(objects_phys)  $(objects_astro)
	gfortran -fopenmp -o c $(objects) $(objects_math) $(objects_phys)  $(objects_astro)


./obj/c.o : c.f90	
	gfortran -c -o ./obj/c.o c.f90
./obj/neutrino.o : neutrino.f90	
	gfortran -c -o ./obj/neutrino.o neutrino.f90
./obj/electron_in_BF.o : electron_in_BF.f90	
	gfortran -c -o ./obj/electron_in_BF.o electron_in_BF.f90
./obj/neutrino_Ann.o : neutrino_Ann.f90	
	gfortran -fopenmp -c -o ./obj/neutrino_Ann.o neutrino_Ann.f90
./obj/neutrino_Syn.o : neutrino_Syn.f90	
	gfortran -c -o ./obj/neutrino_Syn.o neutrino_Syn.f90
./obj/ee_pairs.o : ee_pairs.f90	
	gfortran -c -o ./obj/ee_pairs.o ee_pairs.f90
./obj/pressure.o : pressure.f90	
	gfortran -c -o ./obj/pressure.o pressure.f90


./obj/mf_fnlib_ai.o : mf_fnlib_ai.f 	
	gfortran -c -o ./obj/mf_fnlib_ai.o mf_fnlib_ai.f
./obj/mf_fnlib_aie.o : mf_fnlib_aie.f 	
	gfortran -c -o ./obj/mf_fnlib_aie.o mf_fnlib_aie.f
./obj/mf_fnlib_csevl.o : mf_fnlib_csevl.f 	
	gfortran -c -o ./obj/mf_fnlib_csevl.o mf_fnlib_csevl.f
./obj/mf_fnlib_inits.o : mf_fnlib_inits.f 	
	gfortran -c -o ./obj/mf_fnlib_inits.o mf_fnlib_inits.f
./obj/mf_fnlib_r9aimp.o : mf_fnlib_r9aimp.f 	
	gfortran -c -o ./obj/mf_fnlib_r9aimp.o mf_fnlib_r9aimp.f

./obj/af_ref_frame.o : ../AstroFun/af_ref_frame.f90	
	gfortran -c -o ./obj/af_ref_frame.o ../AstroFun/af_ref_frame.f90


./obj/PhFun_ee_pairs.o : ../PhFun/PhFun_ee_pairs.f90	
	gfortran -c -o ./obj/PhFun_ee_pairs.o ../PhFun/PhFun_ee_pairs.f90
./obj/PhFun_eMaxvell.o : ../PhFun/PhFun_eMaxvell.f90	
	gfortran -c -o ./obj/PhFun_eMaxvell.o ../PhFun/PhFun_eMaxvell.f90
./obj/PhFun_ComptonBF_cons_laws.o : ../PhFun/PhFun_ComptonBF_cons_laws.f90	
	gfortran -c -o ./obj/PhFun_ComptonBF_cons_laws.o ../PhFun/PhFun_ComptonBF_cons_laws.f90
./obj/PhFun_ComptonBF_CrossSection.o : ../PhFun/PhFun_ComptonBF_CrossSection.f90	
	gfortran -c -o ./obj/PhFun_ComptonBF_CrossSection.o ../PhFun/PhFun_ComptonBF_CrossSection.f90
./obj/PhFun_Landau_Level_width.o : ../PhFun/PhFun_Landau_Level_width.f90
	gfortran -c -o ./obj/PhFun_Landau_Level_width.o ../PhFun/PhFun_Landau_Level_width.f90
./obj/PhFun_Cyclotron_absorb.o : ../PhFun/PhFun_Cyclotron_absorb.f90
	gfortran -c -o ./obj/PhFun_Cyclotron_absorb.o ../PhFun/PhFun_Cyclotron_absorb.f90
./obj/PhFun_magnetic_absorb.o : ../PhFun/PhFun_magnetic_absorb.f90
	gfortran -c -o ./obj/PhFun_magnetic_absorb.o ../PhFun/PhFun_magnetic_absorb.f90
./obj/PhFun_ComptonBF_Stokes.o : ../PhFun/PhFun_ComptonBF_Stokes.f90
	gfortran -c -o ./obj/PhFun_ComptonBF_Stokes.o ../PhFun/PhFun_ComptonBF_Stokes.f90
./obj/PhFun_Gaunt_factors.o : ../PhFun/PhFun_Gaunt_factors.f90
	gfortran -c -o ./obj/PhFun_Gaunt_factors.o ../PhFun/PhFun_Gaunt_factors.f90
./obj/PhFun_PlasmaDispFun.o : ../PhFun/PhFun_PlasmaDispFun.f
	gfortran -c -o ./obj/PhFun_PlasmaDispFun.o ../PhFun/PhFun_PlasmaDispFun.f
./obj/PhFun_MagneticDipole.o : ../PhFun/PhFun_MagneticDipole.f90
	gfortran -c -o ./obj/PhFun_MagneticDipole.o ../PhFun/PhFun_MagneticDipole.f90
./obj/PhFun_GR_app.o : ../PhFun/PhFun_GR_app.f90
	gfortran -c -o ./obj/PhFun_GR_app.o ../PhFun/PhFun_GR_app.f90
./obj/PhFun_Polarization.o : ../PhFun/PhFun_Polarization.f90
	gfortran -c -o ./obj/PhFun_Polarization.o ../PhFun/PhFun_Polarization.f90
./obj/PhFun_dielec_tensor.o : ../PhFun/PhFun_dielec_tensor.f90
	gfortran -c -o ./obj/PhFun_dielec_tensor.o ../PhFun/PhFun_dielec_tensor.f90
./obj/PhFun_blackbody.o : ../PhFun/PhFun_blackbody.f90
	gfortran -c -o ./obj/PhFun_blackbody.o ../PhFun/PhFun_blackbody.f90



	
./obj/lib_specfun.o : ../mf/lib_specfun.f90	
	gfortran -c -o ./obj/lib_specfun.o ../mf/lib_specfun.f90
./obj/mf_DIN_djmuz.o : ../mf/mf_DIN_djmuz.f	
	gfortran -c -o ./obj/mf_DIN_djmuz.o ../mf/mf_DIN_djmuz.f
./obj/mf_DIN_besselkn.o : ../mf/mf_DIN_besselkn.f	
	gfortran -c -o ./obj/mf_DIN_besselkn.o ../mf/mf_DIN_besselkn.f
./obj/mf_tcrUL_VF.o : ../mf/mf_tcrUL_VF.f	
	gfortran -c -o ./obj/mf_tcrUL_VF.o ../mf/mf_tcrUL_VF.f
./obj/mf_expmag_VF.o : ../mf/mf_expmag_VF.f	
	gfortran -c -o ./obj/mf_expmag_VF.o ../mf/mf_expmag_VF.f
./obj/mf_bessel_VF.o : ../mf/mf_bessel_VF.f	
	gfortran -c -o ./obj/mf_bessel_VF.o ../mf/mf_bessel_VF.f
./obj/mf_bessel.o : ../mf/mf_bessel.f90	
	gfortran -c -o ./obj/mf_bessel.o ../mf/mf_bessel.f90
./obj/mf_c.o : ../mf/mf_c.f90	
	gfortran -c -o ./obj/mf_c.o ../mf/mf_c.f90
./obj/mf_ermit_pol.o : ../mf/mf_ermit_pol.f90	
	gfortran -c -o ./obj/mf_ermit_pol.o ../mf/mf_ermit_pol.f90
./obj/mf_erm_link.o : ../mf/mf_erm_link.f90	
	gfortran -c -o ./obj/mf_erm_link.o ../mf/mf_erm_link.f90
./obj/mf_fact.o : ../mf/mf_fact.f90	
	gfortran -c -o ./obj/mf_fact.o ../mf/mf_fact.f90
./obj/mf_fi_function.o : ../mf/mf_fi_function.f90	
	gfortran -c -o ./obj/mf_fi_function.o ../mf/mf_fi_function.f90
./obj/mf_gamma-function.o : ../mf/mf_gamma-function.f90	
	gfortran -c -o ./obj/mf_gamma-function.o ../mf/mf_gamma-function.f90
./obj/mf_gener_fun.o : ../mf/mf_gener_fun.f90	
	gfortran -c -o ./obj/mf_gener_fun.o ../mf/mf_gener_fun.f90
./obj/mf_geom_2d.o : ../mf/mf_geom_2d.f90	
	gfortran -c -o ./obj/mf_geom_2d.o ../mf/mf_geom_2d.f90
./obj/mf_geom_3d.o : ../mf/mf_geom_3d.f90	
	gfortran -c -o ./obj/mf_geom_3d.o ../mf/mf_geom_3d.f90
./obj/mf_g.o : ../mf/mf_g.f90	
	gfortran -c -o ./obj/mf_g.o ../mf/mf_g.f90
./obj/mf_integrals.o : ../mf/mf_integrals.f90	
	gfortran -c -o ./obj/mf_integrals.o ../mf/mf_integrals.f90
./obj/mf_integral_2d.o : ../mf/mf_integral_2d.f90	
	gfortran -c -o ./obj/mf_integral_2d.o ../mf/mf_integral_2d.f90
./obj/mf_int_gauss_100.o : ../mf/mf_int_gauss_100.f90	
	gfortran -c -o ./obj/mf_int_gauss_100.o ../mf/mf_int_gauss_100.f90
./obj/mf_int_gauss.o : ../mf/mf_int_gauss.f90	
	gfortran -c -o ./obj/mf_int_gauss.o ../mf/mf_int_gauss.f90
./obj/mf_intGauss.o : ../mf/mf_intGauss.FOR
	gfortran -c -o ./obj/mf_intGauss.o ../mf/mf_intGauss.FOR
./obj/mf_intGauss_p.o : ../mf/mf_intGauss_p.FOR	
	gfortran -c -o ./obj/mf_intGauss_p.o ../mf/mf_intGauss_p.FOR
./obj/mf_int_G.o : ../mf/mf_int_G.FOR	
	gfortran -c -o ./obj/mf_int_G.o ../mf/mf_int_G.FOR
./obj/mf_intSimps.o : ../mf/mf_intSimps.f90	
	gfortran -c -o ./obj/mf_intSimps.o ../mf/mf_intSimps.f90
./obj/mf_intSimp_p.o : ../mf/mf_intSimp_p.FOR	
	gfortran -c -o ./obj/mf_intSimp_p.o ../mf/mf_intSimp_p.FOR
./obj/mf_ksi.o : ../mf/mf_ksi.f90	
	gfortran -c -o ./obj/mf_ksi.o ../mf/mf_ksi.f90
./obj/mf_laguer_exp.o : ../mf/mf_laguer_exp.f90	
	gfortran -c -o ./obj/mf_laguer_exp.o ../mf/mf_laguer_exp.f90
./obj/mf_laguer.o : ../mf/mf_laguer.f90	
	gfortran -c -o ./obj/mf_laguer.o ../mf/mf_laguer.f90
./obj/mf_laguer_gen.o : ../mf/mf_laguer_gen.f90	
	gfortran -c -o ./obj/mf_laguer_gen.o ../mf/mf_laguer_gen.f90
./obj/mf_laguer_n.o : ../mf/mf_laguer_n.f90	
	gfortran -c -o ./obj/mf_laguer_n.o ../mf/mf_laguer_n.f90
./obj/mf_matrix_matrix.o : ../mf/mf_matrix_matrix.f90	
	gfortran -c -o ./obj/mf_matrix_matrix.o ../mf/mf_matrix_matrix.f90
./obj/mf_num_mart.o : ../mf/mf_num_mart.f90	
	gfortran -c -o ./obj/mf_num_mart.o ../mf/mf_num_mart.f90
./obj/mf_numbers.o : ../mf/mf_numbers.f90	
	gfortran -c -o ./obj/mf_numbers.o ../mf/mf_numbers.f90
./obj/mf_vector_sk_mult.o : ../mf/mf_vector_sk_mult.f90	
	gfortran -c -o ./obj/mf_vector_sk_mult.o ../mf/mf_vector_sk_mult.f90
./obj/mf_coord_sys.o : ../mf/mf_coord_sys.f90	
	gfortran -c -o ./obj/mf_coord_sys.o ../mf/mf_coord_sys.f90
./obj/mf_algebra_polinom.o : ../mf/mf_algebra_polinom.f90	
	gfortran -c -o ./obj/mf_algebra_polinom.o ../mf/mf_algebra_polinom.f90
./obj/mf_interpol.o : ../mf/mf_interpol.f90	
	gfortran -c -o ./obj/mf_interpol.o ../mf/mf_interpol.f90
./obj/mf_for_mass.o : ../mf/mf_for_mass.f90	
	gfortran -c -o ./obj/mf_for_mass.o ../mf/mf_for_mass.f90
./obj/mf_random.o : ../mf/mf_random.f90	
	gfortran -c -o ./obj/mf_random.o ../mf/mf_random.f90
./obj/mf_noise.o : ../mf/mf_noise.f90	
	gfortran -c -o ./obj/mf_noise.o ../mf/mf_noise.f90
./obj/mf_equations.o : ../mf/mf_equations.f90	
	gfortran -c -o ./obj/mf_equations.o ../mf/mf_equations.f90
./obj/mf_intExp.o : ../mf/mf_intExp.f90	
	gfortran -c -o ./obj/mf_intExp.o ../mf/mf_intExp.f90




