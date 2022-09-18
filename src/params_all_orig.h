// all the parameters (type, name ,default_value) are defined here
#define PARAMS_ALL()\
PARAM(int, random_seed, 1)\
PARAM(int, number_of_events, 10000)\
PARAM(int, number_of_test_events,10000)\
PARAM(int, user_events,0)\
PARAM(int, beam_type, 0)\
PARAM(line, beam_energy,"1000")\
PARAM(int, beam_particle, 14)\
PARAM(vec, beam_direction, vec(0,0,1))\
/* beam type=1 */ \
PARAM(line, beam_content,"")\
/* beam type=2 */ \
PARAM(string, beam_folder, "flux")\
PARAM(int, beam_file_first, 1)\
PARAM(int, beam_file_limit, 0)\
PARAM(bool, beam_weighted, 0)\
/* beam type=3 */ \
PARAM(string, beam_file, "beam/ND280hist.txt")\
PARAM(vec, beam_offset, vec(0,0,0))\
PARAM(int, beam_placement, 0)\
PARAM(int, beam_test_only, 0)\
PARAM(int, nucleus_p, 0)\
PARAM(int, nucleus_n, 0)\
PARAM(double,nucleus_E_b,0)\
PARAM(double,nucleus_kf,0)\
PARAM(int, nucleus_target,0)\
PARAM(int, nucleus_model,0)\
PARAM(int, target_type,0)\
PARAM(line, target_content,"")\
PARAM(string, geo_file,"")\
PARAM(string, geo_name,"ND280Geometry")\
PARAM(string, geo_volume,"")\
PARAM(vec, geo_o,vec(0,0,0))\
PARAM(vec, geo_d,vec(0,0,0))\
PARAM(bool, dyn_qel_cc,0)\
PARAM(bool, dyn_qel_nc,0)\
PARAM(bool, dyn_res_cc,0)\
PARAM(bool, dyn_res_nc,0)\
PARAM(bool, dyn_dis_cc,0)\
PARAM(bool, dyn_dis_nc,0)\
PARAM(bool, dyn_coh_cc,0)\
PARAM(bool, dyn_coh_nc,0)\
/*       PARAM(int, qel_kinematics,0)*/ \
PARAM(int, qel_vector_ff_set, 2)\
PARAM(int, qel_axial_ff_set, 1)\
PARAM(double, qel_cc_vector_mass, 840)\
PARAM(double, qel_cc_axial_mass,1030)\
PARAM(double, qel_nc_axial_mass,1030)\
PARAM(double, qel_s_axial_mass,1030)\
PARAM(int, qel_strange, 0)\
PARAM(int, qel_strangeEM, 0)\
PARAM(double, delta_s, 0)\
PARAM(bool, flux_correction, 1)\
PARAM(int, delta_FF_set,1)\
PARAM(double, pion_axial_mass, 1)\
PARAM(double, pion_C5A, 1.2)\
PARAM(double, res_dis_cut,0)\
PARAM(int, spp_precision,0)\
/*      PARAM(bool, qel_new,0)\
PARAM(bool, qel_cc_japan,0)\
PARAM(bool, qel_relati,0)\
*/      PARAM(bool, coh_mass_correction,0)\
PARAM(bool, coh_new,1)\
PARAM(bool, kaskada_on,1)\
PARAM(bool, kaskada_newangle,1)\
PARAM(bool, kaskada_redo,0)\
PARAM(bool, pauli_blocking,1)\
PARAM(string, formation_zone, "nofz")\
PARAM(double, tau, 0)\
PARAM(bool, first_step, 1)\
PARAM(double, step, 0.2)\
PARAM(int, xsec, 0)\
PARAM(int, sf_method,0)\
PARAM(bool, cc_smoothing,0)\
PARAM(bool, mixed_order,1)\
PARAM(double,rmin,0)\
PARAM(double,rmax,0)\
PARAM(double,res_dis_blending_start,2800)\
PARAM(double,res_dis_blending_end,3200)

