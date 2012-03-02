#ifndef _data_h_
#define _data_h_
//Atmospheric Neutrino - data from Mereneyi: 0 - 0pi, 1 - pi+, 2 - pi0, 3 - pi-pi+, 4 - 2pi0, 5 - pi+ n*pi0 n>0, 6 - 2pi+ n*pi0 n>=0, 7 - pi- (CC interactions)

const double AtmosphericDeuterium[8] = {0.49, 0.33, 0.09, 0.02, 0.01, 0.01, 0.01, 0};
const double AtmosphericDeuteriumerr[8] = {0.03, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0};

const double AtmosphericNeon[8] = {0.57, 0.22, 0.05, 0.03, 0.04, 0.02, 0.02, 0.01};
const double AtmosphericNeonerr[8] = {0.04, 0.02, 0.02, 0.01, 0.02, 0.01, 0.01, 0.01};

//MiniBooNE data - cross section for 1pi0 production in NC reactions on CH2 as a function of kinetic energy and angle

	//neutrino

const double MBMomxOLD2[11] = {50, 125, 175, 225, 275, 350, 450, 550, 700, 900, 1250 };
const double MBMomyOLD2[11] = {0.17e-42, 1.16e-42, 1.58e-42, 1.53e-42, 1.25e-42, 0.83e-42, 0.44e-42, 0.19e-42, 0.04e-42, 0.027e-42, 0.0033e-42 };
const double MBMomxerrOLD2[11] = {50, 25, 25, 25, 25, 50, 50, 50, 100, 100, 250};
const double MBMomyerrOLD2[11] = {0.05e-42, 0.19e-42, 0.22e-42, 0.2e-42, 0.2e-42, 0.13e-42, 0.083e-42, 0.057e-42, 0.037e-42, 0.016e-42, 0.0032e-42};

const double MBAnglexOLD2[18]= {-0.81,-0.48, -0.23, -0.035, 0.13, 0.26, 0.37, 0.47, 0.56, 0.635, 0.7, 0.76, 0.805, 0.85, 0.89 ,0.93, 0.97,  0.99};
const double MBAngleyOLD2[18] = {0.077e-39, 0.1e-39, 0.13e-39, 0.16e-39, 0.19e-39, 0.215e-39, 0.25e-39, 0.27e-39, 0.3e-39, 0.35e-39, 0.37e-39, 0.42e-39, 0.47e-39, 0.52e-39, 0.6e-39, 0.69e-39, 0.8e-39, 0.9e-39};
const double MBAnglexerrOLD2[18] = {0.18, 0.14, 0.11, 0.095, 0.07, 0.063, 0.051, 0.055, 0.04, 0.037, 0.033, 0.025, 0.022, 0.025, 0.018, 0.018,  0.014, 0.008};
const double MBAngleyerrOLD2[18] = {0.015e-39, 0.019e-39, 0.024e-39, 0.023e-39, 0.037e-39, 0.031e-39, 0.036e-39, 0.04e-39, 0.043e-39, 0.051e-39, 0.058e-39, 0.062e-39, 0.073e-39, 0.083e-39, 0.089e-39, 0.1e-39, 0.14e-39, 0.14e-39};

const double MBMomxOLD1[11] = {50, 125, 175, 225, 275, 350, 450, 550, 700, 900, 1250};
const double MBMomyOLD1[11] = {0.19e-42, 1.18e-42, 1.60e-42, 1.56e-42, 1.27e-42, 0.87e-42, 0.48e-42, 0.23e-42, 0.07e-42, 0.047e-42, 0.019e-42 };
const double MBMomxerrOLD1[11] = {50, 25, 25, 25, 25, 50, 50, 50, 100, 100, 250};
const double MBMomyerrOLD1[11] = {0.05e-42, 0.19e-42, 0.22e-42, 0.2e-42, 0.2e-42, 0.13e-42, 0.085e-42, 0.057e-42, 0.047e-42, 0.019e-42, 0.019e-42};

const double MBAnglexOLD1[18]= {-0.81,-0.48, -0.23, -0.035, 0.13, 0.26, 0.37, 0.47, 0.56, 0.635, 0.7, 0.76, 0.805, 0.85, 0.89 ,0.93, 0.97,  0.99};
const double MBAngleyOLD1[18] = {0.078e-39, 0.11e-39, 0.13e-39, 0.16e-39, 0.194e-39, 0.223e-39, 0.25e-39, 0.28e-39, 0.31e-39, 0.37e-39, 0.39e-39, 0.43e-39, 0.49e-39, 0.55e-39, 0.63e-39, 0.72e-39, 0.84e-39, 0.95e-39};
const double MBAnglexerrOLD1[18] = {0.19, 0.14, 0.1, 0.095, 0.07, 0.06, 0.043, 0.03, 0.042, 0.033, 0.029, 0.021, 0.024, 0.017, 0.021, 0.017,  0.013, 0.017};
const double MBAngleyerrOLD1[18] = {0.023e-39, 0.023e-39, 0.029e-39, 0.029e-39, 0.046e-39, 0.029e-39, 0.046e-39, 0.04e-39, 0.052e-39, 0.058e-39, 0.063e-39, 0.064e-39, 0.082e-39, 0.087e-39, 0.104e-39, 0.122e-39, 0.145e-39, 0.151e-39};

const double MBMomx[11] = {50, 125, 175, 225, 275, 350, 450, 550, 700, 900, 1250};
const double MBMomxerr[11] = {50, 25, 25, 25, 25, 50, 50, 50, 100, 100, 250};
const double MBMomy[11] = {0.18e-42, 1.19e-42, 1.63e-42, 1.58e-42, 1.28e-42, 0.87e-42, 0.47e-42, 0.21e-42, 0.05e-42, 0.03e-42, 0.01e-42};
const double MBMomyerr[11] = {0.06e-42, 0.21e-42, 0.24e-42, 0.21e-42, 0.21e-42, 0.14e-42, 0.09e-42, 0.06e-42, 0.04e-42, 0.02e-42, 0.01e-42};

const double MBAnglex[18] = {-0.81, -0.48, -0.235, -0.035, 0.13, 0.26, 0.37, 0.47, 0.56, 0.635, 0.7, 0.755, 0.805, 0.85, 0.89, 0.93, 0.9625, 0.9875};
const double MBAnglexerr[18] = {0.19, 0.14, 0.105, 0.095, 0.07, 0.06, 0.05, 0.05, 0.04, 0.035, 0.03, 0.025, 0.025, 0.02, 0.02, 0.02, 0.0125, 0.0125};
const double MBAngley[18] = {0.82e-40, 1.1e-40, 1.37e-40, 1.64e-40, 1.99e-40, 2.26e-40, 2.58e-40, 2.82e-40, 3.16e-40, 3.68e-40, 3.94e-40, 4.38e-40, 4.96e-40, 5.49e-40, 6.33e-40, 7.28e-40, 8.42e-40, 9.56e-40};
const double MBAngleyerr[18] = {0.16e-40, 0.22e-40, 0.24e-40, 0.27e-40, 0.43e-40, 0.35e-40, 0.43e-40, 0.46e-40, 0.50e-40, 0.63e-40, 0.65e-40, 0.70e-40, 0.77e-40, 0.91e-40, 1.04e-40, 1.18e-40, 1.45e-40, 1.58e-40};
	
	//anti-neutrino

const double MBantiMomxOLD2[10] = {65, 150, 195, 230, 260, 300, 345, 405, 505, 835};
const double MBantiMomyOLD2[10] = {0.108e-42, 0.49e-42, 0.57e-42, 0.53e-42, 0.46e-42, 0.35e-42, 0.25e-42, 0.15e-42, 0.068e-42, 0.011e-42};
const double MBantiMomxerrOLD2[10] = {61, 23, 19, 11, 15, 19, 27, 30, 61, 262};
const double MBantiMomyerrOLD2[10] = {0.023e-42, 0.083e-42, 0.081e-42, 0.074e-42, 0.061e-42, 0.052e-42, 0.063e-42, 0.063e-42, 0.023e-42, 0.0068e-42};

const double MBantiAnglexOLD2[10]= {-0.804, -0.407, -0.051, 0.26, 0.50, 0.67, 0.80, 0.89, 0.95, 0.99};
const double MBantiAngleyOLD2[10] = {0.037e-39, 0.046e-39, 0.052e-39, 0.045e-39, 0.057e-39, 0.091e-39, 0.13e-39, 0.19e-39, 0.29e-39, 0.39e-39};
const double MBantiAnglexerrOLD2[10] = {0.19, 0.19, 0.17, 0.14, 0.09, 0.07, 0.055, 0.028, 0.028, 0.028};
const double MBantiAngleyerrOLD2[10] = {0.0077e-39, 0.0077e-39, 0.0077e-39, 0.0108e-39, 0.014e-39, 0.015e-39, 0.02e-39, 0.029e-39, 0.040e-39, 0.057e-39};

const double MBantiMomxOLD1[10] = {65, 150, 195, 230, 260, 300, 345, 405, 505, 835};
const double MBantiMomyOLD1[10] = {0.096e-42, 0.47e-42, 0.55e-42, 0.52e-42, 0.45e-42, 0.35e-42, 0.26e-42, 0.17e-42, 0.074e-42, 0.008e-42};
const double MBantiMomxerrOLD1[10] = {61, 23, 18, 16, 18, 20, 25, 34, 65, 178};
const double MBantiMomyerrOLD1[10] = {0.020e-42, 0.079e-42, 0.078e-42, 0.076e-42, 0.066e-42, 0.052e-42, 0.050e-42, 0.030e-42, 0.020e-42, 0.007e-42};

const double MBantiAnglexOLD1[10]= {-0.80, -0.41, -0.05, 0.26, 0.50, 0.67, 0.79, 0.88, 0.93, 0.98};
const double MBantiAngleyOLD1[10] = {0.031e-39, 0.038e-39, 0.045e-39, 0.047e-39, 0.065e-39, 0.097e-39, 0.14e-39, 0.20e-39, 0.29e-39, 0.37e-39};
const double MBantiAnglexerrOLD1[10] = {0.2, 0.19, 0.17, 0.14, 0.10, 0.07, 0.055, 0.033, 0.028, 0.02};
const double MBantiAngleyerrOLD1[10] = {0.0045e-39, 0.009e-39, 0.009e-39, 0.009e-39, 0.012e-39, 0.016e-39, 0.026e-39, 0.036e-39, 0.044e-39, 0.06e-39};

const double MBantiMomx[10] = {65, 150, 190, 225, 260, 300, 345, 405, 505, 835};
const double MBantiMomxerr[10] = {65, 20, 20, 15, 20, 20, 25, 35, 65, 265};
const double MBantiMomy[10] = {1.13e-43, 5.20e-43, 5.86e-43, 5.26e-43, 4.42e-43, 3.68e-43, 2.84e-43, 1.72e-43, 0.71e-43, 0.11e-43};
const double MBantiMomyerr[10] = {0.25e-43, 0.86e-43, 0.86e-43, 0.78e-43, 0.64e-43, 0.55e-43, 0.49e-43, 0.36e-43, 0.19e-43, 0.06e-43};

const double MBantiAnglex[10] = {-0.8, -0.41, -0.05, 0.26, 0.5, 0.67, 0.795, 0.88, 0.935, 0.98};
const double MBantiAnglexerr[10] = {0.2, 0.19, 0.17, 0.14, 0.1, 0.07, 0.055, 0.03, 0.025, 0.02};
const double MBantiAngley[10] = {0.38e-40, 0.40e-40, 0.50e-40, 0.61e-40, 0.69e-40, 1.00e-40, 1.33e-40, 1.94e-40, 2.76e-40, 4.06e-40};
const double MBantiAngleyerr[10] = {0.08e-40, 0.08e-40, 0.10e-40, 0.12e-40, 0.15e-40, 0.20e-40, 0.27e-40, 0.38e-40, 0.50e-40, 0.74e-40};

	//neutrino with antineutrino contamination
	
const double MBMomContx[11]    = {50, 125, 175, 225, 275, 350, 450, 550, 700, 900, 1250};
const double MBMomContxerr[11] = {50, 25, 25, 25, 25, 50, 50, 50, 100, 100, 250};
const double MBMomConty[11]    = {0.17e-42, 1.14e-42, 1.56e-42, 1.51e-42, 1.23e-42, 0.83e-42, 0.45e-42, 0.201e-42, 0.051e-42, 0.033e-42, 0.0058e-42};
const double MBMomContyerr[11] = {0.04e-42, 0.19e-42, 0.21e-42, 0.19e-42, 0.18e-42, 0.13e-42, 0.087e-42, 0.047e-42, 0.031e-42, 0.023e-42, 0};

	
	//antineutrino with neutrino containation
	
const double MBantiMomContx[10]    = {65, 150, 190, 225, 260, 300, 345, 405, 505, 835};
const double MBantiMomContxerr[10] = {65, 20, 20, 15, 20, 20, 25, 35, 65, 265};
const double MBantiMomConty[10]    = {1.260e-43, 5.758e-43, 6.536e-43, 6.005e-43, 5.074e-43, 4.255e-43, 3.371e-43, 2.090e-43, 0.945e-43, 0.178e-43};
const double MBantiMomContyerr[10] = {0.23e-43, 0.81e-43, 0.81e-43, 0.74e-43, 0.58e-43, 0.48e-43, 0.45e-43, 0.36e-43, 0.23e-43, 0.064e-43};

	//cross section
	
const double MBcross = 4.73;
const double MBcrosserr = 0.05;
const double MBcrosserr2 = 0.40;

const double MBcrossanti = 1.42;
const double MBcrossantierr = 0.04;
const double MBcrossantierr2 = 0.14;

	//cross section with contamination
	
const double MBcrossCont = 4.56;
const double MBcrosserrCont = 0.05;
const double MBcrosserr2Cont = 0.38;

const double MBcrossantiCont = 1.75;
const double MBcrossantierrCont = 0.04;
const double MBcrossantierr2Cont = 0.14;

	//CCpi0 total
	
const double MBCCnuen[14]    = {0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.9};
const double MBCCnuenerr[14] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1};

const double MBCCtotal[14]    = {1.76e-39, 3.83e-39, 5.68e-39, 7.31e-39, 9.20e-39, 11.06e-39, 12.42e-39, 13.89e-39, 15.23e-39, 16.38e-39, 18.20e-39, 19.37e-39, 20.80e-39, 21.92e-39};
const double MBCCtotalerr[14] = {0.49e-39, 0.78e-39, 0.97e-39, 1.17e-39, 1.50e-39, 1.85e-39, 2.16e-39, 2.46e-39, 2.82e-39, 3.17e-39, 3.61e-39, 4.15e-39, 4.49e-39, 5.26e-39};

const double MBCCtotalnuance[15] = {7.33e-40, 7.33e-40, 1.83e-39, 3.16e-39, 4.43e-39, 5.76e-39, 7.06e-39, 8.08e-39, 9.22e-39, 1.04e-38, 1.17e-38, 1.31e-38, 1.42e-38, 1.54e-38, 1.68e-38};

	//CCpi0 momentum
	
const double MBCCmom[11]    = {0.05, 0.125, 0.175, 0.225, 0.275, 0.35, 0.45, 0.55, 0.70, 0.90, 1.20};
const double MBCCmomerr[11] = {0.05, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 0.05, 0.10, 0.10, 0.20};

const double MBCCxsec[11]    = {4.92e-39, 26.65e-39, 32.90e-39, 28.99e-39, 19.02e-39, 13.65e-39, 7.41e-39, 4.27e-39, 1.90e-39, 0.87e-39, 0.19e-39};
const double MBCCxsecerr[11] = {1.34e-39, 4.94e-39, 5.00e-39, 4.31e-39, 3.09e-39, 2.49e-39, 2.01e-39, 1.14e-39, 0.63e-39, 0.40e-39, 0.79e-39};

const double MBCCxsecnuance[11] = {2.39264171893e-39, 1.12688576316e-38, 1.54838785738e-38, 1.59580255594e-38, 1.28534559999e-38, 1.0935004137e-38, 6.94826824887e-39, 3.94144327599e-39, 1.87594582892e-39, 5.81473113808e-40, 1.39757610574e-41};


//pi0 Muliplicity as a function of anti-neutrino energy (CC interactions)- generator comparisons (Gallagher Nuclear Physics B (Proc. Suppl.) 139 (2005) 278–285)

	//free proton

const double NeutProton[5] = {0.04, 0.07, 0.17, 0.38, 0.57};
const double NuanceProton[5] = {0.04, 0.07, 0.19, 0.41, 0.68};
const double NuxProton[5] = {0.09, 0.16, 0.31, 0.52, 0.76};
const double NeugenProton[5] = {0.05, 0.10, 0.26, 0.55, 0.88};

	//proton bound in oxygen

const double NeutOxygen[5] = {0.02, 0.04, 0.14, 0.33, 0.54};
const double NuanceOxygen[5] = {0.03, 0.06, 0.19, 0.44, 0.69};
const double NuxOxygen[5] = {0.08, 0.14, 0.26, 0.49, 0.74};
const double NeugenOxygen[5] = {0.04, 0.09, 0.23, 0.48, 0.77};

//pi0 Multiplicity as a function of W2 (CC interactions) - data from SKAT (freon), BEBC (H2 and neon)

const double SKATfreonx[13] = {1.12, 1.42, 1.83, 2.38, 3.02, 3.86, 4.99, 6.32, 8.10, 10.26, 13.34, 17.08, 21.77};
const double SKATfreony[13] = {0.21, 0.25, 0.46, 0.66, 0.78, 0.91, 1.06, 1.20, 1.29, 1.45, 1.82, 1.95, 2.37};
const double SKATfreonxerr[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const double SKATfreonyerr[13] = {0.06, 0.06, 0.06, 0.06, 0.07, 0.06, 0.08, 0.06, 0.07, 0.08, 0.10, 0.15, 0.17};

const double BEBCh2x[8] = {14.71, 18.49, 23.17, 28.92, 37.00, 45.65, 65.25, 96.71};
const double BEBCh2y[8] = {1.38, 1.53, 1.86, 1.82, 1.98, 2.12, 2.27, 2.45};
const double BEBCh2xerr[8] = {0, 0, 0, 0, 0, 0, 0, 0};
const double BEBCh2yerr[8] = {0.15, 0.16, 0.19, 0.17, 0.17, 0.19, 0.15, 0.21};

const double BEBCneonx[10] = {1.26, 1.98, 3.09, 4.90, 7.75, 12.57, 19.65, 30.46, 48.48, 76.76};
const double BEBCneony[10] = {0.22, 0.44, 0.45, 0.94, 1.10, 1.50, 1.80, 2.08, 2.43, 2.51};
const double BEBCneonxerr[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const double BEBCneonyerr[10] = {0.06, 0.15, 0.14, 0.12, 0.10, 0.09, 0.11, 0.12, 0.11, 0.14};

//SciBooNE data - cross section for pi0 production in NC reactions on C8H8 as a funcion of momentum (yerror = 15%)

const double SBxOLD[25] = {20, 60, 100, 140, 180, 220, 260, 300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 700, 740, 780, 820, 860, 900, 940, 980};
const double SByOLD[25] = {0, 13.9, 53.6, 72.6, 93.0, 68.8, 35.4, 34.9, 33.7, 29.4, 11.7, 5.5, 2.5, 5.7, 7.8, 0, 1.5, 0, 0, 0, 0.3, 0, 1.0, 0.5, 0.5};

const double SBx[10] = {40, 120, 200, 280, 360, 440, 520, 600, 680, 760};
const double SBxerr[10] = {40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
const double SBy[10] = {0.021, 0.159, 0.262, 0.153, 0.130, 0.099, 0.067, 0.036, 0.024, 0.052};
const double SByerr[10] = {0.003, 0.013, 0.018, 0.011, 0.010, 0.009, 0.009, 0.008, 0.007, 0.029};

const double SBxang[10] = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9};
const double SBxangerr[10] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
const double SByang[10] = {0.043, 0.043, 0.055, 0.058, 0.063, 0.076, 0.080, 0.110, 0.180, 0.310};
const double SByangerr[10] = {0.008, 0.007, 0.007, 0.007, 0.008, 0.007, 0.009, 0.011, 0.016, 0.029};

//SciBooNE data - NCPi0/alCC

const double SBratio = 0.070;
const double SBratioerr = 0.005;
const double SBratioerrup = 0.004;
const double SBratioerrdown = 0.005;

//K2K data - cross section for 1pi0 production in NC reactions on H2O as a function of momentum

const double K2Kx[9] = {50, 150, 250, 350, 450, 550, 650, 750, 850};
const double K2Ky[9] = {117, 817, 975, 694, 398, 225, 142, 85, 172};

const double K2Kxerr[9] = {50, 50, 50, 50, 50, 50, 50, 50, 50};
const double K2Kyerr[9] = {30, 90, 105, 87, 62, 47, 37, 25, 98};

//K2K data - 1NCPi0/allCC

const double K2Kratio = 0.064;
const double K2Kratioerr = 0.001;
const double K2Kratioerr2 = 0.007;

//NC - single pion production on free nucleon in energy = 2200 MeV (0 - pi0p, 1 - pi0n, 2 - pi+, 3 - pi-, 4 - pi0CC)

const double NCSPP[5] = {297, 177, 180, 237, 526};
const double NCSPPerr[5] = {37, 43, 31, 59, 65};

const double NCSPPnorm[5] = {0.13, 0.08, 0.08, 0.11, 0.24};
const double NCSPPnormerr[5] = {0.02, 0.02, 0.02, 0.03, 0.04};

//pi+Nucleus interaction

const double CarbonEnergyAshery[6] = {85, 125, 165, 205, 245, 315};
const double CarbonAbsAshery[6] = {109, 166, 194, 157, 95, 64};
const double CarbonAbsAsheryErr[6] = {20, 26, 36, 37, 32, 27};
const double CarbonSCXAshery[6] = {35, 38, 46, 45, 47, 45};
const double CarbonSCXAsheryErr[6] = {12, 12, 23, 23, 23, 22};
const double CarbonInelAshery[6] = {143, 213, 207, 210, 224, 200};
const double CarbonInelAsheryErr[6] = {26, 33, 33, 51, 30, 22};
const double CarbonReacAshery[6] = {287, 417, 447, 412, 366, 309};

const double CarbonEnergyNavon[1] = {50};
const double CarbonAbsNavon[1] = {88};
const double CarbonAbsNavonErr[1] = {27};

const double CarbonEnergyJones[4] = {250, 300, 400, 500};
const double CarbonAbsJones[4] = {73, 69, 45, 39};
const double CarbonAbsJonesErr[4] = {15, 14, 9, 19};
const double CarbonSCXJones[4] = {46, 55, 59, 91};

const double CarbonEnergyJones2[2] = {250, 300};
const double CarbonInelJones[2] = {230, 250};
const double CarbonReacJones[2] = {349, 374};

const double IronEnergyAshery[6] = {85, 125, 165, 205, 245, 315};
const double IronAbsAshery[6] = {421, 527, 577, 607, 411, 320};
const double IronAbsAsheryErr[6] = {70, 74, 87, 86, 70, 62};
const double IronSCXAshery[6] = {79, 83, 103, 83, 95, 98};
const double IronSCXAsheryErr[6] = {40, 41, 50, 40, 50, 50};
const double IronInelAshery[6] = {784, 644, 474, 360, 430, 389};
const double IronInelAsheryErr[6] = {115, 123, 130, 125, 110, 100};
const double IronReacAshery[6] = {1284, 1254, 1154, 1050, 936, 807};

const double IronEnergyNavon[1] = {50};
const double IronAbsNavon[1] = {443};
const double IronAbsNavonErr[1] = {58};

const double NickelEnergyJones[3] = {250, 300, 500};
const double NickelAbsJones[3] = {280, 310, 330};
const double NickelAbsJonesErr[3] = {56, 62, 100};
const double NickelSCXJones[3] = {110, 140, 290};

const double NickelEnergyJones2[2] = {250, 300};
const double NickelInelJones[2] = {590, 770};
const double NickelReacJones[2] = {980, 1220};

const double OxygenEnergyIngram[3] = {114, 163, 240};
const double OxygenAbsIngram[3] = {206, 188, 89};
const double OxygenAbsIngramErr[3] = {33, 36, 35};
const double OxygenSCXIngram[3] = {58, 63, 62};
const double OxygenSCXIngramErr[3] = {17, 19, 19};
const double OxygenInelIngram[3] = {191, 259, 249};
const double OxygenInelIngramErr[3] = {12, 17, 16};
const double OxygenReacIngram[3] = {455, 510, 406};

//Proton Transparency

	//NE-18
	
const double CarbonNE18Q2[5]   = {0.33, 1.04, 3.06, 5.00, 6.77};
const double CarbonNE18T[5]    = {0.77, 0.64, 0.63, 0.61, 0.67};
const double CarbonNE18Terr[5] = {0.04, 0.05, 0.06, 0.06, 0.07};

const double IronNE18Q2[5]    = {0.33, 1.04, 3.06, 5.00, 6.77};
const double IronNE18T[5]     = {0.54, 0.50, 0.39, 0.40, 0.43};
const double IronNE18Terr[5]  = {0.03, 0.05, 0.04, 0.05, 0.05};

	//Bates

const double CarbonBatesQ2[4] = {0.66, 1.30, 1.80, 3.26};
const double CarbonBatesT[4]  = {0.61, 0.60, 0.58, 0.58};

const double IronBatesQ2[4] = {0.66, 1.30, 1.80, 3.26};
const double IronBatesT[4]  = {0.46, 0.43, 0.39, 0.41};

//error = 3.2%

//NOMAD (backwards protons and pions)

const double NOMADprotonsQ2[8]      = {0.5, 1.5, 3.0, 5.5, 11.0, 22.5, 40.0, 75.0};
const double NOMADprotonsQ2err[8]   = {0.5, 0.5, 1.0, 1.5, 4.0, 7.5, 10.0, 25.0};
const double NOMADprotonsBack[8]    = {0.055, 0.045, 0.040, 0.038, 0.037, 0.034, 0.028, 0.022};
const double NOMADprotonsBackerr[8] = {0.0017, 0.0010, 0.0010, 0.0013, 0.0014, 0.0013, 0.0018, 0.0};

const double NOMADpionsQ2[6]        = {1.0, 3.0, 7.0, 17.5, 42.5, 80.0};
const double NOMADpionsQ2err[6]        = {1.0, 1.0, 3.0, 7.5, 17.5, 20.0};
const double NOMADpionsBack[6]      = {0.0051, 0.0048, 0.0043, 0.0038, 0.0037, 0.0023};
const double NOMADpionsBackerr[6]   = {0.0002, 0.0002, 0.0002, 0.0002, 0.0005, 0.0009};

const int NOMADmultiProton[4] = {904212, 37634, 2168, 5};
const int NOMADmultiPion[4] = {939617, 4238, 164, 0};

//Pions Transparency

const double PTQ2[5]    = {1.1, 2.1, 3.0, 3.9, 4.7};
const double PTQ2err[5] = {0.08, 0.12, 0.11, 0.15, 0.18};

const double PTppi[5]   = {2.793, 3.187, 3.418, 4.077, 4.412};

const double PTC[5]     = {0.67, 0.65, 0.68, 0.77, 0.70};
const double PTCerr[5]  = {0.03, 0.03, 0.03, 0.04, 0.04};

const double PTAl[5]    = {0.49, 0.52, 0.57, 0.59, 0.71};
const double PTAlerr[5] = {0.02, 0.03, 0.04, 0.06, 0.15};

const double PTCu[5]    = {0.45, 0.45, 0.43, 0.52, 0.53};
const double PTCuerr[5] = {0.02, 0.02, 0.02, 0.02, 0.03};

//Proton Transparency (low energy)

const double PrTrLeEn[7]    = {150.0, 170.0, 190.0, 210.0, 230.0, 250.0, 270.0};
const double PrTrLeAl[7]    = {0.66, 0.62, 0.64, 0.61, 0.64, 0.61, 0.66};
const double PrTrLeAlerr[7] = {0.04, 0.04, 0.04, 0.05, 0.04, 0.04, 0.04};
const double PrTrLeC[7]     = {0.76, 0.72, 0.73, 0.71, 0.73, 0.71, 0.75};
const double PrTrLeCerr[7]  = {0.03, 0.03, 0.04, 0.04, 0.04, 0.04, 0.04};
const double PrTrLeLi[7]    = {0.84, 0.82, 0.83, 0.81, 0.83, 0.81, 0.84};
const double PrTrLeLierr[7] = {0.02, 0.02, 0.02, 0.03, 0.02, 0.03, 0.02}; 

//Pion Transparency (low energy)

	//charged

const double PiTrLeEnp[9]    = {50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0};
const double PiTrLeAlp[9]    = {0.57, 0.44, 0.34, 0.29, 0.30, 0.31, 0.32, 0.34, 0.39};
const double PiTrLeAlperr[9] = {0.04, 0.05, 0.04, 0.05, 0.04, 0.04, 0.06, 0.04, 0.06};
const double PiTrLeCp[9]     = {0.68, 0.57, 0.45, 0.39, 0.41, 0.41, 0.43, 0.46, 0.49};
const double PiTrLeCperr[9]  = {0.03, 0.04, 0.04, 0.05, 0.04, 0.04, 0.06, 0.06, 0.04};
const double PiTrLeLip[9]    = {0.79, 0.69, 0.59, 0.54, 0.55, 0.56, 0.56, 0.60, 0.64};
const double PiTrLeLiperr[9] = {0.03, 0.04, 0.04, 0.07, 0.06, 0.06, 0.04, 0.05, 0.04};

	//neutral
	
const double PiTrLeEn0[9]    = {45.0, 95.0, 145.0, 195.0, 245.0, 295.0, 345.0, 395.0, 445.0};
const double PiTrLeAl0[9]    = {0.61, 0.41, 0.29, 0.30, 0.37, 0.38, 0.39, 0.36, 0.31}; 
const double PiTrLeAl0err[9] = {0.04, 0.07, 0.04, 0.05, 0.04, 0.05, 0.05, 0.04, 0.05}; 
const double PiTrLeC0[9]     = {0.70, 0.53, 0.39, 0.41, 0.48, 0.49, 0.51, 0.48, 0.43};
const double PiTrLeC0err[9]  = {0.05, 0.08, 0.04, 0.04, 0.04, 0.04, 0.03, 0.05, 0.07};
const double PiTrLeLi0[9]    = {0.80, 0.66, 0.53, 0.55, 0.63, 0.63, 0.65, 0.62, 0.57};
const double PiTrLeLi0err[9] = {0.02, 0.05, 0.06, 0.05, 0.03, 0.03, 0.04, 0.03, 0.07};

//SciBooNE CC total xsec

const double sbccneut[6] = {2.76e-39, 5.80e-39, 1.03e-38, 1.38e-38, 1.62e-38, 1.74e-38};
const double sbccnuance[6] = {3.40e-39, 6.39e-39, 1.01e-38, 1.29e-38, 1.56e-38, 1.66e-38};

const double nomadcc[30] = {0.786, 0.763, 0.722, 0.701, 0.716, 0.706, 0.705, 0.697, 0.7, 0.698, 0.698, 0.7, 0.699, 0.694, 0.694, 0.694, 0.677, 0.681, 0.675, 0.682, 0.67, 0.675, 0.684, 0.678, 0.677, 0.674, 0.661, 0.671, 0.667, 0.721};
const double minoscc[13] = {0.748, 0.711, 0.708, 0.722, 0.699, 0.691, 0.708, 0.689, 0.683, 0.686, 0.675, 0.675, 0.676};
const double minosccbar[11] = {0.305, 0.3, 0.303, 0.314, 0.304, 0.316, 0.32, 0.332, 0.325, 0.352, 0.324};

//MiniBooNE CCpi+-like/CCQE-like

const double mbraten[13]    = {0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.50, 1.70, 2.10};
const double mbratenerr[13] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10, 0.10, 0.30};

const double mbrat[13]    = {0.036, 0.100, 0.191, 0.278, 0.371, 0.465, 0.551, 0.607, 0.677, 0.700, 0.777, 0.904, 1.022};
const double mbraterr[13] = {0.005, 0.011, 0.019, 0.028, 0.040, 0.053, 0.066, 0.077, 0.091, 0.097, 0.109, 0.137, 0.161};

const double mbccpien[27]    = {0.550, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.850, 1.950};
const double mbccpienerr[27] = {0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050};

const double mbccpi[27]    = {0.61e-38, 1.15e-38, 1.58e-38, 1.96e-38, 2.41e-38, 2.83e-38, 3.25e-38, 3.73e-38, 4.16e-38, 4.66e-38, 4.97e-38, 5.29e-38, 5.63e-38, 5.91e-38, 6.23e-38, 6.67e-38, 7.03e-38, 7.23e-38, 7.76e-38, 8.08e-38, 8.37e-38, 8.64e-38, 8.85e-38, 9.30e-38, 9.26e-38, 9.71e-38, 9.92e-38};
const double mbccpierr[27] = {0.08e-38, 0.13e-38, 0.17e-38, 0.21e-38, 0.25e-38, 0.29e-38, 0.33e-38, 0.39e-38, 0.44e-38, 0.50e-38, 0.56e-38, 0.61e-38, 0.70e-38, 0.76e-38, 0.84e-38, 0.94e-38, 1.05e-38, 1.14e-38, 1.28e-38, 1.39e-38, 1.49e-38, 1.59e-38, 1.69e-38, 1.85e-38, 1.93e-38, 2.14e-38, 2.37e-38};
	
const double mbccqeen[14]    = {0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.850, 0.950, 1.050, 1.200, 1.400, 1.750};
const double mbccqeenerr[14] = {0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 0.050, 0.100, 0.100, 0.250};

const double mbccqe[14]    = {5.8296e-38, 6.0756e-38, 6.456e-38, 6.9048e-38, 7.2426e-38, 7.5528e-38, 7.728e-38, 7.9332e-38, 8.1288e-38, 8.1804e-38, 8.2602e-38, 8.151e-38, 8.1444e-38, 7.7052e-38};
const double mbccqeerr[14] = {1.1982e-38, 0.9192e-38, 0.798e-38, 0.7254e-38, 0.6744e-38, 0.6534e-38, 0.639e-38, 0.6468e-38, 0.6774e-38, 0.7302e-38, 0.8154e-38, 0.9972e-38, 1.2696e-38, 1.5678e-38};
	
const double sbtoten[6]      = {0.38, 0.62, 0.87, 1.11, 1.43, 2.47};
const double sbtotneut[6]    = {0.276e-38, 0.580e-38, 1.030e-38, 1.380e-38, 1.630e-38, 1.740e-38};
const double sbtotneuterr[6] = {0.075e-38, 0.075e-38, 0.100e-38, 0.170e-38, 0.290e-38, 0.380e-38};
const double sbtotnuan[6]    = {0.340e-38, 0.639e-38, 1.010e-38, 1.290e-38, 1.560e-38, 1.660e-38};
const double sbtotnuanerr[6] = {0.096e-38, 0.081e-38, 0.090e-38, 0.150e-38, 0.280e-38, 0.370e-38};

const double mbccpi0en[14]    = {0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.90};
const double mbccpi0enerr[14] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10};
const double mbccpi0[14]      = {0.176e-38, 0.383e-38, 0.568e-38, 0.731e-38, 0.920e-38, 1.106e-38, 1.242e-38, 1.389e-38, 1.523e-38, 1.638e-38, 1.820e-38, 1.937e-38, 2.080e-38, 2.192e-38};
const double mbccpi0err[14]   = {0.018e-38, 0.023e-38, 0.026e-38, 0.029e-38, 0.034e-38, 0.041e-38, 0.049e-38, 0.060e-38, 0.074e-38, 0.092e-38, 0.124e-38, 0.162e-38, 0.216e-38, 0.224e-38};

//Oset 

const double pq1[4] = {0.90, 0.09, 0.01, 0.00};
const double pq2[4] = {0.69, 0.25, 0.05, 0.01};

const double pa1[5] = {0.81, 0.17, 0.02, 0.00, 0.00};
const double pa2[5] = {0.37, 0.41, 0.17, 0.04, 0.01};

#endif
