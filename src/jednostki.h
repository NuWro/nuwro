//----------------------------------------------------------------------
// Jednostki fizyczne w uk�adzie hbar=c=MeV=1
// (Pozwala nie pisa� hbar ani c we wzorach.
// aby otrzyma� wynik np. w cm  pisz  cout<<wynik/cm <<"[cm]"<<endl;
// idea podpatrzona w geancie, ale tam przyj�to inne jednostki podstawowe.) 
//
// Plik sprawdzony: wszelkie zmiany prosz� opatrywa� komentarzem
// (C. Juszczak 2001)
//----------------------------------------------------------------------
#ifndef _jednostki_h_
#define _jednostki_h_

///
/// Mathematical constants 
///

static const double Pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862;
//const double e  = 2.7182818284590452353602874713526624977572470936999595749669676277240766303535;
static const double Pi2= Pi*Pi;

///
/// natural units
/// 
static const double hbar = 1;
static const double c = 1;
/// Assuming also MeV=1 determines almost everything else

//
/// Physical units
///

// const double MeV = 0.001; // use this instead of MeV=1 to test if MeV is inserted where needed
static const double MeV = 1;
static const double eV = MeV / 1E6;
static const double GeV = 1000 * MeV;
static const double GeV2 = GeV*GeV;

static const double sek = 1 / (6.58211889 * 1E-22 * MeV);
static const double metr = sek / 299792458;
static const double kg = eV / (c * c) / (1.782661731 * 1E-36);


static const double J = kg * metr * metr / (sek * sek);

///
/// the following two constants (should be moved to a separate file one day)
///

const double G = 1.16639 * 1E-5 / (GeV * GeV);
//const double cos2thetac = 0.97418 * 0.97418;	//cos(Cabibbo angle)^2
const double cos2thetac = 0.974213 * 0.974213;	//assuming thetaC=13.04 deg
const double sin2thetaW = 0.2312215;  // weak mixing angle
//----------------------------------------------------------------------


// 
// Length [L]
//
static const double meter  = metr;                  
static const double meter2 = meter*meter;
static const double meter3 = meter*meter*meter;

static const double millimeter  = 0.001*meter;                        
static const double millimeter2 = millimeter*millimeter;
static const double millimeter3 = millimeter*millimeter*millimeter;

static const double centimeter  = 10.*millimeter;   
static const double centimeter2 = centimeter*centimeter;
static const double centimeter3 = centimeter*centimeter*centimeter;

static const double kilometer = 1000.*meter;                   
static const double kilometer2 = kilometer*kilometer;
static const double kilometer3 = kilometer*kilometer*kilometer;

static const double parsec = 3.0856775807e+16*meter;

static const double micrometer = 1.e-6 *meter;             
static const double  nanometer = 1.e-9 *meter;
static const double  angstrom  = 1.e-10*meter;
static const double  fermi     = 1.e-15*meter;
static const double  fermi2     = fermi*fermi;
static const double  fermi3     = fermi2*fermi;
static const double  fm 	    = 1.e-15*meter;

static const double      barn = 1.e-28*meter2;
static const double millibarn = 1.e-3 *barn;
static const double microbarn = 1.e-6 *barn;
static const double  nanobarn = 1.e-9 *barn;
static const double  picobarn = 1.e-12*barn;

// symbols
static const double mm  = millimeter;                        
static const double mm2 = millimeter2;
static const double mm3 = millimeter3;

static const double cm  = centimeter;   
static const double cm2 = centimeter2;
static const double cm3 = centimeter3;

static const double m  = meter;                  
static const double m2 = meter2;
static const double m3 = meter3;

static const double km  = kilometer;                   
static const double km2 = kilometer2;
static const double km3 = kilometer3;

static const double pc = parsec;

//
// Angle
//
static const double radian      = 1.;                  
static const double milliradian = 1.e-3*radian;
static const double degree = (3.14159265358979323846/180.0)*radian;

static const double   steradian = 1.;
	
// symbols
static const double rad  = radian;	
static const double mrad = milliradian;
static const double sr   = steradian;
static const double deg  = degree;

//
// Time [T]
//
static const double second      = sek;
static const double nanosecond  = 1.e-9 *second;
static const double millisecond = 1.e-3 *second;
static const double microsecond = 1.e-6 *second;
static const double  picosecond = 1.e-12*second;

static const double hertz = 1./second;
static const double kilohertz = 1.e+3*hertz;
static const double megahertz = 1.e+6*hertz;

// symbols
static const double ns = nanosecond;			
static const double  s = second;
static const double ms = millisecond;

//
// Mass [E][T^2][L^-2]
//
static const double  kilogram = kg;   
static const double      gram = 1.e-3*kilogram;
static const double milligram = 1.e-3*gram;

// symbols
//static const double  kg = kilogram;
static const double   g = gram;
static const double  mg = milligram;

//
// Electric current [Q][T^-1]
//
static const double      ampere = 1.;
static const double milliampere = 1.e-3*ampere;
static const double microampere = 1.e-6*ampere;
static const double  nanoampere = 1.e-9*ampere;

//
// Electric charge [Q]
//
static const double coulomb = ampere*second;
static const double e_SI  = 1.60217733e-19;	// positron charge in coulomb
static const double eplus = e_SI*coulomb ;		// positron charge

//
// Energy [E]
//
static const double joule = kg*m*m/(s*s);

//static const double     electronvolt = e_SI*joule;

static const double     electronvolt = eV;;
static const double kiloelectronvolt = 1.e+3*electronvolt;
static const double megaelectronvolt = 1.e+6*electronvolt; 
static const double gigaelectronvolt = 1.e+9*electronvolt;
static const double teraelectronvolt = 1.e+12*electronvolt;
static const double petaelectronvolt = 1.e+15*electronvolt;

// symbols
//static const double MeV = megaelectronvolt;
//static const double  eV = electronvolt;
static const double keV = kiloelectronvolt;
//static const double GeV = gigaelectronvolt;
static const double TeV = teraelectronvolt;
static const double PeV = petaelectronvolt;

//
// Power [E][T^-1]
//
static const double watt = joule/second;	// watt = 6.24150 e+3 * MeV/ns

//
// Force [E][L^-1]
//
static const double newton = joule/meter;	// newton = 6.24150 e+9 * MeV/mm

//
// Pressure [E][L^-3]
//
//#define pascal hep_pascal                          // a trick to avoid warnings 
static const double pascal = newton/m2;	   // pascal = 6.24150 e+3 * MeV/mm3
static const double bar        = 100000*pascal; // bar    = 6.24150 e+8 * MeV/mm3
static const double atmosphere = 101325*pascal; // atm    = 6.32420 e+8 * MeV/mm3

//
// Electric potential [E][Q^-1]
//
static const double megavolt = megaelectronvolt/eplus;
static const double kilovolt = 1.e-3*megavolt;
static const double     volt = 1.e-6*megavolt;

//
// Electric resistance [E][T][Q^-2]
//
static const double ohm = volt/ampere;	// ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

//
// Electric capacitance [Q^2][E^-1]
//
static const double farad = coulomb/volt;	// farad = 6.24150e+24 * eplus/Megavolt
static const double millifarad = 1.e-3*farad;
static const double microfarad = 1.e-6*farad;
static const double  nanofarad = 1.e-9*farad;
static const double  picofarad = 1.e-12*farad;

//
// Magnetic Flux [T][E][Q^-1]
//
static const double weber = volt*second;	// weber = 1000*megavolt*ns

//
// Magnetic Field [T][E][Q^-1][L^-2]
//
static const double tesla     = volt*second/meter2;	// tesla =0.001*megavolt*ns/mm2

static const double gauss     = 1.e-4*tesla;
static const double kilogauss = 1.e-1*tesla;

//
// Inductance [T^2][E][Q^-2]
//
static const double henry = weber/ampere;	// henry = 1.60217e-7*MeV*(ns/eplus)**2

//
// Temperature
//
static const double kelvin = 1.;

//
// Amount of substance
//
static const double mole = 1.;

//
// Activity [T^-1]
//
static const double becquerel = 1./second ;
static const double curie = 3.7e+10 * becquerel;

//
// Absorbed dose [L^2][T^-2]
//
static const double gray = joule/kilogram ;

//
// Luminous intensity [I]
//
static const double candela = 1.;

//
// Luminous flux [I]
//
static const double lumen = candela*steradian;

//
// Illuminance [I][L^-2]
//
static const double lux = lumen/meter2;

//
// Miscellaneous
//
static const double perCent     = 0.01 ;
static const double perThousand = 0.001;
static const double perMillion  = 0.000001;



#endif
