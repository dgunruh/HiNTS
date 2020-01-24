package util;

public class Constants {
	
	public static final double pi = 3.14159265358979323846;
	public static final double tpi = 2.0 * pi;
	public static final double fpi = 4.0 * pi;
	public static final double sqrtpi = 1.77245385090551602729; 
	public static final double sqrtpm1 = 1.0 / sqrtpi;
	public static final double sqrt2  = 1.41421356237309504880;
	
	//  ... physical constants, si (nist codata 2006), web version 5.1
	//      http://physics.nist.gov/constants
	public static final double h_planck_si      = 6.62606896e-34  ;// // j s
	public static final double k_boltzmann_si   = 1.3806504e-23   ;// // j k^-1 
	public static final double electron_si      = 1.602176487e-19 ;// // c
	public static final double electronvolt_si  = 1.602176487e-19 ;// // j  
	public static final double electronmass_si  = 9.10938215e-31  ;// // kg
	public static final double hartree_si       = 4.35974394e-18  ;// // j
	public static final double rydberg_si       = hartree_si/2.0  ;// // j
	public static final double bohr_radius_si   = 0.52917720859e-10 ;//// m
	public static final double amu_si           = 1.660538782e-27  ;//// kg
	public static final double c_si             = 2.99792458e+8    ;//// m sec^-1
	 
	// covert voltage in volt to ry atomic units 
	public static final double volt_ry          = electron_si/sqrt2/rydberg_si ;
	public static final double ry_volt          = 1.0/volt_ry ;
	
	//convert rydbergs to meters
	public static final double ry_to_m			= bohr_radius_si; // m
	
	// ... physical constants, atomic units:
	// ... au for "hartree" atomic units (e = m = hbar = 1)
	// ... ry for "rydberg" atomic units (e^2=2, m=1/2, hbar=1)
	public static final double	k_boltzmann_au   = k_boltzmann_si / hartree_si;
	public static final double	k_boltzmann_ry   = k_boltzmann_si / rydberg_si;
		  //
		  // ... unit conversion factors: energy and masses
		  //
	public static final double	autoev           = hartree_si / electronvolt_si;
	public static final double	rytoev           = autoev / 2.0;
	public static final double  rytojoule		 = rytoev*electronvolt_si;
	public static final double	evtory           = 2.0 / autoev;
	public static final double	amu_au           = amu_si / electronmass_si;
	public static final double	amu_ry           = amu_au / 2.0;
		  //
		  // ... unit conversion factors: atomic unit of time, in s and ps
		  //
	public static final double	au_sec           = h_planck_si/tpi/hartree_si;
	public static final double	au_ps            = au_sec * 1.0e+12;
	public static final double	ry_sec           = h_planck_si/tpi/rydberg_si;
	public static final double	ry_ps            = ry_sec * 1.0e+12;
		  //
		  // ... unit conversion factors: pressure (1 pa = 1 j/m^3, 1gpa = 10 kbar )
		  //
	public static final double	au_gpa           = hartree_si / Math.pow(bohr_radius_si, 3) / 1.0e+9 ;
	public static final double	ry_kbar          = 10.0 * au_gpa / 2.0;
		  //
		  // ... unit conversion factors: 1 debye = 10^-18 esu*cm 
		  // ...                                  = 3.3356409519*10^-30 c*m 
		  // ...                                  = 0.208194346 e*a
		  // ... ( 1 esu = (0.1/c) am, c=299792458 m/s)
		  //
	public static final double	debye_si         = 3.3356409519 * 1.0e-30; // c*m 
	public static final double	au_debye         = electron_si * bohr_radius_si / debye_si;
		  //
		  // unit conversion for temperature 
		  //
	public static final double	ev_to_kelvin = electronvolt_si / k_boltzmann_si;
	public static final double	ry_to_kelvin = rydberg_si / k_boltzmann_si;
	public static final double	kelvintory = k_boltzmann_si / rydberg_si;
		  // 
		  // .. unit conversion factors: energy to wavelength
		  //
	public static final double	evtonm = 1e+9 * h_planck_si * c_si / electronvolt_si;
	public static final double	rytonm = 1e+9 * h_planck_si * c_si / rydberg_si;
		  //
		  //  speed of light in atomic units
		  //
	public static final double	c_au = c_si / bohr_radius_si * au_sec;
		  //
		  // ... zero up to a given accuracy
		  //
	public static final double	eps4  = 1.0e-4;
	public static final double	eps6  = 1.0e-6;
	public static final double	eps8  = 1.0e-8;
	public static final double	eps12 = 1.0e-12;
	public static final double	eps14 = 1.0e-14;
	public static final double	eps16 = 1.0e-16;
	public static final double	eps24 = 1.0e-24;
	public static final double	eps32 = 1.0e-32;
		  //
	public static final double	e2 = 2.0  ;    // the square of the electron charge
		  //
		  //////////// compatibiility
		  //
	public static final double	amconv = amu_ry;
	public static final double	bohr_radius_cm = bohr_radius_si * 100.0;
	public static final double	bohr_radius_angs = bohr_radius_cm * 1.0e8;
	public static final double	nmtobohr = 10.0/bohr_radius_angs;
	public static final double	bohrtonm = 1.0/nmtobohr;
	public static final double	angstrom_au = 1.0/bohr_radius_angs;
	public static final double	dip_debye = au_debye;
	public static final double	au_terahertz  = au_ps;
	public static final double	au_to_ohmcmm1 = 46000.0 ; // (ohm cm)^-1
	public static final double	ry_to_thz = 1.0 / au_terahertz / fpi;
	public static final double	ry_to_cmm1 = 1.e+10 * ry_to_thz / c_si;





}
