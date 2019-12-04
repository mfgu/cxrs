implements("const");
provide("const");

static variable 
        h= 6.6255E-27, % ;plank's constant, erg sec.
        hc= 12.398421E3, % ;plank's constant * light speed, eV*A
        hbc= 1.97327053E3, % ; h_bar*c, eV*A
        kb= 8.6174E-5, % ; Boltzmann constant, ev/K
        e= 1.60218e-19, % ; electron charge, coulomb
        c= 2.9979e10, % ; speed of light cm/s
        c10= 2.9979, % ; speed of light 10^10 cm/s
        eV= 1.602e-12, % ; ergs/eV
        Mp= 1.6726e-24, % ; mass of proton, gm
        Me= 9.1094e-28, % ; mass of electron, gm
        Mp_MeV= 938.272, % ; mass of proton MeV
        Me_keV= 511.999, % ; mass of electron keV
        Mp_keV= 938.272E3, % ; mass of proton keV
        Me_eV= 511.999E3, % ; mass of electron eV
        re= 2.81794092, % ;electron classical radius, fm
        RBohr= 5.29177e-9, % ;bohr radius, cm
        Ryd_eV= 13.606, % ;rydberg energy, eV
        Hartree_eV= 27.2113962, % ;Hartree, eV
        Alpha= 7.29735308E-3, % ;Fine structure constant.
        Rate_AU = 4.13413733E16, % ;Atomic Rate unit 1/s
        Rate_AU10 = 4.13413733E6, % ;Atomic Rate unit 10^10/s
        Area_AU20 = 2.80028560859E3, % ;Atomic Area unit 10^-20 cm2
        FWHM = 2.35482005, %  ; conversion from sig to FWHM for gaussian.
        sig_t = 6.6524616, %  ; thompson cross section in 10^{-25} cm2
        G = 6.673E-11, % ; newton constant of gravity in SI.
        G_cgs = 6.673E-8, % ; newton constant of gravity in cgs.
        Maxwellian = 1.12837967, % ; constant before the maxwellian distribution.
        Msun33= 1.989, % ;solar mass, 10^33 gm
        Lsun33= 3.826, % ;solar luminosity, 10^33 ergs/sec
        Rsun= 6.596e10, % ;solar radius, cm
        AU= 1.496e13, % ;sun to earth distance, cm
        pc= 3.086e18, % ;parsec, cm
        barn= 1.e-24, % ;barn, cm2
        AtomicSymbol = ["", "H", "He", "Li", "Be", "B", "C", "N", "O", "F",
			"Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
			"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", 
			"Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
			"Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", 
			"Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
			"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", 
			"Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
			"Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
			"Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
			"Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
			"Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"],
        RomanNumber = ["","I","II","III","IV","V","VI","VII","VIII","IX","X",
		       "XI","XII","XIII","XIV","XV","XVI","XVII","XVIII",
		       "XIX","XX","XXI","XXII","XXIII","XXIV","XXV","XXVI",
		       "XXVII","XXVIII","XXIX","XXX"];
static define E2V(e, m) {
  variable k;

  if (m == 0) {
    k = e/Hartree_eV;
    k = 2.0*k*(1.0 + 0.5*Alpha^2*k);
    k = Alpha^2*k;
    k = sqrt(k/(1.0+k));
    k /= Alpha;
    k *= RBohr*Rate_AU*1E-10;
  } else {
    m = Mp_keV*m/Me_keV;
    k = e/Hartree_eV;
    k = 2.0*m*k*(1.0+0.5*Alpha^2*k/m);
    k = Alpha^2*k;
    k = sqrt(k/(m*m+k));
    k /= Alpha;
    k *= RBohr*Rate_AU*1E-10;
  }

  return k;
}

static define gasdens(p, t) {
  return (p*133.32/(1.38e-23*t))*1e-6;
}
