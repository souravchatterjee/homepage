import scipy.integrate

def inspiral_time_peters(a0,e0,m1,m2):
    	"""
    	Computes the inspiral time, in Gyr, for a binary
	    a0 in Au, and masses in solar masses
    	"""

    	coef = 6.086768e-11 #G^3 / c^5 in au, gigayear, solar mass units
    	beta = (64./5.) * coef * m1 * m2 * (m1+m2)

    	if e0 == 0:
        	return a0**4 / (4*beta)

    	c0 = a0 * (1.-e0**2.) * e0**(-12./19.) * (1.+(121./304.)*e0**2.)**(-870./2299.)

    	time_integrand = lambda e: e**(29./19.)*(1.+(121./304.)*e**2.)**(1181./2299.) / (1.-e**2.)**1.5
    	integral,abserr = scipy.integrate.quad(time_integrand,0,e0)

    	return integral * (12./19.) * c0**4. / beta
