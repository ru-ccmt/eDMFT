"""SciPy spherical harmonics in the polar-angle-first LAPW convention."""
try:
    from scipy.special import sph_harm_y
except ImportError:
    from scipy.special import sph_harm

    def sph_harm_y(ell, m, theta, phi):
        """Return Y_ell^m(theta, phi) with theta polar and phi azimuthal."""
        # SciPy before 1.15 takes m before ell and azimuth before polar angle.
        return sph_harm(m, ell, phi, theta)
