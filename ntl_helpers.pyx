# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

# Cursed typo
# Fixed in https://github.com/sagemath/sage/pull/36071
from sage.rings.polynomial.polynomial_zz_pex cimport Polynomial_ZZ_pX as Polynomial_ZZ_pEX

from cysignals.signals cimport sig_on, sig_off
from sage.libs.ntl.ZZ_pEX cimport ZZ_pEX_reverse

def ntl_reverse(f):    
    """
    SageMath currently doesn't use fast reversing
    from NTL, so we include it as an extension.
    """
    r = f.parent()(0)
    
    # First convert polynomials to c-type
    f_c = (<Polynomial_ZZ_pEX>f)
    r_c = (<Polynomial_ZZ_pEX>r)

    # Reverse within NTL (Skips conversion)
    sig_on()
    ZZ_pEX_reverse(r_c.x, f_c.x)
    sig_off()

    return r_c