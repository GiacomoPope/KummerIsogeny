# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

"""
Additional NTL library bindings missing from Sage.
"""
from cysignals.signals cimport sig_on, sig_off

from sage.libs.ntl.types cimport ZZ_pE_c, vec_ZZ_pE_c
from sage.libs.ntl.ZZ_pEX cimport (
    ZZ_pEX_eval, 
    ZZ_pEX_eval_vec, 
    ZZ_pEX_reverse, 
    ZZ_pEX_reverse_hi, 
    ZZ_pEX_trunc,
    ZZ_pEX_MulTrunc, 
    ZZ_pEX_InvTrunc,
    ZZ_pEX_MulMod,
    ZZ_pEX_RightShift
)

from sage.rings.integer cimport Integer
from sage.libs.ntl.types cimport ZZ_pX_c, ZZ_p_c, ZZ_c
from sage.libs.ntl.ZZ_pE cimport ZZ_pE_to_ZZ_pX
from sage.libs.ntl.ZZ_pX cimport ZZ_pX_deg, ZZ_pX_coeff
from sage.libs.ntl.ZZ_p cimport ZZ_p_rep
from sage.libs.ntl.convert cimport ZZ_to_mpz
from sage.libs.ntl.ntl_ZZ_pE cimport ntl_ZZ_pE
from sage.libs.ntl.ntl_ZZ_pEX cimport ntl_ZZ_pEX

# Cursed typo
# Fixed in https://github.com/sagemath/sage/pull/36071
from sage.rings.polynomial.polynomial_zz_pex cimport Polynomial_ZZ_pX as Polynomial_ZZ_pEX


def ntl_reverse(f):    
    # Create a new r in a stupid way
    # I'm not sure how much of a performance
    # penalty this is...
    r = f.parent()(0)
    
    # Convert to c-type
    f_c = (<Polynomial_ZZ_pEX>f)
    r_c = (<Polynomial_ZZ_pEX>r)

    sig_on()
    ZZ_pEX_reverse(r_c.x, f_c.x)
    sig_off()

    return r_c

def ntl_reverse_self(f):        
    # Convert to c-type
    f_c = (<Polynomial_ZZ_pEX>f)

    sig_on()
    ZZ_pEX_reverse(f_c.x, f_c.x)
    sig_off()

    return f_c

def ntl_mul_and_truncate(f, g, m):
    # Create a new r in a stupid way
    # I'm not sure how much of a performance
    # penalty this is...
    r = f.parent()(0)

    f_c = (<Polynomial_ZZ_pEX>f)
    g_c = (<Polynomial_ZZ_pEX>g)
    r_c = (<Polynomial_ZZ_pEX>r)

    if m > 0:
        sig_on()
        ZZ_pEX_MulTrunc(r_c.x, f_c.x, g_c.x, m)
        sig_off()

    return r_c

def ntl_inv_and_truncate(f, m):
    # Create a new r in a stupid way
    # I'm not sure how much of a performance
    # penalty this is...
    r = f.parent()(0)

    f_c = (<Polynomial_ZZ_pEX>f)
    r_c = (<Polynomial_ZZ_pEX>r)

    if m > 0:
        sig_on()
        ZZ_pEX_InvTrunc(r_c.x, f_c.x, m)
        sig_off()

    return r_c

def ntl_right_shift(f, m):
    r = f.parent()(0)

    f_c = (<Polynomial_ZZ_pEX>f)
    r_c = (<Polynomial_ZZ_pEX>r)

    sig_on()
    ZZ_pEX_RightShift(r_c.x, f_c.x, m)
    sig_off()

    return r_c

def ntl_right_shift_self(f, m):
    f_c = (<Polynomial_ZZ_pEX>f)

    sig_on()
    ZZ_pEX_RightShift(f_c.x, f_c.x, m)
    sig_off()

    return f_c
