/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "fmpz.h"

__GMP_DECLSPEC extern void * (*__gmp_allocate_func) (size_t);
__GMP_DECLSPEC extern void * (*__gmp_reallocate_func) (void *, size_t, size_t);
__GMP_DECLSPEC extern void   (*__gmp_free_func) (void *, size_t);
#define ALLOC(x) (x)->_mp_alloc
#define SIZ(x) (x)->_mp_size
#define PTR(x) (x)->_mp_d

mpz_ptr _fmpz_new_mpz2(mp_size_t limbs)
{
    mpz_ptr mf = (mpz_ptr) flint_malloc(sizeof(__mpz_struct));
    /* mpz_init2(mf, 2*FLINT_BITS); */

    limbs = FLINT_MAX(2, limbs);

    ALLOC(mf) = limbs;
    SIZ(mf) = 0;
    PTR(mf) = __gmp_allocate_func(sizeof(mp_limb_t) * limbs);
    return mf;
}

void _fmpz_clear_mpz(fmpz f)
{
    mpz_ptr mf = COEFF_TO_PTR(f);
    if (ALLOC(mf))
        __gmp_free_func(PTR(mf), sizeof(mp_limb_t) * ALLOC(mf));
    flint_free(mf);
}

void _fmpz_cleanup_mpz_content(void)
{
}

void _fmpz_cleanup(void)
{
}

__mpz_struct * _fmpz_promote(fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))  /* f is small so promote it first */
    {
        __mpz_struct * mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
        return mf;
    }
    else  /* f is large already, just return the pointer */
        return COEFF_TO_PTR(*f);
}

__mpz_struct * _fmpz_promote_val(fmpz_t f)
{
    fmpz c = *f;
    if (!COEFF_IS_MPZ(c))  /* f is small so promote it */
    {
        __mpz_struct * mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
        flint_mpz_set_si(mf, c);
        return mf;
    }
    else  /* f is large already, just return the pointer */
        return COEFF_TO_PTR(*f);
}

void _fmpz_demote_val(fmpz_t f)
{
    __mpz_struct * mf = COEFF_TO_PTR(*f);
    int size = mf->_mp_size;

    if (!(((unsigned int) size + 1U) & ~2U))  /* size +-1 */
    {
        ulong uval = mf->_mp_d[0];

        if (uval <= (ulong) COEFF_MAX)
        {
            _fmpz_clear_mpz(*f);
            *f = size * (fmpz) uval;
        }
    }
    else if (size == 0)  /* value is 0 */
    {
        _fmpz_clear_mpz(*f);
        *f = 0;
    }

    /* don't do anything if value has to be multi precision */
}

void _fmpz_init_readonly_mpz(fmpz_t f, const mpz_t z)
{
   __mpz_struct * mf = (__mpz_struct *) flint_malloc(sizeof(__mpz_struct));
    *f = PTR_TO_COEFF(mf);
    *mf = *z;
}

void _fmpz_clear_readonly_mpz(mpz_t z)
{
    if (((z->_mp_size == 1 || z->_mp_size == -1) && (z->_mp_d[0] <= COEFF_MAX))
        || (z->_mp_size == 0))
    {
        if (ALLOC(z))
            __gmp_free_func(PTR(z), sizeof(mp_limb_t) * ALLOC(z));
    }
}
