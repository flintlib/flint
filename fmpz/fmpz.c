/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

/* Always free larger mpz's to avoid wasting too much heap space */
#define FLINT_MPZ_MAX_CACHE_LIMBS 64

/* The number of new mpz's allocated at a time */
#define MPZ_BLOCK 64 

#if HAVE_TLS
#define TLS_PREFIX __thread
#else
#define TLS_PREFIX
#endif

TLS_PREFIX __mpz_struct ** mpz_free_arr = NULL;
TLS_PREFIX ulong mpz_free_num = 0;
TLS_PREFIX ulong mpz_free_alloc = 0;

__mpz_struct * _fmpz_new_mpz(void)
{
    if (mpz_free_num != 0)
    {
        return mpz_free_arr[--mpz_free_num];
    }
    else
    {
        __mpz_struct * z = flint_malloc(sizeof(__mpz_struct));
        mpz_init(z);
        return z;
    }
}

void _fmpz_clear_mpz(fmpz f)
{
    __mpz_struct * ptr = COEFF_TO_PTR(f);

    if (ptr->_mp_alloc > FLINT_MPZ_MAX_CACHE_LIMBS)
        mpz_realloc2(ptr, 1);

    if (mpz_free_num == mpz_free_alloc)
    {
        mpz_free_alloc = FLINT_MAX(64, mpz_free_alloc * 2);
        mpz_free_arr = flint_realloc(mpz_free_arr, mpz_free_alloc * sizeof(__mpz_struct *));
    }

    mpz_free_arr[mpz_free_num++] = ptr;
}

void _fmpz_cleanup_mpz_content(void)
{
    ulong i;

    for (i = 0; i < mpz_free_num; i++)
    {
        mpz_clear(mpz_free_arr[i]);
        flint_free(mpz_free_arr[i]);
    }
}

void _fmpz_cleanup(void)
{
    _fmpz_cleanup_mpz_content();
    flint_free(mpz_free_arr);
}

__mpz_struct * _fmpz_promote(fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f)) /* f is small so promote it first */
    {
        __mpz_struct * mpz_ptr = _fmpz_new_mpz();
        (*f) = PTR_TO_COEFF(mpz_ptr);
        return mpz_ptr;
    }
    else /* f is large already, just return the pointer */
        return COEFF_TO_PTR(*f);
}

__mpz_struct * _fmpz_promote_val(fmpz_t f)
{
    fmpz c = (*f);
    if (!COEFF_IS_MPZ(c)) /* f is small so promote it */
    {
        __mpz_struct * mpz_ptr = _fmpz_new_mpz();
        (*f) = PTR_TO_COEFF(mpz_ptr);
        mpz_set_si(mpz_ptr, c);
        return mpz_ptr;
    }
    else /* f is large already, just return the pointer */
        return COEFF_TO_PTR(c);
}

void _fmpz_demote_val(fmpz_t f)
{
    __mpz_struct * mpz_ptr = COEFF_TO_PTR(*f);
    int size = mpz_ptr->_mp_size;

    if (!(((unsigned int) size + 1U) & ~2U))  /* size +-1 */
    {
        ulong uval = mpz_ptr->_mp_d[0];

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
   __mpz_struct *ptr;
   *f = 0L;
   ptr = _fmpz_promote(f);

   mpz_clear(ptr);
   *ptr = *z;
}

void _fmpz_clear_readonly_mpz(mpz_t z)
{
    if (((z->_mp_size == 1 || z->_mp_size == -1) && (z->_mp_d[0] <= COEFF_MAX))
        || (z->_mp_size == 0))
    {
        mpz_clear(z);
    }

}
