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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"

/* The number of new mpz's allocated at a time */
#define MPZ_BLOCK 64 

/* The array of mpz's used by the F_mpz type */
__mpz_struct * fmpz_arr;

/* Total number of mpz's initialised in F_mpz_arr */
ulong fmpz_allocated = 0;

/* An array of pointers to mpz's which are not being used presently */
ulong * fmpz_unused_arr;

ulong fmpz_num_unused = 0;

__mpz_struct * _fmpz_new_mpz(void)
{
    if (fmpz_num_unused == 0) /* time to allocate MPZ_BLOCK more mpz_t's */
    {
        ulong i;
        if (fmpz_allocated) /* realloc mpz_t's and unused array */
        {
            fmpz_arr = (__mpz_struct *) flint_realloc(fmpz_arr, (fmpz_allocated + MPZ_BLOCK) * sizeof(__mpz_struct));
            fmpz_unused_arr = (ulong *) flint_realloc(fmpz_unused_arr, (fmpz_allocated + MPZ_BLOCK) * sizeof(ulong));
        } else /* first time alloc of mpz_t's and unused array */
        {
            fmpz_arr = (__mpz_struct *) flint_malloc(MPZ_BLOCK*sizeof(__mpz_struct)); 
            fmpz_unused_arr = (ulong *) flint_malloc(MPZ_BLOCK*sizeof(ulong));
        }
        
        /* initialise the new mpz_t's and unused array */
        for (i = 0; i < MPZ_BLOCK; i++)
        {
            mpz_init(fmpz_arr + fmpz_allocated + i);
            fmpz_unused_arr[fmpz_num_unused] = fmpz_allocated + i;
            fmpz_num_unused++;
        }
        
        fmpz_allocated += MPZ_BLOCK;   
    } 
    
    fmpz_num_unused--;
    
    return fmpz_arr + fmpz_unused_arr[fmpz_num_unused];
}

void _fmpz_clear_mpz(fmpz f)
{
    fmpz_unused_arr[fmpz_num_unused] = (f ^ (1L << (FLINT_BITS - 2)));
    fmpz_num_unused++;   
}

void _fmpz_cleanup_mpz_content(void)
{
    long i;
    for (i = 0; i < fmpz_num_unused; i++)
        mpz_clear(fmpz_arr + fmpz_unused_arr[i]);
}

void _fmpz_cleanup(void)
{
    _fmpz_cleanup_mpz_content();

    if (fmpz_num_unused) flint_free(fmpz_unused_arr);
    if (fmpz_allocated) flint_free(fmpz_arr);
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

void _fmpz_init_readonly_mpz(fmpz_t f, mpz_t z)
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
