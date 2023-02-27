/*
    Copyright (C) 2009, 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

#if FLINT_USES_PTHREAD
#include <pthread.h>

static pthread_once_t fmpz_initialised = PTHREAD_ONCE_INIT;
pthread_mutex_t fmpz_lock;
#endif

/* Always free larger mpz's to avoid wasting too much heap space */
#define FLINT_MPZ_MAX_CACHE_LIMBS 64

/* The number of new mpz's allocated at a time */
#define MPZ_BLOCK 64 

/* there's no point using TLS here as GC doesn't support it */
__mpz_struct ** mpz_free_arr = NULL;
__mpz_struct ** mpz_arr = NULL;
ulong mpz_num = 0;
ulong mpz_alloc = 0;
ulong mpz_free_num = 0;
ulong mpz_free_alloc = 0;

#if FLINT_USES_PTHREAD
void fmpz_lock_init()
{
   pthread_mutex_init(&fmpz_lock, NULL);
}
#endif

__mpz_struct * _fmpz_new_mpz(void)
{
    __mpz_struct * z = NULL;

#if FLINT_USES_PTHREAD
    pthread_once(&fmpz_initialised, fmpz_lock_init);
    pthread_mutex_lock(&fmpz_lock);
#endif

    if (mpz_free_num != 0)
        z = mpz_free_arr[--mpz_free_num];
    else
    {
        z = flint_malloc(sizeof(__mpz_struct));

        if (mpz_num == mpz_alloc) /* store pointer to prevent gc cleanup */
        {
            mpz_alloc = FLINT_MAX(64, mpz_alloc * 2);
            mpz_arr = flint_realloc(mpz_arr, mpz_alloc * sizeof(__mpz_struct *));
        }
        mpz_arr[mpz_num++] = z;

        mpz_init(z);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&fmpz_lock);
#endif

    return z;
}

void _fmpz_clear_mpz(fmpz f)
{
    __mpz_struct * ptr = COEFF_TO_PTR(f);

    if (ptr->_mp_alloc > FLINT_MPZ_MAX_CACHE_LIMBS)
        mpz_realloc2(ptr, 1);

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&fmpz_lock);
#endif

    if (mpz_free_num == mpz_free_alloc)
    {
        mpz_free_alloc = FLINT_MAX(64, mpz_free_alloc * 2);
        mpz_free_arr = flint_realloc(mpz_free_arr, mpz_free_alloc * sizeof(__mpz_struct *));
    }

    mpz_free_arr[mpz_free_num++] = ptr;

#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&fmpz_lock);
#endif
}

void _fmpz_cleanup_mpz_content(void)
{
    ulong i;

    for (i = 0; i < mpz_free_num; i++)
    {
        mpz_clear(mpz_free_arr[i]);
        flint_free(mpz_free_arr[i]);
    }

    /* TODO: remove selected mpz's from mpz_arr too and compact */
    mpz_free_num = mpz_free_alloc = 0;
}

void _fmpz_cleanup(void)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&fmpz_lock);
#endif

    _fmpz_cleanup_mpz_content();
    flint_free(mpz_free_arr);
    mpz_free_arr = NULL;

#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&fmpz_lock);
#endif
}

__mpz_struct * _fmpz_promote(fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f)) /* f is small so promote it first */
    {
        __mpz_struct * mf = _fmpz_new_mpz();
        (*f) = PTR_TO_COEFF(mf);
        return mf;
    }
    else /* f is large already, just return the pointer */
        return COEFF_TO_PTR(*f);
}

__mpz_struct * _fmpz_promote_val(fmpz_t f)
{
    fmpz c = (*f);
    if (!COEFF_IS_MPZ(c)) /* f is small so promote it */
    {
        __mpz_struct * mf = _fmpz_new_mpz();
        (*f) = PTR_TO_COEFF(mf);
        flint_mpz_set_si(mf, c);
        return mf;
    }
    else /* f is large already, just return the pointer */
        return COEFF_TO_PTR(c);
}

void _fmpz_demote_val(fmpz_t f)
{
    __mpz_struct * mf = COEFF_TO_PTR(*f);
    int size = mf->_mp_size;

    if (size == 1 || size == -1)
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
   __mpz_struct *ptr;
   *f = WORD(0);
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
