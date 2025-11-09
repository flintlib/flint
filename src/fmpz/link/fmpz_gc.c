/*
    Copyright (C) 2009, 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"

#if FLINT_USES_PTHREAD
#include <pthread.h>

static pthread_once_t fmpz_initialised = PTHREAD_ONCE_INIT;
pthread_mutex_t fmpz_lock;
#endif

/* Always free larger mpz's to avoid wasting too much heap space */
#define FLINT_MPZ_MAX_CACHE_LIMBS 64

#if FLINT_MPZ_MAX_CACHE_LIMBS < MPZ_MIN_ALLOC
# error
#endif

/* The number of new mpz's allocated at a time */
#define MPZ_BLOCK 64

/* there's no point using TLS here as GC doesn't support it */
mpz_ptr * mpz_free_arr = NULL;
mpz_ptr * mpz_arr = NULL;
ulong mpz_num = 0;
ulong mpz_alloc = 0;
ulong mpz_free_num = 0;
ulong mpz_free_alloc = 0;

#if FLINT_USES_PTHREAD
static void fmpz_lock_init()
{
   pthread_mutex_init(&fmpz_lock, NULL);
}
#endif

mpz_ptr _fmpz_new_mpz(void)
{
    mpz_ptr z = NULL;

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
            mpz_alloc = FLINT_MAX(MPZ_BLOCK, 2 * mpz_alloc);
            mpz_arr = flint_realloc(mpz_arr, mpz_alloc * sizeof(mpz_ptr));
        }
        mpz_arr[mpz_num++] = z;

        mpz_init2(z, MPZ_MIN_ALLOC * FLINT_BITS);
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&fmpz_lock);
#endif

    return z;
}

void _fmpz_clear_mpz(fmpz f)
{
    mpz_ptr ptr = COEFF_TO_PTR(f);

    FLINT_ASSERT(ptr->_mp_alloc >= MPZ_MIN_ALLOC);

    if (ptr->_mp_alloc > FLINT_MPZ_MAX_CACHE_LIMBS)
        mpz_realloc(ptr, MPZ_MIN_ALLOC);

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&fmpz_lock);
#endif

    if (mpz_free_num == mpz_free_alloc)
    {
        mpz_free_alloc = FLINT_MAX(MPZ_BLOCK, 2 * mpz_free_alloc);
        mpz_free_arr = flint_realloc(mpz_free_arr, mpz_free_alloc * sizeof(mpz_ptr));
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

mpz_ptr _fmpz_promote(fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f)) /* f is small so promote it first */
    {
        mpz_ptr mf = _fmpz_new_mpz();
        (*f) = PTR_TO_COEFF(mf);
        return mf;
    }
    else /* f is large already, just return the pointer */
        return COEFF_TO_PTR(*f);
}

mpz_ptr _fmpz_promote_val(fmpz_t f)
{
    fmpz c = (*f);
    if (!COEFF_IS_MPZ(c)) /* f is small so promote it */
    {
        mpz_ptr mf = _fmpz_new_mpz();
        (*f) = PTR_TO_COEFF(mf);
        flint_mpz_set_si(mf, c);
        return mf;
    }
    else /* f is large already, just return the pointer */
        return COEFF_TO_PTR(c);
}

void _fmpz_demote_val(fmpz_t f)
{
    mpz_ptr mf = COEFF_TO_PTR(*f);
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
    mpz_ptr ptr;
    *f = WORD(0);
    ptr = _fmpz_promote(f);

    mpz_clear(ptr);
    *ptr = *z;
}

void _fmpz_clear_readonly_mpz(mpz_t z)
{
    int size = z->_mp_size;

    if (size == 0 || ((size == 1 || size == -1) && (z->_mp_d[0] <= COEFF_MAX)))
        mpz_clear(z);
}
