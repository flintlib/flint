/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef __unix__
# include <unistd.h> /* sysconf */
#endif

#if defined(_WIN32) || defined(WIN32)
# include <windows.h> /* GetSystemInfo */
#endif

#include "gmpcompat.h"
#include "fmpz.h"

#if FLINT_USES_PTHREAD
# include <pthread.h>
# include <stdatomic.h>
#endif

#if FLINT_USES_PTHREAD
typedef struct
{
   _Atomic(int) count;
   pthread_t thread;
   void * address;
} fmpz_block_header_s;
#else
typedef struct
{
   int count;
   void * address;
} fmpz_block_header_s;
#endif

/* Always free larger mpz's to avoid wasting too much heap space */
#define FLINT_MPZ_MAX_CACHE_LIMBS 64

#define PAGES_PER_BLOCK 16

/* The number of new mpz's allocated at a time */
#define MPZ_BLOCK 64

FLINT_TLS_PREFIX __mpz_struct ** mpz_free_arr = NULL;
FLINT_TLS_PREFIX ulong mpz_free_num = 0;
FLINT_TLS_PREFIX ulong mpz_free_alloc = 0;

static slong flint_page_size;
static slong flint_mpz_structs_per_block;
static slong flint_page_mask;

slong flint_get_page_size(void)
{
#if defined(__unix__)
   return sysconf(_SC_PAGESIZE);
#elif defined(_WIN32) || defined(WIN32)
   SYSTEM_INFO si;
   GetSystemInfo(&si);
   return si.dwPageSize;
#else
   return 4096;
#endif
}

void * flint_align_ptr(void * ptr, slong size)
{
    slong mask = ~(size - 1);

    return (void *)((mask & (slong) ptr) + size);
}

__mpz_struct * _fmpz_new_mpz(void)
{
    if (mpz_free_num == 0) /* allocate more mpz's */
    {
        void * aligned_ptr, * ptr;

        slong i, j, num, block_size, skip;

        flint_page_size = flint_get_page_size();
        block_size = PAGES_PER_BLOCK*flint_page_size;
        flint_page_mask = ~(flint_page_size - 1);

        /* get new block */
        ptr = flint_malloc(block_size + flint_page_size);

        /* align to page boundary */
        aligned_ptr = flint_align_ptr(ptr, flint_page_size);

        /* set free count to zero and determine if this is the main thread */
        ((fmpz_block_header_s *) ptr)->count = 0;
#if FLINT_USES_PTHREAD
        ((fmpz_block_header_s *) ptr)->thread = pthread_self();
#endif
        /* how many __mpz_structs worth are dedicated to header, per page */
        skip = (sizeof(fmpz_block_header_s) - 1)/sizeof(__mpz_struct) + 1;

        /* total number of number of __mpz_structs worth per page */
        num = flint_page_size/sizeof(__mpz_struct);

        flint_mpz_structs_per_block = PAGES_PER_BLOCK*(num - skip);

        for (i = 0; i < PAGES_PER_BLOCK; i++)
        {
            __mpz_struct * page_ptr = (__mpz_struct *)((slong) aligned_ptr + i*flint_page_size);

            /* set pointer in each page to start of entire block */
            ((fmpz_block_header_s *) page_ptr)->address = ptr;

            for (j = skip; j < num; j++)
            {
                mpz_init2(page_ptr + j, 2*FLINT_BITS);

                /*
                   Cannot be lifted from loop due to possibility of
                   gc calling _fmpz_clear_mpz during call to mpz_init_2
                */
                if (mpz_free_num >= mpz_free_alloc)
                {
                    mpz_free_alloc = FLINT_MAX(mpz_free_num + 1, mpz_free_alloc * 2);
                    mpz_free_arr = flint_realloc(mpz_free_arr, mpz_free_alloc * sizeof(__mpz_struct *));
                }

                mpz_free_arr[mpz_free_num++] = page_ptr + j;
            }
        }
    }

    return mpz_free_arr[--mpz_free_num];
}

void _fmpz_clear_mpz(fmpz f)
{
    __mpz_struct * ptr = COEFF_TO_PTR(f);

    /* check free count for block is zero, else this mpz came from a thread */
    fmpz_block_header_s * header_ptr = (fmpz_block_header_s *)((slong) ptr & flint_page_mask);

    header_ptr = (fmpz_block_header_s *) header_ptr->address;

    /* clean up if this is left over from another thread */
#if FLINT_USES_PTHREAD
    if (header_ptr->count != 0 || !pthread_equal(header_ptr->thread, pthread_self()))
#else
    if (header_ptr->count != 0)
#endif
    {
        int new_count;

        mpz_clear(ptr);

#if FLINT_USES_PTHREAD
       new_count = atomic_fetch_add(&(header_ptr->count), 1) + 1;
#else
       new_count = ++header_ptr->count;
#endif
        if (new_count == flint_mpz_structs_per_block)

            flint_free(header_ptr);
    } else
    {
        if (ptr->_mp_alloc > FLINT_MPZ_MAX_CACHE_LIMBS)
            mpz_realloc2(ptr, 2*FLINT_BITS);

        if (mpz_free_num == mpz_free_alloc)
        {
            mpz_free_alloc = FLINT_MAX(64, mpz_free_alloc * 2);
            mpz_free_arr = flint_realloc(mpz_free_arr, mpz_free_alloc * sizeof(__mpz_struct *));
        }

        mpz_free_arr[mpz_free_num++] = ptr;
    }
}

void _fmpz_cleanup_mpz_content(void)
{
    ulong i;

    for (i = 0; i < mpz_free_num; i++)
    {
       int new_count;
       fmpz_block_header_s * ptr;

       mpz_clear(mpz_free_arr[i]);

       /* update count of cleared mpz's for block */
       ptr = (fmpz_block_header_s *)((slong) mpz_free_arr[i] & ~(flint_page_size - 1));

       ptr = (fmpz_block_header_s *) ptr->address;

#if FLINT_USES_PTHREAD
       new_count = atomic_fetch_add(&(ptr->count), 1) + 1;
#else
       new_count = ++ptr->count;
#endif
       if (new_count == flint_mpz_structs_per_block)
          flint_free(ptr);
    }

    mpz_free_num = mpz_free_alloc = 0;
}

void _fmpz_cleanup(void)
{
    _fmpz_cleanup_mpz_content();
    flint_free(mpz_free_arr);
    mpz_free_arr = NULL;
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
