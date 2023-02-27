/* 
    Copyright (C) 2009, 2011, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"
#include "thread_support.h"

typedef struct
{
    volatile mp_size_t * i;
    slong num;
    mp_size_t coeff_limbs;
    mp_size_t output_limbs;
    mp_srcptr limbs;
    mp_limb_t ** poly;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
}
split_limbs_arg_t;

void
_split_limbs_worker(void * arg_ptr)
{
    split_limbs_arg_t arg = *((split_limbs_arg_t *) arg_ptr);
    slong num = arg.num;
    mp_size_t skip;
    mp_size_t coeff_limbs = arg.coeff_limbs;
    mp_size_t output_limbs = arg.output_limbs;
    mp_srcptr limbs = arg.limbs;
    mp_limb_t ** poly = arg.poly;
    mp_size_t i, end;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
	i = *arg.i;
        end = *arg.i = FLINT_MIN(i + 16, num);
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (i >= num)
            return;

        for ( ; i < end; i++)
        {
           skip = i*coeff_limbs;

           flint_mpn_zero(poly[i], output_limbs + 1);
           flint_mpn_copyi(poly[i], limbs + skip, coeff_limbs);
        }
    }
}

mp_size_t fft_split_limbs(mp_limb_t ** poly, mp_srcptr limbs, 
          mp_size_t total_limbs, mp_size_t coeff_limbs, mp_size_t output_limbs)
{
    mp_size_t i, shared_i = 0, skip, length = (total_limbs - 1)/coeff_limbs + 1;
    mp_size_t num = total_limbs/coeff_limbs;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong num_threads;
    thread_pool_handle * threads;
    split_limbs_arg_t * args;
    
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    num_threads = flint_request_threads(&threads,
                            FLINT_MIN(flint_get_num_threads(), (num + 15)/16));

    args = (split_limbs_arg_t *)
                     flint_malloc(sizeof(split_limbs_arg_t)*(num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
       args[i].i = &shared_i;
       args[i].num = num;
       args[i].coeff_limbs = coeff_limbs;
       args[i].output_limbs = output_limbs;
       args[i].limbs = limbs;
       args[i].poly = poly;
#if FLINT_USES_PTHREAD
       args[i].mutex = &mutex;       
#endif
    }

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                                                _split_limbs_worker, &args[i]);

    _split_limbs_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    flint_give_back_threads(threads, num_threads);

    flint_free(args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    i = num;
    skip = i*coeff_limbs;
   
    if (i < length) 
        flint_mpn_zero(poly[i], output_limbs + 1);
   
    if (total_limbs > skip) 
        flint_mpn_copyi(poly[i], limbs + skip, total_limbs - skip);
   
    return length;
}

typedef struct
{
    volatile mp_size_t * i;
    slong length;
    mp_size_t coeff_limbs;
    mp_size_t output_limbs;
    mp_srcptr limbs;
    flint_bitcnt_t top_bits;
    mp_limb_t mask;
    mp_limb_t ** poly;
#if FLINT_USES_PTHREAD
    pthread_mutex_t * mutex;
#endif
}
split_bits_arg_t;

void
_split_bits_worker(void * arg_ptr)
{
    split_bits_arg_t arg = *((split_bits_arg_t *) arg_ptr);
    slong length = arg.length;
    mp_size_t coeff_limbs = arg.coeff_limbs;
    mp_size_t output_limbs = arg.output_limbs;
    mp_srcptr limbs = arg.limbs;
    flint_bitcnt_t top_bits = arg.top_bits;
    mp_limb_t mask = arg.mask;
    mp_limb_t ** poly = arg.poly;
    flint_bitcnt_t shift_bits;
    mp_srcptr limb_ptr;
    mp_size_t i, end;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
	i = *arg.i;
        end = *arg.i = FLINT_MIN(i + 16, length - 1);
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (i >= length - 1)
            return;

        for ( ; i < end; i++)
        {
            flint_mpn_zero(poly[i], output_limbs + 1);
      
            limb_ptr = limbs + i*(coeff_limbs - 1) + (i*top_bits)/FLINT_BITS;
            shift_bits = (i*top_bits) % FLINT_BITS;

            if (!shift_bits)
            {
                flint_mpn_copyi(poly[i], limb_ptr, coeff_limbs);
                poly[i][coeff_limbs - 1] &= mask;
                limb_ptr += (coeff_limbs - 1);
                shift_bits += top_bits;
            } else
            {
                mpn_rshift(poly[i], limb_ptr, coeff_limbs, shift_bits);
                limb_ptr += (coeff_limbs - 1);
                shift_bits += top_bits;

                if (shift_bits >= FLINT_BITS)
                {
                   limb_ptr++;
                   poly[i][coeff_limbs - 1] +=
                       (limb_ptr[0] << (FLINT_BITS - (shift_bits - top_bits)));
                   shift_bits -= FLINT_BITS; 
                }
         
                poly[i][coeff_limbs - 1] &= mask;
            } 
        }
    }
}

mp_size_t fft_split_bits(mp_limb_t ** poly, mp_srcptr limbs, 
               mp_size_t total_limbs, flint_bitcnt_t bits, mp_size_t output_limbs)
{
    mp_size_t i, shared_i = 0, coeff_limbs, limbs_left;
    mp_size_t length = (FLINT_BITS*total_limbs - 1)/bits + 1;
    flint_bitcnt_t shift_bits, top_bits = ((FLINT_BITS - 1) & bits);
    mp_srcptr limb_ptr;
    mp_limb_t mask;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong num_threads;
    thread_pool_handle * threads;
    split_bits_arg_t * args;
   
    if (top_bits == 0)
        return fft_split_limbs(poly, limbs, total_limbs, bits/FLINT_BITS, output_limbs);

    coeff_limbs = (bits/FLINT_BITS) + 1;
    mask = (WORD(1)<<top_bits) - WORD(1);
    shift_bits = WORD(0);
    limb_ptr = limbs;                      
    
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    num_threads = flint_request_threads(&threads,
                     FLINT_MIN(flint_get_num_threads(), (length - 1 + 15)/16));

    args = (split_bits_arg_t *)
                      flint_malloc(sizeof(split_bits_arg_t)*(num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
       args[i].i = &shared_i;
       args[i].length = length;
       args[i].coeff_limbs = coeff_limbs;
       args[i].output_limbs = output_limbs;
       args[i].limbs = limbs;
       args[i].top_bits = top_bits;
       args[i].mask = mask;
       args[i].poly = poly;
#if FLINT_USES_PTHREAD
       args[i].mutex = &mutex;       
#endif
    }

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                                                 _split_bits_worker, &args[i]);

    _split_bits_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    flint_give_back_threads(threads, num_threads);

    flint_free(args);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    i = length - 1;
    limb_ptr = limbs + i*(coeff_limbs - 1) + (i*top_bits)/FLINT_BITS;
    shift_bits = (i*top_bits) % FLINT_BITS;

    flint_mpn_zero(poly[i], output_limbs + 1);
   
    limbs_left = total_limbs - (limb_ptr - limbs);
   
    if (!shift_bits)
        flint_mpn_copyi(poly[i], limb_ptr, limbs_left);
    else
        mpn_rshift(poly[i], limb_ptr, limbs_left, shift_bits);                   
     
    return length;
}

