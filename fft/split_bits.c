/* 
    Copyright (C) 2009, 2011, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"
#include "thread_support.h"

typedef struct
{
   volatile slong * i;
   slong num;
   mp_size_t coeff_limbs;
   mp_size_t output_limbs;
   mp_srcptr limbs;
   mp_limb_t ** poly;
   pthread_mutex_t * mutex;
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
    slong i, end;

    while (1)
    {
        pthread_mutex_lock(arg.mutex);
        i = *arg.i;
        end = *arg.i = FLINT_MIN(i + 16, num);
        pthread_mutex_unlock(arg.mutex);

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
    pthread_mutex_t mutex;
    slong num_threads;
    thread_pool_handle * threads;
    split_limbs_arg_t * args;
    
    pthread_mutex_init(&mutex, NULL);

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
       args[i].mutex = &mutex;       
    }

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                                                _split_limbs_worker, &args[i]);

    _split_limbs_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

    flint_give_back_threads(threads, num_threads);

    flint_free(args);

    pthread_mutex_destroy(&mutex);

    i = num;
    skip = i*coeff_limbs;
   
    if (i < length) 
        flint_mpn_zero(poly[i], output_limbs + 1);
   
    if (total_limbs > skip) 
        flint_mpn_copyi(poly[i], limbs + skip, total_limbs - skip);
   
    return length;
}

mp_size_t fft_split_bits(mp_limb_t ** poly, mp_srcptr limbs, 
               mp_size_t total_limbs, flint_bitcnt_t bits, mp_size_t output_limbs)
{
   mp_size_t i, coeff_limbs, limbs_left, length = (FLINT_BITS*total_limbs - 1)/bits + 1;
   flint_bitcnt_t shift_bits, top_bits = ((FLINT_BITS - 1) & bits);
   mp_srcptr limb_ptr;
   mp_limb_t mask;
   
   if (top_bits == 0)
      return fft_split_limbs(poly, limbs, total_limbs, bits/FLINT_BITS, output_limbs);

   coeff_limbs = (bits/FLINT_BITS) + 1;
   mask = (WORD(1)<<top_bits) - WORD(1);
   shift_bits = WORD(0);
   limb_ptr = limbs;                      
    
#pragma omp parallel for private(i, limb_ptr, shift_bits)
   for (i = 0; i < length - 1; i++)
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
            poly[i][coeff_limbs - 1] += (limb_ptr[0] << (FLINT_BITS - (shift_bits - top_bits)));
            shift_bits -= FLINT_BITS; 
         }
         
         poly[i][coeff_limbs - 1] &= mask;
      } 
   }
   
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

