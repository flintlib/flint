/* 

Copyright 2009, 2011 William Hart. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY William Hart ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL William Hart OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of William Hart.

*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fft.h"

/* set p = 2^wn + 1 */
void set_p(mpz_t p, mp_size_t n, mp_bitcnt_t w)
{
   mpz_set_ui(p, 1);
   mpz_mul_2exp(p, p, n*w);
   mpz_add_ui(p, p, 1);
}

void ref_fft_butterfly_twiddle(mpz_t s, mpz_t t, mpz_t i1, mpz_t i2, 
                   mpz_t p, mp_size_t i, mp_size_t w, mp_bitcnt_t b1, mp_bitcnt_t b2)
{
   mpz_add(s, i1, i2);
   mpz_sub(t, i1, i2);
   mpz_mul_2exp(s, s, b1);
   mpz_mul_2exp(t, t, b2);
   mpz_mod(s, s, p);
   mpz_mod(t, t, p);
}

void ref_ifft_butterfly_twiddle(mpz_t s, mpz_t t, mpz_t i1, mpz_t i2, 
      mpz_t p, mp_size_t i, mp_size_t n, mp_size_t w, mp_bitcnt_t b1, mp_bitcnt_t b2)
{
   mpz_mul_2exp(i1, i1, 2*n*w - b1);
   mpz_mul_2exp(i2, i2, 2*n*w - b2);
   mpz_add(s, i1, i2);
   mpz_sub(t, i1, i2);
   mpz_mod(s, s, p);
   mpz_mod(t, t, p);
}

int
main(void)
{
    mp_size_t c, bits, j, k, n, w, limbs;
    mpz_t p, ma, mb, m2a, m2b, mn1, mn2;
    mp_limb_t * nn1, * nn2, * r1, * r2;
    mp_bitcnt_t b1, b2;
   
    flint_rand_t state;

    printf("fft/ifft_butterfly_twiddle....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    mpz_init(p);
    mpz_init(ma);
    mpz_init(mb);
    mpz_init(m2a);
    mpz_init(m2b);
    mpz_init(mn1);
    mpz_init(mn2);
   
    for (bits = FLINT_BITS; bits < 20*FLINT_BITS; bits += FLINT_BITS)
    {
        for (j = 1; j < 10; j++)
        {
            for (k = 1; k <= FLINT_BITS; k <<= 1)
            {
                n = bits/k;
                w = j*k;

                limbs = (n*w)/FLINT_BITS;

                for (c = 0; c < n; c++)
                {
                    b1 = n_randint(state, n*w);
                    b2 = n_randint(state, n*w);
                    nn1 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    nn2 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    r1 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    r2 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    random_fermat(nn1, state, limbs);
                    random_fermat(nn2, state, limbs);
                     
                    fermat_to_mpz(mn1, nn1, limbs);
                    fermat_to_mpz(mn2, nn2, limbs);
                    set_p(p, n, w);
            
                    fft_butterfly_twiddle(r1, r2, nn1, nn2, limbs, b1, b2);
                    fermat_to_mpz(m2a, r1, limbs);
                    fermat_to_mpz(m2b, r2, limbs);
                    
                    mpz_mod(m2a, m2a, p);
                    mpz_mod(m2b, m2b, p);
                    ref_fft_butterfly_twiddle(ma, mb, mn1, mn2, p, c, w, b1, b2);

                    if (mpz_cmp(ma, m2a) != 0)
                    {
                        printf("FAIL:\n");
                        printf("fft_butterfly_twiddle error a\n");
                        printf("limbs = %ld\n", limbs);
                        printf("n = %ld, w = %ld, k = %ld, c = %ld\n", n, w, k, c);
                        gmp_printf("want %Zx\n\n", ma);
                        gmp_printf("got  %Zx\n", m2a);
                        abort();
                    }
                    if (mpz_cmp(mb, m2b) != 0)
                    {
                        printf("FAIL:\n");
                        printf("fft_butterfly_twiddle error b\n");
                        printf("limbs = %ld\n", limbs);
                        printf("n = %ld, w = %ld, k = %ld, c = %ld\n", n, w, k, c);
                        gmp_printf("want %Zx\n\n", mb);
                        gmp_printf("got  %Zx\n", m2b);
                        abort();
                    }
                    
                    flint_free(nn1);
                    flint_free(nn2);
                    flint_free(r1);
                    flint_free(r2);
                }
            }
        }
    }

    for (bits = FLINT_BITS; bits < 20*FLINT_BITS; bits += FLINT_BITS)
    {
        for (j = 1; j < 10; j++)
        {
            for (k = 1; k <= FLINT_BITS; k <<= 1)
            {
                n = bits/k;
                w = j*k;

                limbs = (n*w)/FLINT_BITS;

                for (c = 0; c < n; c++)
                {
                    b1 = n_randint(state, n*w);
                    b2 = n_randint(state, n*w);
                    nn1 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    nn2 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    r1 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    r2 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    random_fermat(nn1, state, limbs);
                    random_fermat(nn2, state, limbs);
                     
                    fermat_to_mpz(mn1, nn1, limbs);
                    fermat_to_mpz(mn2, nn2, limbs);
                    set_p(p, n, w);
            
                    ifft_butterfly_twiddle(r1, r2, nn1, nn2, limbs, b1, b2);
                    fermat_to_mpz(m2a, r1, limbs);
                    fermat_to_mpz(m2b, r2, limbs);
                    
                    mpz_mod(m2a, m2a, p);
                    mpz_mod(m2b, m2b, p);
                    ref_ifft_butterfly_twiddle(ma, mb, mn1, mn2, p, c, n, w, b1, b2);

                    if (mpz_cmp(ma, m2a) != 0)
                    {
                        printf("FAIL:\n");
                        printf("ifft_butterfly_twiddle error a\n");
                        printf("limbs = %ld\n", limbs);
                        printf("n = %ld, w = %ld, k = %ld, c = %ld\n", n, w, k, c);
                        gmp_printf("want %Zx\n\n", ma);
                        gmp_printf("got  %Zx\n", m2a);
                        abort();
                    }
                    if (mpz_cmp(mb, m2b) != 0)
                    {
                        printf("FAIL:\n");
                        printf("ifft_butterfly_twiddle error b\n");
                        printf("limbs = %ld\n", limbs);
                        printf("n = %ld, w = %ld, k = %ld, c = %ld\n", n, w, k, c);
                        gmp_printf("want %Zx\n\n", mb);
                        gmp_printf("got  %Zx\n", m2b);
                        abort();
                    }
                    
                    flint_free(nn1);
                    flint_free(nn2);
                    flint_free(r1);
                    flint_free(r2);
                }
            }
        }
    }

    mpz_clear(p);
    mpz_clear(ma);
    mpz_clear(mb);
    mpz_clear(m2a);
    mpz_clear(m2b);
    mpz_clear(mn1);
    mpz_clear(mn2);

    flint_randclear(state);
    
    printf("PASS\n");
    return 0;
}
