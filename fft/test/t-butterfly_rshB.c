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
#include "mpn_extras.h"
#include "fft.h"

/* set p = 2^wn + 1 */
void set_p(mpz_t p, mp_size_t n, mp_bitcnt_t w)
{
   mpz_set_ui(p, 1);
   mpz_mul_2exp(p, p, n*w);
   mpz_add_ui(p, p, 1);
}

void ref_butterfly_rshB(mpz_t t, mpz_t u, mpz_t i1, mpz_t i2, 
                                 mpz_t p, mp_size_t x, mp_size_t y)
{
   mpz_t mult1, mult2;

   mpz_init(mult1);
   mpz_init(mult2);

   mpz_set_ui(mult1, 1);
   mpz_mul_2exp(mult1, mult1, x*FLINT_BITS);
   mpz_invert(mult1, mult1, p);
   mpz_set_ui(mult2, 1);
   mpz_mul_2exp(mult2, mult2, y*FLINT_BITS);
   mpz_invert(mult2, mult2, p);
   mpz_mul(mult1, mult1, i1);
   mpz_mul(mult2, mult2, i2);
   mpz_add(t, mult1, mult2);
   mpz_sub(u, mult1, mult2);
   mpz_mod(t, t, p);
   mpz_mod(u, u, p);

   mpz_clear(mult1);
   mpz_clear(mult2);
}

int
main(void)
{
    mp_size_t c, bits, j, k, x, y, n, w, limbs;
    mpz_t p, ma, mb, m2a, m2b, mn1, mn2;
    mp_limb_t * nn1, * nn2, * r1, * r2;
   
    flint_rand_t state;

    printf("butterfly_rshB....");
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

                for (c = 0; c < limbs; c++)
                {
                    x = n_randint(state, limbs);
                    y = n_randint(state, limbs);
                    
                    nn1 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    nn2 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    r1 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    r2 = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    flint_mpn_rrandom(nn1, state->gmp_state, limbs);
                    random_fermat(nn1, state, limbs);
                    random_fermat(nn2, state, limbs);
                     
                    fermat_to_mpz(mn1, nn1, limbs);
                    fermat_to_mpz(mn2, nn2, limbs);
                    set_p(p, n, w);
            
                    butterfly_rshB(r1, r2, nn1, nn2, limbs, x, y);
                    fermat_to_mpz(m2a, r1, limbs);
                    fermat_to_mpz(m2b, r2, limbs);
                    
                    mpz_mod(m2a, m2a, p);
                    mpz_mod(m2b, m2b, p);
                    ref_butterfly_rshB(ma, mb, mn1, mn2, p, x, y);

                    if (mpz_cmp(ma, m2a) != 0)
                    {
                        printf("FAIL:\n");
                        printf("butterfly_rshB error a\n");
                        printf("x = %ld, y = %ld, limbs = %ld\n", x, y, limbs);
                        printf("n = %ld, w = %ld, k = %ld, c = %ld\n", n, w, k, c);
                        gmp_printf("want %Zx\n\n", ma);
                        gmp_printf("got  %Zx\n", m2a);
                        abort();
                    }
                    if (mpz_cmp(mb, m2b) != 0)
                    {
                        printf("FAIL:\n");
                        printf("butterfly_rshB error b\n");
                        printf("x = %ld, y = %ld, limbs = %ld\n", x, y, limbs);
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
