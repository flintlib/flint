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

int
main(void)
{
    mp_bitcnt_t bits;
    mp_size_t j, k, n, w, limbs, d;
    mp_limb_t * nn, * r;
    mpz_t p, m1, m2, mn1, mn2;

    flint_rand_t state;

    printf("mul_2expmod_2expp1....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    mpz_init(m1);
    mpz_init(m2);
    mpz_init(mn1);
    mpz_init(mn2);
    mpz_init(p);

    /* normalisation mod p = 2^wn + 1 where B divides nw and n is a power of 2 */
    for (bits = FLINT_BITS; bits < 16*FLINT_BITS; bits += FLINT_BITS)
    {
        for (j = 1; j < 32; j++)
        {
            for (k = 1; k <= FLINT_BITS; k <<= 1)
            {
                for (d = 0; d < FLINT_BITS; d++)
                {
                    n = bits/k;
                    w = j*k;
                    limbs = (n*w)/GMP_LIMB_BITS;
            
                    nn = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    r  = flint_malloc((limbs + 1)*sizeof(mp_limb_t));
                    random_fermat(nn, state, limbs);
                    fermat_to_mpz(mn1, nn, limbs);
                    set_p(p, n, w);
            
                    mpn_mul_2expmod_2expp1(r, nn, limbs, d);
                    fermat_to_mpz(m2, r, limbs);
                    mpz_mod(m2, m2, p);
                    
                    mpz_mod(m1, mn1, p);
                    mpz_mul_2exp(m1, m1, d);
                    mpz_mod(m1, m1, p);
                    
                    if (mpz_cmp(m1, m2) != 0)
                    {
                        printf("FAIL:\n");
                        printf("mpn_mul_2expmod_2expp1 error\n");
                        gmp_printf("want %Zx\n\n", m1);
                        gmp_printf("got  %Zx\n", m2);
                        abort();
                    }

                    flint_free(nn);
                    flint_free(r);
                }
            }
        }
    }

    mpz_clear(mn2);
    mpz_clear(mn1);
    mpz_clear(m2);
    mpz_clear(m1);
    mpz_clear(p);

    flint_randclear(state);
    
    printf("PASS\n");
    return 0;
}
