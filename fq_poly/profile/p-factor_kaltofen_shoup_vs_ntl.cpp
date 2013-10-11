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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "profiler.h"
#include "NTL-interface.h"

#include <NTL/ZZ_pEXFactoring.h>

#define nalgs 2
#define ncases 1
#define cpumin 10

int
main(int argc, char** argv)
{
    slong s[nalgs];

    int c, n, len, ext, reps = 0;
    fmpz_t p, temp;
    ZZ prime;
    ZZ_pX mod;
    fq_poly_t f;
    fq_ctx_t ctx;
    fmpz_mod_poly_t fmod;
    
    flint_rand_t state;
    flint_randinit(state);

    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    fmpz_init(temp);
       
    fmpz_set_str(temp, argv[2], 10);
    ext = fmpz_get_si(temp);

    fmpz_set_str(temp, argv[3], 10);
    len = fmpz_get_si(temp);

    // Set the prime
    fmpz_get_ZZ(prime, p);
    ZZ_p::init(prime);
    
    // Set the modulus
    BuildIrred(mod, ext + 1);
    fmpz_mod_poly_init(fmod, p);
    //fmpz_mod_poly_randtest_not_zero(fmod, state, ext + 1);
    //fmpz_mod_poly_set_coeff_ui(fmod, ext, 1);
    fmpz_mod_poly_set_ZZ_pX(fmod, mod);
      
    ZZ_pE::init(mod);  

    fq_ctx_init_modulus(ctx, p, ext, fmod, "a");
    
    fq_poly_init(f);

    for (c = 0; c < nalgs; c++)
    {
        s[c] = 0L;
    }
       
    for (n = 0; n < ncases; n++)
    {
        timeit_t t[nalgs];
        int l, loops = 1;
        ZZ_pEX nf;
        vec_pair_ZZ_pEX_long v;
        fq_poly_factor_t res;

        /*
           Construct random elements of fq
        */
        {
            fq_poly_randtest(f, state, len, ctx);
            fq_poly_make_monic(f, f, ctx);
            fq_poly_get_ZZ_pEX(nf, f);
        }
        
    loop:
        fflush(stdout);
        timeit_start(t[0]);
        for (l = 0; l < loops; l++)
            fq_poly_factor_init(res, ctx);
            fq_poly_factor_kaltofen_shoup(res, f, ctx);
            fq_poly_factor_clear(res);
        timeit_stop(t[0]);

        timeit_start(t[1]);
        for (l = 0; l < loops; l++)
            v = CanZass(nf);
        timeit_stop(t[1]);

        for (c = 0; c < nalgs; c++)
            if (t[c]->cpu <= cpumin)
            {
                flint_printf("%d ", t[c]->cpu);
                loops *= 2;
                goto loop;
            }
        
        for (c = 0; c < nalgs; c++)
            s[c] += t[c]->cpu;
        reps += loops;
    }
        
    for (c = 0; c < nalgs; c++)
    {
        flint_printf("%20f ", s[c] / (double) reps);
        fflush(stdout);
    }
    printf("\n");
        
    fq_poly_clear(f);
    fmpz_mod_poly_clear(fmod);
    fq_ctx_clear(ctx);
    fmpz_clear(p);
    fmpz_clear(temp);

    flint_randclear(state);
}
