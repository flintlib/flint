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

    Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "profiler.h"

#define NUMEX (10)
#define BITSLOW (50)
#define BITSINC (100)
#define BITSHIGH (2000)
#define MINCPU (10)
#define TOLERANCE (1.05)

#define TIME_INNER(ALG, TIMEVAR) do \
{ \
    /* warm-up */ \
    for (i=0, curseq=seqs; i<NUMEX; ++i, curseq += len) \
    { \
        ALG(polys+i, curseq, len); \
    } \
    \
    loops = 1; \
    do { \
        timeit_start(timer); \
        \
        for (i=0, curseq=seqs; i<NUMEX; ++i, curseq += len) \
        { \
            for (j=0; j<loops; ++j) ALG(polys+i, curseq, len); \
        } \
        \
        timeit_stop(timer); \
        loops *= 2; \
    } \
    while (timer->cpu <= MINCPU); \
    \
    /* cool-down */ \
    for (i=0, curseq=seqs; i<NUMEX; ++i, curseq += len) \
    { \
        ALG(polys+i, curseq, len); \
    } \
    \
    (TIMEVAR) = ((double)timer->cpu) / (NUMEX * loops); \
} while(0)

void time_algs(double* times, flint_rand_t state, const fmpz* primes, slong len)
{
    fmpz_mod_poly_struct polys[NUMEX];
    fmpz *seqs, *curseq, *gen;
    slong genlen = len/2, i, j, loops;
    timeit_t timer;

    seqs = _fmpz_vec_init(len*NUMEX);
    gen = _fmpz_vec_init(genlen);

    /* initialize polys and sequences */
    for (i=0, curseq = seqs; i<NUMEX; ++i, curseq += len)
    {
        fmpz_mod_poly_init2(polys+i, primes+i, len+1);

        for (j=0; j < genlen; ++j) 
        {
            fmpz_randm(gen + j, state, primes+i);
            fmpz_randm(curseq + j, state, primes+i);
        }

        for (; j < len; ++j)
        {
            _fmpz_vec_dot(curseq+j, curseq+(j-genlen), gen, genlen);
            fmpz_mod(curseq+j, curseq+j, primes+i);
        }
    }

    _fmpz_vec_clear(gen, genlen);
    
    TIME_INNER(fmpz_mod_poly_minpoly_bm, times[0]);

    TIME_INNER(fmpz_mod_poly_minpoly_hgcd, times[1]);

    _fmpz_vec_clear(seqs, len*NUMEX);
    for (i=0; i<NUMEX; ++i) fmpz_mod_poly_clear(polys+i);
}
    
int main(void)
{
    mp_bitcnt_t bits;
    fmpz primes[NUMEX];
    double times[2];
    slong i, len, len2, lower, upper;
    
    FLINT_TEST_INIT(state);

    for (i=0; i<NUMEX; ++i) fmpz_init(primes+i);

    for (bits = BITSLOW; bits < BITSHIGH; bits += BITSINC)
    {
        for (i=0; i<NUMEX; ++i) fmpz_randprime(primes+i, state, bits, 0);

        lower = -1;
        len = 2;
        while (1)
        {
            time_algs(times, state, primes, len);
            if (times[1] > TOLERANCE*times[0]) lower = len;
            else if (lower > 0 && times[0] > TOLERANCE*times[1]) break;
            len = len*3/2;
        }
        upper = len;

        while (upper - lower + 1 > 3)
        {
            len = lower + (upper-lower)/3;
            len2 = lower + (upper-lower)*2/3;

            time_algs(times, state, primes, len);
            if (times[1] > TOLERANCE*times[0]) lower = len;
            else if (times[0] > TOLERANCE*times[1])
            {
                upper = len;
                continue;
            }

            time_algs(times, state, primes, len2);
            if (times[1] > TOLERANCE*times[0]) lower = len2;
            else if (times[0] > TOLERANCE*times[1]) upper = len2;
            else if (lower < len) break;
        }

        flint_printf("bits: %wd, lower: %wd, upper: %wd\n", bits, lower, upper);
        fflush(stdout);
    }

    for (i=0; i<NUMEX; ++i) fmpz_clear(primes+i);

    FLINT_TEST_CLEANUP(state);

    return 0;
}
