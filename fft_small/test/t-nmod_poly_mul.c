
#include "fft_small.h"
#include "profiler.h"


int main(void)
{
    flint_bitcnt_t nbits;
    mpn_ctx_t R;
    nmod_t mod;
    FLINT_TEST_INIT(state);

    flint_printf("nmod_poly_mul....");
    fflush(stdout);

    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    for (nbits = 1; nbits <= FLINT_BITS; nbits ++)
    {
        ulong * a, * b, * c, * d;
        ulong an, bn, zn, i, reps;

        nmod_init(&mod, n_randbits(state, 1 + n_randint(state, FLINT_BITS)));

        for (reps = 0; reps < 30; reps++)
        {
            an = 1 + n_randint(state, 10000);
            bn = 1 + n_randint(state, an);
            zn = an + bn - 1;

            a = FLINT_ARRAY_ALLOC(an, ulong);
            b = FLINT_ARRAY_ALLOC(bn, ulong);
            c = FLINT_ARRAY_ALLOC(zn, ulong);
            d = FLINT_ARRAY_ALLOC(zn, ulong);

            for (i = 0; i < an; i++)
                a[i] = n_randint(state, mod.n);

            for (i = 0; i < bn; i++)
                b[i] = n_randint(state, mod.n);

            _nmod_poly_mul(c, a, an, b, bn, mod);
            _mpn_ctx_nmod_poly_mul(R, d, a, an, b, bn, mod);

            for (i = 0; i < zn; i++)
            {
                if (c[i] != d[i])
                {
                    flint_printf("error at index %wu\n", i);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            flint_free(a);
            flint_free(b);
            flint_free(c);
            flint_free(d);
        }
    }

    mpn_ctx_clear(R);

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

