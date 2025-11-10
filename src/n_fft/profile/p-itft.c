#include "profiler.h"
#include "nmod_vec.h"
#include "fft_small.h"
#include "n_fft.h"

/* TODO
 * - add profiling for lazy variants that do not reduce to [0,n)
 **/

#define NUM_PRIMES 7

typedef struct
{
   ulong prime;
   ulong depth;
   ulong iolen;
} info_t;

#define SAMPLE(fun, _variant)                                      \
void sample_##fun##_variant(void * arg, ulong count)               \
{                                                                  \
    info_t * info = (info_t *) arg;                                \
    const ulong p = info->prime;                                   \
    const ulong depth = info->depth;                               \
    const ulong iolen = info->iolen;                               \
                                                                   \
    const ulong len = (UWORD(1) << depth);                         \
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));  \
                                                                   \
    /* modulus, roots of unity */                                  \
    n_fft_ctx_t F;                                                 \
    n_fft_ctx_init2(F, depth, p);                                  \
                                                                   \
    FLINT_TEST_INIT(state);                                        \
                                                                   \
    ulong * coeffs = _nmod_vec_init(len);                          \
    for (ulong k = 0; k < iolen; k++)                              \
        coeffs[k] = n_randint(state, p);                           \
                                                                   \
    for (ulong i = 0; i < count; i++)                              \
    {                                                              \
        prof_start();                                              \
        for (ulong j = 0; j < rep; j++)                            \
            fun##_variant(coeffs, iolen, F);                       \
        prof_stop();                                               \
    }                                                              \
                                                                   \
    _nmod_vec_clear(coeffs);                                       \
    n_fft_ctx_clear(F);                                            \
    FLINT_TEST_CLEAR(state);                                       \
}                                                                  \

SAMPLE(n_fft_itft, )

void sample_sd_ifft(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    const ulong p = info->prime;
    const ulong depth = info->depth;
    const ulong iolen = info->iolen;

    const ulong len = UWORD(1) << depth;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    sd_fft_ctx_t Q;
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_ctx_fit_depth(Q, depth);

    ulong sz = sd_fft_ctx_data_size(depth)*sizeof(double);

    FLINT_TEST_INIT(state);

    nmod_t mod;
    nmod_init(&mod, p);
    ulong * coeffs = _nmod_vec_init(iolen);
    _nmod_vec_rand(coeffs, state, iolen, mod);

    double* data = flint_aligned_alloc(4096, n_round_up(sz, 4096));
    for (ulong i = 0; i < iolen; i++)
        data[i] = coeffs[i];

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
            sd_ifft_trunc(Q, data, depth, iolen);
        prof_stop();
    }

    sd_fft_ctx_clear(Q);
    FLINT_TEST_CLEAR(state);
}

int main()
{
    flint_printf("- depth is log(fft length)\n");
    flint_printf("- timing TFT for several parameters\n");
    flint_printf("depth\t");
    for (ulong ii = 0; ii < 8; ii++)
        flint_printf("sd_ifft\t");
    for (ulong ii = 0; ii < 8; ii++)
        flint_printf("itft\t");
    flint_printf("\n");

    ulong primes[NUM_PRIMES] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        1073479681,          // 30 bits, 1 + 2**30 - 2**18 == 1 + 2**18 * (2**12 - 1)
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
        4611686018427322369  // 62 bits: 1 + 2**62 - 2**16 == 1 + 2**16 * (2**46 - 1)
    };
    ulong max_depths[NUM_PRIMES] = { 18, 18, 25, 25, 25, 25, 16 };

    for (ulong k = 4; k < 6; k++)
    {
        for (ulong depth = 4; depth <= max_depths[k]; depth++)
        {
            const ulong len = UWORD(1) << depth;
            const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

            info_t info;
            info.prime = primes[k];
            info.depth = depth;

            double min_sd[8];
            double min[8];
            double max;

            /* iolen at 8 points well distributed in (len/2, len] */
            ulong iolens[8] = {len/2 + 1*len/16,
                               len/2 + 2*len/16,
                               len/2 + 3*len/16,
                               len/2 + 4*len/16,
                               len/2 + 5*len/16,
                               len/2 + 6*len/16,
                               len/2 + 7*len/16,
                               len/2 + 8*len/16};
            /* make sure multiples of 4 */
            if (len == 16)
            {
                iolens[0] = iolens[1] = iolens[2] = iolens[3] = 12;
                iolens[4] = iolens[5] = iolens[6] = iolens[7] = 16;
            }
            if (len == 32)
            {
                iolens[0] = iolens[1] = 20;
                iolens[2] = iolens[3] = 24;
                iolens[4] = iolens[5] = 28;
                iolens[6] = iolens[7] = 32;
            }

            flint_printf("%ld\t", info.depth);

            for (ulong ii = 0; ii < 8; ii++)
            {
                info.iolen = iolens[ii];
                if (k < 5) prof_repeat(min_sd+ii, &max, sample_sd_ifft, (void *) &info);
                else min_sd[ii] = 0.;
                prof_repeat(min+ii, &max, sample_n_fft_itft, (void *) &info);
            }

            for (ulong ii = 0; ii < 8; ii++)
                flint_printf("%.1e\t", min_sd[ii]/(double)1000000/rep);
            for (ulong ii = 0; ii < 8; ii++)
                flint_printf("%.1e\t", min[ii]/(double)1000000/rep);
            flint_printf("\n");
        }
    }
    return 0;
}

/** 50 bit prime, commit "introduce_nmod_fft ????"
 *
 * Output on zen4 (AMD Ryzen 7 PRO 7840U), 50 bits
 * depth   sd_ifft sd_ifft sd_ifft sd_ifft sd_ifft sd_ifft sd_ifft sd_ifft itft    itft    itft    itft    itft    itft    itft    itft
 * 4       1.8e-08 1.8e-08 1.7e-08 1.7e-08 1.7e-08 1.8e-08 1.8e-08 1.7e-08 6.8e-08 6.5e-08 6.5e-08 6.5e-08 4.8e-08 4.6e-08 4.6e-08 4.6e-08
 * 5       2.5e-08 2.4e-08 2.4e-08 2.4e-08 2.4e-08 2.6e-08 2.4e-08 2.4e-08 1.0e-07 1.0e-07 1.4e-07 1.4e-07 1.6e-07 1.7e-07 1.1e-07 1.1e-07
 * 6       5.5e-08 5.4e-08 5.3e-08 5.4e-08 5.4e-08 5.5e-08 5.7e-08 5.4e-08 1.9e-07 2.2e-07 3.0e-07 2.7e-07 3.0e-07 3.4e-07 4.2e-07 2.4e-07
 * 7       1.0e-07 1.0e-07 1.0e-07 1.0e-07 1.0e-07 1.0e-07 1.0e-07 1.1e-07 3.8e-07 4.4e-07 6.0e-07 5.7e-07 6.2e-07 6.6e-07 8.5e-07 5.8e-07
 * 8       2.6e-07 2.6e-07 2.5e-07 2.5e-07 2.5e-07 2.5e-07 2.6e-07 2.5e-07 8.9e-07 9.9e-07 1.3e-06 1.2e-06 1.3e-06 1.4e-06 1.8e-06 1.2e-06
 * 9       5.2e-07 5.4e-07 5.0e-07 5.0e-07 5.1e-07 5.0e-07 5.0e-07 5.0e-07 1.8e-06 2.0e-06 2.6e-06 2.6e-06 2.8e-06 3.0e-06 3.7e-06 2.8e-06
 * 10      1.0e-06 9.8e-07 9.8e-07 9.7e-07 1.2e-06 1.2e-06 1.2e-06 1.2e-06 4.0e-06 4.6e-06 5.9e-06 5.6e-06 6.0e-06 6.6e-06 7.9e-06 5.5e-06
 * 11      1.8e-06 1.8e-06 2.1e-06 2.1e-06 2.5e-06 2.5e-06 2.7e-06 2.7e-06 8.0e-06 8.9e-06 1.1e-05 1.2e-05 1.2e-05 1.3e-05 1.6e-05 1.3e-05
 * 12      3.6e-06 4.1e-06 4.3e-06 4.6e-06 4.8e-06 5.1e-06 5.5e-06 5.8e-06 1.8e-05 2.0e-05 2.5e-05 2.5e-05 2.7e-05 2.9e-05 3.5e-05 2.6e-05
 * 13      8.1e-06 8.3e-06 9.1e-06 9.7e-06 1.1e-05 1.1e-05 1.2e-05 1.3e-05 3.7e-05 4.1e-05 5.1e-05 5.3e-05 5.6e-05 6.1e-05 7.3e-05 6.0e-05
 * 14      1.7e-05 1.7e-05 1.9e-05 2.0e-05 2.3e-05 2.4e-05 2.6e-05 2.7e-05 8.0e-05 9.1e-05 1.1e-04 1.1e-04 1.2e-04 1.3e-04 1.5e-04 1.2e-04
 * 15      3.5e-05 3.6e-05 4.0e-05 4.2e-05 5.2e-05 5.1e-05 5.3e-05 5.6e-05 1.6e-04 1.8e-04 2.2e-04 2.3e-04 2.5e-04 2.7e-04 3.2e-04 2.7e-04
 * 16      7.3e-05 7.7e-05 8.7e-05 9.1e-05 1.0e-04 1.0e-04 1.1e-04 1.2e-04 3.6e-04 4.0e-04 5.0e-04 4.9e-04 5.3e-04 5.7e-04 6.8e-04 5.4e-04
 * 17      1.6e-04 1.6e-04 1.8e-04 2.0e-04 2.1e-04 2.2e-04 2.4e-04 2.5e-04 7.2e-04 8.1e-04 9.9e-04 1.0e-03 1.1e-03 1.2e-03 1.4e-03 1.2e-03
 * 18      3.3e-04 3.5e-04 4.0e-04 4.1e-04 4.9e-04 4.7e-04 5.0e-04 5.2e-04 1.6e-03 1.8e-03 2.2e-03 2.2e-03 2.4e-03 2.6e-03 3.0e-03 2.4e-03
 * 19      7.2e-04 7.6e-04 8.3e-04 8.6e-04 9.6e-04 1.0e-03 1.1e-03 1.1e-03 3.4e-03 3.7e-03 4.5e-03 4.5e-03 4.9e-03 5.2e-03 6.1e-03 5.5e-03
 * 20      1.7e-03 1.7e-03 1.9e-03 1.9e-03 2.2e-03 2.3e-03 2.4e-03 2.4e-03 7.1e-03 8.1e-03 9.6e-03 9.8e-03 1.0e-02 1.1e-02 1.3e-02 1.1e-02
 * 21      4.2e-03 4.3e-03 4.5e-03 4.6e-03 5.1e-03 5.4e-03 5.6e-03 5.7e-03 1.5e-02 1.7e-02 2.0e-02 2.1e-02 2.2e-02 2.3e-02 2.7e-02 2.4e-02
 * 22      9.7e-03 1.0e-02 1.1e-02 1.1e-02 1.1e-02 1.2e-02 1.4e-02 1.4e-02 3.3e-02 3.7e-02 4.3e-02 4.4e-02 4.8e-02 5.1e-02 5.9e-02 4.7e-02
 * 23      2.1e-02 2.1e-02 2.4e-02 2.5e-02 2.8e-02 2.8e-02 2.8e-02 3.0e-02 6.8e-02 7.5e-02 8.9e-02 9.5e-02 1.0e-01 1.1e-01 1.2e-01 1.1e-01
 * 24      4.4e-02 4.4e-02 5.2e-02 4.9e-02 5.4e-02 6.0e-02 6.1e-02 6.4e-02 1.5e-01 1.7e-01 2.0e-01 2.0e-01 2.2e-01 2.3e-01 2.6e-01 2.2e-01
 * 25      7.8e-02 8.5e-02 9.1e-02 9.4e-02 1.1e-01 1.1e-01 1.2e-01 1.2e-01 2.9e-01 3.2e-01 3.8e-01 3.9e-01 4.2e-01 4.5e-01 5.1e-01 4.5e-01
 *
 *
 *
 * For comparison, NTL (60 bits, no avx):
 * depth   fft -->                                                 <--fft |ifft -->                                                <--ifft
 * 4       5.7e-08 5.7e-08 6.3e-08 6.1e-08 7.0e-08 6.9e-08 7.7e-08 6.4e-08 7.7e-08 7.3e-08 8.5e-08 7.1e-08 9.1e-08 9.7e-08 1.0e-07 7.2e-08
 * 5       9.2e-08 9.0e-08 1.1e-07 1.0e-07 1.2e-07 1.2e-07 1.3e-07 1.2e-07 1.2e-07 1.3e-07 1.4e-07 1.3e-07 1.7e-07 1.5e-07 1.7e-07 1.5e-07
 * 6       1.7e-07 1.9e-07 2.0e-07 2.0e-07 2.3e-07 2.4e-07 2.7e-07 2.6e-07 2.3e-07 2.4e-07 2.6e-07 2.8e-07 2.9e-07 3.0e-07 3.2e-07 3.0e-07
 * 7       3.7e-07 3.8e-07 4.1e-07 4.3e-07 4.8e-07 5.2e-07 5.4e-07 5.7e-07 4.5e-07 4.8e-07 5.5e-07 5.4e-07 6.1e-07 6.2e-07 6.6e-07 6.9e-07
 * 8       7.7e-07 8.1e-07 8.8e-07 9.5e-07 1.1e-06 1.1e-06 1.2e-06 1.3e-06 9.5e-07 1.1e-06 1.1e-06 1.1e-06 1.3e-06 1.3e-06 1.5e-06 1.5e-06
 * 9       1.7e-06 1.8e-06 2.0e-06 2.1e-06 2.3e-06 2.5e-06 2.6e-06 2.8e-06 2.2e-06 2.2e-06 2.4e-06 2.5e-06 2.8e-06 3.0e-06 3.1e-06 3.1e-06
 * 10      3.7e-06 4.0e-06 4.4e-06 4.6e-06 5.1e-06 5.3e-06 5.7e-06 6.5e-06 4.4e-06 4.7e-06 5.3e-06 5.4e-06 6.2e-06 6.5e-06 6.6e-06 6.8e-06
 * 11      8.0e-06 8.9e-06 9.5e-06 1.0e-05 1.1e-05 1.2e-05 1.3e-05 1.3e-05 9.2e-06 9.9e-06 1.1e-05 1.2e-05 1.3e-05 1.3e-05 1.4e-05 1.5e-05
 * 12      1.8e-05 1.8e-05 2.0e-05 2.2e-05 2.4e-05 2.6e-05 2.8e-05 2.9e-05 2.0e-05 2.2e-05 2.4e-05 2.6e-05 2.8e-05 2.9e-05 3.1e-05 3.4e-05
 * 13      3.7e-05 3.9e-05 4.4e-05 4.6e-05 5.3e-05 5.6e-05 5.9e-05 6.2e-05 4.3e-05 4.7e-05 5.1e-05 5.4e-05 5.9e-05 6.2e-05 6.9e-05 6.8e-05
 * 14      7.9e-05 8.6e-05 9.2e-05 1.0e-04 1.1e-04 1.1e-04 1.2e-04 1.3e-04 9.2e-05 9.7e-05 1.1e-04 1.1e-04 1.2e-04 1.4e-04 1.4e-04 1.4e-04
 * 15      1.7e-04 1.8e-04 2.0e-04 2.1e-04 2.3e-04 2.4e-04 2.7e-04 2.9e-04 1.9e-04 2.0e-04 2.2e-04 2.4e-04 2.7e-04 2.9e-04 2.9e-04 3.1e-04
 * 16      3.5e-04 4.0e-04 4.2e-04 4.5e-04 4.9e-04 5.3e-04 5.9e-04 6.0e-04 4.0e-04 4.5e-04 4.8e-04 5.3e-04 5.7e-04 5.9e-04 6.3e-04 6.6e-04
 * 17      7.9e-04 8.2e-04 8.9e-04 9.5e-04 1.1e-03 1.2e-03 1.2e-03 1.3e-03 8.6e-04 9.2e-04 1.1e-03 1.1e-03 1.2e-03 1.3e-03 1.3e-03 1.4e-03
 * 18      1.7e-03 1.8e-03 2.0e-03 2.1e-03 2.5e-03 2.5e-03 2.7e-03 2.8e-03 1.9e-03 2.2e-03 2.3e-03 2.4e-03 2.6e-03 2.8e-03 3.1e-03 3.1e-03
 * 19      3.6e-03 3.8e-03 4.2e-03 4.8e-03 4.9e-03 5.3e-03 5.5e-03 5.8e-03 4.4e-03 4.4e-03 4.7e-03 4.9e-03 5.4e-03 6.1e-03 6.1e-03 6.3e-03
 * 20      7.4e-03 8.5e-03 9.2e-03 9.4e-03 1.0e-02 1.1e-02 1.3e-02 1.3e-02 8.5e-03 9.0e-03 9.8e-03 1.1e-02 1.2e-02 1.2e-02 1.3e-02 1.3e-02
 * 21      1.7e-02 1.7e-02 1.9e-02 2.1e-02 2.4e-02 2.4e-02 2.6e-02 2.7e-02 1.7e-02 2.1e-02 2.1e-02 2.2e-02 2.4e-02 2.7e-02 2.9e-02 2.8e-02
 * 22      3.5e-02 3.9e-02 4.1e-02 4.4e-02 4.7e-02 5.3e-02 5.5e-02 5.7e-02 3.8e-02 4.0e-02 4.6e-02 4.8e-02 5.1e-02 5.4e-02 6.0e-02 6.1e-02
 * 23      7.9e-02 7.7e-02 8.3e-02 9.0e-02 1.0e-01 1.1e-01 1.1e-01 1.2e-01 7.6e-02 8.8e-02 9.1e-02 9.7e-02 1.1e-01 1.2e-01 1.2e-01 1.2e-01
 * 24      1.7e-01 1.7e-01 1.8e-01 2.0e-01 2.2e-01 2.2e-01 2.4e-01 2.6e-01 1.6e-01 1.7e-01 2.0e-01 2.0e-01 2.2e-01 2.4e-01 2.5e-01 2.7e-01
 */
