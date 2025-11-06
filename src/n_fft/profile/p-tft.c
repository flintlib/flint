#include "profiler.h"
#include "nmod_vec.h"
#include "fft_small.h"
#include "n_fft.h"
#include "n_fft/impl.h"

/* TODO
 * - add profiling for lazy variants that do not reduce to [0,n)
 **/

#define NUM_PRIMES 7

typedef struct
{
   ulong prime;
   ulong depth;
   ulong ilen;
   ulong olen;
} info_t;

#define SAMPLE(fun, _variant)                                      \
void sample_##fun##_variant(void * arg, ulong count)               \
{                                                                  \
    info_t * info = (info_t *) arg;                                \
    const ulong p = info->prime;                                   \
    const ulong depth = info->depth;                               \
    const ulong ilen = info->ilen;                                 \
    const ulong olen = info->olen;                                 \
                                                                   \
    const ulong len = (UWORD(1) << depth);                         \
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));  \
                                                                   \
    /* modulus, roots of unity */                                  \
    n_fft_ctx_t F;                                                 \
    n_fft_ctx_init2(F, depth, p);                                  \
    n_fft_args_t Fargs;                                            \
    n_fft_set_args(Fargs, F->mod, F->tab_w);                       \
                                                                   \
    FLINT_TEST_INIT(state);                                        \
                                                                   \
    ulong * coeffs = _nmod_vec_init(FLINT_MAX(ilen, len));         \
    for (ulong k = 0; k < ilen; k++)                               \
        coeffs[k] = n_randint(state, p);                           \
                                                                   \
    for (ulong i = 0; i < count; i++)                              \
    {                                                              \
        prof_start();                                              \
        for (ulong j = 0; j < rep; j++)                            \
            fun##_variant(coeffs, ilen, olen, Fargs);              \
        prof_stop();                                               \
    }                                                              \
                                                                   \
    _nmod_vec_clear(coeffs);                                       \
    n_fft_ctx_clear(F);                                            \
    FLINT_TEST_CLEAR(state);                                       \
}                                                                  \

#define SAMPLE_OLEN(fun, _variant)                                               \
void sample_##fun##_variant(void * arg, ulong count)                             \
{                                                                                \
    info_t * info = (info_t *) arg;                                              \
    const ulong p = info->prime;                                                 \
    const ulong depth = info->depth;                                             \
    const ulong olen = info->olen;                                               \
                                                                                 \
    const ulong len = (UWORD(1) << depth);                                       \
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));                \
                                                                                 \
    /* modulus, roots of unity */                                                \
    n_fft_ctx_t F;                                                               \
    n_fft_ctx_init2(F, depth, p);                                                \
    n_fft_args_t Fargs;                                                          \
    n_fft_set_args(Fargs, F->mod, F->tab_w);                                     \
                                                                                 \
    FLINT_TEST_INIT(state);                                                      \
                                                                                 \
    ulong * coeffs = _nmod_vec_init(len);                                        \
    for (ulong k = 0; k < len; k++)                                              \
        coeffs[k] = n_randint(state, p);                                         \
                                                                                 \
    for (ulong i = 0; i < count; i++)                                            \
    {                                                                            \
        prof_start();                                                            \
        for (ulong j = 0; j < rep; j++)                                          \
            fun##_variant(coeffs, olen, depth, 0, Fargs);                        \
        prof_stop();                                                             \
    }                                                                            \
                                                                                 \
    _nmod_vec_clear(coeffs);                                                     \
    n_fft_ctx_clear(F);                                                          \
    FLINT_TEST_CLEAR(state);                                                     \
}                                                                                \

SAMPLE(tft_lazy_1_4, )

SAMPLE_OLEN(tft_node_lazy_4_4, )

void sample_sd_fft(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    const ulong p = info->prime;
    const ulong depth = info->depth;
    ulong ilen = info->ilen;
    const ulong olen = info->olen;

    const ulong len = UWORD(1) << depth;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    if (ilen > len)
        ilen = len;

    sd_fft_ctx_t Q;
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_ctx_fit_depth(Q, depth);

    ulong sz = sd_fft_ctx_data_size(depth)*sizeof(double);

    FLINT_TEST_INIT(state);

    nmod_t mod;
    nmod_init(&mod, p);
    ulong * coeffs = _nmod_vec_init(ilen);
    _nmod_vec_rand(coeffs, state, ilen, mod);

    double* data = flint_aligned_alloc(4096, n_round_up(sz, 4096));
    for (ulong i = 0; i < ilen; i++)
        data[i] = coeffs[i];

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
            sd_fft_trunc(Q, data, depth, ilen, olen);
        prof_stop();
    }

    sd_fft_ctx_clear(Q);
    FLINT_TEST_CLEAR(state);
}

int main()
{
    flint_printf("- depth is log(fft length)\n");
    flint_printf("- timing TFT for several parameters\n");
    flint_printf("depth\tilen\tsd_fft\tsd_fft\tsd_fft\tsd_fft\ttft\ttft\ttft\ttft\n");

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
        for (ulong depth = 5; depth <= max_depths[k]; depth++)
        {
            const ulong len = UWORD(1) << depth;
            const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

            info_t info;
            info.prime = primes[k];
            info.depth = depth;

            double min_sd[5];
            double min[5];
            double min_olen[5];
            double max;

            /* olen in {len/2 + 8, len/2 + len/4, len} */
            ulong ilens[5] = {len/10, len/2 - 4, len/2 + 4, len, 3*len};
            ulong olens[4] = {len/2 + 4, 3*len/4, len - 4, len};
            for (ulong ili = 0; ili < 5; ili++)
            {
                info.ilen = FLINT_MAX(8, 4 * (ilens[ili] / 4));
                flint_printf("%ld\t%ld\t", info.depth, info.ilen);

                for (ulong oli = 0; oli < 4; oli++)
                {
                    info.olen = olens[oli];
                    if (k < 5) prof_repeat(min_sd+oli, &max, sample_sd_fft, (void *) &info);
                    prof_repeat(min+oli, &max, sample_tft_lazy_1_4, (void *) &info);
                    prof_repeat(min_olen+oli, &max, sample_tft_node_lazy_4_4, (void *) &info);
                }

                flint_printf("%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\n",
                        min_sd[0]/(double)1000000/rep,
                        min_sd[1]/(double)1000000/rep,
                        min_sd[2]/(double)1000000/rep,
                        min_sd[3]/(double)1000000/rep,
                        min[0]/(double)1000000/rep,
                        min[1]/(double)1000000/rep,
                        min[2]/(double)1000000/rep,
                        min[3]/(double)1000000/rep,
                        min_olen[0]/(double)1000000/rep,
                        min_olen[1]/(double)1000000/rep,
                        min_olen[2]/(double)1000000/rep,
                        min_olen[3]/(double)1000000/rep
                        );
            }
        }
    }
    return 0;
}

/** 50 bit prime, commit "introduce_nmod_fft ????"
 *
 * Output on zen4 (AMD Ryzen 7 PRO 7840U)
 * 
 * depth   ilen    sd_fft  sd_fft  sd_fft  sd_fft  tft     tft     tft     tft
 * 5       8       2.2e-08 2.2e-08 2.2e-08 2.2e-08 5.0e-08 5.5e-08 6.7e-08 7.5e-08 6.5e-08 7.0e-08 7.7e-08 8.7e-08
 * 5       12      2.3e-08 2.2e-08 2.2e-08 2.2e-08 4.1e-08 4.6e-08 5.6e-08 6.2e-08 6.5e-08 7.0e-08 7.8e-08 8.2e-08
 * 5       20      2.4e-08 2.3e-08 2.4e-08 2.3e-08 4.8e-08 5.4e-08 6.3e-08 6.6e-08 6.5e-08 7.1e-08 7.7e-08 8.3e-08
 * 5       32      2.7e-08 2.6e-08 2.6e-08 2.8e-08 5.0e-08 5.3e-08 6.2e-08 6.9e-08 6.5e-08 7.1e-08 7.7e-08 8.5e-08
 * 5       96      2.7e-08 2.7e-08 2.6e-08 2.6e-08 7.3e-08 7.7e-08 8.7e-08 9.2e-08 6.5e-08 7.0e-08 7.8e-08 8.4e-08
 * 6       8       5.8e-08 5.9e-08 5.6e-08 5.6e-08 8.7e-08 1.1e-07 1.4e-07 1.4e-07 1.4e-07 1.6e-07 2.0e-07 2.0e-07
 * 6       28      5.8e-08 5.6e-08 5.6e-08 5.8e-08 9.2e-08 1.1e-07 1.5e-07 1.5e-07 1.4e-07 1.6e-07 2.1e-07 2.0e-07
 * 6       36      5.7e-08 5.5e-08 5.5e-08 5.7e-08 9.2e-08 1.1e-07 1.4e-07 1.5e-07 1.4e-07 1.6e-07 2.0e-07 2.0e-07
 * 6       64      6.1e-08 6.0e-08 5.9e-08 6.0e-08 9.7e-08 1.2e-07 1.5e-07 1.5e-07 1.4e-07 1.6e-07 2.1e-07 1.9e-07
 * 6       192     6.1e-08 5.9e-08 6.1e-08 6.0e-08 1.1e-07 1.4e-07 1.7e-07 1.8e-07 1.4e-07 1.6e-07 2.1e-07 1.9e-07
 * 7       12      1.2e-07 1.1e-07 1.1e-07 1.1e-07 1.6e-07 2.2e-07 2.8e-07 2.9e-07 3.4e-07 4.0e-07 4.6e-07 4.4e-07
 * 7       60      1.1e-07 1.1e-07 1.1e-07 1.1e-07 2.1e-07 2.8e-07 3.7e-07 3.6e-07 3.6e-07 4.0e-07 4.7e-07 4.4e-07
 * 7       68      1.1e-07 1.1e-07 1.1e-07 1.2e-07 2.6e-07 3.0e-07 3.8e-07 3.8e-07 3.4e-07 4.0e-07 4.7e-07 4.4e-07
 * 7       128     1.2e-07 1.2e-07 1.1e-07 1.1e-07 2.7e-07 3.1e-07 3.9e-07 3.8e-07 3.4e-07 3.9e-07 4.7e-07 4.4e-07
 * 7       384     1.2e-07 1.1e-07 1.1e-07 1.1e-07 3.1e-07 3.7e-07 4.3e-07 4.2e-07 3.5e-07 4.1e-07 4.7e-07 4.4e-07
 * 8       24      3.0e-07 2.9e-07 2.9e-07 2.9e-07 3.5e-07 4.8e-07 6.5e-07 6.5e-07 7.5e-07 8.4e-07 1.0e-06 1.0e-06
 * 8       124     2.9e-07 2.8e-07 2.8e-07 2.8e-07 5.3e-07 6.4e-07 8.7e-07 8.6e-07 7.6e-07 8.4e-07 1.0e-06 1.0e-06
 * 8       132     2.9e-07 2.8e-07 2.8e-07 2.8e-07 5.6e-07 6.7e-07 8.8e-07 8.5e-07 7.3e-07 8.4e-07 1.1e-06 1.0e-06
 * 8       256     2.8e-07 2.7e-07 2.7e-07 2.7e-07 5.6e-07 6.8e-07 8.8e-07 8.5e-07 7.2e-07 8.4e-07 1.0e-06 1.0e-06
 * 8       768     2.8e-07 2.7e-07 2.7e-07 2.7e-07 6.5e-07 8.1e-07 9.8e-07 9.4e-07 7.3e-07 8.5e-07 1.0e-06 1.0e-06
 * 9       48      6.3e-07 6.1e-07 6.1e-07 6.0e-07 7.8e-07 1.1e-06 1.4e-06 1.5e-06 1.6e-06 1.9e-06 2.3e-06 2.3e-06
 * 9       252     6.2e-07 6.0e-07 5.9e-07 5.9e-07 1.2e-06 1.4e-06 1.9e-06 1.9e-06 1.7e-06 1.9e-06 2.3e-06 2.3e-06
 * 9       260     5.8e-07 5.7e-07 5.5e-07 5.6e-07 1.3e-06 1.5e-06 2.0e-06 2.0e-06 1.6e-06 1.9e-06 2.3e-06 2.3e-06
 * 9       512     5.7e-07 5.6e-07 5.4e-07 5.4e-07 1.3e-06 1.6e-06 2.0e-06 2.0e-06 1.6e-06 1.8e-06 2.3e-06 2.3e-06
 * 9       1536    5.6e-07 5.6e-07 5.4e-07 5.4e-07 1.5e-06 1.7e-06 2.3e-06 2.2e-06 1.6e-06 1.9e-06 2.4e-06 2.3e-06
 * 10      100     9.5e-07 9.4e-07 1.2e-06 1.2e-06 1.9e-06 2.6e-06 3.6e-06 3.6e-06 3.4e-06 4.1e-06 5.1e-06 5.1e-06
 * 10      508     9.8e-07 9.4e-07 1.2e-06 1.3e-06 2.7e-06 3.3e-06 4.4e-06 4.4e-06 3.5e-06 4.1e-06 5.1e-06 5.2e-06
 * 10      516     1.0e-06 9.8e-07 1.4e-06 1.3e-06 2.8e-06 3.4e-06 4.5e-06 4.4e-06 3.4e-06 4.1e-06 5.1e-06 5.1e-06
 * 10      1024    1.0e-06 9.9e-07 1.3e-06 1.3e-06 3.2e-06 3.5e-06 4.5e-06 4.4e-06 3.4e-06 4.1e-06 5.1e-06 5.2e-06
 * 10      3072    1.0e-06 9.9e-07 1.3e-06 1.3e-06 3.1e-06 3.8e-06 4.9e-06 4.7e-06 3.4e-06 4.1e-06 5.1e-06 5.1e-06
 * 11      204     1.6e-06 1.9e-06 2.6e-06 2.5e-06 4.2e-06 5.9e-06 7.9e-06 8.0e-06 7.4e-06 9.4e-06 1.1e-05 1.1e-05
 * 11      1020    1.8e-06 2.1e-06 2.8e-06 2.8e-06 5.6e-06 7.3e-06 9.6e-06 9.7e-06 7.4e-06 9.0e-06 1.1e-05 1.1e-05
 * 11      1028    1.9e-06 2.1e-06 2.8e-06 2.8e-06 6.1e-06 7.6e-06 9.8e-06 9.9e-06 7.3e-06 9.0e-06 1.1e-05 1.1e-05
 * 11      2048    2.0e-06 2.2e-06 2.8e-06 3.0e-06 6.2e-06 7.7e-06 9.9e-06 1.0e-05 7.6e-06 9.0e-06 1.1e-05 1.1e-05
 * 11      6144    2.0e-06 2.2e-06 3.0e-06 2.8e-06 6.6e-06 8.2e-06 1.1e-05 1.1e-05 7.3e-06 9.0e-06 1.1e-05 1.1e-05
 * 12      408     3.0e-06 3.9e-06 5.3e-06 5.3e-06 9.8e-06 1.4e-05 1.8e-05 1.9e-05 1.6e-05 1.9e-05 2.5e-05 2.5e-05
 * 12      2044    3.4e-06 4.4e-06 6.0e-06 5.8e-06 1.3e-05 1.7e-05 2.1e-05 2.2e-05 1.5e-05 2.0e-05 2.4e-05 2.5e-05
 * 12      2052    3.5e-06 4.5e-06 5.8e-06 5.8e-06 1.3e-05 1.7e-05 2.2e-05 2.1e-05 1.6e-05 1.9e-05 2.4e-05 2.5e-05
 * 12      4096    3.8e-06 4.6e-06 6.2e-06 6.2e-06 1.4e-05 1.7e-05 2.2e-05 2.1e-05 1.6e-05 1.9e-05 2.4e-05 2.5e-05
 * 12      12288   3.8e-06 4.6e-06 6.0e-06 6.2e-06 1.5e-05 1.8e-05 2.4e-05 2.3e-05 1.5e-05 1.9e-05 2.4e-05 2.5e-05
 * 13      816     6.4e-06 8.8e-06 1.2e-05 1.2e-05 2.1e-05 3.1e-05 4.0e-05 4.0e-05 3.4e-05 4.2e-05 5.3e-05 5.2e-05
 * 13      4092    7.3e-06 9.4e-06 1.2e-05 1.2e-05 2.6e-05 3.5e-05 4.7e-05 4.7e-05 3.4e-05 4.2e-05 5.3e-05 5.3e-05
 * 13      4100    7.4e-06 9.7e-06 1.2e-05 1.3e-05 2.9e-05 3.7e-05 4.7e-05 4.8e-05 3.7e-05 4.2e-05 5.3e-05 5.2e-05
 * 13      8192    7.6e-06 1.0e-05 1.3e-05 1.3e-05 3.0e-05 3.7e-05 4.8e-05 4.9e-05 3.4e-05 4.4e-05 5.3e-05 5.3e-05
 * 13      24576   7.7e-06 9.8e-06 1.3e-05 1.3e-05 3.2e-05 4.0e-05 5.1e-05 5.1e-05 3.3e-05 4.2e-05 5.4e-05 5.2e-05
 * 14      1636    1.3e-05 1.8e-05 2.5e-05 2.7e-05 4.7e-05 6.7e-05 9.1e-05 9.3e-05 7.0e-05 9.1e-05 1.1e-04 1.1e-04
 * 14      8188    1.6e-05 2.1e-05 2.7e-05 2.7e-05 5.9e-05 7.8e-05 1.0e-04 1.0e-04 7.0e-05 9.0e-05 1.1e-04 1.2e-04
 * 14      8196    1.6e-05 2.1e-05 2.9e-05 2.7e-05 6.0e-05 7.9e-05 1.0e-04 1.0e-04 7.0e-05 9.0e-05 1.1e-04 1.1e-04
 * 14      16384   1.7e-05 2.1e-05 2.7e-05 2.8e-05 6.1e-05 8.0e-05 1.0e-04 1.0e-04 7.1e-05 9.0e-05 1.1e-04 1.1e-04
 * 14      49152   1.7e-05 2.3e-05 2.8e-05 2.7e-05 6.6e-05 8.4e-05 1.1e-04 1.1e-04 7.1e-05 9.2e-05 1.1e-04 1.1e-04
 * 15      3276    2.7e-05 3.8e-05 5.1e-05 5.3e-05 1.0e-04 1.4e-04 1.9e-04 1.9e-04 1.5e-04 1.9e-04 2.4e-04 2.4e-04
 * 15      16380   3.3e-05 4.4e-05 6.1e-05 5.8e-05 1.2e-04 1.6e-04 2.2e-04 2.2e-04 1.5e-04 1.9e-04 2.4e-04 2.4e-04
 * 15      16388   3.3e-05 4.4e-05 5.8e-05 6.0e-05 1.3e-04 1.7e-04 2.2e-04 2.2e-04 1.5e-04 1.9e-04 2.5e-04 2.4e-04
 * 15      32768   3.4e-05 4.5e-05 5.8e-05 5.8e-05 1.3e-04 1.7e-04 2.2e-04 2.2e-04 1.5e-04 2.0e-04 2.5e-04 2.4e-04
 * 15      98304   3.4e-05 4.5e-05 5.8e-05 6.2e-05 1.4e-04 1.8e-04 2.3e-04 2.3e-04 1.5e-04 1.9e-04 2.4e-04 2.4e-04
 * 16      6552    5.8e-05 8.2e-05 1.1e-04 1.2e-04 2.2e-04 3.2e-04 4.9e-04 4.3e-04 3.1e-04 4.0e-04 5.9e-04 5.2e-04
 * 16      32764   6.6e-05 9.0e-05 1.2e-04 1.2e-04 2.6e-04 3.6e-04 4.8e-04 4.8e-04 3.1e-04 4.0e-04 5.1e-04 5.2e-04
 * 16      32772   6.6e-05 9.0e-05 1.2e-04 1.3e-04 2.8e-04 3.6e-04 4.8e-04 4.7e-04 3.2e-04 4.0e-04 5.2e-04 5.2e-04
 * 16      65536   6.7e-05 9.3e-05 1.2e-04 1.2e-04 2.8e-04 3.6e-04 4.8e-04 4.7e-04 3.1e-04 4.1e-04 5.1e-04 5.3e-04
 * 16      196608  7.0e-05 9.3e-05 1.2e-04 1.2e-04 3.0e-04 3.9e-04 5.0e-04 5.2e-04 3.1e-04 4.0e-04 5.1e-04 5.3e-04
 * 17      13104   1.3e-04 1.8e-04 2.4e-04 2.4e-04 4.8e-04 6.8e-04 9.1e-04 9.2e-04 6.6e-04 8.7e-04 1.1e-03 1.1e-03
 * 17      65532   1.4e-04 1.9e-04 2.6e-04 2.6e-04 5.5e-04 7.7e-04 1.0e-03 1.0e-03 6.7e-04 8.6e-04 1.1e-03 1.1e-03
 * 17      65540   1.4e-04 1.9e-04 2.6e-04 2.6e-04 6.0e-04 7.9e-04 1.1e-03 1.1e-03 6.6e-04 8.7e-04 1.2e-03 1.1e-03
 * 17      131072  1.4e-04 2.0e-04 2.6e-04 2.6e-04 6.0e-04 7.9e-04 1.0e-03 1.0e-03 6.6e-04 8.7e-04 1.1e-03 1.1e-03
 * 17      393216  1.5e-04 2.1e-04 2.6e-04 2.8e-04 6.6e-04 8.7e-04 1.1e-03 1.1e-03 6.6e-04 8.6e-04 1.1e-03 1.1e-03
 * 18      26212   2.7e-04 3.8e-04 5.3e-04 5.7e-04 1.0e-03 1.5e-03 2.0e-03 2.0e-03 1.4e-03 1.8e-03 2.3e-03 2.3e-03
 * Output on meteorlake (Intel(R) Core(TM) Ultra 7 165H)
 *
 */
