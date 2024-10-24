#include "profiler.h"
#include "nmod_vec.h"
#include "fft_small.h"
#include "n_fft.h"

#define num_primes 5

typedef struct
{
   ulong prime;
   ulong depth;
   ulong maxdepth;
   ulong stride;
} info_t;

#define SAMPLE(fun, _variant)                                                    \
void sample_##fun##_variant(void * arg, ulong count)                             \
{                                                                                \
    info_t * info = (info_t *) arg;                                              \
    const ulong p = info->prime;                                                 \
    const ulong depth = info->depth;                                             \
    const ulong maxdepth = info->maxdepth;                                       \
    const ulong stride = info->stride;                                           \
                                                                                 \
    const ulong len = stride * (UWORD(1) << depth);                              \
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));                \
                                                                                 \
    /* modulus, roots of unity */                                                \
    nmod_t mod;                                                                  \
    nmod_init(&mod, p);                                                          \
    ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> maxdepth, mod); \
    ulong w = nmod_pow_ui(w0, 1UL<<(maxdepth - depth), mod);                     \
    n_fft_ctx_t F;                                                               \
    n_fft_ctx_init2_root(F, w, depth, depth, p);                                 \
                                                                                 \
    FLINT_TEST_INIT(state);                                                      \
                                                                                 \
    ulong * coeffs = _nmod_vec_init(len);                                        \
    _nmod_vec_randtest(coeffs, state, len, mod);                                 \
                                                                                 \
    for (ulong i = 0; i < count; i++)                                            \
    {                                                                            \
        prof_start();                                                            \
        for (ulong j = 0; j < rep; j++)                                          \
            n_fft_##fun##_variant(coeffs, depth, F);                             \
        prof_stop();                                                             \
    }                                                                            \
                                                                                 \
    n_fft_ctx_clear(F);                                                          \
    FLINT_TEST_CLEAR(state);                                                     \
}                                                                                \

SAMPLE(dft, )
//SAMPLE(n_fft_dft, _stride)

void sample_sd_fft(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    const ulong p = info->prime;
    const ulong depth = info->depth;

    const ulong len = UWORD(1) << depth;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    sd_fft_ctx_t Q;
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_ctx_fit_depth(Q, depth);

    ulong sz = sd_fft_ctx_data_size(depth)*sizeof(double);

    FLINT_TEST_INIT(state);

    nmod_t mod;
    nmod_init(&mod, p);
    ulong * coeffs = _nmod_vec_init(len);
    _nmod_vec_randtest(coeffs, state, len, mod);

    double* data = flint_aligned_alloc(4096, n_round_up(sz, 4096));
    for (ulong i = 0; i < len; i++)
        data[i] = coeffs[i];

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
            sd_fft_trunc(Q, data, depth, len, len);
        prof_stop();
    }

    sd_fft_ctx_clear(Q);
    FLINT_TEST_CLEAR(state);
}

int main()
{
    flint_printf("- depth is log(fft length)\n");
    flint_printf("- timing DFT (length power of 2) for several bit lengths and depths\n");
    flint_printf("depth\tsd_fft\trec4\n");

    // FIXME FLINT_BITS issue
    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
    };
    ulong max_depths[num_primes] = { 18, 25, 25, 25, 25 };

    for (ulong k = 3; k < 4; k++)
    {
        for (ulong depth = 3; depth <= max_depths[k]; depth++)
        {
            printf("%ld\t", depth);

            info_t info;
            info.prime = primes[k];
            info.maxdepth = max_depths[k];
            info.depth = depth;
            info.stride = 1;

            const ulong len = UWORD(1) << depth;
            const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

            double min[15];
            double max;

            prof_repeat(min+0, &max, sample_sd_fft, (void *) &info);
            prof_repeat(min+1, &max, sample_dft, (void *) &info);

            flint_printf("%.1e\t%.1e\t\n",
                    min[0]/(double)1000000/rep,
                    min[1]/(double)1000000/rep
                    );
        }
    }
    return 0;
}
