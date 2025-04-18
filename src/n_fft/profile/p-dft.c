#include "nmod_poly.h"
#include "profiler.h"
#include "nmod_vec.h"
#include "fft_small.h"
#include "n_fft.h"

#define num_primes 7

typedef struct
{
   ulong prime;
   ulong depth;
   ulong stride;
} info_t;

#define SAMPLE(fun, _variant)                                                    \
void sample_##fun##_variant(void * arg, ulong count)                             \
{                                                                                \
    info_t * info = (info_t *) arg;                                              \
    const ulong p = info->prime;                                                 \
    const ulong depth = info->depth;                                             \
    const ulong stride = info->stride;                                           \
                                                                                 \
    const ulong len = stride * (UWORD(1) << depth);                              \
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));                \
                                                                                 \
    /* modulus, roots of unity */                                                \
    n_fft_ctx_t F;                                                               \
    n_fft_ctx_init2(F, depth, p);                                                \
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
            n_fft_##fun##_variant(coeffs, depth, F);                             \
        prof_stop();                                                             \
    }                                                                            \
                                                                                 \
    _nmod_vec_clear(coeffs);                                                     \
    n_fft_ctx_clear(F);                                                          \
    FLINT_TEST_CLEAR(state);                                                     \
}                                                                                \

void n_fft_mul_poly(nmod_poly_t ab, nmod_poly_t a, nmod_poly_t b, n_fft_ctx_t F)
{
    const ulong depth = n_clog2(a->length + b->length - 1);
    const ulong len = UWORD(1) << depth;

    nn_ptr va = _nmod_vec_init(len);
    nn_ptr vb = _nmod_vec_init(len);
    _nmod_vec_set(va, a->coeffs, a->length);
    _nmod_vec_zero(va + a->length, len - a->length);
    _nmod_vec_set(vb, b->coeffs, b->length);
    _nmod_vec_zero(vb + b->length, len - b->length);

    n_fft_dft_lazy_1_1(va, depth, F);
    n_fft_dft_lazy_1_1(vb, depth, F);
    for (ulong k = 0; k < len; k++)
    {
        va[k] = nmod_mul(va[k], vb[k], a->mod);
        // Note: a variant that only needs < 2n, thus allowing to ignore the
        // reductions to [0..n) after DFT, such as using lazy_1_4 and then
        //            if (va[k] >= 2*F->mod) va[k] -= 2*F->mod;
        //            if (vb[k] >= 2*F->mod) vb[k] -= 2*F->mod;
        //            ulong p_hi, p_lo;
        //            umul_ppmm(p_hi, p_lo, va[k], vb[k]);
        //            //NMOD_RED2(va[k], p_hi, p_lo, a->mod);
        //            {
        //                ulong q0xx, q1xx, r1xx;
        //                const ulong u1xx = ((p_hi)<<a->mod.norm)
        //                + (((a->mod).norm == 0) ? UWORD(0) : (p_lo)>>(FLINT_BITS - (a->mod).norm));
        //                const ulong u0xx = (p_lo)<<(a->mod).norm;
        //                const ulong nxx = (a->mod).n<<(a->mod).norm;
        //                umul_ppmm(q1xx, q0xx, (a->mod).ninv, u1xx);
        //                add_ssaaaa(q1xx, q0xx, q1xx, q0xx, u1xx, u0xx);
        //                r1xx = (u0xx - (q1xx + 1)*nxx);
        //                if (r1xx > q0xx) r1xx += nxx;
        //                va[k] = (r1xx>>(a->mod).norm);
        //                //if (r1xx < nxx) r = (r1xx>>(a->mod).norm);
        //                //else r = ((r1xx - nxx)>>(a->mod).norm);
        //            }

        // basically gives no speedup (including when simplifying the end of
        // NMOD_RED2 using the fact that idft accepts [0..2n))
    }
    n_fft_idft(va, depth, F);
    nmod_poly_fit_length(ab, len);
    _nmod_poly_set_length(ab, len);
    _nmod_vec_set(ab->coeffs, va, len);
    _nmod_poly_normalise(ab);
    _nmod_vec_clear(va);
    _nmod_vec_clear(vb);
}


void sample_polymul(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    const ulong p = info->prime;
    const ulong depth = info->depth;
    const ulong stride = info->stride;

    const ulong len = stride * (UWORD(1) << depth);
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    /* modulus, roots of unity */
    n_fft_ctx_t F;
    n_fft_ctx_init2(F, depth, p);
    nmod_t mod;
    nmod_init(&mod, p);

    /* polynomials */
    /* FIXME: vary the lengths */
    /* note: make sure depth of FFT context is large enough! */
    const ulong pdepth = depth-1;
    const ulong plen = stride * (UWORD(1) << pdepth);
    nmod_poly_t a;
    nmod_poly_t b;
    nmod_poly_t ab;
    nmod_poly_init(a, p);
    nmod_poly_init(b, p);
    nmod_poly_init(ab, p);

    FLINT_TEST_INIT(state);

    nmod_poly_fit_length(a, plen);
    _nmod_poly_set_length(a, plen);
    nmod_poly_fit_length(b, plen);
    _nmod_poly_set_length(b, plen);
    for (ulong k = 0; k < plen; k++)
    {
        a->coeffs[k] = n_randint(state, p);
        b->coeffs[k] = n_randint(state, p);
    }
    _nmod_poly_normalise(a);
    _nmod_poly_normalise(b);

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
            n_fft_mul_poly(ab, a, b, F);
        prof_stop();
        //n_fft_mul_poly(ab, a, b, F);
        //nmod_poly_t ab2;
        //nmod_poly_init(ab2, p);
        //nmod_poly_mul(ab2, a, b);
        //if (!nmod_poly_equal(ab, ab2))
        //{
        //    printf("%ld: %ld - %ld - %ld\n", i, nmod_poly_degree(ab), nmod_poly_degree(a), nmod_poly_degree(b));
        //    printf("!!!WRONG!!!\n");
        //    if (pdepth == 2)
        //    {
        //        flint_printf("%{nmod_poly}\n", ab);
        //        flint_printf("%{nmod_poly}\n", ab2);
        //    }
        //}
        //nmod_poly_clear(ab2);
    }

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    n_fft_ctx_clear(F);
    FLINT_TEST_CLEAR(state);
}

SAMPLE(dft, )
SAMPLE(idft, )
SAMPLE(dft_t, )
SAMPLE(idft_t, )
SAMPLE(dft_lazy_1_1, )
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
    flint_printf("depth\tsd_fft\tdft\tidft\tdft_t\tidft_t\n");

    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        1073479681,          // 30 bits, 1 + 2**30 - 2**18 == 1 + 2**18 * (2**12 - 1)
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
        4611686018427322369  // 62 bits: 1 + 2**62 - 2**16 == 1 + 2**16 * (2**46 - 1)
    };
    ulong max_depths[num_primes] = { 18, 18, 25, 25, 25, 25, 16 };

    for (ulong k = 5; k < 7; k++)
    {
        for (ulong depth = 3; depth <= max_depths[k]; depth++)
        {
            printf("%ld\t", depth);

            info_t info;
            info.prime = primes[k];
            info.depth = depth;
            info.stride = 1;

            const ulong len = UWORD(1) << depth;
            const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

            double min[15];
            double max;

            //prof_repeat(min+0, &max, sample_sd_fft, (void *) &info);
            prof_repeat(min+1, &max, sample_dft, (void *) &info);
            prof_repeat(min+2, &max, sample_idft, (void *) &info);
            prof_repeat(min+3, &max, sample_dft_t, (void *) &info);
            prof_repeat(min+4, &max, sample_idft_t, (void *) &info);
            prof_repeat(min+5, &max, sample_dft_lazy_1_1, (void *) &info);
            prof_repeat(min+6, &max, sample_polymul, (void *) &info);

            flint_printf("%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\n",
                    min[0]/(double)1000000/rep,
                    min[1]/(double)1000000/rep,
                    min[2]/(double)1000000/rep,
                    min[3]/(double)1000000/rep,
                    min[4]/(double)1000000/rep,
                    min[5]/(double)1000000/rep,
                    min[6]/(double)1000000/rep
                    );
        }
    }
    return 0;
}
