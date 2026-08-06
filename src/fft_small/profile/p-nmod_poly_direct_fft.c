/*
    Benchmark _nmod_poly_mul_mid_mpn_ctx when the modulus itself
    can serve as the FFT prime.
*/

#include "nmod.h"
#include "nmod_poly.h"
#include "fft_small.h"
#include "profiler.h"
#include "ulong_extras.h"

static ulong
find_fft_prime(ulong bits, ulong depth)
{
    ulong lo, hi, step, k, klo, khi;

    if (bits <= depth || bits > 50)
        return 0;

    lo = UWORD(1) << (bits - 1);
    hi = UWORD(1) << bits;
    step = UWORD(1) << depth;
    klo = (lo <= 1) ? 1 : ((lo - 1 + step - 1) >> depth);
    khi = (hi - 2) >> depth;

    if ((klo & 1) == 0)
        klo++;

    for (k = klo; k <= khi; k += 2)
    {
        ulong p = (k << depth) + 1;

        if (n_is_prime(p) && fft_small_mulmod_satisfies_bounds(p))
            return p;
    }

    return 0;
}

static double
bench_one(ulong p, ulong n, int threads, flint_rand_t state, mpn_ctx_t R)
{
    nmod_t mod;
    ulong * a, * b, * z;
    ulong i, reps;
    timeit_t timer;
    double best;

    nmod_init(&mod, p);

    a = FLINT_ARRAY_ALLOC(n, ulong);
    b = FLINT_ARRAY_ALLOC(n, ulong);
    z = FLINT_ARRAY_ALLOC(2*n - 1, ulong);

    for (i = 0; i < n; i++)
    {
        a[i] = n_randint(state, p);
        b[i] = n_randint(state, p);
    }

    flint_set_num_threads(threads);

    _nmod_poly_mul_mid_mpn_ctx(z, 0, 2*n - 1, a, n, b, n, mod, R);

    reps = 1;
    if (n <= 1500)
        reps = 8;
    else if (n <= 3000)
        reps = 5;
    else if (n <= 6000)
        reps = 3;

    best = 1e100;
    for (i = 0; i < 3; i++)
    {
        double t;

        timeit_start_us(timer);
        for (ulong j = 0; j < reps; j++)
            _nmod_poly_mul_mid_mpn_ctx(z, 0, 2*n - 1, a, n, b, n, mod, R);
        timeit_stop_us(timer);

        t = ((double) timer->wall)/reps;
        best = FLINT_MIN(best, t);
    }

    flint_free(a);
    flint_free(b);
    flint_free(z);

    return best;
}

int main(void)
{
    const ulong ns[] = {1500, 3000, 6000, 12000};
    const int thread_counts[] = {1, 8};
    flint_rand_t state;
    mpn_ctx_t R;

    flint_rand_init(state);
    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    flint_printf("threads,bits,p,n,depth,usec\n");

    for (ulong ti = 0; ti < sizeof(thread_counts)/sizeof(thread_counts[0]); ti++)
    {
        int threads = thread_counts[ti];

        for (ulong ni = 0; ni < sizeof(ns)/sizeof(ns[0]); ni++)
        {
            ulong n = ns[ni];
            ulong zn = 2*n - 1;
            ulong ztrunc = n_round_up(zn, BLK_SZ);
            ulong depth = n_max(LG_BLK_SZ, n_clog2(ztrunc));

            for (ulong bits = depth + 1; bits <= 50; bits++)
            {
                ulong p = find_fft_prime(bits, depth);

                if (p != 0)
                {
                    double t = bench_one(p, n, threads, state, R);
                    flint_printf("%d,%wu,%wu,%wu,%wu,%.3f\n",
                                 threads, bits, p, n, depth, t);
                    fflush(stdout);
                }
            }
        }
    }

    mpn_ctx_clear(R);
    flint_rand_clear(state);
    flint_cleanup();

    return 0;
}
