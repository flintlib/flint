#include "ulong_extras.h"
#include "fft_small.h"
#include "profiler.h"


void flint_print_d_fixed(double x, ulong l)
{
    ulong i;
    TMP_INIT;
    TMP_START;
    char* s = TMP_ARRAY_ALLOC(l + 1, char);

    for (i = 0; i < l; i++)
        s[i] = ' ';
    s[l] = 0;

    ulong y = fabs(rint(x));
    while (l > 0)
    {
        s[--l] = '0' + (y%10);
        y = y/10;
        if (y == 0)
            break;
    }

    printf("%s", s);

    TMP_END;
}

void flint_print_d_fixed_dot(double x, ulong l, ulong r)
{
    ulong i;
    TMP_INIT;
    TMP_START;
    char* s = TMP_ARRAY_ALLOC(l + r + 2, char);

    for (i = 0; i < l + r + 1; i++)
        s[i] = ' ';
    s[l + r + 1] = 0;

    s[l] = '.';
    ulong y = fabs(rint(x*n_pow(10, r)));
    while (r > 0)
    {
        --r;
        s[l+1+r] = '0' + (y%10);
        y = y/10;
    }
    while (l > 0)
    {
        --l;
        s[l] = '0' + (y%10);
        y = y/10;
        if (y == 0)
            break;
    }

    printf("%s", s);

    TMP_END;
}


void record_time(
    int first,
    double len,
    double t1, double t2)
{
    t1 *= 1e6;
    t2 *= 1e6;
    double depth = log2(len);

    flint_printf(first ? "depth" : "     ");
    flint_print_d_fixed_dot(depth, 3,1);
    flint_printf(":");
    flint_print_d_fixed_dot(5*len*depth/t1, 4,2);
    flint_print_d_fixed_dot(5*len*depth/t2, 4,2);
    flint_printf("\n");
}


#define EXP_FAC 16

void profile_v2_fft(sd_fft_ctx_t Q, ulong minL, ulong maxL)
{
    timeit_t timer;

    minL = n_max(minL, 10);
    maxL = n_max(maxL, minL + 1);

    timeit_start(timer);
    sd_fft_ctx_fit_depth(Q, maxL);
    timeit_stop(timer);
    flint_printf("depth %wu setup: %wd ms\n", maxL, timer->wall);

    minL = n_max(minL, LG_BLK_SZ + 1);
    for (ulong L = minL+1; L <= maxL; L++)
    {
        double fft_time, ifft_time;
        ulong  sz = sd_fft_ctx_data_size(L)*sizeof(double);
        double* data = flint_aligned_alloc(4096, n_round_up(sz, 4096));

        // do 1/2*2^L < otrunc <= 2^L
        ulong otrunc = n_pow2(L-1);
        while (otrunc < n_pow2(L))
        {
            otrunc = n_round_up(1+1+(EXP_FAC+1)*otrunc/EXP_FAC, BLK_SZ);
            otrunc = n_min(otrunc, n_pow2(L));
            ulong itrunc = n_round_up(otrunc/2, BLK_SZ);
            for (ulong i = 0; i < n_pow2(L); i++)
                data[i] = 0;
            ulong nreps = 1 + 300000000/(otrunc*n_clog2(otrunc));

            timeit_start(timer);
            for (ulong i = 0; i < nreps; i++)
                sd_fft_trunc(Q, data, L, itrunc, otrunc);
            timeit_stop(timer);
            fft_time = timer->wall;

            timeit_start(timer);
            for (ulong i = 0; i < nreps; i++)
                sd_ifft_trunc(Q, data, L, otrunc);
            timeit_stop(timer);
            ifft_time = timer->wall;

            record_time(otrunc == n_pow2(L), otrunc,
                        fft_time/nreps, ifft_time/nreps);
        }

        flint_free(data);
    }
}


int main(void)
{
    sd_fft_ctx_t Q;
    sd_fft_ctx_init_prime(Q, UWORD(0x0003f00000000001));
    flint_printf("--- v2 fft/ifft gflops ---\n");
    profile_v2_fft(Q, 9, 20);
    sd_fft_ctx_clear(Q);

    return 0;
}
