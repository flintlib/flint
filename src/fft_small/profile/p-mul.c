#include <gmp.h>
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


void profile_mul(
    mpn_ctx_t Q,
    int min_lg_bit_size, int max_lg_bit_size,
    int for_latex,
    int with_precomp)
{
    timeit_t timer;
    int first = 1;
    min_lg_bit_size = FLINT_MAX(min_lg_bit_size, 12);
    for (int lg_bit_size = min_lg_bit_size; lg_bit_size < max_lg_bit_size; lg_bit_size++)
    {
        ulong max_len = n_pow2(1+lg_bit_size-6);
        ulong* data = FLINT_ARRAY_ALLOC(max_len*2, ulong);
        for (ulong i = 0; i < max_len; i++)
            data[i] = -(ulong)(1+i);

        double total_precomp = 0;
        double max_precomp = 0;
        ulong nprecomp_samples = 0;
        for (int ci = 0; ci < 10; ci++)
        {
            ulong cn = max_len/2*pow(2.0, ci*1.0/10);
            ulong cbits = cn*64;
            double lgcbits = log2(cbits);
            ulong nsamples = 0;
            double total_time = 0;
            double max_time = 0;
            double min_time = 1.0e100;
            for (ulong an = (cn+1)/2; an <= 3*cn/4; an += 1+an/12) /* usually 1+an/12 */
            {
                ulong bn = cn - an;
                if (!(bn <= an && an < cn))
                    continue;
                nsamples++;
                ulong* a = data;
                ulong* b = data + an;
                ulong* c = data + max_len;
                double precomp = 0;
                if (with_precomp)
                {
                    timeit_start_us(timer);
                    mpn_ctx_mpn_mul(Q, c, a, an, b, bn);
                    timeit_stop_us(timer);
                    precomp = (double)timer->wall;
                    nprecomp_samples++;
                }
                ulong nreps = 1 + 30000000/(cn*n_clog2(cn));
                timeit_start_us(timer);
                for (ulong rep = 0; rep < nreps; rep++)
                    mpn_ctx_mpn_mul(Q, c, a, an, b, bn);
                timeit_stop_us(timer);
                double time = ((double)timer->wall)/nreps;

                precomp -= time;
                precomp = precomp*1e5/(lgcbits*cbits);
                max_precomp = FLINT_MAX(max_precomp, precomp);
                total_precomp += precomp;

                time = time*1e5/(lgcbits*cbits);
                total_time += time;
                min_time = FLINT_MIN(min_time, time);
                max_time = FLINT_MAX(max_time, time);
            }

            if (nsamples < 1)
                continue;

            if (for_latex)
            {
                if (!first)
                    flint_printf(",");
                flint_printf("(%f,%f),", lgcbits, min_time);
                flint_printf("(%f,%f),", lgcbits, total_time/nsamples);
                flint_printf("(%f,%f)", lgcbits, max_time);
            }
            else
            {
                flint_print_d_fixed_dot(lgcbits, 2, 2);
                flint_printf(":");
                flint_print_d_fixed_dot(total_time/nsamples, 3, 3);
                flint_print_d_fixed_dot(min_time, 3, 3);
                flint_print_d_fixed_dot(max_time, 3, 3);
                flint_printf("\n");
            }
            fflush(stdout);
            first = 0;
        }

        if (!for_latex && nprecomp_samples > 0)
        {
            flint_printf("precomp: avg %f, max %f\n", total_precomp/nprecomp_samples, max_precomp);
        }

        flint_free(data);
    }
}


void profile_mul_gmp(
    int min_lg_bit_size, int max_lg_bit_size,
    int for_latex,
    int with_precomp)
{
    timeit_t timer;
    int first = 1;
    min_lg_bit_size = FLINT_MAX(min_lg_bit_size, 12);
    for (int lg_bit_size = min_lg_bit_size; lg_bit_size < max_lg_bit_size; lg_bit_size++)
    {
        ulong max_len = n_pow2(1+lg_bit_size-6);
        ulong* data = FLINT_ARRAY_ALLOC(max_len*2, ulong);
        for (ulong i = 0; i < max_len; i++)
            data[i] = -(ulong)(1+i);

        double total_precomp = 0;
        double max_precomp = 0;
        ulong nprecomp_samples = 0;
        for (int ci = 0; ci < 10; ci++)
        {
            ulong cn = max_len/2*pow(2.0, ci*1.0/10);
            ulong cbits = cn*64;
            double lgcbits = log2(cbits);
            ulong nsamples = 0;
            double total_time = 0;
            double max_time = 0;
            double min_time = 1.0e100;
            for (ulong an = (cn+1)/2; an <= 3*cn/4; an += 1+an/12)
            {
                ulong bn = cn - an;
                if (!(bn <= an && an < cn))
                    continue;
                nsamples++;
                ulong* a = data;
                ulong* b = data + an;
                ulong* c = data + max_len;
                double precomp = 0;
                if (with_precomp)
                {
                    timeit_start(timer);
                    mpn_mul(c, a, an, b, bn);
                    timeit_stop(timer);
                    precomp = (double)timer->wall;
                    nprecomp_samples++;
                }
                ulong nreps = 1 + 20000000/(cn*n_clog2(cn));
                timeit_start(timer);
                for (ulong rep = 0; rep < nreps; rep++)
                    mpn_mul(c, a, an, b, bn);
                timeit_stop(timer);
                precomp -= ((double)timer->wall)/nreps;
                max_precomp = FLINT_MAX(max_precomp, precomp*1e8/(lgcbits*cbits));
                total_precomp += precomp*1e8/(lgcbits*cbits);
                double time = ((double)timer->wall)*1e8/(lgcbits*cbits*nreps);
                total_time += time;
                min_time = FLINT_MIN(min_time, time);
                max_time = FLINT_MAX(max_time, time);
            }

            if (nsamples < 1)
                continue;

            if (for_latex)
            {
                if (!first)
                    flint_printf(",");
                flint_printf("(%f,%f),", lgcbits, min_time);
                flint_printf("(%f,%f),", lgcbits, total_time/nsamples);
                flint_printf("(%f,%f)", lgcbits, max_time);
            }
            else
            {
                flint_print_d_fixed_dot(lgcbits, 2, 2);
                flint_printf(":");
                flint_print_d_fixed_dot(total_time/nsamples, 3, 3);
                flint_print_d_fixed_dot(min_time, 3, 3);
                flint_print_d_fixed_dot(max_time, 3, 3);
                flint_printf("\n");
            }
            fflush(stdout);
            first = 0;
        }

        if (!for_latex && nprecomp_samples > 0)
        {
            flint_printf("precomp: %f, %f\n", total_precomp/nprecomp_samples, max_precomp);
        }

        flint_free(data);
    }
}

/*
some notes on precomp:

(1) the global twiddle factors need to be precomputed
(2) when the big buffer for temp space needs to be reallocated, the accesses
    to the new space all incur page faults. These occur out of order in the
    beginning of the calculation and contribute noticably to the run time.

Therefore, there is a penalty for the first run of a computation of a certain
size. If the data comes out like

24.00:  2.856  2.831  2.980
24.10:  2.827  2.769  2.885
24.20:  2.984  2.948  3.055
24.30:  2.939  2.889  3.038
24.40:  3.062  3.055  3.101
24.50:  3.049  2.982  3.097
24.60:  3.082  3.037  3.144
24.70:  3.053  2.971  3.169
24.80:  2.999  2.899  3.083
24.90:  2.979  2.950  2.993
precomp: avg 0.073270, max 0.852562

this means that the penalty on average was 0.073270/3 = 2.5%
and at most 0.85/3 = 28.3%
*/

int main(void)
{
    mpn_ctx_t R;
#if 0
    int cpu_affinities[32];
    for (int i = 0; i < 32; i++)
        cpu_affinities[i] = i;
#endif
    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    /* the majority of the precomp is reallocating the temp buffer (2) */
//    mpn_ctx_fit_buffer(R, 1610620928);

flint_printf(" --- fft_small 1 thread  --- \n");
    flint_set_num_threads(1);
    profile_mul(R, 20, 26, 0, 1);

#if 0
flint_printf(" --- fft_small 8 threads --- \n");
    flint_set_num_threads(8);
    profile_mul(R, 20, 28, 1, 1);
    mpn_ctx_clear(R);
#endif

#if 0
flint_printf(" --- gmp --- \n");
    profile_mul_gmp(12, 31, 1, 1);
#endif

    mpn_ctx_clear(R);

    return 0;
}
