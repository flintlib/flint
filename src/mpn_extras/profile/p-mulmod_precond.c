#include "mpn_extras.h"
#include "profiler.h"

slong ns[] = { 2, 3, 4, 6, 8, 9, 10, 11, 12, 16, 20, 24, 32, 48, 64, 96, 128, 192, 256, 320, 384, 0 };

int main()
{
    slong i, ni, n, npre, num, numi;
    flint_rand_t state;
    flint_rand_init(state);

    mp_ptr a, apre, b, d, dnormed, dinv, r1, r2, t, u;
    ulong norm;
    int normed;
    double t1, t2, FLINT_SET_BUT_UNUSED(tcpu);

    flint_printf("        ");
    for (numi = 0; (num = ns[numi]) != 0 && num <= 128; numi++)
        flint_printf("%7wd", num);
    flint_printf("\n");

    for (normed = 0; normed <= 1; normed++)
    {
        flint_printf("normed = %d\n", normed);

        for (ni = 0; (n = ns[ni]) != 0; ni++)
        {
            flint_printf("%8wd", n);

            for (numi = 0; (num = ns[numi]) != 0 && num <= 128; numi++)
            {
                npre = flint_mpn_mulmod_precond_alloc(n);

                a = flint_malloc(sizeof(mp_limb_t) * n);
                apre = flint_malloc(sizeof(mp_limb_t) * npre);
                b = flint_malloc(sizeof(mp_limb_t) * (num * n));
                d = flint_malloc(sizeof(mp_limb_t) * n);
                dnormed = flint_malloc(sizeof(mp_limb_t) * n);
                dinv = flint_malloc(sizeof(mp_limb_t) * n);
                t = flint_malloc(sizeof(mp_limb_t) * 2 * n);
                u = flint_malloc(sizeof(mp_limb_t) * (n + 1));
                r1 = flint_malloc(sizeof(mp_limb_t) * (num * n));
                r2 = flint_malloc(sizeof(mp_limb_t) * (num * n));

                flint_mpn_rrandom(d, state, n);
                do
                {
                    flint_mpn_rrandom(d + n - 1, state, 1);
                    if (normed == 0)
                        d[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));
                    else
                        d[n - 1] >>= 1;
                }
                while (d[n - 1] == 0);

                norm = flint_clz(d[n - 1]);
                if (norm == 0)
                    mpn_copyi(dnormed, d, n);
                else
                    mpn_lshift(dnormed, d, n, norm);
                flint_mpn_preinvn(dinv, dnormed, n);

                /* reduce a, b mod d */
                flint_mpn_rrandom(a, state, n);
                mpn_tdiv_qr(t, a, 0, a, n, d, n);

                for (i = 0; i < num; i++)
                {
                    flint_mpn_rrandom(b + i * n, state, n);
                    mpn_tdiv_qr(t, b + i * n, 0, b + i * n, n, d, n);
                }

                TIMEIT_START
                if (norm == 0)
                {
                    for (i = 0; i < num; i++)
                        flint_mpn_mulmod_preinvn(r1 + i * n, a, b + i * n, n, dnormed, dinv, norm);
                }
                else
                {
                    mpn_lshift(t, a, n, norm);
                    for (i = 0; i < num; i++)
                    {
                        flint_mpn_mulmod_preinvn(r1 + i * n, t, b + i * n, n, dnormed, dinv, 0);
                        mpn_rshift(r1 + i * n, r1 + i * n, n, norm);
                    }
                }
                TIMEIT_STOP_VALUES(tcpu, t1)

                TIMEIT_START
                flint_mpn_mulmod_precond_precompute(apre, a, n, dnormed, dinv, norm);
                for (i = 0; i < num; i++)
                    flint_mpn_mulmod_precond(r2 + i * n, apre, b + i * n, n, dnormed, dinv, norm);
                TIMEIT_STOP_VALUES(tcpu, t2)

                flint_printf("  %.3f", n, num, t1 / t2);
                fflush(stdout);

                if (mpn_cmp(r1, r2, n * num))
                {
                    flint_printf("FAIL\n");
                    flint_abort();
                }

                flint_free(a);
                flint_free(apre);
                flint_free(b);
                flint_free(d);
                flint_free(dnormed);
                flint_free(dinv);
                flint_free(t);
                flint_free(u);
                flint_free(r1);
                flint_free(r2);
            }

            flint_printf("\n");
        }
    }

    flint_rand_clear(state);
}

