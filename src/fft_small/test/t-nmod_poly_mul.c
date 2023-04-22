#if defined(__AVX2__)

#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_poly.h"
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
        ulong * a, *b, * c, * d;
        ulong an, zn, zl, zh, sz, i, reps;

        flint_printf("sqr nbits: %02wu", nbits);
        fflush(stdout);

        for (reps = 0; reps < 100; reps++)
        {
            flint_set_num_threads(1 + n_randint(state, 10));

            nmod_init(&mod, n_randbits(state, nbits));

            an = 1 + n_randint(state, 7000);
            zn = an + an - 1;
            zl = n_randint(state, zn+10);
            zh = n_randint(state, zn+20);

            sz = FLINT_MAX(zl, zh);
            sz = FLINT_MAX(sz, zn);

            a = FLINT_ARRAY_ALLOC(an, ulong);
            b = FLINT_ARRAY_ALLOC(an, ulong);
            c = FLINT_ARRAY_ALLOC(sz, ulong);
            d = FLINT_ARRAY_ALLOC(sz, ulong);

            for (i = 0; i < an; i++)
                b[i] = a[i] = n_randint(state, mod.n);

            flint_mpn_zero(c, sz);
            _nmod_poly_mul(c, a, an, a, an, mod);
            _nmod_poly_mul_mid_mpn_ctx(d, zl, zh, a, an, b, an, mod, R);

            for (i = zl; i < zh; i++)
            {
                if (c[i] != d[i-zl])
                {
                    flint_printf("(squaring) mulmid error at index %wu\n", i);
                    flint_printf("zl=%wu, zh=%wu, an=%wu\n", zl, zh, an);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            flint_free(a);
            flint_free(b);
            flint_free(c);
            flint_free(d);
        }

        for (i = 0; i < 13; i++)
            flint_printf("%c", '\b');
        fflush(stdout);
    }

    for (nbits = 1; nbits <= FLINT_BITS; nbits ++)
    {
        ulong * a, * b, * c, * d;
        ulong an, bn, zn, zl, zh, sz, i, reps;

        flint_printf("mul nbits: %02wu", nbits);
        fflush(stdout);

        for (reps = 0; reps < 100; reps++)
        {
            flint_set_num_threads(1 + n_randint(state, 10));

            nmod_init(&mod, n_randbits(state, nbits));

            an = 1 + n_randint(state, 7000);
            bn = 1 + n_randint(state, an);
            zn = an + bn - 1;
            zl = n_randint(state, zn+10);
            zh = n_randint(state, zn+20);

            sz = FLINT_MAX(zl, zh);
            sz = FLINT_MAX(sz, zn);

            a = FLINT_ARRAY_ALLOC(an, ulong);
            b = FLINT_ARRAY_ALLOC(bn, ulong);
            c = FLINT_ARRAY_ALLOC(sz, ulong);
            d = FLINT_ARRAY_ALLOC(sz, ulong);

            for (i = 0; i < an; i++)
                a[i] = n_randint(state, mod.n);

            for (i = 0; i < bn; i++)
                b[i] = n_randint(state, mod.n);

            flint_mpn_zero(c, sz);
            _nmod_poly_mul(c, a, an, b, bn, mod);
            _nmod_poly_mul_mid_mpn_ctx(d, zl, zh, a, an, b, bn, mod, R);

            for (i = zl; i < zh; i++)
            {
                if (c[i] != d[i-zl])
                {
                    flint_printf("mulmid error at index %wu\n", i);
                    flint_printf("zl=%wu, zh=%wu, an=%wu, bn=%wu\n", zl, zh, an, bn);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            flint_free(a);
            flint_free(b);
            flint_free(c);
            flint_free(d);
        }

        for (i = 0; i < 13; i++)
            flint_printf("%c", '\b');
        fflush(stdout);
    }


    for (nbits = 2; nbits <= FLINT_BITS; nbits++)
    {
        ulong * a, * b, * q1, * q2, * q3, * r1, * r2, * r3;
        ulong an, bn, qn, i, reps;
        nmod_poly_divrem_precomp_struct M[1];

        flint_printf("divrem nbits: %02wu", nbits);
        fflush(stdout);

        for (reps = 0; reps < 100; reps++)
        {
            flint_set_num_threads(1 + n_randint(state, 10));

            nmod_init(&mod, n_randbits(state, nbits));

            bn = 2 + n_randint(state, 5000);
            qn = 1 + n_randint(state, 5000);
            an = bn + qn - 1;

            a = FLINT_ARRAY_ALLOC(an, ulong);
            b = FLINT_ARRAY_ALLOC(bn, ulong);
            q1 = FLINT_ARRAY_ALLOC(qn, ulong);
            q2 = FLINT_ARRAY_ALLOC(qn, ulong);
            q3 = FLINT_ARRAY_ALLOC(qn, ulong);
            r1 = FLINT_ARRAY_ALLOC(bn, ulong);
            r2 = FLINT_ARRAY_ALLOC(bn, ulong);
            r3 = FLINT_ARRAY_ALLOC(bn, ulong);

            for (i = 0; i < an; i++)
                a[i] = n_randint(state, mod.n);

            for (i = 0; i < bn; i++)
                b[i] = n_randint(state, mod.n);

            while (n_gcd(b[bn-1], mod.n) != 1)
                b[bn-1] = n_randint(state, mod.n);

            _nmod_poly_divrem(q1, r1, a, an, b, bn, mod);

            _nmod_poly_divrem_mpn_ctx(q2, r2, a, an, b, bn, mod, R);

            ulong prec = qn + n_randint(state, 200);
            _nmod_poly_divrem_precomp_init(M, b, bn, prec, mod, R);
            _nmod_poly_divrem_precomp(q3, r3, a, an, M, mod, R);
            _nmod_poly_divrem_precomp_clear(M);            

            for (i = qn; i > 0; i--)
            {
                if (q1[i-1] != q2[i-1] || q1[i-1] != q3[i-1])
                {
                    flint_printf("quotient error at index %wu\n", i-1);
                    flint_printf("qn=%wu, an=%wu, bn=%wu\n", qn, an, bn);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            for (i = bn-1; i > 0; i--)
            {
                if (r1[i-1] != r2[i-1] || r1[i-1] != r3[i-1])
                {
                    flint_printf("remainder error at index %wu\n", i-1);
                    flint_printf("r1[i]=%wu, r2[i]=%wu, bn=%wu\n", r1[i-1], r2[i-1]);
                    flint_printf("qn=%wu, an=%wu, bn=%wu\n", qn, an, bn);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            flint_free(a);
            flint_free(b);
            flint_free(q1);
            flint_free(q2);
            flint_free(r1);
            flint_free(r2);
        }

        for (i = 0; i < 16; i++)
            flint_printf("%c", '\b');
        fflush(stdout);
    }

    mpn_ctx_clear(R);

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
#else

int main()
{
    return 0;
}

#endif
