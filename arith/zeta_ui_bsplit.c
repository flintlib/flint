#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "arith.h"

typedef struct
{
    fmpz_t A;
    fmpz_t B;
    fmpz_t C;
    fmpz_t D;
    fmpz_t E;
    fmpz_t Q1;
    fmpz_t Q2;
    fmpz_t Q3;
}
zeta_bsplit_state;

static __inline__ void
zeta_bsplit_init(zeta_bsplit_state * S)
{
    fmpz_init(S->A);
    fmpz_init(S->B);
    fmpz_init(S->C);
    fmpz_init(S->D);
    fmpz_init(S->E);
    fmpz_init(S->Q1);
    fmpz_init(S->Q2);
    fmpz_init(S->Q3);
}

static __inline__ void
zeta_bsplit_clear(zeta_bsplit_state * S)
{
    fmpz_clear(S->A);
    fmpz_clear(S->B);
    fmpz_clear(S->C);
    fmpz_clear(S->D);
    fmpz_clear(S->E);
    fmpz_clear(S->Q1);
    fmpz_clear(S->Q2);
    fmpz_clear(S->Q3);
}


static __inline__ void
zeta_coeff_k(zeta_bsplit_state * S, long k, long n, long s)
{
    if (k + 1 < 0)
    {
        fmpz_set_si(S->D, 1);
        fmpz_set_si(S->Q1, 1);
    }
    else if (k + 1 > n)
    {
        fmpz_zero(S->D);
        fmpz_set_si(S->Q1, 1);
    }
    else
    {
        fmpz_set_si(S->D, 2 * (n + (k + 1) - 1));
        fmpz_mul_si(S->D, S->D, n + 1 - (k + 1));
        fmpz_set_si(S->Q1, k + 1);
        fmpz_mul_si(S->Q1, S->Q1, 2*(k + 1) - 1);
    }

    if (k - 1 < 0)
    {
        fmpz_zero(S->E);
        fmpz_set_si(S->Q2, 1);
    }
    else if (k - 1 >= n)
    {
        fmpz_set_si(S->E, 1);
        fmpz_set_si(S->Q2, 1);
    }
    else
    {
        fmpz_set_si(S->E, ((k - 1) % 2) ? -1 : 1);
        fmpz_set_si(S->Q2, k);
        fmpz_pow_ui(S->Q2, S->Q2, s);
    }

    fmpz_mul(S->Q3, S->Q1, S->Q2);
    fmpz_mul(S->A, S->E, S->Q1);
    fmpz_zero(S->B);
    fmpz_set(S->C, S->Q1);
}

void
zeta_bsplit(zeta_bsplit_state * L, long a, long b, long n, long s, int resumable)
{
    if (a + 1 == b)
    {
        zeta_coeff_k(L, a, n, s);
    }
    else
    {
        zeta_bsplit_state Rs;
        zeta_bsplit_state * R = &Rs;

        long m = (a + b) / 2;

        zeta_bsplit(L, m, b, n, s, 1);

        zeta_bsplit_init(R);
        zeta_bsplit(R, a, m, n, s, 1);

        fmpz_mul(L->E, L->E, R->Q2);
        fmpz_addmul(L->E, R->E, L->Q2);

        fmpz_mul(L->B, L->B, R->D);
        fmpz_addmul(L->B, L->A, R->C);
        fmpz_mul(L->B, L->B, R->Q2);
        fmpz_addmul(L->B, R->B, L->Q3);

        if (resumable)
        {
            fmpz_mul(L->A, L->A, R->Q3);
            fmpz_addmul(L->A, R->A, L->Q3);
        }

        fmpz_mul(L->C, L->C, R->D);
        fmpz_addmul(L->C, R->C, L->Q1);
        fmpz_mul(L->Q2, L->Q2, R->Q2);

        if (resumable)
        {
            fmpz_mul(L->D, L->D, R->D);
            fmpz_mul(L->Q1, L->Q1, R->Q1);
            fmpz_mul(L->Q3, L->Q3, R->Q3);
        }

        zeta_bsplit_clear(R);
    }
}

static __inline__ void
mpfr_set_fmpz(mpfr_t c, const fmpz_t b, mpfr_rnd_t rnd)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_set_z(c, COEFF_TO_PTR(*b), rnd);
    else
        mpfr_set_si(c, *b, rnd);
}

static __inline__ void
mpfr_mul_fmpz(mpfr_t c, mpfr_srcptr a, const fmpz_t b, mpfr_rnd_t rnd)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_mul_z(c, a, COEFF_TO_PTR(*b), rnd);
    else
        mpfr_mul_si(c, a, *b, rnd);
}

static __inline__ void
mpfr_div_fmpz(mpfr_t c, mpfr_srcptr a, const fmpz_t b, mpfr_rnd_t rnd)
{
    if (COEFF_IS_MPZ(*b))
        mpfr_div_z(c, a, COEFF_TO_PTR(*b), rnd);
    else
        mpfr_div_si(c, a, *b, rnd);
}

void
mpfr_zeta_ui_bsplit(mpfr_t x, ulong s, mpfr_rnd_t rnd)
{
    zeta_bsplit_state sm;
    zeta_bsplit_state * sum = &sm;
    long prec, wp, n;
    mpfr_t t, u;
    mpz_t w;
    long i;

    if (s < 2)
    {
        printf("Exception (mpfr_zeta_ui_bsplit). s < 2.\n");
        abort();
    }

    prec = mpfr_get_prec(x);
    /* 1/log(3+sqrt(8)) = 0.39321985... */
    n = 0.39322 * prec + 10;
    wp = prec + 30;

    zeta_bsplit_init(sum);
    zeta_bsplit(sum, 0, n + 1, n, s, 0);

    mpfr_init2(t, wp);
    mpfr_init2(u, wp);

    /* (W*C - B) / (Q2 * C) */
    mpfr_set_fmpz(t, sum->E, MPFR_RNDD);
    mpfr_set_fmpz(u, sum->C, MPFR_RNDD);
    mpfr_mul(t, t, u, MPFR_RNDD);
    mpfr_set_fmpz(u, sum->B, MPFR_RNDD);
    mpfr_sub(t, t, u, MPFR_RNDD);

    mpfr_set_fmpz(u, sum->Q2, MPFR_RNDD);
    mpfr_mul_fmpz(u, u, sum->C, MPFR_RNDD);
    mpfr_div(t, t, u, MPFR_RNDD);

    /* Multiply by 1/(1 - 2^(1-s)) */
    mpz_init(w);

    for (i = wp; i >= 0; i -= (s - 1))
        mpz_setbit(w, i);

    mpfr_set_z_2exp(u, w, -wp, MPFR_RNDD);
    mpfr_mul(x, t, u, rnd);

    mpz_clear(w);
    mpfr_clear(t);
    mpfr_clear(u);
    zeta_bsplit_clear(sum);
}
