/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "arb.h"
#include "radix.h"
#include "double_extras.h"
#include "mag.h"
#include "profiler.h"
#include "hypgeom.h"
#include "fmpz_poly.h"

/* Don't call arb_const_e because it caches the result. */
static void
arb_const_e_eval(arb_t s, slong prec)
{
    hypgeom_t series;
    arb_t t;

    arb_init(t);
    hypgeom_init(series);

    fmpz_poly_set_str(series->A, "1  1");
    fmpz_poly_set_str(series->B, "1  1");
    fmpz_poly_set_str(series->P, "1  1");
    fmpz_poly_set_str(series->Q, "2  0 1");

    prec += FLINT_CLOG2(prec);
    arb_hypgeom_infsum(s, t, series, prec, prec);
    arb_div(s, s, t, prec);

    hypgeom_clear(series);
    arb_clear(t);
}

char * radix_integer_get_str_decimal(char * s, const radix_integer_t x, const radix_t radix)
{
    return radix_get_str_decimal(s, x->d, FLINT_ABS(x->size), x->size < 0, radix);
}

slong radix_digits_to_limbs(slong n, const radix_t radix)
{
    return (n + radix->exp - 1) / radix->exp;
}

void
print_condensed(const char * s, slong ndigits, slong n)
{
    char * tmp = flint_malloc(n + 1);

    strncpy(tmp, s, 1);
    tmp[1] = '\0';
    flint_printf("%s.", tmp);

    strncpy(tmp, s + 1, n);
    tmp[n] = '\0';

    flint_printf("%s{...%wd digits...}%s\n", tmp, ndigits - 2 * n - 1, s + ndigits - n);

    flint_free(tmp);
}

/* Calculate sqrt(a) as a fixed-point number with n fraction limbs. */
void
radix_integer_fixed_sqrt_1(radix_integer_t res, ulong a, slong n, const radix_t radix)
{
    nn_ptr rd;
    rd = radix_integer_fit_limbs(res, n + 1, radix);
    radix_rsqrt_1_approx(rd, a, n, radix);
    if (a == 2)
        rd[n] = radix_mul_two(rd, rd, n, radix);
    else
        rd[n] = radix_mul_1(rd, rd, n, a, radix);
    res->size = n + 1;
}

#define BSPLIT_E_MPN_BASECASE 2

static void
bsplit_e(radix_integer_t P, radix_integer_t Q, slong a, slong b, const radix_t radix)
{
    if (b - a == 1)
    {
        radix_integer_one(P, radix);
        radix_integer_set_ui(Q, b, radix);
    }
    else if (b - a == 2 && b < (UWORD(1) << (FLINT_BITS / 2)))
    {
        radix_integer_set_ui(P, a + 3, radix);
        radix_integer_set_ui(Q, (a + 1) * (a + 2), radix);
    }
    else if (b - a <= BSPLIT_E_MPN_BASECASE && b < (UWORD(1) << (FLINT_BITS / 2)))
    {
        ulong p[BSPLIT_E_MPN_BASECASE];
        ulong q[BSPLIT_E_MPN_BASECASE];
        ulong p2, q2;
        slong pn, qn, k;

        p[0] = a + 3;
        q[0] = (a + 1) * (a + 2);
        pn = qn = 1;

        for (k = a + 2; k < b; k++)
        {
            if (k + 1 < b) { p2 = k + 3; q2 = (k + 1) * (k + 2); k++; }
            else           { p2 = 1; q2 = k + 1; }

            p[pn] = mpn_mul_1(p, p, pn, q2);
            p[pn] += mpn_add_1(p, p, pn, p2);
            pn += (p[pn] != 0);
            q[qn] = mpn_mul_1(q, q, qn, q2);
            qn += (q[qn] != 0);
        }

        P->size = radix_set_mpn(radix_integer_fit_limbs(P,
            radix_set_mpn_need_alloc(pn, radix), radix), p, pn, radix);
        Q->size = radix_set_mpn(radix_integer_fit_limbs(Q,
            radix_set_mpn_need_alloc(qn, radix), radix), q, qn, radix);
    }
    else
    {
        slong m = a + (b - a) / 2;

        radix_integer_t P1, Q1, P2, Q2;

        radix_integer_init(P1, radix);
        radix_integer_init(Q1, radix);
        radix_integer_init(P2, radix);
        radix_integer_init(Q2, radix);

        bsplit_e(P1, Q1, a, m, radix);
        bsplit_e(P2, Q2, m, b, radix);

        radix_integer_mul(P, P1, Q2, radix);
        radix_integer_add(P, P, P2, radix);
        radix_integer_mul(Q, Q1, Q2, radix);

        radix_integer_clear(P1, radix);
        radix_integer_clear(Q1, radix);
        radix_integer_clear(P2, radix);
        radix_integer_clear(Q2, radix);
    }
}

void
radix_integer_fixed_e(radix_integer_t res, slong n, const radix_t radix)
{
    slong N;
    radix_integer_t P, Q;

    /* Select truncation N */
    {
        mag_t one, eps, err;

        double logB = log(LIMB_RADIX(radix));
        N = n * logB / d_lambertw(n * logB * exp(-1.0)) + 1;

        mag_init(one);
        mag_init(eps);
        mag_init(err);

        mag_one(one);
        mag_set_ui(eps, LIMB_RADIX(radix));
        mag_inv_lower(eps, eps);
        mag_pow_ui_lower(eps, eps, n);

        while (1)
        {
            mag_exp_tail(err, one, N + 1);
            if (mag_cmp(err, eps) <= 0)
                break;
            else
                N++;
        }

        mag_clear(one);
        mag_clear(err);
        mag_clear(eps);
    }

    radix_integer_init(P, radix);
    radix_integer_init(Q, radix);

    bsplit_e(P, Q, 1, N, radix);
    radix_integer_lshift_limbs(P, P, n, radix);
    radix_integer_tdiv_q(res, P, Q, radix);
    radix_integer_set_limb(res, res, n, 2, radix);

    radix_integer_clear(P, radix);
    radix_integer_clear(Q, radix);
}

int main()
{
    slong ndigits;
    radix_t radix;
    char * s;

    radix_init(radix, 10, 0);

    flint_set_num_threads(1);

    for (ndigits = 100000; ndigits <= 1000000000; ndigits *= 10)
    {
        radix_integer_t x;
        arb_t y;

        slong prec_bits = 3.32192809488736 * ndigits + 1;
        slong prec_radix = radix_digits_to_limbs(ndigits, radix);

        flint_printf("\n10^%wd digits (%wd bits; %wd decimal limbs = %wd digits)\n",
            n_clog(ndigits, 10),
            prec_bits, prec_radix, prec_radix * radix->exp);

        arb_init(y);
        flint_printf("Arb compute sqrt(2):     ");
        TIMEIT_START;
        arb_sqrt_ui(y, 2, prec_bits);
        TIMEIT_STOP;

        flint_printf("Arb to string:           ");
        s = NULL;
        TIMEIT_START;
        flint_free(s);
        s = arb_get_str(y, ndigits, ARB_STR_CONDENSE * 20);
        TIMEIT_STOP;

        flint_printf("%s\n", s);
        arb_clear(y);

        radix_integer_init(x, radix);
        flint_printf("Radix compute sqrt(2):   ");
        TIMEIT_START;
        radix_integer_fixed_sqrt_1(x, 2, prec_radix, radix);
        TIMEIT_STOP;

        flint_printf("Radix to string:         ");
        s = NULL;
        TIMEIT_START;
        flint_free(s);
        s = radix_integer_get_str_decimal(NULL, x, radix);
        TIMEIT_STOP;

        print_condensed(s, ndigits, 20);
        radix_integer_clear(x, radix);
    }

    for (ndigits = 100000; ndigits <= 1000000000; ndigits *= 10)
    {
        radix_integer_t x;
        arb_t y;

        slong prec_bits = 3.32192809488736 * ndigits + 1;
        slong prec_radix = radix_digits_to_limbs(ndigits, radix);

        flint_printf("\n10^%wd digits (%wd bits; %wd decimal limbs = %wd digits)\n",
            n_clog(ndigits, 10),
            prec_bits, prec_radix, prec_radix * radix->exp);

        arb_init(y);
        flint_printf("Arb compute e:     ");
        TIMEIT_START;
        arb_const_e_eval(y, prec_bits);
        TIMEIT_STOP;

        flint_printf("Arb to string:     ");
        s = NULL;
        TIMEIT_START;
        flint_free(s);
        s = arb_get_str(y, ndigits, ARB_STR_CONDENSE * 20);
        TIMEIT_STOP;

        flint_printf("%s\n", s);
        arb_clear(y);

        radix_integer_init(x, radix);
        flint_printf("Radix compute e:   ");
        TIMEIT_START;
        radix_integer_fixed_e(x, prec_radix, radix);
        TIMEIT_STOP;

        flint_printf("Radix to string:   ");
        s = NULL;
        TIMEIT_START;
        flint_free(s);
        s = radix_integer_get_str_decimal(NULL, x, radix);
        TIMEIT_STOP;

        print_condensed(s, ndigits, 20);
        radix_integer_clear(x, radix);
    }

    radix_clear(radix);
    flint_cleanup_master();
    return 0;
}

