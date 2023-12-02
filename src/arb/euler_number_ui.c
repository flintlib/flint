/*
    Copyright (C) 2016, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "arb.h"
#include "acb_dirichlet.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

#if FLINT64
#define ARB_EULER_NUMBER_TAB_SIZE 25
#else
#define ARB_EULER_NUMBER_TAB_SIZE 15
#endif

static const ulong arb_euler_number_tab[] =
{
    1, 1, 5, 61, 1385, 50521, 2702765, 199360981,
#if FLINT64
    UWORD(19391512145), UWORD(2404879675441), UWORD(370371188237525),
    UWORD(69348874393137901), UWORD(15514534163557086905)
#endif
};

static double
arb_euler_number_mag(double n)
{
    double x;
    x = n + 2;
    x += ((n + 1) * log(n + 1) - n) * 1.44269504088897;  /* 1/log(2) */
    x -= 1.6514961294723*(n+1);  /* log2(pi) */
    return x;
}

void
arb_euler_number_ui_beta(arb_t res, ulong n, slong prec)
{
    slong pi_prec;
    arb_t t;
    const signed char chi[4] = {0, 1, 0, -1};

    pi_prec = prec + 2 * FLINT_BIT_COUNT(n);
    arb_init(t);

    /* |E_n| = 2 n! beta(n+1) / (pi/2)^(n+1) */
    arb_const_pi(t, pi_prec);
    arb_mul_2exp_si(t, t, -1);
    arb_pow_ui(t, t, n + 1, pi_prec);

    _acb_dirichlet_euler_product_real_ui(res, n + 1, chi, 4, 1, prec);

    arb_mul(res, res, t, prec);
    arb_fac_ui(t, n, pi_prec);  /* todo: prec should be enough */
    arb_div(res, t, res, prec);
    arb_mul_2exp_si(res, res, 1);

    if (n % 4 == 2)
        arb_neg(res, res);

    arb_clear(t);
}

#define LOW_MASK ((UWORD(1) << (FLINT_BITS / 2)) - 1)

typedef struct
{
    ulong n;
    ulong ninv;
    ulong F;
}
nmod_redc_t;

static void
nmod_redc_init(nmod_redc_t * rmod, nmod_t mod)
{
    ulong n, ninv2;
    int bits;

    n = mod.n;
    rmod->n = n;
    rmod->F = (UWORD(1) << (FLINT_BITS / 2));
    NMOD_RED(rmod->F, rmod->F, mod);

    /* Newton's method for 2-adic inversion */
    ninv2 = -n;   /* already correct mod 8 */
    for (bits = 3; bits < FLINT_BITS / 2; bits *= 2)
        ninv2 = 2*ninv2 + n * ninv2 * ninv2;

    rmod->ninv = ninv2 & LOW_MASK;
}

static inline ulong
nmod_redc_fast(ulong x, ulong n, ulong ninv2)
{
    ulong y = (x * ninv2) & LOW_MASK;
    ulong z = x + (n * y);
    return z >> (FLINT_BITS / 2);
}

static inline ulong
nmod_redc(ulong x, ulong n, ulong ninv2)
{
    ulong y = nmod_redc_fast(x, n, ninv2);
    if (y >= n)
       y -= n;
    return y;
}

static ulong
nmod_redc_mul_fast(ulong a, ulong b, nmod_redc_t rmod)
{
    return nmod_redc_fast(a * b, rmod.n, rmod.ninv);
}

static ulong
nmod_redc_mul(ulong a, ulong b, nmod_redc_t rmod)
{
    return nmod_redc(a * b, rmod.n, rmod.ninv);
}


static ulong
nmod_to_redc(ulong x, nmod_t mod, nmod_redc_t rmod)
{
    return nmod_mul(x, rmod.F, mod);
}

static ulong
nmod_from_redc(ulong x, nmod_redc_t rmod)
{
    return nmod_redc(x, rmod.n, rmod.ninv);
}

ulong
nmod_redc_pow_ui(ulong a, ulong exp, nmod_redc_t rmod)
{
    ulong x;

    while ((exp & 1) == 0)
    {
        a = nmod_redc_mul(a, a, rmod);
        exp >>= 1;
    }

    x = a;

    while (exp >>= 1)
    {
        a = nmod_redc_mul(a, a, rmod);

        if (exp & 1)
            x = nmod_redc_mul(x, a, rmod);
    }

    return x;
}

ulong
euler_mod_p_powsum_1(ulong n, ulong p)
{
    slong j;
    ulong s, t;
    nmod_t mod;

    if (n % 2 == 1)
        return 0;

    n = n % (p - 1);

    if (n == 0)
        return UWORD_MAX;

    nmod_init(&mod, p);
    s = 1;

    for (j = 3; j <= p - 2; j += 2)
    {
        t = nmod_pow_ui(j, n, mod);
        s = nmod_sub(t, s, mod);
    }

    if (p % 4 == 1)
        s = nmod_neg(s, mod);

    s = nmod_add(s, s, mod);
    return s;
}

ulong
euler_mod_p_powsum_noredc(ulong n, ulong p, const unsigned int * divtab)
{
    unsigned int * pows;
    slong i, N, horner_point;
    ulong s, t, z;
    ulong v2n, power_of_two;
    nmod_t mod;
    TMP_INIT;

    if (n % 2 == 1)
        return 0;

    n = n % (p - 1);

    if (n == 0)
        return UWORD_MAX;

    N = p / 4;

    nmod_init(&mod, p);

    TMP_START;
    pows = TMP_ALLOC(sizeof(unsigned int) * (N / 3 + 1));

    s = z = 0;

    /* Evaluate as a polynomial in 2^n */
    power_of_two = 1;
    while (power_of_two * 2 <= N)
        power_of_two *= 2;

    horner_point = 1;
    v2n = nmod_pow_ui(2, n, mod);

    for (i = 1; i <= N / 3; i += 2)
    {
        if (divtab[i] == 1)
            t = nmod_pow_ui(i, n, mod);
        else
            t = nmod_mul(pows[divtab[i]], pows[divtab[i + 1]], mod);

        pows[i] = t;
        s = nmod_add(s, t, mod);

        if (i == horner_point)
        {
            while (i == horner_point && power_of_two != 1)
            {
                z = nmod_add(s, nmod_mul(v2n, z, mod), mod);
                power_of_two /= 2;
                horner_point = N / power_of_two;
                if (horner_point % 2 == 0)
                    horner_point--;
            }
        }
    }

    /* Same as above, but here we don't need to write the powers. */
    for ( ; i <= N; i += 2)
    {
        if (divtab[i] == 1)
            t = nmod_pow_ui(i, n, mod);
        else
            t = nmod_mul(pows[divtab[i]], pows[divtab[i + 1]], mod);

        s = nmod_add(s, t, mod);

        if (i == horner_point)
        {
            while (i == horner_point && power_of_two != 1)
            {
                z = nmod_add(s, nmod_mul(v2n, z, mod), mod);
                power_of_two /= 2;
                horner_point = N / power_of_two;
                if (horner_point % 2 == 0)
                    horner_point--;
            }
        }
    }

    s = nmod_add(s, nmod_mul(v2n, z, mod), mod);

    if (p % 4 == 3)
        s = nmod_neg(s, mod);

    t = nmod_inv(nmod_pow_ui(4, p - n - 2, mod), mod);
    s = nmod_mul(s, t, mod);

    TMP_END;

    return s;
}

ulong
euler_mod_p_powsum_redc(ulong n, ulong p, const unsigned int * divtab)
{
    unsigned int * pows;
    slong i, N, horner_point;
    ulong s, t, z;
    ulong v2n, power_of_two;
    nmod_t mod;
    nmod_redc_t rmod;
    TMP_INIT;

    if (n % 2 == 1)
        return 0;

    n = n % (p - 1);

    if (n == 0)
        return UWORD_MAX;

    N = p / 4;

    nmod_init(&mod, p);
    nmod_redc_init(&rmod, mod);

    TMP_START;
    pows = TMP_ALLOC(sizeof(unsigned int) * (N / 3 + 1));

    s = z = 0;

    /* Evaluate as a polynomial in 2^n */
    power_of_two = 1;
    while (power_of_two * 2 <= N)
        power_of_two *= 2;

    horner_point = 1;
    v2n = nmod_redc_pow_ui(nmod_to_redc(2, mod, rmod), n, rmod);

    for (i = 1; i <= N / 3; i += 2)
    {
        if (divtab[i] == 1)
            t = nmod_redc_pow_ui(nmod_to_redc(i, mod, rmod), n, rmod);
        else
            t = nmod_redc_mul(pows[divtab[i]], pows[divtab[i + 1]], rmod);

        pows[i] = t;
        s += t;

        if (i == horner_point)
        {
            while (i == horner_point && power_of_two != 1)
            {
                NMOD_RED(s, s, mod);
                z = nmod_add(s, nmod_redc_mul(v2n, z, rmod), mod);
                power_of_two /= 2;
                horner_point = N / power_of_two;
                if (horner_point % 2 == 0)
                    horner_point--;
            }
        }
    }

    /* Same as above, but here we don't need to write the powers. */
    for ( ; i <= N; i += 2)
    {
        if (divtab[i] == 1)
            t = nmod_redc_pow_ui(nmod_to_redc(i, mod, rmod), n, rmod);
        else
            t = nmod_redc_mul_fast(pows[divtab[i]], pows[divtab[i + 1]], rmod);

        s += t;

        if (i == horner_point)
        {
            while (i == horner_point && power_of_two != 1)
            {
                NMOD_RED(s, s, mod);
                z = nmod_add(s, nmod_redc_mul(v2n, z, rmod), mod);
                power_of_two /= 2;
                horner_point = N / power_of_two;
                if (horner_point % 2 == 0)
                    horner_point--;
            }
        }
    }

    NMOD_RED(s, s, mod);
    s = nmod_add(s, nmod_redc_mul(v2n, z, rmod), mod);
    s = nmod_from_redc(s, rmod);

    if (p % 4 == 3)
        s = nmod_neg(s, mod);

    t = nmod_inv(nmod_pow_ui(4, p - n - 2, mod), mod);
    s = nmod_mul(s, t, mod);

    TMP_END;

    return s;
}

ulong
euler_mod_p_powsum(ulong n, ulong p, const unsigned int * divtab)
{
    if (p < (UWORD(1) << (FLINT_BITS / 2 - 1)))
        return euler_mod_p_powsum_redc(n, p, divtab);
    else
        return euler_mod_p_powsum_noredc(n, p, divtab);
}

static void
divisor_table_odd(unsigned int * tab, slong len)
{
    slong i, j;

    tab[0] = 0;

    for (i = 1; i < len; i += 2)
    {
        tab[i] = 1;
        tab[i + 1] = i;
    }

    for (i = 3; i < len; i += 2)
    {
        for (j = 3; j <= i && i * j < len; j += 2)
        {
            tab[i * j]     = j;
            tab[i * j + 1] = i;
        }
    }
}

typedef struct
{
    ulong n;
    const unsigned int * divtab;
    mp_ptr primes;
    mp_ptr residues;
}
mod_p_param_t;

static void
mod_p_worker(slong i, void * param)
{
    mod_p_param_t * p = (mod_p_param_t *) param;

    p->residues[i] = euler_mod_p_powsum(p->n, p->primes[i], p->divtab);
}

#define TIMING 0
#define DEBUG 0

/* todo: optimize basecase and move to flint */
void
_arb_tree_crt(fmpz_t r, fmpz_t m, mp_srcptr residues, mp_srcptr primes, slong len);

void
arb_fmpz_euler_number_ui_multi_mod(fmpz_t num, ulong n, double alpha)
{
    n_primes_t prime_iter;
    slong i, bits, mod_bits, zeta_bits, num_primes;
    ulong p;
    mp_ptr primes, residues;
    mag_t primes_product;
    unsigned int * divtab_odd;
    fmpz_t M;
#if TIMING
    double t1, t2;
#endif

    if (n % 2 != 0)
    {
        fmpz_zero(num);
        return;
    }

    if (n < ARB_EULER_NUMBER_TAB_SIZE)
    {
        if (n % 4 == 0)
            fmpz_set_ui(num, arb_euler_number_tab[n / 2]);
        else
            fmpz_neg_ui(num, arb_euler_number_tab[n / 2]);
        return;
    }

    if (alpha < 0)
    {
        if (n < 2000)
            alpha = 0.0;
        else if (n < 6000)
            alpha = 0.002 + 1.0e-5 * (n - 2000);
        else if (n < 90000)
            alpha = -0.132 + 0.02 * log(n);
        else
            alpha = FLINT_MIN(0.0085 * log(n), 0.11);
    }

#if TIMING
    t1 = clock();
#endif

    bits = arb_euler_number_mag(n) + 2;
    mod_bits = bits * alpha;
    zeta_bits = bits - mod_bits;

    num_primes = 0;
    mag_init(primes_product);
    mag_one(primes_product);

    n_primes_init(prime_iter);

    p = 5;
    n_primes_jump_after(prime_iter, 5);

    for ( ; mag_cmp_2exp_si(primes_product, mod_bits) < 0; p = n_primes_next(prime_iter))
    {
        if (n % (p - 1) != 0)
        {
            mag_mul_ui_lower(primes_product, primes_product, p);
            num_primes++;
        }
    }

#if DEBUG
    printf("\nn = %lu, bits = %lu, num_primes = %ld\n", n, bits, num_primes);
#endif

    primes = flint_malloc(sizeof(mp_limb_t) * num_primes);
    residues = flint_malloc(sizeof(mp_limb_t) * num_primes);

    p = 5;
    n_primes_jump_after(prime_iter, 5);

    for (i = 0; i < num_primes; p = n_primes_next(prime_iter))
    {
        if (n % (p - 1) != 0)
        {
            primes[i] = p;
            i++;
        }
    }

    n_primes_clear(prime_iter);

    if (num_primes == 0)
    {
        divtab_odd = NULL;
    }
    else
    {
        p = primes[num_primes - 1];
        divtab_odd = flint_malloc(sizeof(unsigned int) * (p / 4 + 2));
        divisor_table_odd(divtab_odd, p / 4 + 1);
    }

#if TIMING
    t2 = clock();
    printf("init time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
    printf("num_primes = %ld\n", num_primes);
#endif

    {
        mod_p_param_t param;
        param.n = n;
        param.primes = primes;
        param.residues = residues;
        param.divtab = divtab_odd;

        flint_parallel_do(mod_p_worker, &param, num_primes, 0, FLINT_PARALLEL_STRIDED /* | FLINT_PARALLEL_VERBOSE */);
    }


#if TIMING
    t2 = clock();
    printf("mod time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
    printf("start CRT\n");
    t1 = clock();
#endif

    fmpz_init(M);
    _arb_tree_crt(num, M, residues, primes, num_primes);
    fmpz_mod(num, num, M);

    if (n % 4 == 2)
    {
        fmpz_sub(num, M, num);
        fmpz_neg(num, num);
    }

#if TIMING
    printf("end CRT\n");
    t2 = clock();
    printf("CRT time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
    t1 = clock();
#endif

    if (zeta_bits > 0)
    {
        slong prec;
        arb_t b;
        fmpz_t t;

        arb_init(b);
        fmpz_init(t);

        for (prec = zeta_bits + 10; ; prec += 32)
        {
            arb_euler_number_ui_beta(b, n, prec);
            arb_sub_fmpz(b, b, num, prec);
            arb_div_fmpz(b, b, M, prec);

            if (arb_get_unique_fmpz(t, b))
            {
                fmpz_addmul(num, t, M);
                break;
            }

            flint_printf("euler: n = %wu, bits = %wd, mod = %wd, zeta = %wd: get_unique_fmpz failed!\n", n, bits, mod_bits, zeta_bits);
        }

        arb_clear(b);
        fmpz_clear(t);
    }

#if TIMING
    printf("end zeta\n");
    t2 = clock();
    printf("zeta time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
#endif

    flint_free(primes);
    flint_free(residues);
    flint_free(divtab_odd);
    fmpz_clear(M);
    mag_clear(primes_product);
}

void
arb_fmpz_euler_number_ui(fmpz_t res, ulong n)
{
    arb_t x;
    double mag;

    if (n % 2 != 0)
    {
        fmpz_zero(res);
        return;
    }

    if (n < ARB_EULER_NUMBER_TAB_SIZE)
    {
        if (n % 4 == 0)
            fmpz_set_ui(res, arb_euler_number_tab[n / 2]);
        else
            fmpz_neg_ui(res, arb_euler_number_tab[n / 2]);
        return;
    }

    if (n < 2000)
    {
        mag = arb_euler_number_mag(n);

        arb_init(x);
        arb_euler_number_ui_beta(x, n, mag + 5);
        if (!arb_get_unique_fmpz(res, x))
        {
            flint_throw(FLINT_ERROR, "arb_fmpz_euler_number_ui: unexpected inaccuracy\n");
        }
        arb_clear(x);
    }
    else
    {
        arb_fmpz_euler_number_ui_multi_mod(res, n, -1.0);
    }
}

void
arb_euler_number_ui(arb_t res, ulong n, slong prec)
{
    double mag;

    if (n % 2 != 0)
    {
        arb_zero(res);
        return;
    }

    if (n < ARB_EULER_NUMBER_TAB_SIZE)
    {
        arb_set_ui(res, arb_euler_number_tab[n / 2]);
        if (n % 4 == 2)
            arb_neg(res, res);
        arb_set_round(res, res, prec);
        return;
    }

    mag = arb_euler_number_mag(n);

    if (prec > 0.9 * mag)
    {
        fmpz_t t;
        fmpz_init(t);
        arb_fmpz_euler_number_ui(t, n);
        arb_set_round_fmpz(res, t, prec);
        fmpz_clear(t);
    }
    else
    {
        arb_euler_number_ui_beta(res, n, prec + 5);
        arb_set_round(res, res, prec);
    }
}
