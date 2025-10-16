/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

static ulong
n_multiplicative_order(ulong x, ulong p, ulong pinv, n_factor_t * p1fac)
{
    ulong m, q, mm;
    slong i;
    m = p - 1;

    for (i = 0; i < p1fac->num; i++)
    {
        q = p1fac->p[i];

        while (m % q == 0)
        {
            mm = m / q;
            if (n_powmod2_preinv(x, mm, p, pinv) != 1)
                break;
            m = mm;
        }
    }

    return m;
}

int
nmod_poly_irreducible_binomial(nmod_poly_t res, ulong n)
{
    ulong q = res->mod.n;
    ulong e, a, qinv;
    n_factor_t nfac, q1fac;
    slong i;

    if (n == 0)
        return 0;

    if (n == 1)
    {
        nmod_poly_one(res);
        nmod_poly_set_coeff_ui(res, 1, 1);
        return 1;
    }

    /* There is an irreducible binomial x^n + a mod q iff
       p | q - 1 for each prime p | n and
       if 4 | n then 4 | q - 1. */
    if (n % 4 == 0)
        if (q % 4 != 1)
            return 0;
    n_factor_init(&nfac);
    n_factor(&nfac, n, 1);
    for (i = 0; i < nfac.num; i++)
        if ((q - 1) % nfac.p[i] != 0)
            return 0;

    nmod_poly_zero(res);
    nmod_poly_set_coeff_ui(res, n, 1);

    /* Let e = ord_q(-a). Then we have irreducibility iff
       (p divides e  but  p does not divide (q-1)/e)
       for each prime factor p of n. */
    n_factor_init(&q1fac);
    n_factor(&q1fac, q - 1, 1);
    qinv = n_preinvert_limb(q);

    for (a = 1; a < q; a++)
    {
        e = n_multiplicative_order(q - a, q, qinv, &q1fac);

        for (i = 0; i < nfac.num; i++)
            if (!(e % nfac.p[i] == 0 && ((q - 1) / e) % nfac.p[i] != 0))
                goto next;

        nmod_poly_set_coeff_ui(res, 0, a);
        return 1;

        next:
            continue;
    }

    flint_throw(FLINT_ERROR, "Failed to construct irreducible binomial (n = %wu, p = %wu)", n, q);
}



/* Todo: refine and extend this table */
static const unsigned int sieve_limit_tab[][10] = {
     /* p     1    2     3     4     5      6      7     8     */
    {   2,    1,   1,    1,    1,    24,    24,    256,  256,  0 },
    {   3,    1,   1,    1,    1,    128,   512,   512,  3072, 0 },
    {   5,    1,   1,    1,    48,   512,   1536,  3072, 0 },
    {   7,    1,   1,    1,    64,   1536,  3072,  0 },
    {  11,    1,   1,    64,   768,  2048,  0 },
    {  13,    1,   1,    48,   768,  0 },
    {  17,    1,   1,    384,  768,  0 },
    {  19,    1,   16,   96,   1536, 0 },
    {  23,    1,   16,   256,  2048, 0 },
    {  29,    1,   12,   512,  0 },
    {  31,    1,   12,   384,  0 },
    {  37,    1,   12,   384,  0 },
    {  41,    1,   24,   1024, 0 },
    {  43,    1,   32,   512,  0 },
    {  47,    1,   32,   768,  0 },
    {  53,    1,   48,   1024, 0 },
    {  59,    1,   32,   1024, 0 },
    {  61,    1,   64,   1024, 0 },
    {  67,    1,   64,   1536, 0 },
    {  71,    1,   64,   1536, 0 },
    {  73,    1,   64,   1536, 0 },
    {  79,    1,   64,   1536, 0 },
    {  83,    1,   64,   2048, 0 },
    {  89,    1,   64,   2048, 0 },
    {  97,    4,   128,  2048, 0 },
    { 101,    8,   128,  3072, 0 },
    { 103,    4,   96,   3072, 0 },
    { 107,    4,   96,   0 },
    { 127,    4,   128,  0 },
    { 179,    4,   192,  0 },
    { 269,    8,   256,  0 },
    { 283,    8,   384,  0 },
    { 431,    8,   512,  0 },
    { 487,    8,   768,  0 },
    { 653,    12,  768,  0 },
    { 727,    16,  1024, 0 },
    { 953,    16,  1536, 0 },
    { 1009,   24,  2048, 0 },
    { 1549,   32, 0 },
    { 2333,   48, 0 },
    { 3877,   64, 0 },
    { 7789,   96, 0 },
    { 14107,  128, 0 },
    { 28069,  192, 0 },
    { 35831,  256, 0 },
    { 82457,  384, 0 },
    { 121949, 768, 0 },
    { 180247, 1024, 0 },
    { 0 },
};

/* Return b such that we want to do trial division by all irreducible
   factors up to degree b. */
static ulong trinomial_sieve_limit(ulong p, ulong n)
{
    slong i;
    ulong ptab, b, limit;

    for (i = 0; ; i++)
    {
        ptab = sieve_limit_tab[i][0];
        if (ptab == 0)
            return 0;

        if (ptab >= p)
        {
            limit = 0;
            for (b = 1; sieve_limit_tab[i][b] != 0; b++)
            {
                if (n >= sieve_limit_tab[i][b])
                    limit = b;
            }

            return limit;
        }
    }
}


/* Return the product of all irreducible polynomial of degree r. */
/* Well-known fact: x^(p^r) - x contains all factors of degree dividing r,
   so we just need to divide out the smaller products. */
void
nmod_poly_product_all_irreducibles_deg(nmod_poly_t res, ulong r)
{
    ulong k, d, step, p = res->mod.n;

    FLINT_ASSERT(r >= 1 && r <= 8);

    nmod_poly_zero(res);

    if (r == 1)
    {
        nmod_poly_set_coeff_ui(res, 1, 1);
        nmod_poly_set_coeff_ui(res, p, 1);
        return;
    }

    if (r == 4 || r == 6)
    {
        d = n_pow(p, r) - n_pow(p, 2);
        step = p * p - 1;
    }
    else if (r == 8)
    {
        d = n_pow(p, r) - n_pow(p, 4);
        step = n_pow(p, 4) - 1;
    }
    else
    {
        d = n_pow(p, r) - p;
        step = p - 1;
    }

    for (k = 0; k <= d; k += step)
        nmod_poly_set_coeff_ui(res, k, 1);

    if (r == 6)
    {
        nmod_poly_t t;
        nmod_poly_init(t, p);
        nmod_poly_product_all_irreducibles_deg(t, 3);
        nmod_poly_divexact(res, res, t);
        nmod_poly_clear(t);
    }
}

/* Hack: nmod_poly_gcd is not optimised for unbalanced GCD with the smaller
   operand sparse, so we deal with this here. */
void
_nmod_poly_inplace_rem_sparse_monic(nn_ptr R,
        slong lenA, nn_srcptr Bcoeffs, const slong * Bexps, slong nzB, slong lenB, nmod_t mod)
{
    slong i, j, k;
    slong n = lenB - 1;
    ulong c;

    for (i = lenA - 1; i >= n; i--)
    {
        c = R[i];

        for (k = nzB - 2; k >= 0; k--)
        {
            j = Bexps[k];
            R[j + i - n] = nmod_sub(R[j + i - n], nmod_mul(c, Bcoeffs[k], mod), mod);
        }
    }
}

void
nmod_poly_gcd_with_sparse(nmod_poly_t res, const nmod_poly_t A, const nmod_poly_t B,
        nn_srcptr Bcoeffs, const slong * Bexps, slong nzB)
{
    if (A->length <= B->length)
    {
        nmod_poly_gcd(res, A, B);
    }
    else
    {
        nn_ptr R;
        slong Rlen = B->length - 1;
        nmod_poly_t T;
        TMP_INIT;
        TMP_START;
        R = TMP_ALLOC(sizeof(ulong) * A->length);
        _nmod_vec_set(R, A->coeffs, A->length);
        _nmod_poly_inplace_rem_sparse_monic(R, A->length, Bcoeffs, Bexps, nzB, B->length, res->mod);
        NMOD_VEC_NORM(R, Rlen);
        T->coeffs = R;
        T->length = Rlen;
        T->alloc = Rlen;
        nmod_poly_gcd(res, T, B);
        TMP_END;
    }
}

/* Discriminant of x^n + ax^k + b  mod p > 2 (Swan's formula). */
ulong
_nmod_poly_trinomial_discriminant(ulong n, ulong k, ulong a, ulong b, nmod_t mod)
{
    ulong d, n1, k1;
    ulong D1, D2, D3;

    d = n_gcd(n, k);
    n1 = n / d;
    k1 = k / d;

    D1 = nmod_pow_ui(b, k - 1, mod);
    D2 = nmod_mul(nmod_pow_ui(n, n1, mod), nmod_pow_ui(b, n1 - k1, mod), mod);
    D3 = nmod_mul(nmod_pow_ui(n - k, n1 - k1, mod), nmod_pow_ui(k, k1, mod), mod);
    D3 = nmod_mul(D3, nmod_pow_ui(a, n1, mod), mod);
    if (n1 % 2)
        D3 = nmod_neg(D3, mod);
    D2 = nmod_sub(D2, D3, mod);
    D2 = nmod_pow_ui(D2, d, mod);
    D1 = nmod_mul(D1, D2, mod);
    if ((n * (n - 1) / 2) % 2)
        D1 = nmod_neg(D1, mod);

    return D1;
}

int
nmod_poly_irreducible_trinomial(nmod_poly_t res, ulong n)
{
    ulong p = res->mod.n;
    nmod_t mod = res->mod;
    ulong a, b, c, k, l;
    ulong a2, b2, e2n, cn, ck, cn2, ck2, D;
    slong i;
    nn_ptr invctab = NULL;
    int want_pruning;
    int found = 0;
    nmod_poly_t sieve[10], h;
    slong num_sieve = 0;
    ulong Bcoeffs[3];
    slong Bexps[3];
    ulong sieve_limit;

    if (n <= 1)
        return 0;

    if (p == 2 && (n % 8 == 0))
        return 0;

    /* Ciet, Quisquater and Sica, A Short Note on Irreducible Trinomials in Binary Fields */
    if (p == 2 && ((n % 24 == 13) || (n % 24 == 19)) && n_is_prime(n))
        return 0;

    sieve_limit = trinomial_sieve_limit(p, n);
    want_pruning = (sieve_limit > 0);

    if (want_pruning)
    {
        invctab = _nmod_vec_init(p);
        for (c = 1; c < p; c++)
            invctab[c] = nmod_inv(c, mod);

        nmod_poly_init_mod(h, mod);

        for (a = 2; a <= sieve_limit && a < n; a++)
        {
            nmod_poly_init_mod(sieve[num_sieve], mod);
            nmod_poly_product_all_irreducibles_deg(sieve[num_sieve], a);
            num_sieve++;
        }
    }

#define ALREADY_CHECKED(a,b,a2,b2) ((a2) < (a) || ((a2) == (a) && (b2) < (b)))

    /* Iterate over x^n + a x^k + b */
    for (a = 1; a < p; a++)
    {
        for (b = 1; b < p; b++)
        {
            /* Mirror symmetry: f(x) is irreducible iff x^n f(1/x) is irreducible */
            for (k = 1; k * 2 <= n; k++)
            {
                if (p > 2)
                {
                    /* See if we have already checked f(-x). */
                    /* n even, k even   x^n + ax^k + b ->  x^n + ax^k + b */
                    /* n even, k odd    x^n + ax^k + b ->  x^n - ax^k + b */
                    /* n odd,  k even   x^n + ax^k + b ->  x^n - ax^k - b */
                    /* n odd,  k odd    x^n + ax^k + b ->  x^n + ax^k - b */
                    b2 = (n % 2) ? p - b : b;
                    a2 = ((n ^ k) % 2) ? p - a : a;
                    if (ALREADY_CHECKED(a, b, a2, b2))
                        goto next_trinomial;
                }

                /* Swan-Stickelberger */
                if (p == 2)
                {
                    if (n % 2 == k % 2)
                        l = n - k;
                    else
                        l = k;

                    if (n % 2 == 0 && l % 2 == 1 && n != 2 * l && ((n * l) / 2) % 4 <= 1)
                        goto next_trinomial;
                    if (n % 2 == 1 && l % 2 == 0 && ((2*n) % l != 0) && (n % 8 == 3 || n % 8 == 5))
                        goto next_trinomial;
                    if (n % 2 == 1 && l % 2 == 0 && ((2*n) % l == 0) && (n % 8 == 1 || n % 8 == 7))
                        goto next_trinomial;
                }
                else
                {
                    D = _nmod_poly_trinomial_discriminant(n, k, a, b, mod);

                    if (nmod_pow_ui(D, (p - 1) / 2, mod) == ((n % 2 == 1) ? p - 1 : 1))
                        goto next_trinomial;
                }

                if (want_pruning && p > 2)
                {
                    for (c = 1; 2 * c < p; c++)
                    {
                        ck = nmod_pow_ui(c, k, mod);
                        cn = nmod_mul(ck, nmod_pow_ui(c, n - k, mod), mod);

                        if (c >= 2)
                        {
                            /* See if we have already checked f(c*x). */
                            e2n = invctab[cn];
                            b2 = nmod_mul(b, e2n, mod);
                            a2 = nmod_mul(nmod_mul(a, ck, mod), e2n, mod);
                            if (ALREADY_CHECKED(a, b, a2, b2))
                                goto next_trinomial;

                            /* See if we have already checked f(-c*x). */
                            b2 = (n % 2) ? nmod_neg(b2, mod) : b2;
                            a2 = ((n ^ k) % 2) ? nmod_neg(a2, mod) : a2;
                            if (ALREADY_CHECKED(a, b, a2, b2))
                                goto next_trinomial;
                        }

                        /* Trial division: check if c^n + a c^k + b = 0. */
                        if (nmod_add(cn, nmod_mul(a, ck, mod), mod) == p - b)
                            goto next_trinomial;

                        /* Trial division: check if (-c)^n + a (-c)^k + b = 0. */
                        cn2 = (n % 2) ? nmod_neg(cn, mod) : cn;
                        ck2 = (k % 2) ? nmod_neg(ck, mod) : ck;
                        if (nmod_add(cn2, nmod_mul(a, ck2, mod), mod) == p - b)
                            goto next_trinomial;
                    }
                }

                nmod_poly_zero(res);
                nmod_poly_set_coeff_ui(res, n, 1);
                nmod_poly_set_coeff_ui(res, 0, b);
                nmod_poly_set_coeff_ui(res, k, a);

                for (i = 0; i < num_sieve; i++)
                {
                    Bcoeffs[0] = b;
                    Bcoeffs[1] = a;
                    Bcoeffs[2] = 1;
                    Bexps[0] = 0;
                    Bexps[1] = k;
                    Bexps[2] = n;

                    nmod_poly_gcd_with_sparse(h, sieve[i], res, Bcoeffs, Bexps, 3);
                    if (!nmod_poly_is_one(h))
                        goto next_trinomial;
                }

                if (nmod_poly_is_irreducible_ddf(res))
                {
                    found = 1;
                    goto cleanup;
                }

                next_trinomial:
                    continue;
            }
        }
    }

cleanup:
    if (want_pruning)
    {
        _nmod_vec_clear(invctab);
        for (k = 0; k < num_sieve; k++)
            nmod_poly_clear(sieve[k]);
        nmod_poly_clear(h);
    }

    return found;
}

int
nmod_poly_irreducible_tetranomial(nmod_poly_t res, ulong n)
{
    ulong p = res->mod.n;
    nmod_t mod = res->mod;
    ulong a, b, c, d, k, l;
    ulong dk, dl, dn, dk2, dl2, dn2, s;
    slong i;
    int want_pruning;
    int found = 0;
    nmod_poly_t sieve[10], h;
    slong num_sieve = 0;
    ulong Bcoeffs[4];
    slong Bexps[4];
    ulong sieve_limit;

    if (n <= 2 || p == 2)
        return 0;

    sieve_limit = trinomial_sieve_limit(p, n);
    want_pruning = (sieve_limit > 0);

    if (want_pruning)
    {
        nmod_poly_init_mod(h, mod);

        for (a = 2; a <= sieve_limit && a < n; a++)
        {
            nmod_poly_init_mod(sieve[num_sieve], mod);
            nmod_poly_product_all_irreducibles_deg(sieve[num_sieve], a);
            num_sieve++;
        }
    }

    /* Iterate over all x^n + a x^k + b x^l + c */
    for (a = 1; a < p; a++)
    {
        for (b = 1; b < p; b++)
        {
            for (c = 1; c < p; c++)
            {
                for (k = 2; k < n; k++)
                {
                    for (l = 1; l < k; l++)
                    {
                        if (want_pruning)
                        {
                            for (d = 1; 2 * d < p; d++)
                            {
                                dk = nmod_pow_ui(d, k, mod);
                                dl = nmod_pow_ui(d, l, mod);
                                dn = nmod_pow_ui(d, n, mod);

                                /* Trial division: check if d^n + a d^k + b d^l + c = 0. */
                                s = nmod_add(dn, nmod_mul(a, dk, mod), mod);
                                s = nmod_add(s, nmod_mul(b, dl, mod), mod);
                                if (s == p - c)
                                    goto next_tetranomial;

                                dn2 = (n % 2) ? nmod_neg(dn, mod) : dn;
                                dk2 = (k % 2) ? nmod_neg(dk, mod) : dk;
                                dl2 = (l % 2) ? nmod_neg(dl, mod) : dl;
                                /* Trial division: check if (-d)^n + a (-d)^k + b (-d)^l + c = 0. */
                                s = nmod_add(dn2, nmod_mul(a, dk2, mod), mod);
                                s = nmod_add(s, nmod_mul(b, dl2, mod), mod);
                                if (s == p - c)
                                    goto next_tetranomial;
                            }
                        }

                        nmod_poly_zero(res);
                        nmod_poly_set_coeff_ui(res, n, 1);
                        nmod_poly_set_coeff_ui(res, k, a);
                        nmod_poly_set_coeff_ui(res, l, b);
                        nmod_poly_set_coeff_ui(res, 0, c);

                        for (i = 0; i < num_sieve; i++)
                        {
                            Bcoeffs[0] = c;
                            Bcoeffs[1] = b;
                            Bcoeffs[2] = a;
                            Bcoeffs[3] = 1;
                            Bexps[0] = 0;
                            Bexps[1] = l;
                            Bexps[2] = k;
                            Bexps[3] = n;

                            nmod_poly_gcd_with_sparse(h, sieve[i], res, Bcoeffs, Bexps, 4);
                            if (!nmod_poly_is_one(h))
                                goto next_tetranomial;
                        }

                        {
                            ulong D, D1, D2;
                            D = nmod_poly_discriminant(res);
                            D1 = nmod_pow_ui(D, (p - 1) / 2, mod);
                            D2 = (n % 2 == 1) ? p - 1 : 1; 
                            if (D1 == D2)
                                goto next_tetranomial;
                        }

                        if (nmod_poly_is_irreducible_ddf(res))
                        {
                            found = 1;
                            goto cleanup;
                        }

                        next_tetranomial:
                            continue;
                    }
                }
            }
        }
    }

cleanup:
    if (want_pruning)
    {
        for (k = 0; k < num_sieve; k++)
            nmod_poly_clear(sieve[k]);
        nmod_poly_clear(h);
    }

    return found;
}

/* Only intended for GF(2), so we only bother with coefficients 1 */
int
nmod_poly_irreducible_pentanomial(nmod_poly_t res, ulong n)
{
    ulong p = res->mod.n;
    nmod_t mod = res->mod;
    ulong k, l, m, a;
    slong i;
    int want_pruning;
    int found = 0;
    nmod_poly_t sieve[10], h;
    slong num_sieve = 0;
    ulong Bcoeffs[5];
    slong Bexps[5];
    ulong sieve_limit;

    if (n <= 3)
        return 0;

    sieve_limit = trinomial_sieve_limit(p, n);
    want_pruning = (sieve_limit > 0);

    if (want_pruning)
    {
        nmod_poly_init_mod(h, mod);

        for (a = 2; a <= sieve_limit && a < n; a++)
        {
            nmod_poly_init_mod(sieve[num_sieve], mod);
            nmod_poly_product_all_irreducibles_deg(sieve[num_sieve], a);
            num_sieve++;
        }
    }

    Bcoeffs[0] = 1;
    Bcoeffs[1] = 1;
    Bcoeffs[2] = 1;
    Bcoeffs[3] = 1;
    Bcoeffs[4] = 1;

    /* Iterate over all x^n + x^k + x^l + x^m + 1 */
    for (k = 3; k < n; k++)
    {
        for (l = 2; l < k; l++)
        {
            for (m = 1; m < l; m++)
            {
                nmod_poly_zero(res);
                nmod_poly_set_coeff_ui(res, n, 1);
                nmod_poly_set_coeff_ui(res, k, 1);
                nmod_poly_set_coeff_ui(res, l, 1);
                nmod_poly_set_coeff_ui(res, m, 1);
                nmod_poly_set_coeff_ui(res, 0, 1);

                for (i = 0; i < num_sieve; i++)
                {
                    Bexps[0] = 0;
                    Bexps[1] = m;
                    Bexps[2] = l;
                    Bexps[3] = k;
                    Bexps[4] = n;

                    nmod_poly_gcd_with_sparse(h, sieve[i], res, Bcoeffs, Bexps, 5);
                    if (!nmod_poly_is_one(h))
                        goto next_pentanomial;
                }

                if (nmod_poly_is_irreducible_ddf(res))
                {
                    found = 1;
                    goto cleanup;
                }

                next_pentanomial:
                    continue;
            }
        }
    }

cleanup:
    if (want_pruning)
    {
        for (k = 0; k < num_sieve; k++)
            nmod_poly_clear(sieve[k]);
        nmod_poly_clear(h);
    }

    return found;
}

void
nmod_poly_minimal_irreducible(nmod_poly_t res, ulong n)
{
    FLINT_ASSERT(n != 0);

    if (n == 1)
    {
        nmod_poly_zero(res);
        nmod_poly_set_coeff_ui(res, 1, 1);
        return;
    }

    if (nmod_poly_irreducible_binomial(res, n))
        return;

    if (nmod_poly_irreducible_trinomial(res, n))
        return;

    if (nmod_poly_irreducible_tetranomial(res, n))
        return;

    if (nmod_poly_irreducible_pentanomial(res, n))
        return;

    /* The conjecture is that a tetranomial (p != 2) or pentanomial (p = 2)
       always works. */
    flint_throw(FLINT_ERROR, "Failed to construct minimal irreducible polynomial (n = %wu, p = %wu)", n, res->mod.n);
}

