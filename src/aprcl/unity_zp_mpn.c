/*
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "mpn_mod.h"
#include "aprcl.h"

/*
    mpn_mod based arithmetic in Z[zeta_{p^k}]/(n), used by the APRCL
    Jacobi sum primality test when the modulus n is in range for the
    mpn_mod module (currently 2 to 16 limbs). This is the bottleneck
    of aprcl_is_prime() and is performance-critical.

    Elements are stored as a dense coefficient vector of length
    d = deg Phi_{p^k} = (p - 1) p^(k-1), each coefficient being a
    fully reduced mpn_mod residue of nlimbs limbs. In contrast with
    the fmpz_mod based unity_zp code, elements are always reduced
    modulo the cyclotomic polynomial Phi_{p^k}.

    Multiplication and squaring first compute the unreduced integer
    polynomial product with delayed reduction (each product coefficient
    is an exact integer of s limbs, with s up to 2 nlimbs + 1, computed
    by _flint_mpn_poly_mulmid with the classical algorithm or Karatsuba), then
    perform the cyclotomic reduction on the unreduced accumulators, and
    finally do a single reduction mod n per output coefficient. The cyclotomic
    reduction uses

        c_i' = c_i + c_{i + p^k}                     (x^{p^k} = 1)
        r_i  = c_i' - c_{d + (i mod p^(k-1))}        (reduction by Phi)

    since x^e for e in [d, p^k) satisfies
    x^e = -(x^{e-d} + x^{e-d+p^(k-1)} + ... + x^{e-d+(p-2)p^(k-1)}),
    so each output coefficient i receives a contribution from exactly
    one exponent e = d + (i mod p^(k-1)) in [d, p^k).

    To keep all mpn operations unsigned we add the padding constant
    n * B^(nlimbs + 1) (B = 2^FLINT_BITS) before subtracting; this is
    a multiple of n exceeding any accumulator value, and the padded
    value still fits in 2 nlimbs + 2 limbs.
*/

#define UNITY_ZP_MPN_NLIMBS(f) MPN_MOD_CTX_NLIMBS((f)->ctx)
#define COEFF(f, i, n) ((f)->coeffs + (i) * (n))

/* Memory management *********************************************************/

void
unity_zp_mpn_init(unity_zp_mpn f, ulong p, ulong exp, gr_ctx_t ctx)
{
    f->p = p;
    f->exp = exp;
    f->pk1 = n_pow(p, exp - 1);
    f->pk = f->pk1 * p;
    f->d = (p - 1) * f->pk1;
    f->ctx = (gr_ctx_struct *) ctx;
    f->coeffs = (nn_ptr) flint_calloc(f->d * MPN_MOD_CTX_NLIMBS(ctx),
                                      sizeof(ulong));

    /* the products are always d x d with the residues as unsigned
       coefficients: fix the algorithm and the output size once (the
       squaring cutoff is not below the multiplication one, so the
       multiplication's choice and size are valid for both) */
    {
        slong n = MPN_MOD_CTX_NLIMBS(ctx);
        int norm = MPN_MOD_CTX_NORM(ctx);
        slong cutoff = _flint_mpn_poly_mulmid_cutoff(n, n, 0);
        f->use_karatsuba = _flint_mpn_poly_mulmid_use_karatsuba(f->d, norm, FLINT_MPN_POLY_MUL_UNSIGNED, cutoff);
        f->slimbs = _flint_mpn_poly_mulmid_slimbs(f->d, n, norm, n, norm, FLINT_MPN_POLY_MUL_UNSIGNED, cutoff);
    }
}

void
unity_zp_mpn_clear(unity_zp_mpn f)
{
    flint_free(f->coeffs);
}

void
unity_zp_mpn_set_zero(unity_zp_mpn f)
{
    flint_mpn_zero(f->coeffs, f->d * UNITY_ZP_MPN_NLIMBS(f));
}

void
unity_zp_mpn_swap(unity_zp_mpn f, unity_zp_mpn g)
{
    FLINT_SWAP(nn_ptr, f->coeffs, g->coeffs);
}

void
unity_zp_mpn_copy(unity_zp_mpn f, const unity_zp_mpn g)
{
    flint_mpn_copyi(f->coeffs, g->coeffs, f->d * UNITY_ZP_MPN_NLIMBS(f));
}

void
unity_zp_mpn_coeff_set_ui(unity_zp_mpn f, ulong ind, ulong x)
{
    FLINT_ASSERT(ind < f->d);
    GR_MUST_SUCCEED(mpn_mod_set_ui(COEFF(f, ind, UNITY_ZP_MPN_NLIMBS(f)),
                                   x, f->ctx));
}

/*
    Sets f to the value of the fmpz_mod based element g, which need not be
    reduced modulo Phi_{p^k} (any degree < p^k is accepted).
*/
void
unity_zp_mpn_set_unity_zp(unity_zp_mpn f, const unity_zp g)
{
    slong i, j, len;
    slong n = UNITY_ZP_MPN_NLIMBS(f);
    ulong t[MPN_MOD_MAX_LIMBS];

    len = g->poly->length;
    FLINT_ASSERT((ulong) len <= f->pk);

    unity_zp_mpn_set_zero(f);

    for (i = 0; i < len; i++)
    {
        if (fmpz_is_zero(g->poly->coeffs + i))
            continue;

        GR_MUST_SUCCEED(mpn_mod_set_fmpz(t, g->poly->coeffs + i, f->ctx));

        if ((ulong) i < f->d)
        {
            GR_MUST_SUCCEED(mpn_mod_add(COEFF(f, i, n), COEFF(f, i, n),
                                        t, f->ctx));
        }
        else
        {
            /* x^i = -(x^{i-d} + x^{i-d+p^(k-1)} + ...) */
            for (j = 0; j < (slong) f->p - 1; j++)
            {
                ulong ind = i - f->d + j * f->pk1;
                GR_MUST_SUCCEED(mpn_mod_sub(COEFF(f, ind, n),
                                            COEFF(f, ind, n), t, f->ctx));
            }
        }
    }
}

/* Multiplication and squaring ***********************************************/

/*
    Given the unreduced product coefficients c[0], ..., c[2 d - 2] computed
    by _flint_mpn_poly_mulmid, each an s-limb nonnegative integer < B^s with
    B = 2^FLINT_BITS, performs the cyclotomic
    reduction

        f_i = c_i + c_{i + p^k} - c_{d + (i mod p^(k-1))}

    with a single mod n reduction per output coefficient. To keep the
    intermediate value nonnegative we add the padding constant
    pad = n 2^t with t = FLINT_BITS s - bits(n) + 1, a multiple of n
    satisfying B^s <= pad < 2 B^s, so that the total fits s + 1 limbs.
*/
static void
_unity_zp_mpn_cyclotomic_reduce(unity_zp_mpn f, nn_srcptr c, slong s)
{
    slong i;
    slong n = UNITY_ZP_MPN_NLIMBS(f);
    slong d = f->d;
    slong lastc = 2 * d - 2;
    flint_bitcnt_t t = FLINT_BITS * s - MPN_MOD_CTX_MODULUS_BITS(f->ctx) + 1;
    slong q = t / FLINT_BITS, r = t % FLINT_BITS;
    slong pn;
    ulong pm[MPN_MOD_MAX_LIMBS + 1];
    ulong u[2 * MPN_MOD_MAX_LIMBS + 2];

    /* pad = pm B^q with pm = n 2^r, so that B^s <= pad < 2 B^s */
    if (r != 0)
    {
        pm[n] = mpn_lshift(pm, MPN_MOD_CTX_MODULUS(f->ctx), n, r);
        pn = n + 1;
    }
    else
    {
        flint_mpn_copyi(pm, MPN_MOD_CTX_MODULUS(f->ctx), n);
        pn = n;
    }

    FLINT_ASSERT(q + pn <= s + 1);

    for (i = 0; i < d; i++)
    {
        slong e = d + ((ulong) i % f->pk1);
        slong un;
        ulong cy;

        /* u = c_i */
        flint_mpn_copyi(u, c + i * s, s);
        u[s] = 0;

        /* u += c_{i + p^k} */
        if (i + (slong) f->pk <= lastc)
        {
            cy = mpn_add(u, u, s + 1, c + (i + f->pk) * s, s);
            FLINT_ASSERT(cy == 0);
        }

        /* u += pad, u -= c_{d + (i mod p^(k-1))} */
        if (e <= lastc)
        {
            cy = mpn_add(u + q, u + q, s + 1 - q, pm, pn);
            FLINT_ASSERT(cy == 0);
            cy = mpn_sub(u, u, s + 1, c + e * s, s);
            FLINT_ASSERT(cy == 0);
        }

        (void) cy;  /* used for asserts only */

        un = s + 1;
        MPN_NORM(u, un);
        GR_MUST_SUCCEED(mpn_mod_set_mpn(COEFF(f, i, n), u, un, f->ctx));
    }
}

/* f = g h (or g^2 when h == g), computing the exact integer product with
   delayed reduction via _flint_mpn_poly_mulmid, which chooses between the
   classical algorithm and Karatsuba. */
static void
_unity_zp_mpn_mul(unity_zp_mpn f, const unity_zp_mpn g, const unity_zp_mpn h)
{
    slong d = f->d;
    slong n = UNITY_ZP_MPN_NLIMBS(f);
    int norm = MPN_MOD_CTX_NORM(f->ctx);
    slong slimbs = f->slimbs;
    nn_ptr c;
    TMP_INIT;

    TMP_START;
    c = TMP_ALLOC((2 * d - 1) * slimbs * sizeof(ulong));

    /* the products are typically tiny (d is usually below 16), so the
       fixed costs of the general entry point are avoided in the classical
       range; the algorithm and output size were fixed at initialization */
    if (!f->use_karatsuba)
        _flint_mpn_poly_mulmid_classical(c, g->coeffs, d, n, h->coeffs, d, n, 0, 2 * d - 1, slimbs, FLINT_MPN_POLY_MUL_UNSIGNED);
    else
        _flint_mpn_poly_mulmid(c, g->coeffs, d, n, norm, h->coeffs, d, n, norm, 0, 2 * d - 1, slimbs, FLINT_MPN_POLY_MUL_UNSIGNED);

    _unity_zp_mpn_cyclotomic_reduce(f, c, slimbs);

    TMP_END;
}

void
unity_zp_mpn_mul(unity_zp_mpn f, const unity_zp_mpn g, const unity_zp_mpn h)
{
    _unity_zp_mpn_mul(f, g, h);
}

void
unity_zp_mpn_sqr(unity_zp_mpn f, const unity_zp_mpn g)
{
    _unity_zp_mpn_mul(f, g, g);
}

void
unity_zp_mpn_mul_scalar_ui(unity_zp_mpn f, const unity_zp_mpn g, ulong s)
{
    slong i;
    slong n = UNITY_ZP_MPN_NLIMBS(f);

    for (i = 0; i < (slong) f->d; i++)
        GR_MUST_SUCCEED(mpn_mod_mul_ui(COEFF(f, i, n), COEFF(g, i, n),
                                       s, f->ctx));
}

/* Powering ******************************************************************/

void
unity_zp_mpn_pow_ui(unity_zp_mpn f, const unity_zp_mpn g, ulong pow)
{
    slong i;
    unity_zp_mpn base, temp;

    if (pow == 0)
    {
        unity_zp_mpn_set_zero(f);
        unity_zp_mpn_coeff_set_ui(f, 0, 1);
        return;
    }

    unity_zp_mpn_init(base, f->p, f->exp, f->ctx);
    unity_zp_mpn_init(temp, f->p, f->exp, f->ctx);

    unity_zp_mpn_copy(base, g);
    unity_zp_mpn_copy(f, g);

    for (i = FLINT_BIT_COUNT(pow) - 2; i >= 0; i--)
    {
        unity_zp_mpn_sqr(temp, f);
        unity_zp_mpn_swap(temp, f);

        if ((pow >> i) & 1)
        {
            unity_zp_mpn_mul(temp, f, base);
            unity_zp_mpn_swap(temp, f);
        }
    }

    unity_zp_mpn_clear(base);
    unity_zp_mpn_clear(temp);
}

/*
    Sliding window powering, f = g^pow; same algorithm as
    unity_zp_pow_sliding_fmpz.
*/
void
unity_zp_mpn_pow_sliding_fmpz(unity_zp_mpn f, const unity_zp_mpn g,
                              const fmpz_t pow)
{
    ulong h, k, value;
    slong i, j, m;
    unity_zp_mpn temp;
    unity_zp_mpn * g_powers;

    unity_zp_mpn_init(temp, f->p, f->exp, f->ctx);

    /* selects optimal k value for pow */
    k = _unity_zp_pow_select_k(pow);
    m = n_pow(2, k - 1);

    /*
        g_powers store odd powers of g up to 2^k - 1;
        g_powers[(i + 1) / 2] = g^i
    */
    g_powers = (unity_zp_mpn *) flint_malloc(sizeof(unity_zp_mpn) * (m + 1));

    unity_zp_mpn_init(g_powers[0], f->p, f->exp, f->ctx);
    unity_zp_mpn_coeff_set_ui(g_powers[0], 0, 1);

    unity_zp_mpn_init(g_powers[1], f->p, f->exp, f->ctx);
    unity_zp_mpn_copy(g_powers[1], g);

    /* temp = g^2 */
    unity_zp_mpn_sqr(temp, g);

    for (i = 2; i <= m; i++)
    {
        unity_zp_mpn_init(g_powers[i], f->p, f->exp, f->ctx);
        unity_zp_mpn_mul(g_powers[i], g_powers[i - 1], temp);
    }

    unity_zp_mpn_set_zero(f);
    unity_zp_mpn_coeff_set_ui(f, 0, 1);

    i = fmpz_bits(pow) - 1;
    while (i >= 0)
    {
        if (fmpz_tstbit(pow, i) == 0)
        {
            unity_zp_mpn_sqr(temp, f);
            unity_zp_mpn_swap(temp, f);
            i--;
        }
        else
        {
            /*
                find length of chain; chain is the longest bitstring
                of length at most k ending in 1
            */
            j = FLINT_MAX(i - (slong) k + 1, 0);
            while (fmpz_tstbit(pow, j) == 0 && j <= i)
                j++;

            /* f = f^(2^(i - j + 1)) */
            for (h = 0; h < (ulong) (i - j + 1); h++)
            {
                unity_zp_mpn_sqr(temp, f);
                unity_zp_mpn_swap(temp, f);
            }

            /* value = binary number (e_i, ... , e_j) in decimal base */
            value = 0;
            for (h = 0; h < (ulong) (i - j + 1); h++)
                value += fmpz_tstbit(pow, j + h) << h;

            /* f = f * g^value */
            unity_zp_mpn_mul(temp, f, g_powers[(value + 1) / 2]);
            unity_zp_mpn_swap(temp, f);

            i = j - 1;
        }
    }

    for (i = 0; i <= m; i++)
        unity_zp_mpn_clear(g_powers[i]);
    flint_free(g_powers);

    unity_zp_mpn_clear(temp);
}

/* Automorphisms *************************************************************/

/*
    Computes f such that \sigma_x(f) = g; f and g must not alias.
    Same algorithm as unity_zp_aut_inv.
*/
void
unity_zp_mpn_aut_inv(unity_zp_mpn f, const unity_zp_mpn g, ulong x)
{
    ulong i, j, g_ind, p_pow_preinv;
    slong n = UNITY_ZP_MPN_NLIMBS(f);

    FLINT_ASSERT(f->coeffs != g->coeffs);

    p_pow_preinv = n_preinvert_limb(f->pk);

    unity_zp_mpn_set_zero(f);

    /* for i = 0, 1, ..., d - 1 set f[i] = g[x * i mod p^k] */
    for (i = 0; i < f->d; i++)
    {
        g_ind = n_mulmod2_preinv(x, i, f->pk, p_pow_preinv);

        if (g_ind < f->d)
            GR_MUST_SUCCEED(mpn_mod_set(COEFF(f, i, n), COEFF(g, g_ind, n),
                                        f->ctx));
    }

    /*
        for i = d, ..., p^k - 1 and j = 1, ..., p - 1 set
        f[i - j * p^(k-1)] -= g[x * i mod p^k]
    */
    for (i = f->d; i < f->pk; i++)
    {
        g_ind = n_mulmod2_preinv(x, i, f->pk, p_pow_preinv);

        if (g_ind >= f->d)
            continue;

        for (j = 1; j < f->p; j++)
        {
            ulong f_ind = i - j * f->pk1;
            GR_MUST_SUCCEED(mpn_mod_sub(COEFF(f, f_ind, n),
                                        COEFF(f, f_ind, n),
                                        COEFF(g, g_ind, n), f->ctx));
        }
    }
}

/* Comparison ****************************************************************/

/*
    If f = \zeta_{p^k}^h for some h, returns h; otherwise returns -1.
    Uses the fact that for h < d, \zeta^h is a basis vector, while for
    d <= h < p^k the reduced representation of \zeta^h has the p - 1
    coefficients (h - d) + j p^(k-1), 0 <= j < p - 1, all equal to -1
    and all other coefficients zero.
*/
slong
unity_zp_mpn_is_unity(const unity_zp_mpn f)
{
    slong i;
    slong n = UNITY_ZP_MPN_NLIMBS(f);
    slong ones = 0, negones = 0;
    slong one_ind = -1, negone_ind = -1;
    ulong negone[MPN_MOD_MAX_LIMBS];

    GR_MUST_SUCCEED(mpn_mod_neg_one(negone, f->ctx));

    for (i = 0; i < (slong) f->d; i++)
    {
        nn_srcptr ci = COEFF(f, i, n);

        if (flint_mpn_zero_p(ci, n))
            continue;

        if (ci[0] == 1 && flint_mpn_zero_p(ci + 1, n - 1))
        {
            ones++;
            one_ind = i;
        }
        else if (flint_mpn_equal_p(ci, negone, n))
        {
            negones++;
            if (negone_ind == -1)
                negone_ind = i;
        }
        else
        {
            return -1;
        }
    }

    if (ones == 1 && negones == 0)
        return one_ind;

    if (ones == 0 && negones == (slong) f->p - 1)
    {
        ulong e0 = (ulong) negone_ind % f->pk1;
        ulong j;

        if ((ulong) negone_ind != e0)
            return -1;

        for (j = 0; j < f->p - 1; j++)
        {
            nn_srcptr ci = COEFF(f, e0 + j * f->pk1, n);
            if (!flint_mpn_equal_p(ci, negone, n))
                return -1;
        }

        return f->d + e0;
    }

    return -1;
}

/* Jacobi test check functions ***********************************************/

/*
    The following are mpn_mod based versions of
    _aprcl_is_prime_jacobi_check_pk / _check_22 / _check_2k. Each takes the
    fmpz_mod based Jacobi sums as input, converts them, and runs the whole
    check using mpn_mod arithmetic. They return 1 and set *h if the modulus
    is in range for mpn_mod, and return 0 (leaving *h unchanged) if the
    caller must fall back to the generic code.
*/

static int
_aprcl_jacobi_mpn_start(gr_ctx_t ctx, const unity_zp j)
{
    const fmpz * nmod = fmpz_mod_ctx_modulus(j->ctx);

    if (fmpz_sgn(nmod) <= 0)
        return 0;

    if (gr_ctx_init_mpn_mod(ctx, nmod) != GR_SUCCESS)
        return 0;

    return 1;
}

int
_aprcl_is_prime_jacobi_check_pk_mpn(slong * h, const unity_zp j,
                                    const fmpz_t u, ulong v)
{
    ulong i, r;
    gr_ctx_t ctx;
    unity_zp_mpn J, j0, jv, temp, aut;

    if (!_aprcl_jacobi_mpn_start(ctx, j))
        return 0;

    r = n_pow(j->p, j->exp);

    unity_zp_mpn_init(J, j->p, j->exp, ctx);
    unity_zp_mpn_init(j0, j->p, j->exp, ctx);
    unity_zp_mpn_init(jv, j->p, j->exp, ctx);
    unity_zp_mpn_init(temp, j->p, j->exp, ctx);
    unity_zp_mpn_init(aut, j->p, j->exp, ctx);

    unity_zp_mpn_set_unity_zp(J, j);

    unity_zp_mpn_coeff_set_ui(j0, 0, 1);
    unity_zp_mpn_coeff_set_ui(jv, 0, 1);

    for (i = 1; i <= r; i++)
    {
        if (i % j->p == 0)
            continue;

        /* j0 *= \sigma_i^{-1}(J^i) */
        unity_zp_mpn_pow_ui(temp, J, i);
        unity_zp_mpn_aut_inv(aut, temp, i);
        unity_zp_mpn_mul(temp, j0, aut);
        unity_zp_mpn_swap(temp, j0);

        /* jv *= \sigma_i^{-1}(J^{(v * i) / r}) */
        unity_zp_mpn_pow_ui(temp, J, (v * i) / r);
        unity_zp_mpn_aut_inv(aut, temp, i);
        unity_zp_mpn_mul(temp, jv, aut);
        unity_zp_mpn_swap(temp, jv);
    }

    /* j0 = j0^u * jv */
    unity_zp_mpn_pow_sliding_fmpz(temp, j0, u);
    unity_zp_mpn_mul(j0, jv, temp);

    *h = unity_zp_mpn_is_unity(j0);

    unity_zp_mpn_clear(J);
    unity_zp_mpn_clear(j0);
    unity_zp_mpn_clear(jv);
    unity_zp_mpn_clear(temp);
    unity_zp_mpn_clear(aut);

    gr_ctx_clear(ctx);

    return 1;
}

int
_aprcl_is_prime_jacobi_check_22_mpn(slong * h, const unity_zp j,
                                    const fmpz_t u, ulong v, ulong q)
{
    gr_ctx_t ctx;
    unity_zp_mpn J, j0, jv, j_pow;

    if (!_aprcl_jacobi_mpn_start(ctx, j))
        return 0;

    unity_zp_mpn_init(J, 2, 2, ctx);
    unity_zp_mpn_init(j0, 2, 2, ctx);
    unity_zp_mpn_init(jv, 2, 2, ctx);
    unity_zp_mpn_init(j_pow, 2, 2, ctx);

    unity_zp_mpn_set_unity_zp(J, j);

    /* j0 = q * J^2 */
    unity_zp_mpn_sqr(j_pow, J);
    unity_zp_mpn_mul_scalar_ui(j0, j_pow, q);

    /* jv = 1 if v == 1; jv = J^2 if v == 3 */
    if (v == 1)
        unity_zp_mpn_coeff_set_ui(jv, 0, 1);
    else if (v == 3)
        unity_zp_mpn_swap(jv, j_pow);

    /* j0 = j0^u * jv */
    unity_zp_mpn_pow_sliding_fmpz(j_pow, j0, u);
    unity_zp_mpn_mul(j0, jv, j_pow);

    *h = unity_zp_mpn_is_unity(j0);

    unity_zp_mpn_clear(J);
    unity_zp_mpn_clear(j0);
    unity_zp_mpn_clear(jv);
    unity_zp_mpn_clear(j_pow);

    gr_ctx_clear(ctx);

    return 1;
}

int
_aprcl_is_prime_jacobi_check_2k_mpn(slong * h, const unity_zp j,
                                    const unity_zp j2_1, const unity_zp j2_2,
                                    const fmpz_t u, ulong v)
{
    ulong i, r;
    gr_ctx_t ctx;
    unity_zp_mpn J, J2_1, J2_2, j_j0, j0, jv, temp, aut;

    if (!_aprcl_jacobi_mpn_start(ctx, j))
        return 0;

    r = n_pow(j->p, j->exp);

    unity_zp_mpn_init(J, 2, j->exp, ctx);
    unity_zp_mpn_init(J2_1, 2, j->exp, ctx);
    unity_zp_mpn_init(J2_2, 2, j->exp, ctx);
    unity_zp_mpn_init(j_j0, 2, j->exp, ctx);
    unity_zp_mpn_init(j0, 2, j->exp, ctx);
    unity_zp_mpn_init(jv, 2, j->exp, ctx);
    unity_zp_mpn_init(temp, 2, j->exp, ctx);
    unity_zp_mpn_init(aut, 2, j->exp, ctx);

    unity_zp_mpn_set_unity_zp(J, j);
    unity_zp_mpn_set_unity_zp(J2_1, j2_1);
    unity_zp_mpn_set_unity_zp(J2_2, j2_2);

    unity_zp_mpn_coeff_set_ui(j0, 0, 1);
    unity_zp_mpn_coeff_set_ui(jv, 0, 1);

    /* j_j0 = J(2, q) * J_3(q) */
    unity_zp_mpn_mul(j_j0, J, J2_1);

    /* for i in 1..r and (i == 1 or i == 3 mod 8) */
    for (i = 1; i < r;)
    {
        ulong step;

        for (step = 0; step < 2; step++)
        {
            /* j0 *= \sigma_i^{-1}(j_j0^i) */
            unity_zp_mpn_pow_ui(temp, j_j0, i);
            unity_zp_mpn_aut_inv(aut, temp, i);
            unity_zp_mpn_mul(temp, j0, aut);
            unity_zp_mpn_swap(temp, j0);

            /* jv *= \sigma_i^{-1}(j_j0^{(v * i) / r}) */
            unity_zp_mpn_pow_ui(temp, j_j0, (v * i) / r);
            unity_zp_mpn_aut_inv(aut, temp, i);
            unity_zp_mpn_mul(temp, jv, aut);
            unity_zp_mpn_swap(temp, jv);

            /* from i == 1 mod 8 to i == 3 mod 8, then wrap to 1 mod 8 */
            i += (step == 0) ? 2 : 6;
        }
    }

    /* if v not congruent to 1 or 3 modulo 8 then jv *= J_2(q)^2 */
    if (v % 8 != 1 && v % 8 != 3)
    {
        unity_zp_mpn_sqr(temp, J2_2);
        unity_zp_mpn_mul(j_j0, jv, temp);
        unity_zp_mpn_swap(j_j0, jv);
    }

    /* j0 = j0^u * jv */
    unity_zp_mpn_pow_sliding_fmpz(temp, j0, u);
    unity_zp_mpn_mul(j0, jv, temp);

    *h = unity_zp_mpn_is_unity(j0);

    unity_zp_mpn_clear(J);
    unity_zp_mpn_clear(J2_1);
    unity_zp_mpn_clear(J2_2);
    unity_zp_mpn_clear(j_j0);
    unity_zp_mpn_clear(j0);
    unity_zp_mpn_clear(jv);
    unity_zp_mpn_clear(temp);
    unity_zp_mpn_clear(aut);

    gr_ctx_clear(ctx);

    return 1;
}
