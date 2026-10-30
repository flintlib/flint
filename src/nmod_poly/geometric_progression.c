/*
    Copyright (C) 2025, Vincent Neiger, Éric Schost, Mael Hostettler
    Copyright (C) 2026, Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly/impl.h"
#if FLINT_HAVE_FFT_SMALL
#include "fft_small.h"
#endif

/* specialized for moduli that support n_mulmod_shoup */
static
void _nmod_geometric_progression_evaluate_init_nonfullword(nmod_geometric_progression_t G,
                                                           ulong r, slong len, nmod_t mod,
                                                           ulong q, ulong q_pr_quo, ulong q_pr_rem,
                                                           ulong inv_r,
                                                           ulong inv_q, ulong inv_q_pr_quo, ulong inv_q_pr_rem)
{
    /* G->ev_f = sum_{0 <= i < 2*len - 1} q**(i*i/2) * x**i */
    /* G->ev_s[i] = 1 / q**(i*i/2) = 1 / r**(i*i)           */
    nmod_poly_init2_preinv(G->ev_f, mod.n, mod.ninv, 2*len - 1);
    G->ev_s = _nmod_vec_init(len);

    G->ev_f->length = 2*len - 1;

    /* precomputations for Shoup multiplication */
    ulong r_pow_i2 = r;
    ulong r_pr_quo = _nmod_mulmod_precomp_shoup(r_pow_i2, mod);

    G->ev_f->coeffs[0] = 1;
    for (slong i = 1; i < 2*len - 1; i++)
    {
        G->ev_f->coeffs[i] = n_mulmod_shoup(r_pow_i2, G->ev_f->coeffs[i - 1], r_pr_quo, mod.n);
        n_mulmod_and_precomp_shoup(&r_pow_i2, &r_pr_quo, q, r_pow_i2, q_pr_quo, q_pr_rem, r_pr_quo, mod.n);
    }

    r_pr_quo = _nmod_mulmod_precomp_shoup(inv_r, mod);
    G->ev_s[0] = 1;
    for (slong i = 1; i < len; i++)
    {
        G->ev_s[i] = n_mulmod_shoup(inv_r, G->ev_s[i - 1], r_pr_quo, mod.n);
        n_mulmod_and_precomp_shoup(&inv_r, &r_pr_quo, inv_q, inv_r, inv_q_pr_quo, inv_q_pr_rem, r_pr_quo, mod.n);
    }
}

/* general variant */
static
void _nmod_geometric_progression_evaluate_init(nmod_geometric_progression_t G,
                                               ulong r, slong len, nmod_t mod,
                                               ulong q, ulong inv_r, ulong inv_q)
{
    /* G->ev_f = sum_{0 <= i < 2*len - 1} q**(i*i/2) * x**i */
    /* G->ev_s[i] = 1 / q**(i*i/2)                          */
    nmod_poly_init2_preinv(G->ev_f, mod.n, mod.ninv, 2*len - 1);
    G->ev_s = _nmod_vec_init(len);

    G->ev_f->length = 2*len - 1;
    G->ev_f->coeffs[0] = 1;
    for (slong i = 1; i < 2*len - 1; i++)
    {
        G->ev_f->coeffs[i] = nmod_mul(G->ev_f->coeffs[i - 1], r, mod);
        r = nmod_mul(r, q, mod);  /* r**(2*i+1) */
    }

    G->ev_s[0] = 1;
    for (slong i = 1; i < len; i++)
    {
        G->ev_s[i] = nmod_mul(G->ev_s[i - 1], inv_r, mod);
        inv_r = nmod_mul(inv_r, inv_q, mod);
    }
}

static
void _nmod_geometric_progression_interpolate_init_nonfullword(nmod_geometric_progression_t G,
                                                              slong len, nmod_t mod,
                                                              ulong q, ulong q_pr_quo, ulong q_pr_rem,
                                                              ulong inv_q, ulong inv_q_pr_quo, ulong inv_q_pr_rem)
{
    /* quantities for Newton interpolation/evaluation/change-of-basis */
    /* see [Bostan - Schost, J.Complexity 2005, Section 5.1]          */
    /* write u_i for prod_{1 <= k <= i} (q**k - 1),                   */
    /*   and q_i for q**(i*(i-1)/2)                                   */

    /* coeff(G->int_f, i) = (-1)**i * q_i / u_i */
    nmod_poly_init2_preinv(G->int_f, mod.n, mod.ninv, len);
    G->int_f->length = len;

    /* G->int_s1[i] = 1 / u_i */
    /* G->int_s2[i] = u_i / q_i */
    G->int_s1 = _nmod_vec_init(len);
    G->int_s2 = _nmod_vec_init(len);

    G->int_f->coeffs[0] = 1;
    G->int_s2[0] = 1;

    const ulong one_precomp = _nmod_mulmod_precomp_shoup(UWORD(1), mod);
    ulong q_pow_i = 1;
    ulong q_pow_i_pr = one_precomp;
    ulong inv_q_pow_i = 1;
    ulong inv_q_pow_i_pr = one_precomp;
    ulong inv_q_i = 1;
    ulong inv_q_i_pr = one_precomp;
    ulong prod_diff = 1;

    for (slong i = 1; i < len; i++)
    {
        ulong inv_q_pow_i_pr_rem = n_mulmod_precomp_shoup_rem_from_quo(inv_q_pow_i_pr, mod.n);
        n_mulmod_and_precomp_shoup(&inv_q_i, &inv_q_i_pr, inv_q_pow_i, inv_q_i, inv_q_pow_i_pr, inv_q_pow_i_pr_rem, inv_q_i_pr, mod.n);  /* 1 / q_i */
        n_mulmod_and_precomp_shoup(&inv_q_pow_i, &inv_q_pow_i_pr, inv_q, inv_q_pow_i, inv_q_pr_quo, inv_q_pr_rem, inv_q_pow_i_pr, mod.n);  /* 1 / q**i */
        G->int_f->coeffs[i] = n_mulmod_shoup(q_pow_i, G->int_f->coeffs[i-1], q_pow_i_pr, mod.n);  /* q_i */
        n_mulmod_and_precomp_shoup(&q_pow_i, &q_pow_i_pr, q, q_pow_i, q_pr_quo, q_pr_rem, q_pow_i_pr, mod.n);  /* q**i */
        G->int_s1[i-1] = q_pow_i - 1;                       /* temporarily, q**i - 1 */
        prod_diff = nmod_mul(q_pow_i - 1, prod_diff, mod);  /* u_i */
        G->int_s2[i] = n_mulmod_shoup(inv_q_i, prod_diff, inv_q_i_pr, mod.n);  /* u_i / q_i */
    }

    G->int_s1[len-1] = nmod_inv(prod_diff, mod);  /* 1 / u_{len-1} */
    for (slong i = len - 1; i > 0; i--)
    {
        ulong w_i = G->int_s1[i];   /* 1 / u_i */
        if (i % 2)                  /* i odd, -q_i / u_i */
            G->int_f->coeffs[i] = mod.n - nmod_mul(G->int_f->coeffs[i], w_i, mod);
        else                        /* i even,  q_i / u_i */
            G->int_f->coeffs[i] = nmod_mul(G->int_f->coeffs[i], w_i, mod);
        G->int_s1[i-1] = nmod_mul(G->int_s1[i-1], w_i, mod);  /* 1 / u_{i-1} */
    }
}

static
void _nmod_geometric_progression_interpolate_init(nmod_geometric_progression_t G,
                                                  slong len, nmod_t mod,
                                                  ulong q, ulong inv_q)
{
    /* quantities for Newton interpolation/evaluation/change-of-basis */
    /* see [Bostan - Schost, J.Complexity 2005, Section 5.1]          */
    /* write u_i for prod_{1 <= k <= i} (q**k - 1),                   */
    /*   and q_i for q**(i*(i-1)/2)                                   */

    /* coeff(G->int_f, i) = (-1)**i * q_i / u_i */
    nmod_poly_init2_preinv(G->int_f, mod.n, mod.ninv, len);
    G->int_f->length = len;

    /* G->int_s1[i] = 1 / u_i */
    /* G->int_s2[i] = u_i / q_i */
    G->int_s1 = _nmod_vec_init(len);
    G->int_s2 = _nmod_vec_init(len);

    G->int_f->coeffs[0] = 1;
    G->int_s2[0] = 1;

    ulong q_pow_i = 1;
    ulong inv_q_pow_i = 1;
    ulong inv_q_i = 1;
    ulong prod_diff = 1;

    for (slong i = 1; i < len; i++)
    {
        inv_q_i = nmod_mul(inv_q_i, inv_q_pow_i, mod);      /* 1 / q_i */
        inv_q_pow_i = nmod_mul(inv_q_pow_i, inv_q, mod);    /* 1 / q**i */
        G->int_f->coeffs[i] = nmod_mul(G->int_f->coeffs[i-1], q_pow_i, mod);  /* q_i */
        q_pow_i = nmod_mul(q_pow_i, q, mod);                /* q**i */
        G->int_s1[i-1] = q_pow_i - 1;                       /* temporarily, q**i - 1 */
        prod_diff = nmod_mul(q_pow_i - 1, prod_diff, mod);  /* u_i */
        G->int_s2[i] = nmod_mul(prod_diff, inv_q_i, mod);   /* u_i / q_i */
    }

    G->int_s1[len-1] = nmod_inv(prod_diff, mod);  /* 1 / u_{len-1} */
    for (slong i = len - 1; i > 0; i--)
    {
        ulong w_i = G->int_s1[i];   /* 1 / u_i */
        if (i % 2)                  /* i odd,  - q_i / u_i */
            G->int_f->coeffs[i] = mod.n - nmod_mul(G->int_f->coeffs[i], w_i, mod);
        else                        /* i even, q_i / u_i */
            G->int_f->coeffs[i] = nmod_mul(G->int_f->coeffs[i], w_i, mod);
        G->int_s1[i-1] = nmod_mul(G->int_s1[i-1], w_i, mod);  /* 1 / u_{i-1} */
    }
}

static
void _nmod_geometric_progression_extrapolate_init(nmod_geometric_progression_t G,
                                                  slong len, nmod_t mod,
                                                  ulong q, ulong inv_q)
{
    /* ext_ff: sum_{i=0}^{len-1} x**i / (q**(i+1) - 1)    [for forward]  */
    /* ext_fb: sum_{i=0}^{len-1} x**i / (q**(-i-1) - 1)   [for backward] */
    nmod_poly_init2_preinv(G->ext_ff, mod.n, mod.ninv, len-1);
    nmod_poly_init2_preinv(G->ext_fb, mod.n, mod.ninv, len-1);
    G->ext_ff->length = len-1;
    G->ext_fb->length = len-1;

    /* ext_s1f[i] = prod_{k=1}^{i} (q**k - 1)    [for forward]  */
    /* ext_s1b[i] = prod_{k=1}^{i} (q**(-k) - 1) [for backward] */
    /* ext_s2[i] = 1 / ext_s1f[i]                [for both]     */
    /* ext_s3[i] = 1 / ext_s1b[i]                [for both]     */
    G->ext_s1f = _nmod_vec_init(len);
    G->ext_s1b = _nmod_vec_init(len);
    G->ext_s2 = _nmod_vec_init(len);
    G->ext_s3 = _nmod_vec_init(len);

    G->ext_s1f[0] = 1;
    G->ext_s1b[0] = 1;

    ulong q_pow_i = q;
    ulong inv_q_pow_i = inv_q;

    for (slong i = 1; i < len; i++)
    {
        G->ext_ff->coeffs[i-1] = q_pow_i - 1;       /* temporarily, q**i - 1 */
        G->ext_fb->coeffs[i-1] = inv_q_pow_i - 1;   /* temporarily, q**(-i) - 1 */
        G->ext_s1f[i] = nmod_mul(G->ext_s1f[i-1], q_pow_i - 1, mod);
        G->ext_s1b[i] = nmod_mul(G->ext_s1b[i-1], inv_q_pow_i - 1, mod);
        q_pow_i = nmod_mul(q_pow_i, q, mod);
        inv_q_pow_i = nmod_mul(inv_q_pow_i, inv_q, mod);
    }

    G->ext_s2[len-1] = nmod_inv(G->ext_s1f[len-1], mod);
    G->ext_s3[len-1] = nmod_inv(G->ext_s1b[len-1], mod);

    for (slong i = len-1; i > 0; i--)
    {
        G->ext_s2[i-1] = nmod_mul(G->ext_s2[i], G->ext_ff->coeffs[i-1], mod);
        G->ext_s3[i-1] = nmod_mul(G->ext_s3[i], G->ext_fb->coeffs[i-1], mod);
        G->ext_ff->coeffs[i-1] = nmod_mul(G->ext_s2[i], G->ext_s1f[i-1], mod);
        G->ext_fb->coeffs[i-1] = nmod_mul(G->ext_s3[i], G->ext_s1b[i-1], mod);
    }
}

static
void _nmod_geometric_progression_interpolate_extrapolate_init(nmod_geometric_progression_t G,
                                                              slong len, nmod_t mod,
                                                              ulong q, ulong inv_q)
{
    /* ext_ff: sum_{i=0}^{len-1} x**i / (q**(i+1) - 1)    [for forward]  */
    /* ext_fb: sum_{i=0}^{len-1} x**i / (q**(-i-1) - 1)   [for backward] */
    nmod_poly_init2_preinv(G->ext_ff, mod.n, mod.ninv, len-1);
    nmod_poly_init2_preinv(G->ext_fb, mod.n, mod.ninv, len-1);
    G->ext_ff->length = len-1;
    G->ext_fb->length = len-1;

    /* ext_s1f[i] = prod_{k=1}^{i} (q**k - 1)    [for forward]  */
    /* ext_s1b[i] = prod_{k=1}^{i} (q**(-k) - 1) [for backward] */
    /* ext_s2[i] = 1 / ext_s1f[i]                [for both]     */
    /* ext_s3[i] = 1 / ext_s1b[i]                [for both]     */
    G->ext_s1f = _nmod_vec_init(len);
    G->ext_s1b = _nmod_vec_init(len);
    G->ext_s2 = _nmod_vec_init(len);
    G->ext_s3 = _nmod_vec_init(len);

    /* coeff(G->int_f, i)  = (-1)**i * q**(i*(i-1)/2) / (prod_{k=1}^{i} (q**k - 1)) */
    /*                    == (-1)**i * q**(i*(i-1)/2) * ext_s2[i] */
    nmod_poly_init2_preinv(G->int_f, mod.n, mod.ninv, len);
    G->int_f->length = len;

    /* G->int_s1[i] = ext_s2[i] */
    /* G->int_s2[i] = ext_s1f[i] / q**(i*(i-1)/2) */
    G->int_s1 = G->ext_s2;
    G->int_s2 = _nmod_vec_init(len);

    G->ext_s1f[0] = 1;
    G->ext_s1b[0] = 1;
    G->int_f->coeffs[0] = 1;
    G->int_s2[0] = 1;

    ulong q_pow_i = 1;
    ulong inv_q_pow_i = 1;

    for (slong i = 1; i < len; i++)
    {
        G->int_f->coeffs[i] = nmod_mul(G->int_f->coeffs[i-1],
                                       mod.n - q_pow_i, mod);       /* temporarily, (-1)**i * q**(i*(i-1)/2) */
        G->int_s2[i] = nmod_mul(G->int_s2[i-1], inv_q_pow_i, mod);  /* temporarily, 1 / q**(i*(i-1)/2) */
        G->int_s2[i-1] = nmod_mul(G->int_s2[i-1], G->ext_s1f[i-1], mod);  /* int_s2[i-1] is now ok */
        q_pow_i = nmod_mul(q_pow_i, q, mod);
        inv_q_pow_i = nmod_mul(inv_q_pow_i, inv_q, mod);
        G->ext_ff->coeffs[i-1] = q_pow_i - 1;       /* temporarily, q**i - 1 */
        G->ext_fb->coeffs[i-1] = inv_q_pow_i - 1;   /* temporarily, q**(-i) - 1 */
        G->ext_s1f[i] = nmod_mul(G->ext_s1f[i-1], q_pow_i - 1, mod);
        G->ext_s1b[i] = nmod_mul(G->ext_s1b[i-1], inv_q_pow_i - 1, mod);
    }

    G->int_s2[len-1] = nmod_mul(G->int_s2[len-1], G->ext_s1f[len-1], mod);  /* int_s2[len-1] is now ok */
    G->ext_s2[len-1] = nmod_inv(G->ext_s1f[len-1], mod);
    G->ext_s3[len-1] = nmod_inv(G->ext_s1b[len-1], mod);

    for (slong i = len-1; i > 0; i--)
    {
        G->ext_s2[i-1] = nmod_mul(G->ext_s2[i], G->ext_ff->coeffs[i-1], mod);
        G->int_f->coeffs[i] = nmod_mul(G->ext_s2[i], G->int_f->coeffs[i], mod);
        G->ext_s3[i-1] = nmod_mul(G->ext_s3[i], G->ext_fb->coeffs[i-1], mod);
        G->ext_ff->coeffs[i-1] = nmod_mul(G->ext_s2[i], G->ext_s1f[i-1], mod);
        G->ext_fb->coeffs[i-1] = nmod_mul(G->ext_s3[i], G->ext_s1b[i-1], mod);
    }
}

/* Stored transforms (bit 3 of function) ***********************************/

#if FLINT_HAVE_FFT_SMALL

/* transform of (b, bn) for cyclic convolutions of length a power of two
   >= n >= bn, by operands of length at most an_max; NULL if no plan
   exists */
static mul_precomp_struct *
_fft_init(nn_srcptr b, slong bn, slong n, slong an_max, nmod_t mod)
{
    mul_precomp_struct * M = flint_malloc(sizeof(mul_precomp_struct));

    if (!_mul_precomp_init(M, b, bn, bn, FLINT_MAX(4, n_clog2(n)),
                           FLINT_MIN(an_max, bn), bn, mod, get_default_mpn_ctx()))
    {
        flint_free(M);
        return NULL;
    }

    return M;
}

/* the precomputation of transposed products with (v, vn), multipliers
   with at most an_max nonzero coefficients (see _mul_precomp_init_tr);
   NULL if no plan exists */
static mul_precomp_struct *
_fft_init_tr(nn_srcptr v, slong vn, slong an_max, nmod_t mod)
{
    mul_precomp_struct * M = flint_malloc(sizeof(mul_precomp_struct));

    if (!_mul_precomp_init_tr(M, v, vn, FLINT_MIN(an_max, vn), vn, mod,
                              get_default_mpn_ctx()))
    {
        flint_free(M);
        return NULL;
    }

    return M;
}

static void
_fft_clear(mul_precomp_struct * F)
{
    if (F != NULL)
    {
        _mul_precomp_clear(F);
        flint_free(F);
    }
}

/* The operations below compute products with a fixed polynomial b, in the
   form of the coefficients [zl - s, zh - s) of a * (b + s) with b + s
   truncated to length bn >= an: through _nmod_poly_mullow if zl = s = 0, or
   _nmod_poly_mulmid if zl - s = an - 1. These are also the coefficients
   [zl, zh) of a * b. Whether M, the stored transform of b, should be used:
   if it contains b + s (s + bn <= M->bn), was made for operands of length
   an, and is cheaper than the transforms that the alternative would make
   (it saves one of three, but may be longer and need more primes). For a
   stored transposed transform, see _fft_use_tr. */

/*
    With a stored transposed transform of b (M->transposed: evaluation and
    extrapolation, where zl >= an - 1), the coefficients [zl, zh) of a * b
    are the transposed product of r, a reversed with zl - (an - 1) zeros in
    front, of length rn = zl + 1, and of b: two transforms of the length of
    M, if zh <= M->bn. The alternative, the middle product of (b + s, bn)
    and a, takes three transforms, of the length of the transposed product
    of a and the window or of the window of the full product, whichever is
    shorter, as _nmod_poly_mulmid_fft_small chooses, and converts both
    operands. With np * ztrunc * depth for the cost of a transform, the
    stored transform pays while its two transforms cost at most 3.5 times
    one of the alternative
    TODO only measured on a Xeon for 30- and 60-bit moduli and for a prime of
    the fft_small context; about 4 for a prime transformed directly, hence the
    ratio 3.4 = 17/5 below (this is not expected to change much on other machines,
    but some check might be in order).
*/
static int
_fft_use_tr(const mul_precomp_struct * M, slong zl, slong zh, slong an,
            slong s, slong bn, nmod_t mod)
{
    fft_small_plan_t Dw, Dt;
    fft_small_plan_struct * D;
    slong m = zh - zl;
    ulong e, np;

    if (zl >= zh || zl < an - 1 || (ulong) zh > M->bn
            || (ulong) an > M->P->bound_c)
        return 0;

    /* cheap necessary condition for the cost test below: the alternative
       has depth at most e and at most as many primes */
    e = FLINT_MAX(4, n_clog2(an + bn + 3));
    if (10 * M->P->ztrunc * M->P->depth > 17 * e * n_pow2(e))
        return 0;

    /* only where the alternative would use fft_small */
    if (!_nmod_poly_mulmid_want_fft_small(an, m, mod))
        return 0;

    _fft_small_plan_set_window_nmod(Dw, zl - s, zh - s, an + bn - 1,
        FLINT_MAX(_fft_small_len_trunc(an), _fft_small_len_trunc(bn)));
    _fft_small_plan_set_window_nmod(Dt, 0, m, an + m - 1,
        FLINT_MAX(_fft_small_len_trunc(an), _fft_small_len_trunc(m)));
    D = (m < FLINT_MAX(an, bn) && Dt->ztrunc < Dw->ztrunc) ? Dt : Dw;

    D->R = get_default_mpn_ctx();
    if (M->P->np == 1 || !_fft_small_plan_set_bound(D, FLINT_MIN(an, bn),
                                         2 * NMOD_BITS(mod), MPN_CTX_NCRTS))
        np = M->P->np;
    else
        np = D->np;

    return 10 * M->P->np * M->P->ztrunc * M->P->depth
                <= 17 * np * D->ztrunc * D->depth;
}

static int
_fft_use(const mul_precomp_struct * M, slong zl, slong zh, slong an,
         slong s, slong bn, nmod_t mod)
{
    fft_small_plan_t P, D;
    ulong e;

    if (M->transposed)
        return _fft_use_tr(M, zl, zh, an, s, bn, mod);

    if (zl >= zh || (ulong) (s + bn) > M->bn || (ulong) an > M->P->bound_c)
        return 0;

    /* cheap necessary condition for the cost test below: with M, a
       transform has length at least M->bn and the depth of M; the
       alternative has depth at most e and at most as many primes */
    e = FLINT_MAX(4, n_clog2(an + bn + 3));
    if (2 * M->bn * M->P->depth > 3 * e * n_pow2(e))
        return 0;

    /* only where the alternative would use fft_small */
    if (zl == s ? !_nmod_poly_mullow_want_fft_small(bn, an, zh - zl, 0, mod)
                : !_nmod_poly_mulmid_want_fft_small(an, zh - zl, mod))
        return 0;

    /* the transform length of F for this window, and the plan of the
       alternative; a transform costs about np * ztrunc * depth */
    *P = *M->P;
    if (!_fft_small_plan_fit_window(P, zl, zh, an + M->bn - 1, an, 4))
        return 0;

    D->R = get_default_mpn_ctx();
    _fft_small_plan_set_window_nmod(D, zl - s, zh - s, an + bn - 1, FLINT_MAX(an, bn));
    if (P->np == 1 || !_fft_small_plan_set_bound(D, FLINT_MIN(an, bn),
                                         2 * NMOD_BITS(mod), MPN_CTX_NCRTS))
        D->np = P->np;

    return 2 * P->np * P->ztrunc * P->depth <= 3 * D->np * D->ztrunc * D->depth;
}

/* M, if the largest call, which is the most favourable one, would use it
   (see _fft_use); otherwise M is cleared and NULL is returned */
static mul_precomp_struct *
_fft_if_used(mul_precomp_struct * M, slong zl, slong zh, slong an,
             slong s, slong bn, nmod_t mod)
{
    if (M != NULL && !_fft_use(M, zl, zh, an, s, bn, mod))
    {
        _fft_clear(M);
        return NULL;
    }

    return M;
}

/* the product above, computed with M (which is only read) if _fft_use
   says so, returning 1; otherwise returns 0 */
int
_nmod_geometric_progression_mulmid_fft(nn_ptr z, slong zl, slong zh,
        nn_srcptr a, slong an, const mul_precomp_struct * M,
        slong s, slong bn, nmod_t mod)
{
    if (!_fft_use(M, zl, zh, an, s, bn, mod))
        return 0;

    if (M->transposed)
    {
        slong rn = zl + 1, t = rn - an, i;
        nn_ptr r;
        int done;
        TMP_INIT;

        TMP_START;
        r = TMP_ARRAY_ALLOC(rn, ulong);
        _nmod_vec_zero(r, t);
        for (i = 0; i < an; i++)
            r[t + i] = a[an - 1 - i];
        done = _nmod_poly_mul_tr_precomp(z, zh - zl, r, rn, M, mod,
                                         get_default_mpn_ctx());
        TMP_END;

        return done;
    }

    return _nmod_poly_mul_mid_precomp(z, zl, zh, a, an, M, mod, get_default_mpn_ctx());
}

#else

static struct mul_precomp_struct *
_fft_init(nn_srcptr FLINT_UNUSED(b), slong FLINT_UNUSED(bn), slong FLINT_UNUSED(n),
          slong FLINT_UNUSED(an_max), nmod_t FLINT_UNUSED(mod))
{
    return NULL;
}

static struct mul_precomp_struct *
_fft_init_tr(nn_srcptr FLINT_UNUSED(v), slong FLINT_UNUSED(vn),
             slong FLINT_UNUSED(an_max), nmod_t FLINT_UNUSED(mod))
{
    return NULL;
}

static void
_fft_clear(struct mul_precomp_struct * FLINT_UNUSED(F))
{
}

static struct mul_precomp_struct *
_fft_if_used(struct mul_precomp_struct * M, slong FLINT_UNUSED(zl), slong FLINT_UNUSED(zh),
             slong FLINT_UNUSED(an), slong FLINT_UNUSED(s), slong FLINT_UNUSED(bn),
             nmod_t FLINT_UNUSED(mod))
{
    return M;
}

int
_nmod_geometric_progression_mulmid_fft(nn_ptr FLINT_UNUSED(z), slong FLINT_UNUSED(zl),
        slong FLINT_UNUSED(zh), nn_srcptr FLINT_UNUSED(a), slong FLINT_UNUSED(an),
        const struct mul_precomp_struct * FLINT_UNUSED(M), slong FLINT_UNUSED(s), slong FLINT_UNUSED(bn),
        nmod_t FLINT_UNUSED(mod))
{
    return 0;
}

#endif

/* Stores the transforms of the polynomials which the three functionalities
   multiply by, for those whose largest product would use fft_small and
   would use the stored transform (the products using a stored transform
   are a subset of these, see _fft_use). Each functionality has its own
   bound (0: no transform; negative: the stored transforms of this
   functionality, if any, are kept), and the transform lengths allow any
   call within it, with no wraparound reaching the coefficients sought:
   - evaluate polynomials of length at most ev_len: coefficients
     [plen - 1, plen - 1 + olen) of a * ev_f, where a has length at most
     plen <= ev_len and olen <= len; this involves ev_f[0, ev_len + len - 1),
     whose images through the transpose of the inverse transform are stored
     for transposed products of length ev_len + len - 1 (see
     _mul_precomp_init_tr);
   - interpolate from the first m <= int_len points: coefficients [0, m) of
     f * int_f, f of length at most m; this involves int_f[0, int_len), and
     the convolution has length >= 2*int_len - 1;
   - extrapolate from at most ext_len values: coefficients of tmp * ext_ff
     (or ext_fb) of index at least ilen - 1, tmp of length
     ilen <= min(ext_len, len - 1), transposed products of length len - 1.
     Since ilen + olen <= len, the largest product has
     ilen = min(ext_len, len/2). */
void
_nmod_geometric_progression_init_fft(nmod_geometric_progression_t G,
                                     slong ev_len, slong int_len, slong ext_len)
{
    slong len = G->len;
    nmod_t mod = G->mod;

    FLINT_ASSERT(ev_len <= len && int_len <= len && ext_len <= len);

    G->function |= UWORD(8);

    if (ev_len >= 0)
    {
        _fft_clear(G->ev_f_fft);
        G->ev_f_fft = NULL;
    }

    if (int_len >= 0)
    {
        _fft_clear(G->int_f_fft);
        G->int_f_fft = NULL;
    }

    if (ext_len >= 0)
    {
        _fft_clear(G->ext_ff_fft);
        _fft_clear(G->ext_fb_fft);
        G->ext_ff_fft = G->ext_fb_fft = NULL;
    }

    if ((G->function & UWORD(1)) && ev_len > 0
            && _nmod_poly_mulmid_want_fft_small(ev_len, len, mod))
    {
        slong n = ev_len + len - 1;
        G->ev_f_fft = _fft_if_used(_fft_init_tr(G->ev_f->coeffs, n, ev_len, mod),
                                   ev_len - 1, n, ev_len, 0, n, mod);
    }

    if ((G->function & UWORD(2)) && int_len > 0
            && _nmod_poly_mullow_want_fft_small(int_len, int_len, int_len, 0, mod))
    {
        G->int_f_fft = _fft_if_used(_fft_init(G->int_f->coeffs, int_len, 2*int_len - 1, int_len, mod),
                                    0, int_len, int_len, 0, int_len, mod);
    }

    /* the largest products, forward and backward, have the same shape:
       elen values to the next (or previous) len - elen points */
    if ((G->function & UWORD(4)) && ext_len > 0 && len > 1)
    {
        slong elen = FLINT_MIN(ext_len, len / 2);
        slong an_max = FLINT_MIN(ext_len, len - 1);

        if (elen > 0 && _nmod_poly_mulmid_want_fft_small(elen, len - elen, mod))
        {
            G->ext_ff_fft = _fft_if_used(_fft_init_tr(G->ext_ff->coeffs, len - 1, an_max, mod),
                                         elen - 1, len - 1, elen, 0, len - 1, mod);
            if (G->ext_ff_fft != NULL)
                G->ext_fb_fft = _fft_init_tr(G->ext_fb->coeffs, len - 1, an_max, mod);
        }
    }
}

/* initialize for selection of functionalities:                        */
/*    the lowest 3 bits of `function` act as a mask for the            */
/*    three functionalities, in this order:                            */
/*    evaluate(bit0)+interpolate(bit1)+extrapolate(bit2)               */
/*    and bit3 asks for stored transforms of the polynomials involved  */
/*                                                                     */
/* TODO possible improvements, if precomputation efficiency matters:   */
/*   - extrapolate may benefit from n_mulmod_shoup                     */
/*   - interpolate version with mulmod_shoup is not fully "Shoupified" */
/*   - unrolling may help                                              */
void _nmod_geometric_progression_init_function(nmod_geometric_progression_t G,
                                               ulong r, slong len, nmod_t mod,
                                               ulong function)
{
    G->len = len;
    G->mod = mod;
    G->function = function;
    G->ev_f_fft = NULL;
    G->int_f_fft = NULL;
    G->ext_ff_fft = NULL;
    G->ext_fb_fft = NULL;

    if (NMOD_CAN_USE_SHOUP(mod))
    {
        const ulong q = nmod_mul(r, r, mod);
        const ulong inv_r = nmod_inv(r, mod);
        const ulong inv_q = nmod_mul(inv_r, inv_r, mod);

        ulong q_pr_rem;
        ulong q_pr_quo;
        ulong inv_q_pr_rem;
        ulong inv_q_pr_quo;
        _nmod_mulmod_precomp_shoup_quo_rem(&q_pr_quo, &q_pr_rem, q, mod);
        _nmod_mulmod_precomp_shoup_quo_rem(&inv_q_pr_quo, &inv_q_pr_rem, inv_q, mod);

        if (function & UWORD(1))  /* evaluate */
            _nmod_geometric_progression_evaluate_init_nonfullword(G, r, len, mod, q, q_pr_quo, q_pr_rem, inv_r, inv_q, inv_q_pr_quo, inv_q_pr_rem);

        if ((function & UWORD(6)) == UWORD(6))  /* interpolate+extrapolate */
            _nmod_geometric_progression_interpolate_extrapolate_init(G, len, mod, q, inv_q);
        else if ((function & UWORD(2)) == UWORD(2))  /* interpolate */
            _nmod_geometric_progression_interpolate_init_nonfullword(G, len, mod, q, q_pr_quo, q_pr_rem, inv_q, inv_q_pr_quo, inv_q_pr_rem);
        else if ((function & UWORD(4)) == UWORD(4))  /* extrapolate */
            _nmod_geometric_progression_extrapolate_init(G, len, mod, q, inv_q);
    }
    else
    {
        const ulong q = nmod_mul(r, r, mod);
        const ulong inv_r = nmod_inv(r, mod);
        const ulong inv_q = nmod_mul(inv_r, inv_r, mod);

        if (function & UWORD(1))  /* evaluate */
            _nmod_geometric_progression_evaluate_init(G, r, len, mod, q, inv_r, inv_q);

        if ((function & UWORD(6)) == UWORD(6))  /* interpolate+extrapolate */
            _nmod_geometric_progression_interpolate_extrapolate_init(G, len, mod, q, inv_q);
        else if ((function & UWORD(2)) == UWORD(2))  /* interpolate */
                _nmod_geometric_progression_interpolate_init(G, len, mod, q, inv_q);
        else if ((function & UWORD(4)) == UWORD(4))  /* extrapolate */
                _nmod_geometric_progression_extrapolate_init(G, len, mod, q, inv_q);
    }

    if ((function & UWORD(8)) && len > 0)
        _nmod_geometric_progression_init_fft(G, len, len, len);
}

/* clear for selection of functionalities (see init)  */
static
void _nmod_geometric_progression_clear_function(nmod_geometric_progression_t G, ulong function)
{
    if (function & UWORD(1))  /* evaluate */
    {
        _nmod_vec_clear(G->ev_s);
        nmod_poly_clear(G->ev_f);
    }
    if ((function & UWORD(6)) == UWORD(6))  /* interpolate+extrapolate */
    {
        /* G->int_s1 is ext_s2: no alloc -> not cleared */
        _nmod_vec_clear(G->int_s2);
        nmod_poly_clear(G->int_f);
        _nmod_vec_clear(G->ext_s1f);
        _nmod_vec_clear(G->ext_s1b);
        _nmod_vec_clear(G->ext_s2);
        _nmod_vec_clear(G->ext_s3);
        nmod_poly_clear(G->ext_ff);
        nmod_poly_clear(G->ext_fb);
    }
    else if ((function & UWORD(2)) == UWORD(2))  /* interpolate */
    {
        _nmod_vec_clear(G->int_s1);
        _nmod_vec_clear(G->int_s2);
        nmod_poly_clear(G->int_f);
    }
    else if ((function & UWORD(4)) == UWORD(4))  /* extrapolate */
    {
        _nmod_vec_clear(G->ext_s1f);
        _nmod_vec_clear(G->ext_s1b);
        _nmod_vec_clear(G->ext_s2);
        _nmod_vec_clear(G->ext_s3);
        nmod_poly_clear(G->ext_ff);
        nmod_poly_clear(G->ext_fb);
    }

    _fft_clear(G->ev_f_fft);
    _fft_clear(G->int_f_fft);
    _fft_clear(G->ext_ff_fft);
    _fft_clear(G->ext_fb_fft);
}

/* initialize for all: evaluate+interpolate+extrapolate */
void nmod_geometric_progression_init(nmod_geometric_progression_t G,
                                     ulong r, slong len, nmod_t mod)
{
    _nmod_geometric_progression_init_function(G, r, len, mod, UWORD(7));
}

/* clear for all: evaluate+interpolate+extrapolate */
void nmod_geometric_progression_clear(nmod_geometric_progression_t G)
{
    _nmod_geometric_progression_clear_function(G, G->function);
}
