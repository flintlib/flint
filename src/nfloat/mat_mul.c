/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "gr.h"
#include "gr_mat.h"
#include "gr_generic.h"
#include "acf.h"
#include "acb.h"
#include "nfloat.h"

#include "gr.h"
#include "nfloat.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_special.h"
#include "fmpz_mat.h"


int
nfloat_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);


/* For printing */
#include "arf.h"

/* Arithmetic on fixed-point numbers in (-1,1) */
/* x[0] stores the sign bit, x[1], ..., x[n] store the absolute value */

void
nfixed_print(nn_srcptr x, slong nlimbs, slong exp)
{
    arf_t t;
    arf_init(t);
    _arf_set_mpn_fixed(t, x + 1, nlimbs, nlimbs, x[0], nlimbs * FLINT_BITS, ARF_RND_DOWN);
    arf_mul_2exp_si(t, t, exp);
    arf_printd(t, nlimbs * FLINT_BITS / 3.321928 + 1);
    arf_clear(t);
}


/* todo: don't do this */
#define NFIXED_MAX_NLIMBS (2 * NFLOAT_MAX_LIMBS)

FLINT_FORCE_INLINE
void nfixed_add(nn_ptr res, nn_srcptr a, nn_srcptr b, slong nlimbs)
{
    int asgn, bsgn;
    asgn = a[0];
    bsgn = b[0];

    if (asgn == bsgn)
    {
        res[0] = asgn;
        mpn_add_n(res + 1, a + 1, b + 1, nlimbs);
    }
    else
    {
        res[0] = asgn ^ flint_mpn_signed_sub_n(res + 1, a + 1, b + 1, nlimbs);
    }
}

FLINT_FORCE_INLINE
void nfixed_sub(nn_ptr res, nn_srcptr a, nn_srcptr b, slong nlimbs)
{
    int asgn, bsgn;
    asgn = a[0];
    bsgn = b[0];

    if (asgn != bsgn)
    {
        res[0] = asgn;
        mpn_add_n(res + 1, a + 1, b + 1, nlimbs);
    }
    else
    {
        res[0] = asgn ^ flint_mpn_signed_sub_n(res + 1, a + 1, b + 1, nlimbs);
    }
}

FLINT_FORCE_INLINE
void _nfixed_vec_add(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs)
{
    slong i;

    for (i = 0; i < len; i++)
        nfixed_add(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1), nlimbs);
}

FLINT_FORCE_INLINE
void _nfixed_vec_sub(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs)
{
    slong i;

    for (i = 0; i < len; i++)
        nfixed_sub(res + i * (nlimbs + 1), a + i * (nlimbs + 1), b + i * (nlimbs + 1), nlimbs);
}

FLINT_FORCE_INLINE
void nfixed_mul(nn_ptr res, nn_srcptr a, nn_srcptr b, slong nlimbs)
{
    int asgn, bsgn;
    asgn = a[0];
    bsgn = b[0];
    res[0] = asgn ^ bsgn;
    flint_mpn_mulhigh_n(res + 1, a + 1, b + 1, nlimbs);
}

FLINT_FORCE_INLINE
void nfixed_sqr(nn_ptr res, nn_srcptr a, slong nlimbs)
{
    res[0] = 0;
    flint_mpn_sqrhigh(res + 1, a + 1, nlimbs);
}

FLINT_FORCE_INLINE
void nfixed_div2(nn_ptr res, nn_srcptr a, slong nlimbs)
{
    res[0] = a[0];
    mpn_rshift(res + 1, a + 1, nlimbs, 1);
}

/* A is (m x n), B is (n x p), C is (m x p) */
void
_nfixed_mat_mul_classical(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
{
    slong i, j, k;
    nn_ptr t;
    TMP_INIT;

    TMP_START;

    t = TMP_ALLOC((nlimbs + 1) * sizeof(ulong));

#define A_ENTRY(i, j) ((A) + ((i) * n + (j)) * (nlimbs + 1))
#define B_ENTRY(i, j) ((B) + ((i) * p + (j)) * (nlimbs + 1))
#define C_ENTRY(i, j) ((C) + ((i) * p + (j)) * (nlimbs + 1))

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < p; j++)
        {
            nfixed_mul(C_ENTRY(i, j), A_ENTRY(i, 0), B_ENTRY(0, j), nlimbs);

            for (k = 1; k < n; k++)
            {
                nfixed_mul(t, A_ENTRY(i, k), B_ENTRY(k, j), nlimbs);
                nfixed_add(C_ENTRY(i, j), C_ENTRY(i, j), t, nlimbs);
            }
        }
    }

    TMP_END;

#undef A_ENTRY
#undef B_ENTRY
#undef C_ENTRY
}

/* compute c += (a1 + b1) * (a2 + b2) */
/* val0, val1, val2 are scratch space */
FLINT_FORCE_INLINE void
addmul_addadd(nn_ptr val0, nn_ptr val1, nn_ptr val2, nn_ptr c, nn_srcptr a1, nn_srcptr b1, nn_srcptr a2, nn_srcptr b2, slong nlimbs)
{
    nfixed_add(val1, a1, b1, nlimbs);
    nfixed_add(val2, a2, b2, nlimbs);
    nfixed_mul(val0, val1, val2, nlimbs);
    nfixed_add(c, c, val0, nlimbs);
}

/* compute c += (a1 - b1) * (a2 - b2) */
/* val0, val1, val2 are scratch space */
FLINT_FORCE_INLINE void
addmul_subsub(nn_ptr val0, nn_ptr val1, nn_ptr val2, nn_ptr c, nn_srcptr a1, nn_srcptr b1, nn_srcptr a2, nn_srcptr b2, slong nlimbs)
{
    nfixed_sub(val1, a1, b1, nlimbs);
    nfixed_sub(val2, a2, b2, nlimbs);
    nfixed_mul(val0, val1, val2, nlimbs);
    nfixed_add(c, c, val0, nlimbs);
}

void
_nfixed_mat_mul_waksman(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
{
    slong l, j, k;

    nn_ptr Ctmp = flint_calloc((nlimbs + 1) * ((p + m) + 5), sizeof(ulong));
                                            /* Ctmp itself has m * p entries */
    nn_ptr Crow = Ctmp;                     /* Crow has p entries */
    nn_ptr Ccol = Crow + (nlimbs + 1) * p;  /* Ccol has m entries */
    nn_ptr val0 = Ccol + (nlimbs + 1) * m;  /* val0 has room for 2 sums */
    nn_ptr val1 = val0 + (nlimbs + 1) * 2;  /* val1 has room for 1 sum   */
    nn_ptr val2 = val1 + (nlimbs + 1);      /* val2 has room for 1 sum   */
    nn_ptr crow = val2 + (nlimbs + 1);      /* crow has room for 1 sum   */

#define A_ENTRY(i, j) ((A) + ((i) * n + (j)) * (nlimbs + 1))
#define B_ENTRY(i, j) ((B) + ((i) * p + (j)) * (nlimbs + 1))
#define C_ENTRY(i, j) ((C) + ((i) * p + (j)) * (nlimbs + 1))

#define Crow_ENTRY(ii) (Crow + (ii) * (nlimbs + 1))
#define Ccol_ENTRY(ii) (Ccol + (ii) * (nlimbs + 1))

    slong np = n >> 1;

    for (j = 1; j <= np; j++)
    {
        slong j2 = (j << 1) - 1;

        for (k = 0; k < p; k++)
        {
            addmul_addadd(val0, val1, val2, C_ENTRY(0, k), A_ENTRY(0, j2-1), B_ENTRY(j2, k), A_ENTRY(0, j2), B_ENTRY(j2-1, k), nlimbs);
            addmul_subsub(val0, val1, val2, Crow_ENTRY(k), A_ENTRY(0, j2-1), B_ENTRY(j2, k), A_ENTRY(0, j2), B_ENTRY(j2-1, k), nlimbs);
        }

        for (l = 1; l < m; l++)
        {
            addmul_addadd(val0, val1, val2, C_ENTRY(l, 0), A_ENTRY(l, j2-1), B_ENTRY(j2, 0), A_ENTRY(l, j2), B_ENTRY(j2-1, 0), nlimbs);
            addmul_subsub(val0, val1, val2, Ccol_ENTRY(l), A_ENTRY(l, j2-1), B_ENTRY(j2, 0), A_ENTRY(l, j2), B_ENTRY(j2-1, 0), nlimbs);
        }

        for (k = 1; k < p; k++)
        {
            for (l = 1; l < m; l++)
            {
                addmul_addadd(val0, val1, val2, C_ENTRY(l, k), A_ENTRY(l, j2-1), B_ENTRY(j2, k), A_ENTRY(l, j2), B_ENTRY(j2-1, k), nlimbs);
            }
        }
    }

    for (l = 1; l < m; l++)
    {
        nfixed_add(val1, Ccol_ENTRY(l), C_ENTRY(l, 0), nlimbs);
        nfixed_div2(Ccol_ENTRY(l), val1, nlimbs);
        nfixed_sub(C_ENTRY(l, 0), C_ENTRY(l, 0), Ccol_ENTRY(l), nlimbs);
    }

    nfixed_add(val1, Crow, C_ENTRY(0, 0), nlimbs);
    nfixed_div2(val0, val1, nlimbs);
    nfixed_sub(C_ENTRY(0, 0), C_ENTRY(0, 0), val0, nlimbs);

    for (k = 1; k < p; k++)
    {
        nfixed_add(crow, Crow_ENTRY(k), C_ENTRY(0, k), nlimbs);
        nfixed_div2(val1, crow, nlimbs);
        nfixed_sub(C_ENTRY(0, k), C_ENTRY(0, k), val1, nlimbs);
        nfixed_sub(crow, val1, val0, nlimbs);

        for (l = 1; l < m; l++)
        {
            nfixed_sub(val2, C_ENTRY(l, k), crow, nlimbs);
            nfixed_sub(C_ENTRY(l, k), val2, Ccol_ENTRY(l), nlimbs);
        }
    }

    if ((n & 1) == 1)
    {
        for (l = 0; l < m; l++)
        {
            for (k = 0; k < p; k++)
            {
                nfixed_mul(val0, A_ENTRY(l, n-1), B_ENTRY(n-1, k), nlimbs);
                nfixed_add(C_ENTRY(l, k), C_ENTRY(l, k), val0, nlimbs);
            }
        }
    }

    flint_free(Ctmp);

#undef A_ENTRY
#undef B_ENTRY
#undef C_ENTRY
}

FLINT_FORCE_INLINE void
_nfloat_get_nfixed(nn_ptr res, nn_srcptr x, slong exp, slong fix_nlimbs, gr_ctx_t ctx)
{
    slong rel_exp;

    /* assumes res is already zeroed */
    if (NFLOAT_IS_ZERO(x))
        return;

    rel_exp = NFLOAT_EXP(x) - exp;
    if (rel_exp >= 0)
        flint_abort();

    res[0] = NFLOAT_SGNBIT(x);
    _arf_get_integer_mpn(res + 1, NFLOAT_D(x), NFLOAT_CTX_NLIMBS(ctx), fix_nlimbs * FLINT_BITS + rel_exp);
}

FLINT_FORCE_INLINE int
_nfloat_set_nfixed(nn_ptr res, nn_srcptr x, slong exp, slong fix_nlimbs, gr_ctx_t ctx)
{
    return nfloat_set_mpn_2exp(res, x + 1, fix_nlimbs, exp, x[0], ctx);
}

static void
_nfloat_mat_exp_range(slong * _Amin, slong * _Amax, const gr_mat_t A, gr_ctx_t ctx)
{
    slong Amax, Amin;
    slong m = A->r;
    slong n = A->c;
    slong exp, i, j;
    slong sz = ctx->sizeof_elem;

    Amax = WORD_MIN;
    Amin = WORD_MAX;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            exp = NFLOAT_EXP(GR_MAT_ENTRY(A, i, j, sz));
            Amax = FLINT_MAX(Amax, exp);
            Amin = FLINT_MIN(Amin, exp);
        }
    }

    _Amin[0] = Amin;
    _Amax[0] = Amax;
}

int
_nfloat_mat_mul_fixed_given_exp(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong Aexp, slong Bexp, slong fnlimbs, int waksman, gr_ctx_t ctx)
{
    nn_ptr T, TA, TB, TC;
    slong i, j;
    slong sz = ctx->sizeof_elem;
    slong fdnlimbs;

    slong m = A->r;
    slong n = A->c;
    slong p = B->c;

    /* limbs including sign limb */
    fdnlimbs = fnlimbs + 1;

    T = flint_calloc(fdnlimbs * (m * n + n * p + m * p), sizeof(ulong));

    TA = T;
    TB = TA + fdnlimbs * (m * n);
    TC = TB + fdnlimbs * (n * p);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            _nfloat_get_nfixed(TA + i * fdnlimbs * n + j * fdnlimbs, GR_MAT_ENTRY(A, i, j, sz), Aexp, fnlimbs, ctx);

    for (i = 0; i < n; i++)
        for (j = 0; j < p; j++)
            _nfloat_get_nfixed(TB + i * fdnlimbs * p + j * fdnlimbs, GR_MAT_ENTRY(B, i, j, sz), Bexp, fnlimbs, ctx);

    if (waksman)
        _nfixed_mat_mul_waksman(TC, TA, TB, m, n, p, fnlimbs);
    else
        _nfixed_mat_mul_classical(TC, TA, TB, m, n, p, fnlimbs);

    for (i = 0; i < m; i++)
        for (j = 0; j < p; j++)
            _nfloat_set_nfixed(GR_MAT_ENTRY(C, i, j, sz), TC + i * fdnlimbs * p + j * fdnlimbs, Aexp + Bexp, fnlimbs, ctx);

    flint_free(T);

    return GR_SUCCESS;
}

int
_nfloat_mat_mul_fixed(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int waksman, slong max_extra_bits, gr_ctx_t ctx)
{
    slong Amax, Amin, Bmax, Bmin, Adelta, Bdelta, Aexp, Bexp;
    slong prec;
    slong pad_top, pad_bot, extra_bits, fbits, fnlimbs;
    slong n = A->c;

    if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        return GR_UNABLE;

    prec = NFLOAT_CTX_PREC(ctx);

    _nfloat_mat_exp_range(&Amin, &Amax, A, ctx);
    _nfloat_mat_exp_range(&Bmin, &Bmax, B, ctx);

    if (Amax < NFLOAT_MIN_EXP || Bmax < NFLOAT_MIN_EXP)
        return gr_mat_zero(C, ctx);

    /* Currently, we don't handle zeros. (They pose no problem, but zero entries in
       the output may not be exact. To be done.) */
    if (Amin < NFLOAT_MIN_EXP || Bmin < NFLOAT_MIN_EXP)
        return gr_mat_mul_classical(C, A, B, ctx);

    Adelta = Amax - Amin;
    Bdelta = Bmax - Bmin;

    /* sanity check */
    if (Adelta > 10 * prec || Bdelta > 10 * prec)
        return gr_mat_mul_classical(C, A, B, ctx);

    /*
    To double check: for Waksman,
        * The intermediate entries are bounded by 8n max(|A|,|B|)^2.
        * The error, including error from converting
          the input matrices, is bounded by 8n ulps.
    */

    pad_top = 3 + FLINT_BIT_COUNT(n);
    pad_bot = 3 + FLINT_BIT_COUNT(n);

    extra_bits = Adelta + Bdelta + pad_top + pad_bot;

    if (extra_bits >= max_extra_bits)
        return gr_mat_mul_classical(C, A, B, ctx);

    Aexp = Amax + pad_top;
    Bexp = Bmax + pad_top;
    fbits = prec + extra_bits;
    fnlimbs = (fbits + FLINT_BITS - 1) / FLINT_BITS;

    return _nfloat_mat_mul_fixed_given_exp(C, A, B, Aexp, Bexp, fnlimbs, waksman, ctx);
}

int
nfloat_mat_mul_fixed_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    return _nfloat_mat_mul_fixed(C, A, B, 0, 100000, ctx);
}

int
nfloat_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    return _nfloat_mat_mul_fixed(C, A, B, 1, 100000, ctx);
}


static void
_nfloat_2exp_get_fmpz(fmpz_t res, nfloat_srcptr x, slong fixexp, gr_ctx_t ctx)
{
    slong exp, zn;
    mpz_ptr zz;
    nn_ptr zp;
    int negative;
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);

    if (NFLOAT_IS_SPECIAL(x))
    {
        fmpz_zero(res);
        return;
    }

    exp = NFLOAT_EXP(x) - fixexp;

    if (exp <= 0)
    {
        fmpz_zero(res);
        return;
    }

    /* todo: small case */

    negative = NFLOAT_SGNBIT(x);

    zn = (exp + FLINT_BITS - 1) / FLINT_BITS;
    zz = _fmpz_promote(res);
    zp = FLINT_MPZ_REALLOC(zz, zn);
    _arf_get_integer_mpn(zp, NFLOAT_D(x), nlimbs, exp);
    zz->_mp_size = negative ? -zn : zn;
    _fmpz_demote_val(res);
}

int
nfloat_mat_addmul_block_fallback(gr_mat_t C,
    const gr_mat_t A, const gr_mat_t B,
    slong block_start,
    slong block_end,
    gr_ctx_t ctx)
{
    slong M, P, n;
    slong i, j, sz;
    nn_ptr tmpB;
    slong ndlimbs = NFLOAT_CTX_DATA_NLIMBS(ctx);
    sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    M = A->r;
    P = B->c;

    n = block_end - block_start;

    tmpB = flint_malloc(sizeof(ulong) * ndlimbs * (P * n));

#define AA(ii, jj) GR_MAT_ENTRY(A, ii, block_start + (jj), sz)

    for (i = 0; i < P; i++)
        for (j = 0; j < n; j++)
            flint_mpn_copyi(tmpB + (i * n + j) * ndlimbs, GR_MAT_ENTRY(B, block_start + j, i, sz), ndlimbs);

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < P; j++)
        {
            status |= _nfloat_vec_dot(GR_MAT_ENTRY(C, i, j, sz),
                (block_start == 0) ? NULL : GR_MAT_ENTRY(C, i, j, sz), 0,
                GR_MAT_ENTRY(A, i, block_start, sz),
                tmpB + j * n * ndlimbs, n, ctx);
        }
    }

    flint_free(tmpB);

    return status;
}

int
nfloat_mat_addmul_block_prescaled(gr_mat_t C,
    const gr_mat_t A, const gr_mat_t B,
    slong block_start,
    slong block_end,
    const slong * A_min,  /* A per-row bottom exponent */
    const slong * B_min,  /* B per-row bottom exponent */
    gr_ctx_t ctx)
{
    slong M, P, n;
    slong i, j;
    slong M0, M1, P0, P1, Mstep, Pstep;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong t[NFLOAT_MAX_ALLOC];
    slong e;

    M = A->r;
    P = B->c;

    n = block_end - block_start;

    /* Create sub-blocks to keep matrices nearly square. Necessary? */
#if 1
    Mstep = (M < 2 * n) ? M : n;
    Pstep = (P < 2 * n) ? P : n;
#else
    Mstep = M;
    Pstep = P;
#endif

    for (M0 = 0; M0 < M; M0 += Mstep)
    {
        for (P0 = 0; P0 < P; P0 += Pstep)
        {
            fmpz_mat_t AA, BB, CC;

            M1 = FLINT_MIN(M0 + Mstep, M);
            P1 = FLINT_MIN(P0 + Pstep, P);

            fmpz_mat_init(AA, M1 - M0, n);
            fmpz_mat_init(BB, n, P1 - P0);
            fmpz_mat_init(CC, M1 - M0, P1 - P0);

            /* Convert to fixed-point matrices. */
            for (i = M0; i < M1; i++)
            {
                if (A_min[i] == WORD_MIN)  /* only zeros in this row */
                    continue;

                for (j = 0; j < n; j++)
                    _nfloat_2exp_get_fmpz(fmpz_mat_entry(AA, i - M0, j), GR_MAT_ENTRY(A, i, block_start + j, sz), A_min[i], ctx);
            }

            for (i = P0; i < P1; i++)
            {
                if (B_min[i] == WORD_MIN)  /* only zeros in this column */
                    continue;

                for (j = 0; j < n; j++)
                    _nfloat_2exp_get_fmpz(fmpz_mat_entry(BB, j, i - P0), GR_MAT_ENTRY(B, block_start + j, i, sz), B_min[i], ctx);
            }

            /* The main multiplication */
            fmpz_mat_mul(CC, AA, BB);
            fmpz_mat_clear(AA);
            fmpz_mat_clear(BB);

            /* Add to the result matrix */
            for (i = M0; i < M1; i++)
            {
                for (j = P0; j < P1; j++)
                {
                    e = A_min[i] + B_min[j];

                    /* The first time we write this Cij */
                    if (block_start == 0)
                    {
                        status |= nfloat_set_fmpz(GR_MAT_ENTRY(C, i, j, sz), fmpz_mat_entry(CC, i - M0, j - P0), ctx);
                        status |= nfloat_mul_2exp_si(GR_MAT_ENTRY(C, i, j, sz), GR_MAT_ENTRY(C, i, j, sz), e, ctx);
                    }
                    else
                    {
                        status |= nfloat_set_fmpz(t, fmpz_mat_entry(CC, i - M0, j - P0), ctx);
                        status |= nfloat_mul_2exp_si(t, t, e, ctx);
                        status |= nfloat_add(GR_MAT_ENTRY(C, i, j, sz), GR_MAT_ENTRY(C, i, j, sz), t, ctx);
                    }
                }
            }

            fmpz_mat_clear(CC);
        }
    }

    return status;
}

FLINT_FORCE_INLINE slong
_nfloat_nbits(nfloat_srcptr x, slong nlimbs)
{
    nn_srcptr ad;
    slong bits;

    ad = NFLOAT_D(x);
    bits = FLINT_BITS * nlimbs;

    while (ad[0] == 0)
    {
        bits -= FLINT_BITS;
        ad++;
    }

    bits -= flint_ctz(ad[0]);

    return bits;
}

int
nfloat_complex_mat_addmul_block_fallback(gr_mat_t C,
    const gr_mat_t A, const gr_mat_t B,
    slong block_start,
    slong block_end,
    gr_ctx_t ctx)
{
    slong M, P, n;
    slong i, j, sz;
    nn_ptr tmpB;
    slong ndlimbs = NFLOAT_COMPLEX_CTX_DATA_NLIMBS(ctx);
    sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    M = A->r;
    P = B->c;

    n = block_end - block_start;

    tmpB = flint_malloc(sizeof(ulong) * ndlimbs * (P * n));

#define AA(ii, jj) GR_MAT_ENTRY(A, ii, block_start + (jj), sz)

    for (i = 0; i < P; i++)
        for (j = 0; j < n; j++)
            flint_mpn_copyi(tmpB + (i * n + j) * ndlimbs, GR_MAT_ENTRY(B, block_start + j, i, sz), ndlimbs);

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < P; j++)
        {
            status |= _nfloat_complex_vec_dot(GR_MAT_ENTRY(C, i, j, sz),
                (block_start == 0) ? NULL : GR_MAT_ENTRY(C, i, j, sz), 0,
                GR_MAT_ENTRY(A, i, block_start, sz),
                tmpB + j * n * ndlimbs, n, ctx);
        }
    }

    flint_free(tmpB);

    return status;
}

static void
_nfloat_complex_2exp_get_fmpz_fmpz(fmpz_t res1, fmpz_t res2, nfloat_complex_srcptr x, slong fixexp, gr_ctx_t ctx)
{
    _nfloat_2exp_get_fmpz(res1, NFLOAT_COMPLEX_RE(x, ctx), fixexp, ctx);
    _nfloat_2exp_get_fmpz(res2, NFLOAT_COMPLEX_IM(x, ctx), fixexp, ctx);
}

int
nfloat_complex_mat_addmul_block_prescaled(gr_mat_t C,
    const gr_mat_t A, const gr_mat_t B,
    slong block_start,
    slong block_end,
    const slong * A_min,  /* A per-row bottom exponent */
    const slong * B_min,  /* B per-row bottom exponent */
    gr_ctx_t ctx)
{
    slong M, P, n;
    slong i, j;
    slong M0, M1, P0, P1, Mstep, Pstep;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong t[NFLOAT_MAX_ALLOC];
    slong e;

    M = A->r;
    P = B->c;

    n = block_end - block_start;

    /* Create sub-blocks to keep matrices nearly square. Necessary? */
#if 1
    Mstep = (M < 2 * n) ? M : n;
    Pstep = (P < 2 * n) ? P : n;
#else
    Mstep = M;
    Pstep = P;
#endif

    for (M0 = 0; M0 < M; M0 += Mstep)
    {
        for (P0 = 0; P0 < P; P0 += Pstep)
        {
            fmpz_mat_t Aa, Ab, Ba, Bb, AaBa, AbBb, Cb;

            M1 = FLINT_MIN(M0 + Mstep, M);
            P1 = FLINT_MIN(P0 + Pstep, P);

            fmpz_mat_init(Aa, M1 - M0, n);
            fmpz_mat_init(Ab, M1 - M0, n);

            fmpz_mat_init(Ba, n, P1 - P0);
            fmpz_mat_init(Bb, n, P1 - P0);

            fmpz_mat_init(AaBa, M1 - M0, P1 - P0);
            fmpz_mat_init(AbBb, M1 - M0, P1 - P0);
            fmpz_mat_init(Cb, M1 - M0, P1 - P0);

            /* Convert to fixed-point matrices. */
            for (i = M0; i < M1; i++)
            {
                if (A_min[i] == WORD_MIN)  /* only zeros in this row */
                    continue;

                for (j = 0; j < n; j++)
                    _nfloat_complex_2exp_get_fmpz_fmpz(fmpz_mat_entry(Aa, i - M0, j),
                        fmpz_mat_entry(Ab, i - M0, j), GR_MAT_ENTRY(A, i, block_start + j, sz), A_min[i], ctx);
            }

            for (i = P0; i < P1; i++)
            {
                if (B_min[i] == WORD_MIN)  /* only zeros in this column */
                    continue;

                for (j = 0; j < n; j++)
                    _nfloat_complex_2exp_get_fmpz_fmpz(fmpz_mat_entry(Ba, j, i - P0),
                        fmpz_mat_entry(Bb, j, i - P0), GR_MAT_ENTRY(B, block_start + j, i, sz), B_min[i], ctx);
            }

            fmpz_mat_mul(AaBa, Aa, Ba);
            fmpz_mat_mul(AbBb, Ab, Bb);
            fmpz_mat_add(Aa, Aa, Ab);
            fmpz_mat_add(Ba, Ba, Bb);
            fmpz_mat_mul(Cb, Aa, Ba);
            fmpz_mat_sub(Cb, Cb, AaBa);
            fmpz_mat_sub(Cb, Cb, AbBb);

            /* Add to the result matrix */
            for (i = M0; i < M1; i++)
            {
                for (j = P0; j < P1; j++)
                {
                    nn_ptr cc = NFLOAT_COMPLEX_IM(GR_MAT_ENTRY(C, i, j, sz), ctx);

                    e = A_min[i] + B_min[j];

                    /* The first time we write this Cij */
                    if (block_start == 0)
                    {
                        status |= nfloat_set_fmpz(cc, fmpz_mat_entry(Cb, i - M0, j - P0), ctx);
                        status |= nfloat_mul_2exp_si(cc, cc, e, ctx);
                    }
                    else
                    {
                        status |= nfloat_set_fmpz(t, fmpz_mat_entry(Cb, i - M0, j - P0), ctx);
                        status |= nfloat_mul_2exp_si(t, t, e, ctx);
                        status |= nfloat_add(cc, cc, t, ctx);
                    }
                }
            }

            fmpz_mat_sub(Cb, AaBa, AbBb);

            /* Add to the result matrix */
            for (i = M0; i < M1; i++)
            {
                for (j = P0; j < P1; j++)
                {
                    nn_ptr cc = NFLOAT_COMPLEX_RE(GR_MAT_ENTRY(C, i, j, sz), ctx);

                    e = A_min[i] + B_min[j];

                    /* The first time we write this Cij */
                    if (block_start == 0)
                    {
                        status |= nfloat_set_fmpz(cc, fmpz_mat_entry(Cb, i - M0, j - P0), ctx);
                        status |= nfloat_mul_2exp_si(cc, cc, e, ctx);
                    }
                    else
                    {
                        status |= nfloat_set_fmpz(t, fmpz_mat_entry(Cb, i - M0, j - P0), ctx);
                        status |= nfloat_mul_2exp_si(t, t, e, ctx);
                        status |= nfloat_add(cc, cc, t, ctx);
                    }
                }
            }

            fmpz_mat_clear(Aa);
            fmpz_mat_clear(Ab);

            fmpz_mat_clear(Ba);
            fmpz_mat_clear(Bb);

            fmpz_mat_clear(AaBa);
            fmpz_mat_clear(AbBb);
            fmpz_mat_clear(Cb);
        }
    }

    return status;
}


/* todo: squaring optimizations */
/* note: also supports complex contexts */
int
nfloat_mat_mul_block(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong min_block_size, gr_ctx_t ctx)
{
    slong M, N, P;
    slong *A_min, *A_max, *B_min, *B_max;
    short *A_bits, *B_bits;
    slong *A_bot, *B_bot;
    slong block_start, block_end, i, j, bot, top, max_height;
    slong b, A_max_bits, B_max_bits;
    nfloat_srcptr t;
    double A_density, B_density;
    slong sz = ctx->sizeof_elem;
    slong nlimbs = NFLOAT_CTX_NLIMBS(ctx);
    slong prec = NFLOAT_CTX_PREC(ctx);
    int status = GR_SUCCESS;
    int complex_ctx = (ctx->which_ring == GR_CTX_NFLOAT_COMPLEX);

    M = A->r;
    N = A->c;
    P = B->c;

    if (N != B->r || M != C->r || P != C->c)
        return GR_DOMAIN;

    if (M == 0 || N == 0 || P == 0)
        return gr_mat_zero(C, ctx);

    if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        return GR_UNABLE;

    if (A == C || B == C)
    {
        gr_mat_t T;
        gr_mat_init(T, M, P, ctx);
        status = nfloat_mat_mul_block(T, A, B, min_block_size, ctx);
        status |= gr_mat_swap_entrywise(T, C, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

    /* bottom exponents of A */
    A_bot = flint_malloc(sizeof(slong) * M * N);
    /* minimum bottom exponent in current row */
    A_min = flint_malloc(sizeof(slong) * M);
    /* maximum top exponent in current row */
    A_max = flint_malloc(sizeof(slong) * M);

    B_bot = flint_malloc(sizeof(slong) * N * P);
    B_min = flint_malloc(sizeof(slong) * P);
    B_max = flint_malloc(sizeof(slong) * P);

    /* save space using shorts to store the bit sizes temporarily;
       the block algorithm will not be used at extremely high precision */
    A_bits = flint_malloc(sizeof(short) * M * N);
    B_bits = flint_malloc(sizeof(short) * N * P);

    A_max_bits = B_max_bits = 0;
    A_density = B_density = 0;

    /* Build table of bottom exponents (WORD_MIN signifies a zero),
       and also collect some statistics. */
    if (complex_ctx)
    {
        for (i = 0; i < M; i++)
        {
            for (j = 0; j < N; j++)
            {
                nfloat_srcptr re, im;
                slong h, b1, b2;

                t = GR_MAT_ENTRY(A, i, j, sz);
                re = NFLOAT_COMPLEX_RE(t, ctx);
                im = NFLOAT_COMPLEX_IM(t, ctx);

                if (NFLOAT_IS_ZERO(re) && NFLOAT_IS_ZERO(im))
                {
                    A_bot[i * N + j] = WORD_MIN;
                    A_bits[i * N + j] = 0;
                }
                else if (NFLOAT_IS_ZERO(im))
                {
                    b = _nfloat_nbits(re, nlimbs);
                    A_bot[i * N + j] = NFLOAT_EXP(re) - b;
                    A_bits[i * N + j] = b;
                    A_max_bits = FLINT_MAX(A_max_bits, b);
                    A_density++;
                }
                else if (NFLOAT_IS_ZERO(re))
                {
                    b = _nfloat_nbits(im, nlimbs);
                    A_bot[i * N + j] = NFLOAT_EXP(im) - b;
                    A_bits[i * N + j] = b;
                    A_max_bits = FLINT_MAX(A_max_bits, b);
                    A_density++;
                }
                else
                {
                    b1 = _nfloat_nbits(re, nlimbs);
                    b2 = _nfloat_nbits(im, nlimbs);

                    A_bot[i * N + j] = b = FLINT_MIN(NFLOAT_EXP(re) - b1, NFLOAT_EXP(im) - b2);
                    A_bits[i * N + j] = h = FLINT_MAX(NFLOAT_EXP(re) - b, NFLOAT_EXP(im) - b);
                    A_max_bits = FLINT_MAX(A_max_bits, h);
                    A_density++;
                }
            }
        }

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < P; j++)
            {
                nfloat_srcptr re, im;
                slong h, b1, b2;

                t = GR_MAT_ENTRY(B, i, j, sz);
                re = NFLOAT_COMPLEX_RE(t, ctx);
                im = NFLOAT_COMPLEX_IM(t, ctx);

                if (NFLOAT_IS_ZERO(re) && NFLOAT_IS_ZERO(im))
                {
                    B_bot[i * P + j] = WORD_MIN;
                    B_bits[i * P + j] = 0;
                }
                else if (NFLOAT_IS_ZERO(im))
                {
                    b = _nfloat_nbits(re, nlimbs);
                    B_bot[i * P + j] = NFLOAT_EXP(re) - b;
                    B_bits[i * P + j] = b;
                    B_max_bits = FLINT_MAX(B_max_bits, b);
                    B_density++;
                }
                else if (NFLOAT_IS_ZERO(re))
                {
                    b = _nfloat_nbits(im, nlimbs);
                    B_bot[i * P + j] = NFLOAT_EXP(im) - b;
                    B_bits[i * P + j] = b;
                    B_max_bits = FLINT_MAX(B_max_bits, b);
                    B_density++;
                }
                else
                {
                    b1 = _nfloat_nbits(re, nlimbs);
                    b2 = _nfloat_nbits(im, nlimbs);

                    B_bot[i * P + j] = b = FLINT_MIN(NFLOAT_EXP(re) - b1, NFLOAT_EXP(im) - b2);
                    B_bits[i * P + j] = h = FLINT_MAX(NFLOAT_EXP(re) - b, NFLOAT_EXP(im) - b);
                    B_max_bits = FLINT_MAX(B_max_bits, h);
                    B_density++;
                }
            }
        }
    }
    else
    {
        for (i = 0; i < M; i++)
        {
            for (j = 0; j < N; j++)
            {
                t = GR_MAT_ENTRY(A, i, j, sz);
                if (NFLOAT_IS_ZERO(t))
                {
                    A_bot[i * N + j] = WORD_MIN;
                    A_bits[i * N + j] = 0;
                }
                else
                {
                    b = _nfloat_nbits(t, nlimbs);
                    A_bot[i * N + j] = NFLOAT_EXP(t) - b;
                    A_bits[i * N + j] = b;
                    A_max_bits = FLINT_MAX(A_max_bits, b);
                    A_density++;
                }
            }
        }

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < P; j++)
            {
                t = GR_MAT_ENTRY(B, i, j, sz);
                if (NFLOAT_IS_ZERO(t))
                {
                    B_bot[i * P + j] = WORD_MIN;
                    B_bits[i * P + j] = 0;
                }
                else
                {
                    b = _nfloat_nbits(t, nlimbs);
                    B_bot[i * P + j] = NFLOAT_EXP(t) - b;
                    B_bits[i * P + j] = b;
                    B_max_bits = FLINT_MAX(B_max_bits, b);
                    B_density++;
                }
            }
        }
    }

    A_density = A_density / (M * N);
    B_density = B_density / (N * P);

    /* Don't shift too far when creating integer block matrices. */
    max_height = 1.25 * FLINT_MIN(prec, FLINT_MAX(A_max_bits, B_max_bits)) + 192;

    /* FIXME: this condition is bogus */
    if (A_density < 0.1 && B_density < 0.1 && max_height > 1024)
    {
        status = gr_mat_mul_classical(C, A, B, ctx);
        goto cleanup;
    }

    block_start = 0;
    while (block_start < N)
    {
        /* Find a run of columns of A and rows of B such that the
           bottom exponents differ by at most max_height. */

        block_end = block_start + 1;  /* index is exclusive block_end */

        /* begin with this column of A and row of B */
        for (i = 0; i < M; i++)
        {
            A_max[i] = A_min[i] = A_bot[i * N + block_start];
            A_max[i] += (slong) A_bits[i * N + block_start];
        }

        for (i = 0; i < P; i++)
        {
            B_max[i] = B_min[i] = B_bot[block_start * P + i];
            B_max[i] += (slong) B_bits[block_start * P + i];
        }

        while (block_end < N)
        {
            double size;

            /* End block if memory would be excessive. */
            /* Necessary? */
            /* Should also do initial check above, if C alone is too large. */
            size = (block_end - block_start) * M * (double) A_max_bits;
            size += (block_end - block_start) * P * (double) B_max_bits;
            size += (M * P) * (double) (A_max_bits + B_max_bits);
            size /= 8.0;
            if (size > 2e9)
                goto blocks_built;

            /* check if we can extend with column [block_end] of A */
            for (i = 0; i < M; i++)
            {
                bot = A_bot[i * N + block_end];
                /* zeros are irrelevant */
                if (bot == WORD_MIN || A_max[i] == WORD_MIN)
                    continue;
                top = bot + (slong) A_bits[i * N + block_end];
                /* jump will be too big */
                if (top > A_min[i] + max_height || bot < A_max[i] - max_height)
                    goto blocks_built;
            }

            /* check if we can extend with row [block_end] of B */
            for (i = 0; i < P; i++)
            {
                bot = B_bot[block_end * P + i];
                if (bot == WORD_MIN || B_max[i] == WORD_MIN)
                    continue;
                top = bot + (slong) B_bits[block_end * P + i];
                if (top > B_min[i] + max_height || bot < B_max[i] - max_height)
                    goto blocks_built;
            }

            /* second pass to update the extreme values */
            for (i = 0; i < M; i++)
            {
                bot = A_bot[i * N + block_end];
                top = bot + (slong) A_bits[i * N + block_end];
                if (A_max[i] == WORD_MIN)
                {
                    A_max[i] = top;
                    A_min[i] = bot;
                }
                else if (bot != WORD_MIN)
                {
                    if (bot < A_min[i]) A_min[i] = bot;
                    if (top > A_max[i]) A_max[i] = top;
                }
            }

            for (i = 0; i < P; i++)
            {
                bot = B_bot[block_end * P + i];
                top = bot + (slong) B_bits[block_end * P + i];
                if (B_max[i] == WORD_MIN)
                {
                    B_max[i] = top;
                    B_min[i] = bot;
                }
                else if (bot != WORD_MIN)
                {
                    if (bot < B_min[i]) B_min[i] = bot;
                    if (top > B_max[i]) B_max[i] = top;
                }
            }

            block_end++;
        }

    blocks_built:
        if (block_end - block_start < min_block_size)
        {
            block_end = FLINT_MIN(N, block_start + min_block_size);

            if (complex_ctx)
                status |= nfloat_complex_mat_addmul_block_fallback(C, A, B, block_start, block_end, ctx);
            else
                status |= nfloat_mat_addmul_block_fallback(C, A, B, block_start, block_end, ctx);
        }
        else
        {
            if (complex_ctx)
                status |= nfloat_complex_mat_addmul_block_prescaled(C, A, B, block_start, block_end, A_min, B_min, ctx);
            else
                status |= nfloat_mat_addmul_block_prescaled(C, A, B, block_start, block_end, A_min, B_min, ctx);
        }

        block_start = block_end;
    }

cleanup:
    flint_free(A_bot);
    flint_free(A_max);
    flint_free(A_min);
    flint_free(B_bot);
    flint_free(B_max);
    flint_free(B_min);
    flint_free(A_bits);
    flint_free(B_bits);

    return status;
}

/* Minimum precision for using fixed-point arithmetic */

/* TODO: for *unsigned* matrices, there is a speedup already for
   prec = 192. Consider inlining fixed-point additions/subtractions for
   4 and 5 limbs to extend this to general matrices. */
/* #define NFLOAT_MAT_MUL_FIXED_CUTOFF 192 */
#define NFLOAT_MAT_MUL_FIXED_CUTOFF 320

/* first cutoff:  classical -> fixed_classical
   second cutoff: fixed_classical -> waksman */
static const int nfloat_mat_mul_cutoff_tab[][2] = {
    {0,  0},    /* prec = 0   */
    {0,  0},    /* prec = 64  */
    {0,  0},    /* prec = 128 */
    {32, 32},   /* prec = 192 */
    {8, 20},    /* prec = 256 */
    {4, 15},    /* prec = 320 */
    {3, 10},    /* prec = 384 */
    {3, 10},    /* prec = 448 */
    {3, 8},     /* prec = 512 */
    {10, 10},   /* prec = 576 */
    {4, 5},     /* prec = 640 */
};

/* {4, 4} from this point */
#define NFLOAT_MAT_MUL_CUTOFF_4 704
/* {3, 3} from this point */
#define NFLOAT_MAT_MUL_CUTOFF_3 1600

int
nfloat_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong cutoff1, cutoff2, dim;
    int use_waksman = 0;
    slong prec;
    slong max_extra_prec;

    slong m = A->r;
    slong n = A->c;
    slong p = B->c;

    dim = FLINT_MIN(n, FLINT_MIN(m, p));

    if (dim <= 2 || NFLOAT_CTX_HAS_INF_NAN(ctx))
        return gr_mat_mul_classical(C, A, B, ctx);

    if (dim <= 80)
    {
        prec = NFLOAT_CTX_PREC(ctx);

        if (prec < NFLOAT_MAT_MUL_FIXED_CUTOFF)
            return gr_mat_mul_classical(C, A, B, ctx);

        if (prec >= NFLOAT_MAT_MUL_CUTOFF_3)
            cutoff1 = cutoff2 = 3;
        else if (prec >= NFLOAT_MAT_MUL_CUTOFF_4)
            cutoff1 = cutoff2 = 4;
        else
        {
            cutoff1 = nfloat_mat_mul_cutoff_tab[prec / 64][0];
            cutoff2 = nfloat_mat_mul_cutoff_tab[prec / 64][1];
        }

        if (dim < cutoff1)
            return gr_mat_mul_classical(C, A, B, ctx);

        use_waksman = (dim >= cutoff2);
        max_extra_prec = (prec < 768) ? 64 : prec / 4;

        return _nfloat_mat_mul_fixed(C, A, B, use_waksman, max_extra_prec, ctx);
    }
    else
    {
        return nfloat_mat_mul_block(C, A, B, 70, ctx);
    }
}

static void
_nfloat_complex_mat_exp_range(slong * _Amin, slong * _Amax, const gr_mat_t A, gr_ctx_t ctx)
{
    slong Amax, Amin;
    slong m = A->r;
    slong n = A->c;
    slong exp, i, j;
    slong sz = ctx->sizeof_elem;
    nfloat_srcptr a;

    Amax = WORD_MIN;
    Amin = WORD_MAX;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            a = NFLOAT_COMPLEX_RE(GR_MAT_ENTRY(A, i, j, sz), ctx);
            exp = NFLOAT_EXP(a);
            Amax = FLINT_MAX(Amax, exp);
            Amin = FLINT_MIN(Amin, exp);

            a = NFLOAT_COMPLEX_IM(GR_MAT_ENTRY(A, i, j, sz), ctx);
            exp = NFLOAT_EXP(a);
            Amax = FLINT_MAX(Amax, exp);
            Amin = FLINT_MIN(Amin, exp);
        }
    }

    _Amin[0] = Amin;
    _Amax[0] = Amax;
}

int
_nfloat_complex_mat_mul_fixed_given_exp(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong Aexp, slong Bexp, slong fnlimbs, int waksman, gr_ctx_t ctx)
{
    nn_ptr T, Aa, Ab, Ba, Bb, AaBa, AbBb, Cb;
    slong i, j;
    slong sz = ctx->sizeof_elem;
    slong fdnlimbs;

    slong m = A->r;
    slong n = A->c;
    slong p = B->c;

    /* limbs including sign limb */
    fdnlimbs = fnlimbs + 1;

    T = flint_calloc(fdnlimbs * (2 * m * n + 2 * n * p + 3 * m * p), sizeof(ulong));

    Aa = T;
    Ab = Aa + fdnlimbs * (m * n);
    Ba = Ab + fdnlimbs * (m * n);
    Bb = Ba + fdnlimbs * (n * p);
    Cb = Bb + fdnlimbs * (n * p);
    AaBa = Cb + fdnlimbs * (m * p);
    AbBb = AaBa + fdnlimbs * (m * p);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            _nfloat_get_nfixed(Aa + i * fdnlimbs * n + j * fdnlimbs, NFLOAT_COMPLEX_RE(GR_MAT_ENTRY(A, i, j, sz), ctx), Aexp, fnlimbs, ctx);
            _nfloat_get_nfixed(Ab + i * fdnlimbs * n + j * fdnlimbs, NFLOAT_COMPLEX_IM(GR_MAT_ENTRY(A, i, j, sz), ctx), Aexp, fnlimbs, ctx);
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < p; j++)
        {
            _nfloat_get_nfixed(Ba + i * fdnlimbs * p + j * fdnlimbs, NFLOAT_COMPLEX_RE(GR_MAT_ENTRY(B, i, j, sz), ctx), Bexp, fnlimbs, ctx);
            _nfloat_get_nfixed(Bb + i * fdnlimbs * p + j * fdnlimbs, NFLOAT_COMPLEX_IM(GR_MAT_ENTRY(B, i, j, sz), ctx), Bexp, fnlimbs, ctx);
        }
    }

    /* (Aa+Ab i)(Ba+Bb i) = (Aa Ba - Ab Bb) + (Aa Bb + Ab Ba i) */

    /* (Aa Ba - Ab Bb) + ((Aa + Ab)(Ba + Bb) - Aa Ba - Ab Bb) i */

    if (waksman)
    {
        _nfixed_mat_mul_waksman(AaBa, Aa, Ba, m, n, p, fnlimbs);
        _nfixed_mat_mul_waksman(AbBb, Ab, Bb, m, n, p, fnlimbs);
        _nfixed_vec_add(Aa, Aa, Ab, m * n, fnlimbs);
        _nfixed_vec_add(Ba, Ba, Bb, n * p, fnlimbs);
        _nfixed_mat_mul_waksman(Cb, Aa, Ba, m, n, p, fnlimbs);
        _nfixed_vec_sub(Cb, Cb, AaBa, m * p, fnlimbs);
        _nfixed_vec_sub(Cb, Cb, AbBb, m * p, fnlimbs);

        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfloat_set_nfixed(NFLOAT_COMPLEX_IM(GR_MAT_ENTRY(C, i, j, sz), ctx), Cb + i * fdnlimbs * p + j * fdnlimbs, Aexp + Bexp, fnlimbs, ctx);

        _nfixed_vec_sub(Cb, AaBa, AbBb, m * p, fnlimbs);

        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfloat_set_nfixed(NFLOAT_COMPLEX_RE(GR_MAT_ENTRY(C, i, j, sz), ctx), Cb + i * fdnlimbs * p + j * fdnlimbs, Aexp + Bexp, fnlimbs, ctx);
    }
    else
    {
        _nfixed_mat_mul_classical(AaBa, Aa, Ba, m, n, p, fnlimbs);
        _nfixed_mat_mul_classical(AbBb, Ab, Bb, m, n, p, fnlimbs);
        _nfixed_vec_add(Aa, Aa, Ab, m * n, fnlimbs);
        _nfixed_vec_add(Ba, Ba, Bb, n * p, fnlimbs);
        _nfixed_mat_mul_classical(Cb, Aa, Ba, m, n, p, fnlimbs);
        _nfixed_vec_sub(Cb, Cb, AaBa, m * p, fnlimbs);
        _nfixed_vec_sub(Cb, Cb, AbBb, m * p, fnlimbs);

        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfloat_set_nfixed(NFLOAT_COMPLEX_IM(GR_MAT_ENTRY(C, i, j, sz), ctx), Cb + i * fdnlimbs * p + j * fdnlimbs, Aexp + Bexp, fnlimbs, ctx);

        _nfixed_vec_sub(Cb, AaBa, AbBb, m * p, fnlimbs);

        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++)
                _nfloat_set_nfixed(NFLOAT_COMPLEX_RE(GR_MAT_ENTRY(C, i, j, sz), ctx), Cb + i * fdnlimbs * p + j * fdnlimbs, Aexp + Bexp, fnlimbs, ctx);
    }

    flint_free(T);

    return GR_SUCCESS;
}

int
_nfloat_complex_mat_mul_fixed(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int waksman, slong max_extra_bits, gr_ctx_t ctx)
{
    slong Amax, Amin, Bmax, Bmin, Adelta, Bdelta, Aexp, Bexp;
    slong prec;
    slong pad_top, pad_bot, extra_bits, fbits, fnlimbs;
    slong n = A->c;

    if (NFLOAT_CTX_HAS_INF_NAN(ctx))
        return GR_UNABLE;

    prec = NFLOAT_CTX_PREC(ctx);

    _nfloat_complex_mat_exp_range(&Amin, &Amax, A, ctx);
    _nfloat_complex_mat_exp_range(&Bmin, &Bmax, B, ctx);

    if (Amax < NFLOAT_MIN_EXP || Bmax < NFLOAT_MIN_EXP)
        return gr_mat_zero(C, ctx);

    /* Currently, we don't handle zeros. (They pose no problem, but zero entries in
       the output may not be exact. To be done.) */
    if (Amin < NFLOAT_MIN_EXP || Bmin < NFLOAT_MIN_EXP)
        return gr_mat_mul_classical(C, A, B, ctx);

    Adelta = Amax - Amin;
    Bdelta = Bmax - Bmin;

    /* sanity check */
    if (Adelta > 10 * prec || Bdelta > 10 * prec)
        return gr_mat_mul_classical(C, A, B, ctx);

    /*
    To double check: for Waksman,
        * The intermediate entries are bounded by 8n max(|A|,|B|)^2.
        * The error, including error from converting
          the input matrices, is bounded by 8n ulps.
    */

    pad_top = 3 + FLINT_BIT_COUNT(n);
    pad_bot = 3 + FLINT_BIT_COUNT(n);

    extra_bits = Adelta + Bdelta + pad_top + pad_bot;

    if (extra_bits >= max_extra_bits)
        return gr_mat_mul_classical(C, A, B, ctx);

    Aexp = Amax + pad_top;
    Bexp = Bmax + pad_top;
    fbits = prec + extra_bits;
    fnlimbs = (fbits + FLINT_BITS - 1) / FLINT_BITS;

    return _nfloat_complex_mat_mul_fixed_given_exp(C, A, B, Aexp, Bexp, fnlimbs, waksman, ctx);
}

int
nfloat_complex_mat_mul_fixed_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    return _nfloat_complex_mat_mul_fixed(C, A, B, 0, 100000, ctx);
}

int
nfloat_complex_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    return _nfloat_complex_mat_mul_fixed(C, A, B, 1, 100000, ctx);
}

FLINT_FORCE_INLINE slong
_nfloat_complex_nbits(nfloat_srcptr x, slong nlimbs)
{
    nn_srcptr ad;
    slong bits;

    ad = NFLOAT_D(x);
    bits = FLINT_BITS * nlimbs;

    while (ad[0] == 0)
    {
        bits -= FLINT_BITS;
        ad++;
    }

    bits -= flint_ctz(ad[0]);

    return bits;
}

int
nfloat_complex_mat_mul_block(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong min_block_size, gr_ctx_t ctx)
{
    return nfloat_mat_mul_block(C, A, B, min_block_size, ctx);
}

static int
_nfloat_complex_mat_is_real(const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j;
    slong m = A->r;
    slong n = A->c;
    slong sz = ctx->sizeof_elem;

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            if (!NFLOAT_IS_ZERO(NFLOAT_COMPLEX_IM(GR_MAT_ENTRY(A, i, j, sz), ctx)))
                return 0;

    return 1;
}

/* check if for all entries with nonzero real and imaginary parts, the
   components don't differ too much in magnitude */
static int
_nfloat_complex_mat_parts_are_well_scaled(const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j;
    slong m = A->r;
    slong n = A->c;
    slong sz = ctx->sizeof_elem;
    nn_srcptr a, b;
    slong aexp, bexp, max;

    max = NFLOAT_CTX_PREC(ctx) / 8 + 64;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            a = NFLOAT_COMPLEX_RE(GR_MAT_ENTRY(A, i, j, sz), ctx);
            if (NFLOAT_IS_ZERO(a))
                continue;

            b = NFLOAT_COMPLEX_IM(GR_MAT_ENTRY(A, i, j, sz), ctx);
            if (NFLOAT_IS_ZERO(b))
                continue;

            aexp = NFLOAT_EXP(a);
            bexp = NFLOAT_EXP(b);

            if (FLINT_ABS(aexp - bexp) > max)
                return 0;
        }
    }

    return 1;
}

static void
_real_part(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j;
    slong m = A->r;
    slong n = A->c;
    slong ndlimbs = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            flint_mpn_copyi(((nn_ptr) res->rows[i]) + j * ndlimbs, ((nn_srcptr) A->rows[i]) + j * 2 * ndlimbs, ndlimbs);
}

static void
_imag_part(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j;
    slong m = A->r;
    slong n = A->c;
    slong ndlimbs = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            flint_mpn_copyi(((nn_ptr) res->rows[i]) + j * ndlimbs, ((nn_srcptr) A->rows[i]) + (j * 2 + 1) * ndlimbs, ndlimbs);
}

static void
_set_real_part(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j;
    slong m = A->r;
    slong n = A->c;
    slong ndlimbs = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            flint_mpn_copyi(((nn_ptr) res->rows[i]) + j * 2 * ndlimbs, ((nn_srcptr) A->rows[i]) + j * ndlimbs, ndlimbs);
}

static void
_set_imag_part(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j;
    slong m = A->r;
    slong n = A->c;
    slong ndlimbs = NFLOAT_CTX_DATA_NLIMBS(ctx);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            flint_mpn_copyi(((nn_ptr) res->rows[i]) + (j * 2 + 1) * ndlimbs, ((nn_srcptr) A->rows[i]) + j * ndlimbs, ndlimbs);
}

int
_nfloat_complex_mat_mul_reorder(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int Areal, int Breal, gr_ctx_t ctx)
{
    gr_mat_t X, Y, Z, W;
    gr_ctx_t ctx2;
    slong i, j, ndlimbs;
    int status = GR_SUCCESS;

    slong M = A->r;
    slong N = A->c;
    slong P = B->c;

    nfloat_ctx_init(ctx2, NFLOAT_CTX_PREC(ctx), NFLOAT_CTX_FLAGS(ctx));
    ndlimbs = NFLOAT_CTX_DATA_NLIMBS(ctx2);

    gr_mat_init(X, M, N, ctx2);
    gr_mat_init(Y, N, P, ctx2);
    gr_mat_init(Z, M, P, ctx2);

    if (Areal && Breal)
    {
        _real_part(X, A, ctx2);
        _real_part(Y, B, ctx2);
        status |= gr_mat_mul(Z, X, Y, ctx2);
        _set_real_part(C, Z, ctx2);

        for (i = 0; i < M; i++)
            for (j = 0; j < P; j++)
                status |= nfloat_zero(((nn_ptr) C->rows[i]) + (j * 2 + 1) * ndlimbs, ctx2);
    }
    else if (Areal)
    {
        _real_part(X, A, ctx2);
        _real_part(Y, B, ctx2);
        status |= gr_mat_mul(Z, X, Y, ctx2);
        _set_real_part(C, Z, ctx2);
        _imag_part(Y, B, ctx2);
        status |= gr_mat_mul(Z, X, Y, ctx2);
        _set_imag_part(C, Z, ctx2);
    }
    else if (Breal)
    {
        _real_part(X, A, ctx2);
        _real_part(Y, B, ctx2);
        status |= gr_mat_mul(Z, X, Y, ctx2);
        _set_real_part(C, Z, ctx2);
        _imag_part(X, A, ctx2);
        status |= gr_mat_mul(Z, X, Y, ctx2);
        _set_imag_part(C, Z, ctx2);
    }
    else
    {
        gr_mat_init(W, M, P, ctx2);

        /* Z = re(A) * re(B) */
        _real_part(X, A, ctx2);
        _real_part(Y, B, ctx2);
        status |= gr_mat_mul(Z, X, Y, ctx2);

        /* W = im(A) * im(B) */
        _imag_part(X, A, ctx2);
        _imag_part(Y, B, ctx2);
        status |= gr_mat_mul(W, X, Y, ctx2);

        if (A == C || B == C)
        {
            gr_mat_t T;
            gr_mat_init(T, M, P, ctx2);

            /* T = Z - W */
            status |= gr_mat_sub(T, Z, W, ctx2);

            /* Z = re(A) * im(B) */
            /* W = im(A) * re(B) */
            _real_part(X, A, ctx2);
            status |= gr_mat_mul(Z, X, Y, ctx2);
            _imag_part(X, A, ctx2);
            _real_part(Y, B, ctx2);
            status |= gr_mat_mul(W, X, Y, ctx2);

            /* re(C) = T */
            _set_real_part(C, T, ctx2);

            gr_mat_clear(T, ctx2);

            /* im(C) = Z + W */
            /* todo: add directly to the output */
            status |= gr_mat_add(Z, Z, W, ctx2);
            _set_imag_part(C, Z, ctx2);
        }
        else
        {
            /* re(C) = Z - W */
            status |= gr_mat_sub(Z, Z, W, ctx2);
            _set_real_part(C, Z, ctx2);

            _real_part(X, A, ctx2);
            status |= gr_mat_mul(Z, X, Y, ctx2);
            _imag_part(X, A, ctx2);
            _real_part(Y, B, ctx2);
            status |= gr_mat_mul(W, X, Y, ctx2);

            /* im(C) = Z + W */
            /* todo: add directly to the output */
            status |= gr_mat_add(Z, Z, W, ctx2);
            _set_imag_part(C, Z, ctx2);
        }

        gr_mat_clear(W, ctx2);
    }

    gr_mat_clear(X, ctx2);
    gr_mat_clear(Y, ctx2);
    gr_mat_clear(Z, ctx2);

    gr_ctx_clear(ctx2);

    return status;
}

int
nfloat_complex_mat_mul_reorder(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    return _nfloat_complex_mat_mul_reorder(C, A, B, _nfloat_complex_mat_is_real(A, ctx),
                                                    _nfloat_complex_mat_is_real(B, ctx), ctx);
}

int
nfloat_complex_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong dim;
    slong block_cutoff;
    int use_waksman = 0;
    slong prec;
    slong max_extra_prec;
    int A_real = 0, B_real = 0;

    slong m = A->r;
    slong n = A->c;
    slong p = B->c;

    dim = FLINT_MIN(n, FLINT_MIN(m, p));

    if (dim <= 2 || NFLOAT_CTX_HAS_INF_NAN(ctx))
        return gr_mat_mul_classical(C, A, B, ctx);

    if (dim >= 6)
    {
        A_real = _nfloat_complex_mat_is_real(A, ctx);
        B_real = _nfloat_complex_mat_is_real(B, ctx);

        if (A_real || B_real)
            return _nfloat_complex_mat_mul_reorder(C, A, B, A_real, B_real, ctx);
    }

    prec = NFLOAT_CTX_PREC(ctx);

    if (prec <= 256)
        block_cutoff = 80;
    else if (prec <= 512)
        block_cutoff = 160;
    else if (prec <= 3072)
        block_cutoff = 100;
    else
        block_cutoff = 80;

    if (dim < block_cutoff)
    {
        if (prec <= 128 || (prec <= 256 && n <= 4) || (prec == 576 && n <= 6))
            return gr_mat_mul_classical(C, A, B, ctx);

        if (prec <= 320)
            use_waksman = 0;
        else if (prec <= 448)
            use_waksman = (dim >= 20);
        else if (prec <= 640)
            use_waksman = (dim >= 6);
        else if (prec <= 1536)
            use_waksman = (dim >= 4);
        else
            use_waksman = (dim >= 3);

        max_extra_prec = (prec < 768) ? 64 : prec / 4;

        return _nfloat_complex_mat_mul_fixed(C, A, B, use_waksman, max_extra_prec, ctx);
    }
    else
    {
        if (_nfloat_complex_mat_parts_are_well_scaled(A, ctx) &&
            _nfloat_complex_mat_parts_are_well_scaled(B, ctx))
        {
            return nfloat_complex_mat_mul_block(C, A, B, block_cutoff - 10, ctx);
        }
        else
        {
            return _nfloat_complex_mat_mul_reorder(C, A, B, A_real, B_real, ctx);
        }
    }
}
