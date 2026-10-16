/*
    Copyright (C) 2024, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"

/*
    Classical and Karatsuba multiplication of polynomials with fixed-size
    coefficients (unsigned or two's complement), computing coefficients
    of the product mod 2^(FLINT_BITS slimbs).

    Small limb sizes dispatch to generated loop instances with the dot
    products inlined (see poly_mulmid_kernels.c); larger sizes use the
    generic loops below with flint_mpn_mul and mpn_add_n.
*/

/* scratch space: n1 + n2 limbs for the generic dot product, 2 n1 + 1 limbs
   for squares */
#define SCRATCH_LIMBS(n1, n2) (3 * ((n1) + (n2)) + 3)
#define SCRATCH_STACK_LIMBS 64

/* sdiff = s - (n1 + n2 - 1) in {0, 1, 2} should be a compile-time constant */
FLINT_FORCE_INLINE void
_dot(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, slong n1, slong n2, slong s, int sdiff, int sgn, nn_ptr scratch)
{
    if (sgn == 2)
        _flint_mpn_dot_rev_generic_signmag(r, a, n1, b, n2, len, s, scratch);
    else if (sgn)
        _flint_mpn_dot_rev_generic_signed(r, a, n1, b, n2, len, s, scratch);
    else
        _flint_mpn_dot_rev_generic_sdiff(r, a, n1, b, n2, len, sdiff, scratch);
}

/* {r, s} = x^2 mod 2^(FLINT_BITS s) */
FLINT_FORCE_INLINE void
_sqr(nn_ptr r, nn_srcptr x, slong n, slong s, int sdiff,
    int sgn, nn_ptr scratch)
{
    if (sgn == 2)
    {
        x += 1;
        sgn = 0;
    }

    if (sgn)
    {
        _dot(r, x, x, 1, n, n, s, sdiff, sgn, scratch);
    }
    else if (sdiff == 0)
    {
        nn_ptr t = scratch + 2 * n;
        flint_mpn_sqr(t, x, n);
        flint_mpn_copyi(r, t, s);
    }
    else
    {
        flint_mpn_sqr(r, x, n);
        if (sdiff == 2)
            r[2 * n] = 0;
    }
}

/* {r, s} += x^2 mod 2^(FLINT_BITS s) */
FLINT_FORCE_INLINE void
_addsqr(nn_ptr r, nn_srcptr x, slong n, slong s, int sdiff,
    int sgn, nn_ptr scratch)
{
    nn_ptr t = scratch + 2 * n;

    if (sgn == 2)
    {
        x += 1;
        sgn = 0;
    }

    if (sgn)
    {
        _dot(t, x, x, 1, n, n, s, sdiff, sgn, scratch);
        mpn_add_n(r, r, t, s);
    }
    else
    {
        flint_mpn_sqr(t, x, n);

        if (sdiff == 2)
            r[2 * n] += mpn_add_n(r, r, t, 2 * n);
        else
            mpn_add_n(r, r, t, s);
    }
}

FLINT_FORCE_INLINE void
_mul_loop(nn_ptr res, nn_srcptr f, slong flen, slong n1, nn_srcptr g, slong glen, slong n2,
    slong nlo, slong nhi, slong s, int sdiff, int sgn, nn_ptr scratch)
{
    slong i, top1, top2;
    slong st1 = n1 + (sgn == 2), st2 = n2 + (sgn == 2);

    for (i = nlo; i < nhi; i++)
    {
        top1 = FLINT_MIN(flen - 1, i);
        top2 = FLINT_MIN(glen - 1, i);

        _dot(res + (i - nlo) * s, f + (i - top2) * st1, g + (i - top1) * st2,
            top1 + top2 - i + 1, n1, n2, s, sdiff, sgn, scratch);
    }
}

FLINT_FORCE_INLINE void
_sqr_loop(nn_ptr res, nn_srcptr f, slong flen, slong n,
    slong nlo, slong nhi, slong s, int sdiff, int sgn, nn_ptr scratch)
{
    slong i, start, stop, len = 2 * flen - 1, ilo, ihi;
    slong st = n + (sgn == 2);
    nn_ptr rp;

    /* the first and last coefficients are plain squares */
    ilo = nlo;
    ihi = nhi;

    if (nlo == 0)
    {
        _sqr(res, f, n, s, sdiff, sgn, scratch);
        ilo = 1;
    }

    if (nhi == len && flen >= 2)
    {
        _sqr(res + (len - 1 - nlo) * s, f + (flen - 1) * st, n, s, sdiff, sgn, scratch);
        ihi = len - 1;
    }

    for (i = ilo; i < ihi; i++)
    {
        rp = res + (i - nlo) * s;
        start = FLINT_MAX(0, i - flen + 1);
        stop = FLINT_MIN(flen - 1, (i + 1) / 2 - 1);

        _dot(rp, f + start * st, f + (i - stop) * st, stop - start + 1, n, n, s, sdiff, sgn, scratch);
        mpn_lshift(rp, rp, s, 1);

        if (i % 2 == 0)
            _addsqr(rp, f + (i / 2) * st, n, s, sdiff, sgn, scratch);
    }
}

/* Shapes with a _short instance (assembly products + register
   accumulation), for 4 limbs where there are no column kernels (for 3
   limbs the column kernels are faster at all lengths). */

/*
    Each loop instance is a separate small function: a single function
    with all instances inlined has a very expensive prologue (large frame,
    many saved registers, spills), which is significant for the short
    base cases arising in Karatsuba multiplication.
*/

typedef struct
{
    nn_ptr res;
    nn_srcptr f;
    nn_srcptr g;
    slong flen, glen, nlo, nhi;
    /* only used by the generic instances */
    slong n1, n2, s;
    int sgn;
    nn_ptr scratch;
} args_struct;

/* Generic instances (limb sizes and sign at runtime). */
#define DEF_GENERIC_INSTANCES(SDIFF) \
    FLINT_STATIC_NOINLINE void _mul_generic_##SDIFF(const args_struct * A) \
    { \
        _mul_loop(A->res, A->f, A->flen, A->n1, A->g, A->glen, A->n2, A->nlo, A->nhi, A->s, SDIFF, A->sgn, A->scratch); \
    } \
    FLINT_STATIC_NOINLINE void _sqr_generic_##SDIFF(const args_struct * A) \
    { \
        _sqr_loop(A->res, A->f, A->flen, A->n1, A->nlo, A->nhi, A->s, SDIFF, A->sgn, A->scratch); \
    }

DEF_GENERIC_INSTANCES(0)
DEF_GENERIC_INSTANCES(1)
DEF_GENERIC_INSTANCES(2)

typedef void (* instance_func_t)(const args_struct *);

/* Tables of the generated loops with their natural signatures, indexed
   [sgn][n1][n2][sdiff], for direct dispatch with a single indirect call. */

static const instance_func_t mul_generic_tab[3] = { _mul_generic_0, _mul_generic_1, _mul_generic_2 };
static const instance_func_t sqr_generic_tab[3] = { _sqr_generic_0, _sqr_generic_1, _sqr_generic_2 };

/* The generated loops for given limb sizes and output size, if any:
   the column kernels for up to 3 limbs, the _short (assembly product)
   kernels for 4 limbs. */
FLINT_FORCE_INLINE flint_mpn_poly_mulmid_func_t
_mul_kernel(slong n1, slong n2, int sdiff, int sgn)
{
    if (sgn == 2)
        return (n1 <= 4) ? flint_mpn_poly_mulmid_signmag_tab[n1][n2][sdiff] : NULL;
    else if (n1 <= 3)
        return flint_mpn_poly_mulmid_tab[sgn != 0][n1][n2][sdiff];
    else if (n1 == 4)
        return flint_mpn_poly_mulmid_short_tab[sgn != 0][n1][n2][sdiff];
    else
        return NULL;
}

FLINT_FORCE_INLINE flint_mpn_poly_sqrmid_func_t
_sqr_kernel(slong n1, slong n2, int sdiff, int sgn)
{
    if (n1 != n2)
        return NULL;
    else if (sgn == 2)
        return (n1 <= 4) ? flint_mpn_poly_sqrmid_signmag_tab[n1][n2][sdiff] : NULL;
    else if (n1 <= 3)
        return flint_mpn_poly_sqrmid_tab[sgn != 0][n1][n2][sdiff];
    else if (n1 == 4)
        return flint_mpn_poly_sqrmid_short_tab[sgn != 0][n1][n2][sdiff];
    else
        return NULL;
}

/* Resolved loops for given limb sizes, output size and sign: the generated
   kernels with their natural signatures where available (called directly,
   which matters for the small base cases of Karatsuba), otherwise the
   generic instances, with the scratch space of SCRATCH_LIMBS(n1, n2)
   limbs. */
typedef struct
{
    flint_mpn_poly_mulmid_func_t mul;
    flint_mpn_poly_sqrmid_func_t sqr;
    instance_func_t mul_generic, sqr_generic;
    args_struct A;
} prepared_struct;

static void
_prepare(prepared_struct * P, slong n1, slong n2, slong s, int sgn, nn_ptr scratch)
{
    int sdiff = (int) (s - (n1 + n2 - 1));

    FLINT_ASSERT(n1 >= n2 && n2 >= 1);
    FLINT_ASSERT(s >= n1 + n2 - 1 && s <= n1 + n2 + 1);

    P->mul = _mul_kernel(n1, n2, sdiff, sgn);
    P->sqr = _sqr_kernel(n1, n2, sdiff, sgn);
    P->mul_generic = P->sqr_generic = NULL;

    if (P->mul == NULL)
    {
        P->mul_generic = mul_generic_tab[sdiff];
        if (n1 == n2)
            P->sqr_generic = sqr_generic_tab[sdiff];
        P->A.n1 = n1;
        P->A.n2 = n2;
        P->A.s = s;
        P->A.sgn = sgn;
        P->A.scratch = scratch;
    }
}

FLINT_FORCE_INLINE void
_prepared_mulmid(prepared_struct * P, nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlo, slong nhi)
{
    int squaring = (f == g && flen == glen);

    if (P->mul != NULL)
    {
        if (squaring && P->sqr != NULL)
            P->sqr(res, f, flen, nlo, nhi);
        else
            P->mul(res, f, flen, g, glen, nlo, nhi);
        return;
    }

    P->A.res = res;
    P->A.f = f;
    P->A.g = g;
    P->A.flen = flen;
    P->A.glen = glen;
    P->A.nlo = nlo;
    P->A.nhi = nhi;

    if (squaring && P->sqr_generic != NULL)
        P->sqr_generic(&P->A);
    else
        P->mul_generic(&P->A);
}

/* Separate function so that the fast path below has a minimal frame. */
FLINT_STATIC_NOINLINE void
_flint_mpn_poly_mulmid_classical_generic(nn_ptr res, nn_srcptr f, slong flen, slong n1, nn_srcptr g, slong glen, slong n2, slong nlo, slong nhi, slong s, int sgn)
{
    prepared_struct P;
    ulong stack_scratch[SCRATCH_STACK_LIMBS];
    nn_ptr scratch;

    if (SCRATCH_LIMBS(n1, n2) <= SCRATCH_STACK_LIMBS)
        scratch = stack_scratch;
    else
        scratch = flint_malloc(sizeof(ulong) * SCRATCH_LIMBS(n1, n2));

    _prepare(&P, n1, n2, s, sgn, scratch);
    _prepared_mulmid(&P, res, f, flen, g, glen, nlo, nhi);

    if (scratch != stack_scratch)
        flint_free(scratch);
}

void
_flint_mpn_poly_mulmid_classical(nn_ptr res, nn_srcptr f, slong flen, slong n1, nn_srcptr g, slong glen, slong n2, slong nlo, slong nhi, slong s, int sgn)
{
    int sdiff = (int) (s - (n1 + n2 - 1));

    FLINT_ASSERT(n1 >= n2 && n2 >= 1);
    FLINT_ASSERT(s >= n1 + n2 - 1 && s <= n1 + n2 + 1);
    FLINT_ASSERT(0 <= nlo && nlo < nhi && nhi <= flen + glen - 1);

    /* fast path: a single indirect call to a generated loop */
    if (n1 <= 4)
    {
        if (f == g && flen == glen)
        {
            flint_mpn_poly_sqrmid_func_t func = _sqr_kernel(n1, n2, sdiff, sgn);

            if (func != NULL)
            {
                func(res, f, flen, nlo, nhi);
                return;
            }
        }
        else
        {
            flint_mpn_poly_mulmid_func_t func = _mul_kernel(n1, n2, sdiff, sgn);

            if (func != NULL)
            {
                func(res, f, flen, g, glen, nlo, nhi);
                return;
            }
        }
    }

    _flint_mpn_poly_mulmid_classical_generic(res, f, flen, n1, g, glen, n2, nlo, nhi, s, sgn);
}

/* Karatsuba multiplication ************************************************/

/* res = f + g (coefficientwise) with nlimbs limbs, assuming that the sums
   fit in nlimbs limbs (which for two's complement values does not exclude
   carries between limb groups, so we then need to add coefficientwise) */
FLINT_FORCE_INLINE void
_flint_mpn_poly_add_n(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlimbs, int sgn)
{
    slong m = FLINT_MIN(flen, glen);
    slong i;

    if (!sgn)
        mpn_add_n(res, f, g, nlimbs * m);
    else
        for (i = 0; i < m; i++)
            mpn_add_n(res + i * nlimbs, f + i * nlimbs, g + i * nlimbs, nlimbs);

    if (flen > glen)
        flint_mpn_copyi(res + nlimbs * m, f + nlimbs * m, nlimbs * (flen - m));
    if (glen > flen)
        flint_mpn_copyi(res + nlimbs * m, g + nlimbs * m, nlimbs * (glen - m));
}

/* res = f + g where the inputs have nlimbs limbs and the output nlimbs + 1
   limbs (zero- or sign-extended) */
static void
_flint_mpn_poly_add_n_extend(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlimbs, int sgn)
{
    slong m = FLINT_MIN(flen, glen);
    slong i;
    slong nlimbs2 = nlimbs + 1;

    /* inline additions for small limb counts (a function call per
       coefficient is a significant part of a small Karatsuba product);
       the sign-extended top limbs make the two's complement case uniform */
    if (nlimbs == 1)
    {
        for (i = 0; i < m; i++)
        {
            ulong sf = sgn ? FLINT_SIGN_EXT(f[i]) : 0;
            ulong sg = sgn ? FLINT_SIGN_EXT(g[i]) : 0;
            add_ssaaaa(res[2 * i + 1], res[2 * i], sf, f[i], sg, g[i]);
        }
    }
    else if (nlimbs == 2)
    {
        for (i = 0; i < m; i++)
        {
            ulong sf = sgn ? FLINT_SIGN_EXT(f[2 * i + 1]) : 0;
            ulong sg = sgn ? FLINT_SIGN_EXT(g[2 * i + 1]) : 0;
            add_sssaaaaaa(res[3 * i + 2], res[3 * i + 1], res[3 * i],
                          sf, f[2 * i + 1], f[2 * i], sg, g[2 * i + 1], g[2 * i]);
        }
    }
    else if (nlimbs == 3)
    {
        for (i = 0; i < m; i++)
        {
            ulong sf = sgn ? FLINT_SIGN_EXT(f[3 * i + 2]) : 0;
            ulong sg = sgn ? FLINT_SIGN_EXT(g[3 * i + 2]) : 0;
            add_ssssaaaaaaaa(res[4 * i + 3], res[4 * i + 2], res[4 * i + 1], res[4 * i],
                             sf, f[3 * i + 2], f[3 * i + 1], f[3 * i],
                             sg, g[3 * i + 2], g[3 * i + 1], g[3 * i]);
        }
    }
    else
    {
        for (i = 0; i < m; i++)
        {
            ulong cy = mpn_add_n(res + i * nlimbs2, f + i * nlimbs, g + i * nlimbs, nlimbs);
            if (sgn)
                cy += FLINT_SIGN_EXT(f[(i + 1) * nlimbs - 1]) + FLINT_SIGN_EXT(g[(i + 1) * nlimbs - 1]);
            res[(i + 1) * nlimbs2 - 1] = cy;
        }
    }

    for (i = m; i < flen; i++)
    {
        flint_mpn_copyi(res + i * nlimbs2, f + i * nlimbs, nlimbs);
        res[(i + 1) * nlimbs2 - 1] = sgn ? FLINT_SIGN_EXT(f[(i + 1) * nlimbs - 1]) : 0;
    }

    for (i = m; i < glen; i++)
    {
        flint_mpn_copyi(res + i * nlimbs2, g + i * nlimbs, nlimbs);
        res[(i + 1) * nlimbs2 - 1] = sgn ? FLINT_SIGN_EXT(g[(i + 1) * nlimbs - 1]) : 0;
    }
}

/* res = f - g (coefficientwise, mod 2^(FLINT_BITS nlimbs)); the result is
   meaningful as a two's complement value when it fits in nlimbs limbs */
static void
_flint_mpn_poly_sub_n(nn_ptr res, nn_srcptr f, nn_srcptr g, slong len, slong nlimbs)
{
    slong i;
    for (i = 0; i < len; i++)
        mpn_sub_n(res + i * nlimbs, f + i * nlimbs, g + i * nlimbs, nlimbs);
}

/* res = f - g where the inputs have nlimbs limbs and the output nlimbs + 1
   limbs (two's complement, sign-extended) */
static void
_flint_mpn_poly_sub_n_extend(nn_ptr res, nn_srcptr f, nn_srcptr g, slong len, slong nlimbs, int sgn)
{
    slong i;
    slong nlimbs2 = nlimbs + 1;

    for (i = 0; i < len; i++)
    {
        ulong cy = mpn_sub_n(res + i * nlimbs2, f + i * nlimbs, g + i * nlimbs, nlimbs);
        ulong top = -cy;
        if (sgn)
            top += FLINT_SIGN_EXT(f[(i + 1) * nlimbs - 1]) - FLINT_SIGN_EXT(g[(i + 1) * nlimbs - 1]);
        res[(i + 1) * nlimbs2 - 1] = top;
    }
}

/* res = f extended from nlimbs to nlimbs + 1 limbs (zero- or sign-extended) */
static void
_flint_mpn_poly_copy_extend(nn_ptr res, nn_srcptr f, slong len, slong nlimbs, int sgn)
{
    slong i;
    slong nlimbs2 = nlimbs + 1;

    for (i = 0; i < len; i++)
    {
        flint_mpn_copyi(res + i * nlimbs2, f + i * nlimbs, nlimbs);
        res[(i + 1) * nlimbs2 - 1] = sgn ? FLINT_SIGN_EXT(f[(i + 1) * nlimbs - 1]) : 0;
    }
}

/* coefficientwise res -= x mod 2^(FLINT_BITS slimbs) */
FLINT_FORCE_INLINE void
_flint_mpn_poly_sub_mod(nn_ptr res, nn_srcptr x, slong len, slong slimbs, int sgn)
{
    slong i;

    /* unsigned values never overflow, so we can do a single subtraction */
    if (!sgn)
        mpn_sub_n(res, res, x, len * slimbs);
    else
        for (i = 0; i < len; i++)
            mpn_sub_n(res + i * slimbs, res + i * slimbs, x + i * slimbs, slimbs);
}

FLINT_FORCE_INLINE void
_flint_mpn_poly_add_mod(nn_ptr res, nn_srcptr x, slong len, slong slimbs, int sgn)
{
    slong i;

    if (!sgn)
        mpn_add_n(res, res, x, len * slimbs);
    else
        for (i = 0; i < len; i++)
            mpn_add_n(res + i * slimbs, res + i * slimbs, x + i * slimbs, slimbs);
}

/*
    Operands with per-operand limb counts and headroom. The limb count of
    an operand is either its base value or the base value plus one (after
    an extension), and the headroom (norm) is the number of leading unused
    bits in the top limb, excluding the sign bit in two's complement mode.
*/
typedef struct
{
    nn_srcptr x;
    slong len;
    slong nlimbs;
    slong norm;
    slong base;         /* limb count before any extension */
} operand_struct;

/*
    Context for the Karatsuba-based routines: base limb counts for the
    two operands, and classical instances (resolved lazily) for the
    combinations of sign mode and extension of either operand.
*/
typedef struct
{
    prepared_struct P[2][2][2];     /* [sgn][f extended][g extended] */
    slong nlimbs1, nlimbs2, slimbs, cutoff;
    slong short_cutoff;     /* for Mulders' short products and the middle product */
    nn_ptr scratch;
} kar_ctx_struct;

#define SHORT_CUTOFF_FACTOR 4

/* Short products (Mulders, middle product) cost roughly 0.75-1.5 times a
   full Karatsuba product while the classical versions cost about half of
   a full classical product, so their crossovers are around 4 times the
   Karatsuba cutoff. */

#define KAR_SCRATCH_LIMBS(n1, n2) SCRATCH_LIMBS((n1) + 1, (n2) + 1)

static void
_kar_ctx_init(kar_ctx_struct * K, slong nlimbs1, slong nlimbs2, slong slimbs, slong cutoff, nn_ptr scratch)
{
    int i, j, k;

    K->nlimbs1 = nlimbs1;
    K->nlimbs2 = nlimbs2;
    K->slimbs = slimbs;
    K->cutoff = FLINT_MAX(cutoff, 2);
    K->short_cutoff = SHORT_CUTOFF_FACTOR * K->cutoff;
    K->scratch = scratch;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
            {
                K->P[i][j][k].mul = NULL;
                K->P[i][j][k].mul_generic = NULL;
                K->P[i][j][k].A.scratch = scratch;
            }
}

/* Window [nlo, nhi) of the product of the operands F and G by the classical
   algorithm, ordering the operands as required by the instances. */
static void
_ctx_mulmid(kar_ctx_struct * K, nn_ptr res, const operand_struct * F, const operand_struct * G, slong nlo, slong nhi, int sgn)
{
    prepared_struct * P;
    int ef = (int) (F->nlimbs - F->base);
    int eg = (int) (G->nlimbs - G->base);
    int e1, e2;

    FLINT_ASSERT(ef == 0 || ef == 1);
    FLINT_ASSERT(eg == 0 || eg == 1);

    /* the instances are indexed by the extensions of the operands with base
       limb counts nlimbs1 and nlimbs2 (which may appear in either order) */
    if (F->base == K->nlimbs1 && G->base == K->nlimbs2)
    {
        e1 = ef;
        e2 = eg;
        if (K->nlimbs1 == K->nlimbs2 && e1 < e2)
            FLINT_SWAP(int, e1, e2);
    }
    else
    {
        FLINT_ASSERT(F->base == K->nlimbs2 && G->base == K->nlimbs1);
        e1 = eg;
        e2 = ef;
    }

    P = &K->P[sgn != 0][e1][e2];

    if (P->mul == NULL && P->mul_generic == NULL)
        _prepare(P, FLINT_MAX(F->nlimbs, G->nlimbs), FLINT_MIN(F->nlimbs, G->nlimbs), K->slimbs, sgn, K->scratch);

    if (F->nlimbs >= G->nlimbs)
        _prepared_mulmid(P, res, F->x, F->len, G->x, G->len, nlo, nhi);
    else
        _prepared_mulmid(P, res, G->x, G->len, F->x, F->len, nlo, nhi);
}

FLINT_FORCE_INLINE operand_struct
_op(nn_srcptr x, slong len, const operand_struct * like)
{
    operand_struct r;
    r.x = x;
    r.len = len;
    r.nlimbs = like->nlimbs;
    r.norm = like->norm;
    r.base = like->base;
    return r;
}

FLINT_FORCE_INLINE operand_struct
_op_part(const operand_struct * F, slong start, slong len)
{
    return _op(F->x + start * F->nlimbs, len, F);
}

/* T = F0 + F1 (parts of an operand), extending by a limb when there is
   no headroom left. The result is written to t (space for
   max(len0, len1) * (nlimbs + 1) limbs). */
static operand_struct
_op_sum(nn_ptr t, const operand_struct * F0, const operand_struct * F1, int sgn)
{
    operand_struct T;

    T.x = t;
    T.len = FLINT_MAX(F0->len, F1->len);
    T.base = F0->base;

    if (F0->norm == 0)
    {
        _flint_mpn_poly_add_n_extend(t, F0->x, F0->len, F1->x, F1->len, F0->nlimbs, sgn);
        T.nlimbs = F0->nlimbs + 1;
        T.norm = FLINT_BITS - 1;
    }
    else
    {
        _flint_mpn_poly_add_n(t, F0->x, F0->len, F1->x, F1->len, F0->nlimbs, sgn);
        T.nlimbs = F0->nlimbs;
        T.norm = F0->norm - 1;
    }

    return T;
}

/*
    Full product of F and G by Karatsuba multiplication. The parts and
    their sums are handled per operand, so unbalanced limb counts (and
    extensions of only one operand) are supported.
*/
static void
_karatsuba(nn_ptr res, const operand_struct * F, const operand_struct * G, int sgn, kar_ctx_struct * K)
{
    slong slimbs = K->slimbs;
    slong flen = F->len, glen = G->len;
    slong m, f1len, g1len, tlen, ulen, vlen;
    operand_struct F0, F1, G0, G1, T, U;
    nn_ptr t, u, v;
    int squaring = (F->x == G->x) && (flen == glen) && (F->nlimbs == G->nlimbs);
    TMP_INIT;

    if (FLINT_MIN(flen, glen) < K->cutoff)
    {
        _ctx_mulmid(K, res, F, G, 0, flen + glen - 1, sgn);
        return;
    }

    /* Unbalanced lengths: chunk the longer operand into pieces of the length
       of the shorter one (plain splitting would double the imbalance at each
       level of recursion). */
    if (flen > glen + glen / 2 || glen > flen + flen / 2)
    {
        const operand_struct * L = F, * S = G;
        slong i, m2, slen;
        operand_struct piece;

        if (flen < glen)
        {
            L = G;
            S = F;
        }

        slen = S->len;

        TMP_START;
        t = TMP_ALLOC(sizeof(ulong) * (2 * slen - 1) * slimbs);

        piece = _op_part(L, 0, slen);
        _karatsuba(res, &piece, S, sgn, K);

        for (i = slen; i < L->len; i += slen)
        {
            m2 = FLINT_MIN(slen, L->len - i);
            piece = _op_part(L, i, m2);
            _karatsuba(t, &piece, S, sgn, K);
            /* the first slen - 1 coefficients overlap with the previous chunk */
            _flint_mpn_poly_add_mod(res + i * slimbs, t, slen - 1, slimbs, sgn);
            flint_mpn_copyi(res + (i + slen - 1) * slimbs, t + (slen - 1) * slimbs, m2 * slimbs);
        }

        TMP_END;
        return;
    }

    /* split at X = x^m */
    /* res = f0 g0 + (f0 g1 + f1 g0) X + f1 g1 X^2
           = f0 g0 + ((f0 + f1) (g0 + g1) - f0 g0 - f1 g1) X + f1 g1 X^2 */
    m = (FLINT_MIN(flen, glen) + 1) / 2;
    f1len = flen - m;
    g1len = glen - m;

    F0 = _op_part(F, 0, m);
    F1 = _op_part(F, m, f1len);
    G0 = _op_part(G, 0, m);
    G1 = _op_part(G, m, g1len);

    /* low part: res[0, ..., 2m-2] = f0 g0 */
    _karatsuba(res, &F0, &G0, sgn, K);

    /* res[2m-1] = 0 */
    flint_mpn_zero(res + (2 * m - 1) * slimbs, slimbs);

    /* high part: res[2m, ..., flen+glen-2] = f1 g1 */
    _karatsuba(res + (2 * m) * slimbs, &F1, &G1, sgn, K);

    /* temporary space for the middle part */
    tlen = FLINT_MAX(m, f1len);
    ulen = FLINT_MAX(m, g1len);
    vlen = tlen + ulen - 1;

    TMP_START;
    t = TMP_ALLOC(sizeof(ulong) * (tlen * (F->nlimbs + 1) + ulen * (G->nlimbs + 1) + vlen * slimbs));
    u = t + tlen * (F->nlimbs + 1);
    v = u + ulen * (G->nlimbs + 1);

    T = _op_sum(t, &F0, &F1, sgn);

    if (squaring)
    {
        _karatsuba(v, &T, &T, sgn, K);
    }
    else
    {
        U = _op_sum(u, &G0, &G1, sgn);
        _karatsuba(v, &T, &U, sgn, K);
    }

    /* v -= f0 g0 */
    _flint_mpn_poly_sub_mod(v, res, 2 * m - 1, slimbs, sgn);
    /* v -= f1 g1 */
    _flint_mpn_poly_sub_mod(v, res + 2 * m * slimbs, f1len + g1len - 1, slimbs, sgn);
    /* res += v X^m */
    _flint_mpn_poly_add_mod(res + m * slimbs, v, vlen, slimbs, sgn);

    TMP_END;
}

/* Mulders' splitting point for short products of length n (the measured
   optimum wobbles between 0.6 and 0.7, consistent with the theoretical
   0.694 for the Karatsuba exponent) */
#define MULDERS_SPLIT(n) (((n) * 7) / 10)

/*
    Low product: coefficients [0, n) of F G, assuming flen, glen <= n.
    Uses Mulders' algorithm: a full product of the low parts f0 g0 (split
    at k > n/2), plus recursive low products of the cross terms (which
    coincide when squaring).
*/
static void
_mullow(nn_ptr res, const operand_struct * F, const operand_struct * G, slong n, int sgn, kar_ctx_struct * K)
{
    slong slimbs = K->slimbs;
    slong flen = F->len, glen = G->len;
    slong L = flen + glen - 1;
    slong k, k1, k2, m, tlen;
    operand_struct F0, F1, G0, G1;
    nn_ptr t;
    int squaring = (F->x == G->x) && (flen == glen) && (F->nlimbs == G->nlimbs);
    TMP_INIT;

    FLINT_ASSERT(flen <= n && glen <= n && flen >= 1 && glen >= 1);

    if (L <= n)
    {
        _karatsuba(res, F, G, sgn, K);
        flint_mpn_zero(res + L * slimbs, (n - L) * slimbs);
        return;
    }

    if (n < K->short_cutoff || FLINT_MIN(flen, glen) < K->cutoff)
    {
        _ctx_mulmid(K, res, F, G, 0, n, sgn);
        return;
    }

    TMP_START;

    /* nearly full: compute the full product */
    if (L - n <= n / 8)
    {
        t = TMP_ALLOC(sizeof(ulong) * L * slimbs);
        _karatsuba(t, F, G, sgn, K);
        flint_mpn_copyi(res, t, n * slimbs);
        TMP_END;
        return;
    }

    k = MULDERS_SPLIT(n);
    k = FLINT_MAX(k, (n + 1) / 2);
    k1 = FLINT_MIN(k, flen);
    k2 = FLINT_MIN(k, glen);
    m = n - k;

    F0 = _op_part(F, 0, k1);
    G0 = _op_part(G, 0, k2);

    /* f0 g0, keeping n coefficients */
    tlen = k1 + k2 - 1;
    if (tlen <= n)
    {
        _karatsuba(res, &F0, &G0, sgn, K);
        flint_mpn_zero(res + tlen * slimbs, (n - tlen) * slimbs);
    }
    else
    {
        t = TMP_ALLOC(sizeof(ulong) * tlen * slimbs);
        _karatsuba(t, &F0, &G0, sgn, K);
        flint_mpn_copyi(res, t, n * slimbs);
    }

    t = TMP_ALLOC(sizeof(ulong) * m * slimbs);

    /* x^k f1 g0, low part: g0 truncated to m coefficients */
    if (flen > k)
    {
        F1 = _op_part(F, k, flen - k);
        G0 = _op_part(G, 0, FLINT_MIN(glen, m));
        _mullow(t, &F1, &G0, m, sgn, K);
        _flint_mpn_poly_add_mod(res + k * slimbs, t, m, slimbs, sgn);
        /* when squaring, the other cross term is the same */
        if (squaring)
            _flint_mpn_poly_add_mod(res + k * slimbs, t, m, slimbs, sgn);
    }

    /* x^k f0 g1, low part */
    if (glen > k && !squaring)
    {
        F0 = _op_part(F, 0, FLINT_MIN(flen, m));
        G1 = _op_part(G, k, glen - k);
        _mullow(t, &F0, &G1, m, sgn, K);
        _flint_mpn_poly_add_mod(res + k * slimbs, t, m, slimbs, sgn);
    }

    TMP_END;
}

FLINT_FORCE_INLINE void
_reverse_coeffs(nn_ptr res, nn_srcptr f, slong len, slong nlimbs)
{
    slong i;
    for (i = 0; i < len; i++)
        flint_mpn_copyi(res + i * nlimbs, f + (len - 1 - i) * nlimbs, nlimbs);
}

/* High product: coefficients [nlo, flen + glen - 1) of F G, via reversal. */
static void
_mulhigh(nn_ptr res, const operand_struct * F, const operand_struct * G, slong nlo, int sgn, kar_ctx_struct * K)
{
    slong slimbs = K->slimbs;
    slong flen = F->len, glen = G->len;
    slong L = flen + glen - 1;
    slong n = L - nlo;
    operand_struct FR, GR;
    nn_ptr fr, gr, t;
    int squaring = (F->x == G->x) && (flen == glen) && (F->nlimbs == G->nlimbs);
    TMP_INIT;

    TMP_START;
    fr = TMP_ALLOC(sizeof(ulong) * (flen * F->nlimbs + (squaring ? 0 : glen * G->nlimbs)));
    t = TMP_ALLOC(sizeof(ulong) * n * slimbs);

    _reverse_coeffs(fr, F->x, flen, F->nlimbs);
    FR = _op(fr, FLINT_MIN(flen, n), F);

    if (squaring)
    {
        GR = FR;
    }
    else
    {
        gr = fr + flen * F->nlimbs;
        _reverse_coeffs(gr, G->x, glen, G->nlimbs);
        GR = _op(gr, FLINT_MIN(glen, n), G);
    }

    _mullow(t, &FR, &GR, n, sgn, K);
    _reverse_coeffs(res, t, n, slimbs);

    TMP_END;
}

/* extended copy of an operand (used when switching the middle product to
   two's complement without headroom) */
static operand_struct
_op_extend(nn_ptr t, const operand_struct * F, int sgn)
{
    operand_struct T;
    _flint_mpn_poly_copy_extend(t, F->x, F->len, F->nlimbs, sgn);
    T.x = t;
    T.len = F->len;
    T.nlimbs = F->nlimbs + 1;
    T.norm = FLINT_BITS - 1;
    T.base = F->base;
    return T;
}

/* T = F0 - F1 in two's complement, extending when there is no headroom */
static operand_struct
_op_diff(nn_ptr t, const operand_struct * F0, const operand_struct * F1, slong len)
{
    operand_struct T;

    T.x = t;
    T.len = len;
    T.base = F0->base;

    if (F0->norm == 0)
    {
        _flint_mpn_poly_sub_n_extend(t, F0->x, F1->x, len, F0->nlimbs, 1);
        T.nlimbs = F0->nlimbs + 1;
        T.norm = FLINT_BITS - 1;
    }
    else
    {
        _flint_mpn_poly_sub_n(t, F0->x, F1->x, len, F0->nlimbs);
        T.nlimbs = F0->nlimbs;
        T.norm = F0->norm - 1;
    }

    return T;
}

/*
    Karatsuba middle product (Hanrot, Quercia, Zimmermann): sets res to
    the n coefficients [n - 1, 2n - 1) of A B where A has 2n - 1 and B
    has n coefficients.

    The recursion involves differences of parts of A, so it switches to
    two's complement arithmetic: if the inputs are unsigned, they are
    reinterpreted as signed values when the top limb has a spare bit,
    and extended by a limb otherwise.
*/
static void
_mulmid_kar(nn_ptr res, const operand_struct * A, const operand_struct * B, slong n, int sgn, kar_ctx_struct * K)
{
    slong slimbs = K->slimbs;
    slong k;
    operand_struct A0, A1, A2, B0, B1, T1, T2, T3, AE, BE;
    nn_ptr t, ta, tb;
    TMP_INIT;

    FLINT_ASSERT(A->len == 2 * n - 1 && B->len == n);

    if (n < K->cutoff)
    {
        _ctx_mulmid(K, res, A, B, n - 1, 2 * n - 1, sgn);
        return;
    }

    TMP_START;

    if (n % 2 == 1)
    {
        /* peel off the last coefficient of b and the last output */
        operand_struct A2, Bt, Al;
        A2 = _op_part(A, 1, 2 * n - 3);
        Bt = _op_part(B, 0, n - 1);
        t = TMP_ALLOC(sizeof(ulong) * (n - 1) * slimbs);
        _mulmid_kar(res, &A2, &Bt, n - 1, sgn, K);
        Al = _op_part(A, 0, n - 1);
        Bt = _op_part(B, n - 1, 1);
        _ctx_mulmid(K, t, &Al, &Bt, 0, n - 1, sgn);
        _flint_mpn_poly_add_mod(res, t, n - 1, slimbs, sgn);
        _ctx_mulmid(K, res + (n - 1) * slimbs, A, B, 2 * n - 2, 2 * n - 1, sgn);
        TMP_END;
        return;
    }

    k = n / 2;

    /* switch to two's complement */
    if (!sgn)
    {
        int ea = (A->norm == 0), eb = (B->norm == 0);

        if (ea || eb)
        {
            ta = TMP_ALLOC(sizeof(ulong) * ((ea ? (2 * n - 1) * (A->nlimbs + 1) : 0) + (eb ? n * (B->nlimbs + 1) : 0)));
            tb = ta + (ea ? (2 * n - 1) * (A->nlimbs + 1) : 0);
            AE = ea ? _op_extend(ta, A, 0) : *A;
            BE = eb ? _op_extend(tb, B, 0) : *B;
            if (!ea) AE.norm--;
            if (!eb) BE.norm--;
            _mulmid_kar(res, &AE, &BE, n, 1, K);
            TMP_END;
            return;
        }

        AE = *A;
        BE = *B;
        AE.norm--;
        BE.norm--;
        A = &AE;
        B = &BE;
        sgn = 1;
    }

    A0 = _op_part(A, 0, 2 * k - 1);
    A1 = _op_part(A, k, 2 * k - 1);
    A2 = _op_part(A, 2 * k, 2 * k - 1);
    B0 = _op_part(B, 0, k);
    B1 = _op_part(B, k, k);

    t = TMP_ALLOC(sizeof(ulong) * (k * (B->nlimbs + 1) + 2 * (2 * k - 1) * (A->nlimbs + 1) + k * slimbs));
    T1 = _op_sum(t, &B0, &B1, 1);
    T2 = _op_diff(t + k * (B->nlimbs + 1), &A0, &A1, 2 * k - 1);
    T3 = _op_diff(t + k * (B->nlimbs + 1) + (2 * k - 1) * (A->nlimbs + 1), &A2, &A1, 2 * k - 1);
    t = t + k * (B->nlimbs + 1) + 2 * (2 * k - 1) * (A->nlimbs + 1);

    /* the untouched parts need the same headroom bookkeeping (they are
       unchanged, so their headroom is not reduced, but the instance
       selection depends only on the limb counts) */

    /* P1 = MP(a1, b0 + b1) */
    _mulmid_kar(t, &A1, &T1, k, 1, K);
    /* res_lo = P1 + MP(a0 - a1, b1) */
    _mulmid_kar(res, &T2, &B1, k, 1, K);
    _flint_mpn_poly_add_mod(res, t, k, slimbs, 1);
    /* res_hi = P1 + MP(a2 - a1, b0) */
    _mulmid_kar(res + k * slimbs, &T3, &B0, k, 1, K);
    _flint_mpn_poly_add_mod(res + k * slimbs, t, k, slimbs, 1);

    TMP_END;
}

void
_flint_mpn_poly_mulmid_karatsuba(nn_ptr res, nn_srcptr f, slong flen, slong nlimbs1, slong norm1, nn_srcptr g, slong glen, slong nlimbs2, slong norm2, slong nlo, slong nhi, slong slimbs, slong cutoff, int sgn)
{
    kar_ctx_struct K;
    operand_struct F, G;
    slong L, len, n, mlo, mhi, i;
    TMP_INIT;

    FLINT_ASSERT(flen >= 1 && glen >= 1);
    FLINT_ASSERT(0 <= nlo && nlo < nhi && nhi <= flen + glen - 1);

    /* trim the inputs to the coefficients that contribute to the window */
    flen = FLINT_MIN(flen, nhi);
    glen = FLINT_MIN(glen, nhi);
    {
        slong nlo2 = (flen + glen - 1) - nlo;

        if (flen > nlo2)
        {
            slong trunc = flen - nlo2;
            f += trunc * nlimbs1;
            flen -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }

        if (glen > nlo2)
        {
            slong trunc = glen - nlo2;
            g += trunc * nlimbs2;
            glen -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }
    }

    L = flen + glen - 1;
    len = nhi - nlo;

    F.x = f; F.len = flen; F.nlimbs = nlimbs1; F.norm = norm1; F.base = nlimbs1;
    G.x = g; G.len = glen; G.nlimbs = nlimbs2; G.norm = norm2; G.base = nlimbs2;

    if (flen < glen)
        FLINT_SWAP(operand_struct, F, G);

    TMP_START;
    _kar_ctx_init(&K, F.nlimbs, G.nlimbs, slimbs, cutoff, TMP_ALLOC(sizeof(ulong) * KAR_SCRATCH_LIMBS(nlimbs1, nlimbs2)));

    if (nlo == 0 && nhi == L)
    {
        _karatsuba(res, &F, &G, sgn, &K);
    }
    else if (G.len < K.cutoff || len < K.short_cutoff)
    {
        _ctx_mulmid(&K, res, &F, &G, nlo, nhi, sgn);
    }
    else if (nlo == 0)
    {
        _mullow(res, &F, &G, nhi, sgn, &K);
    }
    else if (nhi == L)
    {
        _mulhigh(res, &F, &G, nlo, sgn, &K);
    }
    else
    {
        /* Middle product region: coefficients [n - 1, flen) of f g (n = glen)
           involve all n coefficients of g and can be computed in blocks of
           n by the Karatsuba middle product. The edges are done classically.
           The middle product consumes a bit of headroom on every branch at
           each level of recursion; when the input limbs are too tight for
           that, the coefficients get extended early, which inflates the
           kernel costs by about ((m+1)/m)^2 for m-limb coefficients and
           moves the crossover with the classical algorithm accordingly
           (measured: about 4x for 1 limb, 3x for 2 limbs, 2x for 3 and 4,
           none from 5 limbs). */
        slong mp_cutoff = K.short_cutoff;

        n = G.len;
        mlo = FLINT_MAX(nlo, n - 1);
        mhi = FLINT_MIN(nhi, F.len);

        if (FLINT_MIN(F.norm, G.norm) < (slong) FLINT_BIT_COUNT(n / K.cutoff))
        {
            slong nmin = FLINT_MIN(F.nlimbs, G.nlimbs);
            mp_cutoff *= (nmin == 1) ? 4 : (nmin == 2) ? 3 : (nmin <= 4) ? 2 : 1;
        }

        /* A block of n outputs starting at i needs the coefficients
           [i - (n - 1), i + n) of f. The last block may extend past the
           window (it is then truncated); if it also extends past f
           (whose coefficients above the window have been trimmed), f is
           zero-padded, since the missing coefficients only contribute to
           outputs beyond the window. */
        if (mhi - mlo < n / 2)
        {
            /* not decomposable into middle-product blocks:
               full product or classical */
            if (len >= L / 4)
            {
                nn_ptr t = TMP_ALLOC(sizeof(ulong) * L * slimbs);
                _karatsuba(t, &F, &G, sgn, &K);
                flint_mpn_copyi(res, t + nlo * slimbs, len * slimbs);
            }
            else
            {
                _ctx_mulmid(&K, res, &F, &G, nlo, nhi, sgn);
            }
        }
        else if (n < mp_cutoff)
        {
            /* decomposable, but below the middle product crossover the
               classical window wins (also against the full product,
               whose Karatsuba sums consume headroom in the same way) */
            _ctx_mulmid(&K, res, &F, &G, nlo, nhi, sgn);
        }
        else
        {
            if (nlo < mlo)
                _ctx_mulmid(&K, res, &F, &G, nlo, mlo, sgn);
            if (nhi > mhi)
                _ctx_mulmid(&K, res + (mhi - nlo) * slimbs, &F, &G, mhi, nhi, sgn);

            for (i = mlo; i < mhi; i += n)
            {
                if (mhi - i >= n)
                {
                    operand_struct Ablock = _op_part(&F, i - (n - 1), 2 * n - 1);
                    _mulmid_kar(res + (i - nlo) * slimbs, &Ablock, &G, n, sgn, &K);
                }
                else
                {
                    /* last partial block: compute the full block starting
                       at i (zero-padding f if needed) and keep the outputs
                       in the window */
                    slong w = mhi - i;
                    slong avail = F.len - (i - (n - 1));
                    slong st = F.nlimbs + (sgn == 2);
                    operand_struct Ablock;
                    nn_ptr t = TMP_ALLOC(sizeof(ulong) * n * slimbs);

                    if (avail >= 2 * n - 1)
                    {
                        Ablock = _op_part(&F, i - (n - 1), 2 * n - 1);
                    }
                    else
                    {
                        nn_ptr fa = TMP_ALLOC(sizeof(ulong) * (2 * n - 1) * st);
                        flint_mpn_copyi(fa, F.x + (i - (n - 1)) * st, avail * st);
                        flint_mpn_zero(fa + avail * st, (2 * n - 1 - avail) * st);
                        Ablock = _op(fa, 2 * n - 1, &F);
                    }

                    _mulmid_kar(t, &Ablock, &G, n, sgn, &K);
                    flint_mpn_copyi(res + (i - nlo) * slimbs, t, w * slimbs);
                }
            }
        }
    }

    TMP_END;
}

void
_flint_mpn_poly_mul_karatsuba(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlimbs, slong slimbs, slong cutoff, int norm, int sgn)
{
    _flint_mpn_poly_mulmid_karatsuba(res, f, flen, nlimbs, norm, g, glen, nlimbs, norm, 0, flen + glen - 1, slimbs, cutoff, sgn);
}

/* for tuning: -1 = automatic */
int _flint_mpn_poly_mulmid_force_method = -1;
slong _flint_mpn_poly_mulmid_force_cutoff = -1;

slong
_flint_mpn_poly_mulmid_cutoff(slong nlimbs1, slong nlimbs2, int squaring)
{
    if (_flint_mpn_poly_mulmid_force_cutoff != -1)
        return _flint_mpn_poly_mulmid_force_cutoff;
    else
        return _flint_mpn_poly_karatsuba_cutoff(nlimbs1, nlimbs2, squaring);
}

/* Lengths of the operands after discarding the coefficients that cannot
   contribute to the window [nlo, nhi), as done by the multiplication
   routines; the shorter one bounds the number of terms in the dot
   products and determines the choice of algorithm. */
static void
_flint_mpn_poly_mulmid_trunc(slong * len1, slong * len2, slong nlo, slong nhi)
{
    slong l1 = FLINT_MIN(*len1, nhi);
    slong l2 = FLINT_MIN(*len2, nhi);
    slong nlo2 = (l1 + l2 - 1) - nlo;

    if (l1 > nlo2)
        l1 = nlo2;
    if (l2 > nlo2)
        l2 = nlo2;

    *len1 = l1;
    *len2 = l2;
}

void
_flint_mpn_poly_mulmid_params(mpn_poly_mul_params_t P, slong len1, slong bits1, slong len2, slong bits2, slong nlo, slong nhi, int squaring)
{
    slong minlen, ub1, ub2, eb1, eb2, n1, n2;
    int sgn1, sgn2, method;
    int force = _flint_mpn_poly_mulmid_force_method;

    _flint_mpn_poly_mulmid_trunc(&len1, &len2, nlo, nhi);
    minlen = FLINT_MIN(len1, len2);

    sgn1 = bits1 < 0;
    sgn2 = bits2 < 0;
    ub1 = FLINT_MAX(FLINT_ABS(bits1), 1);
    ub2 = FLINT_MAX(FLINT_ABS(bits2), 1);

    /* Sign handling. The two's complement kernels are 15-20% slower than
       the unsigned ones, while the bias method costs O(s) per output
       coefficient, so the latter wins for longer dot products, and
       clearly with Karatsuba where the correction is amortized. Both need
       a sign bit; when that costs a limb (coefficients of exactly 64k
       bits), sign-magnitude coefficients with conditionally negated
       products are much cheaper, but are only available with the
       classical algorithm, which then remains competitive with
       Karatsuba (using the larger limb counts) up to several times the
       usual cutoff. */
    if (!sgn1 && !sgn2)
    {
        method = FLINT_MPN_POLY_MUL_UNSIGNED;
    }
    else
    {
        slong n1b = (ub1 + sgn1 + FLINT_BITS - 1) / FLINT_BITS;
        slong n2b = (ub2 + sgn2 + FLINT_BITS - 1) / FLINT_BITS;
        slong n1s = (ub1 + 1 + FLINT_BITS - 1) / FLINT_BITS;
        slong n2s = (ub2 + 1 + FLINT_BITS - 1) / FLINT_BITS;
        slong n1m = (ub1 + FLINT_BITS - 1) / FLINT_BITS;
        slong n2m = (ub2 + FLINT_BITS - 1) / FLINT_BITS;
        slong nn = FLINT_MAX(n1s, n2s);
        int karatsuba = (minlen >= _flint_mpn_poly_mulmid_cutoff(n1b, n2b, squaring));

        method = -1;

        if (n1m * n2m < n1s * n2s && force != FLINT_MPN_POLY_MUL_SIGNED && force != FLINT_MPN_POLY_MUL_BIAS)
        {
            /* limb penalty factor of the two's complement representation */
            slong penalty = (n1s * n2s + n1m * n2m / 2) / (n1m * n2m);
            slong sm_cutoff = _flint_mpn_poly_mulmid_cutoff(n1m, n2m, squaring);

            /* Measured on full-magnitude data, the classical sign-magnitude
               kernels stay ahead of the signed Karatsuba (which pays the
               extra limb throughout: the linear combinations force the
               growth limb in any representation) until about 12x the
               Karatsuba cutoff for penalty >= 3 and 9x for penalty 2.
               The generic sign-magnitude loop is slower than the kernels,
               so only a small factor applies beyond the kernel range. */
            if (n1m + n2m <= 7)
                sm_cutoff *= (penalty >= 3) ? 12 : (penalty >= 2) ? 9 : 2;
            else
                sm_cutoff *= 2 + (penalty >= 2);

            if (minlen < sm_cutoff || force == FLINT_MPN_POLY_MUL_SIGNMAG)
                method = FLINT_MPN_POLY_MUL_SIGNMAG;
        }

        if (method == -1)
        {
            if (!karatsuba)
            {
                /* classical: the generic (branchy) signed dot product is never
                   competitive, so the signed method needs the kernels */
                method = (nn <= FLINT_MPN_DOT_TAB_N && n1s + n2s + 1 <= 8 && (nn <= 3 || minlen <= 12))
                            ? FLINT_MPN_POLY_MUL_SIGNED : FLINT_MPN_POLY_MUL_BIAS;
            }
            else
            {
                method = ((nn <= 2 && minlen <= 64) || (nn == 3 && minlen <= 24))
                            ? FLINT_MPN_POLY_MUL_SIGNED : FLINT_MPN_POLY_MUL_BIAS;
            }

            if (force == FLINT_MPN_POLY_MUL_SIGNED)
                method = FLINT_MPN_POLY_MUL_SIGNED;
            else if (force == FLINT_MPN_POLY_MUL_BIAS)
                method = FLINT_MPN_POLY_MUL_BIAS;
        }
    }

    if (method == FLINT_MPN_POLY_MUL_SIGNED)
    {
        /* two's complement: both operands need a sign bit */
        eb1 = ub1 + 1;
        eb2 = ub2 + 1;
    }
    else if (method == FLINT_MPN_POLY_MUL_SIGNMAG)
    {
        /* magnitudes only; the sign bit is needed in the output */
        eb1 = ub1;
        eb2 = ub2;
    }
    else
    {
        /* unsigned, or bias: only signed operands need an extra bit */
        eb1 = ub1 + sgn1;
        eb2 = ub2 + sgn2;
    }

    n1 = (eb1 + FLINT_BITS - 1) / FLINT_BITS;
    n2 = (eb2 + FLINT_BITS - 1) / FLINT_BITS;

    P->method = method;
    P->nlimbs1 = n1;
    P->norm1 = FLINT_BITS * n1 - eb1;
    P->nlimbs2 = n2;
    P->norm2 = FLINT_BITS * n2 - eb2;
    P->slimbs = _flint_mpn_poly_mulmid_slimbs(minlen, n1, P->norm1, n2, P->norm2, method,
                    _flint_mpn_poly_mulmid_cutoff(n1, n2, squaring));
}

void
_flint_mpn_poly_mulmid(nn_ptr res, nn_srcptr f, slong flen, slong nlimbs1, slong norm1, nn_srcptr g, slong glen, slong nlimbs2, slong norm2, slong nlo, slong nhi, slong slimbs, int method)
{
    int squaring = (f == g && flen == glen && nlimbs1 == nlimbs2);
    slong len1 = flen, len2 = glen, cutoff;
    int sgn = (method == FLINT_MPN_POLY_MUL_SIGNED) ? 1 : (method == FLINT_MPN_POLY_MUL_SIGNMAG) ? 2 : 0;

    _flint_mpn_poly_mulmid_trunc(&len1, &len2, nlo, nhi);
    cutoff = _flint_mpn_poly_mulmid_cutoff(nlimbs1, nlimbs2, squaring);

    FLINT_ASSERT(slimbs >= _flint_mpn_poly_mulmid_slimbs(FLINT_MIN(len1, len2), nlimbs1, norm1, nlimbs2, norm2, method, cutoff));

    if (!_flint_mpn_poly_mulmid_use_karatsuba(FLINT_MIN(len1, len2), FLINT_MIN(norm1, norm2), method, cutoff))
    {
        if (nlimbs1 >= nlimbs2)
            _flint_mpn_poly_mulmid_classical(res, f, flen, nlimbs1, g, glen, nlimbs2, nlo, nhi, slimbs, sgn);
        else
            _flint_mpn_poly_mulmid_classical(res, g, glen, nlimbs2, f, flen, nlimbs1, nlo, nhi, slimbs, sgn);
    }
    else
    {
        _flint_mpn_poly_mulmid_karatsuba(res, f, flen, nlimbs1, norm1, g, glen, nlimbs2, norm2, nlo, nhi, slimbs, cutoff, sgn);
    }
}

void
_flint_mpn_poly_mul(nn_ptr res, nn_srcptr f, slong flen, slong nlimbs1, slong norm1, nn_srcptr g, slong glen, slong nlimbs2, slong norm2, slong slimbs, int method)
{
    _flint_mpn_poly_mulmid(res, f, flen, nlimbs1, norm1, g, glen, nlimbs2, norm2, 0, flen + glen - 1, slimbs, method);
}
