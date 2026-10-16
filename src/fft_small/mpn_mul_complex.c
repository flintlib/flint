/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fft_small.h"

/*
    Complex products on transformed operands (see fft_small.h for the
    output contracts). One engine serves the full, low, high and windowed
    products and squares of flint_mpn_*_complex:

      - each output part is a sum of at most two products of parts,
        re = ar br - ai bi and im = ar bi + ai br (for squares
        re = ar^2 - ai^2 and im = 2 ar ai), keeping only the products of
        nonzero parts;
      - the plan accumulates as many products per slot as the largest
        output part has (a doubled product counts two), signed when some
        output must be exported signed: one with two products, or any
        window not starting at limb 0; otherwise the plan is unsigned,
        packing more bits per chunk;
      - only the nonzero parts are transformed; the products are formed
        on magnitudes with their signs applied in transform space (by
        choosing additions or subtractions and at most one negation per
        output part) for signed exports, and attached afterwards for
        unsigned ones -- a single product with lo = 0, whose sign is
        known;
      - each output part is written over an operand that is dead after
        its first product when there is one (five transforms' storage
        for a general product, three for a square), and the dead
        operands' storage serves the conversions' staging.

    Windows use the window plan (cyclic transforms when the window
    allows); the sign-magnitude outputs use the ordinary plan and the
    staged exports, exact for lo = 0 and truncated otherwise, each
    output reconstructing only the slots of its own products.

    A full product (lo = 0) of a real or imaginary operand by a complex
    one -- each output part a single product, the two sharing a part --
    is exactly the task of _fft_small_mpn_mul_shared, and goes there:
    its fused per-product pipeline was measured 5-30% faster than this
    engine on such shapes (the engine's separate passes over all
    operands cost some 5%, and forming the shorter product at the
    length of the longer up to a further 25% for unbalanced parts).
*/

typedef struct
{
    nn_srcptr d;
    slong n;        /* normalized length, 0 for a zero part */
    int sgn;
}
_part_struct;

typedef struct
{
    int a, b;       /* parts multiplied (a == b: a square) */
    int e;          /* the product's sign: (-1)^e |a| |b| */
}
_term_struct;

typedef struct
{
    _term_struct t[2];
    int nt;         /* number of products; 0 for a zero output */
    int twice;      /* the single product is doubled */
    int is_signed;  /* exported signed (else unsigned, sign attached) */
    int dest;       /* operand slot receiving the output */
}
_out_struct;

#define MODE_WINDOW 0
#define MODE_SIGNED 1

static void
_part_set(_part_struct * p, nn_srcptr d, slong n, int sgn, slong hi)
{
    n = FLINT_MIN(n, hi);
    while (n > 0 && d[n - 1] == 0)
        n--;
    p->d = d;
    p->n = n;
    p->sgn = sgn;
}

static void
_add_term(_out_struct * o, int a, int b, int e)
{
    o->t[o->nt].a = a;
    o->t[o->nt].b = b;
    o->t[o->nt].e = e;
    o->nt++;
}

/* the terms of both output parts; parts 0 ar, 1 ai, 2 br, 3 bi */
static void
_complex_terms(_out_struct * out, const _part_struct * A, int squaring)
{
    static const int pa[2][2] = { { 0, 1 }, { 0, 1 } };
    static const int pb[2][2] = { { 2, 3 }, { 3, 2 } };
    static const int pneg[2][2] = { { 0, 1 }, { 0, 0 } };
    int k, j;

    for (k = 0; k < 2; k++)
    {
        out[k].nt = 0;
        out[k].twice = 0;
    }

    if (squaring)
    {
        /* re = ar^2 - ai^2, im = 2 ar ai */
        if (A[0].n > 0)
            _add_term(out + 0, 0, 0, 0);
        if (A[1].n > 0)
            _add_term(out + 0, 1, 1, 1);
        if (A[0].n > 0 && A[1].n > 0)
        {
            _add_term(out + 1, 0, 1, A[0].sgn ^ A[1].sgn);
            out[1].twice = 1;
        }
        return;
    }

    for (k = 0; k < 2; k++)
        for (j = 0; j < 2; j++)
        {
            int a = pa[k][j], b = pb[k][j];

            if (A[a].n > 0 && A[b].n > 0)
                _add_term(out + k, a, b, A[a].sgn ^ A[b].sgn ^ pneg[k][j]);
        }
}

/* whether operand x is read by a term after term j of output k */
static int
_used_after(const _out_struct * out, int k, int j, int x)
{
    int kk, jj;

    for (kk = k; kk < 2; kk++)
        for (jj = (kk == k) ? j + 1 : 0; jj < out[kk].nt; jj++)
            if (out[kk].t[jj].a == x || out[kk].t[jj].b == x)
                return 1;

    return 0;
}

/* the output's value (signed export) or magnitude (unsigned export)
   into Z, which may be an operand of the first product only */
static void
_form_output(fft_small_op_struct * Z, const _out_struct * o,
             fft_small_op_struct * X, const fft_small_plan_struct * P)
{
    const _term_struct * t = o->t;

    if (t[0].a == t[0].b)
        fft_small_op_sqr(Z, X + t[0].a, P);
    else
        fft_small_op_mul(Z, X + t[0].a, X + t[0].b, P);

    if (o->nt == 2)
    {
        int sub = (t[1].e != t[0].e);

        if (t[1].a == t[1].b)
        {
            if (sub)
                fft_small_op_subsqr(Z, X + t[1].a, P);
            else
                fft_small_op_addsqr(Z, X + t[1].a, P);
        }
        else
        {
            if (sub)
                fft_small_op_submul(Z, X + t[1].a, X + t[1].b, P);
            else
                fft_small_op_addmul(Z, X + t[1].a, X + t[1].b, P);
        }
    }

    if (o->twice)
        fft_small_op_add(Z, Z, Z, P);

    if (o->is_signed && t[0].e)
        fft_small_op_neg(Z, Z, P);
}

/* exports of one output part */
static void
_export_window(nn_ptr z, slong lo, slong hi, const _out_struct * o,
               fft_small_op_struct * Z, ulong terms, ulong nslots,
               const fft_small_plan_struct * P)
{
    slong zn = hi - lo;

    if (o->is_signed)
    {
        fft_small_ifft_trunc(Z, nslots, 0, NULL, P);
        fft_small_export_mpn_signed_window(z, zn, Z, nslots, lo, P);
    }
    else
    {
        /* lo = 0 and a known sign: the residue mod B^hi of the
           magnitude, by the unsigned reconstruction */
        ulong need = fft_small_export_mpn_limbs(P, terms, nslots);
        ulong en = FLINT_MIN((ulong) zn, need);

        fft_small_ifft_trunc(Z, nslots, n_cdiv(en * FLINT_BITS, P->bits),
                             NULL, P);
        fft_small_export_mpn(z, en, Z, P);
        if ((ulong) zn > en)
            flint_mpn_zero(z + en, zn - en);
        if (o->t[0].e)
            mpn_neg(z, z, zn);
    }
}

static slong
_export_signed(nn_ptr z, slong lo, slong hi, const _out_struct * o,
               fft_small_op_struct * Z, ulong terms, ulong nslots,
               const fft_small_plan_struct * P, fft_small_op_pool_struct * pool)
{
    ulong zn = hi - lo, need, bn, len;
    ulong * buf;
    int sgn, fits;

    if (o->is_signed)
    {
        ulong needf = fft_small_export_mpn_signed_limbs(P, nslots);

        need = ((ulong) lo < needf) ? needf - lo : 1;
        bn = FLINT_MAX(zn, need);
        buf = fft_small_stage_begin(z, zn, need, pool, 1);
        fft_small_ifft_trunc(Z, nslots, 0, NULL, P);
        fft_small_export_mpn_signed_trunc(buf, bn, &sgn, Z, nslots, lo, P);
    }
    else
    {
        need = fft_small_export_mpn_limbs(P, terms, nslots);
        bn = FLINT_MAX(zn, need);
        buf = fft_small_stage_begin(z, zn, need, pool, 1);
        fft_small_ifft_trunc(Z, nslots, n_cdiv(need * FLINT_BITS, P->bits),
                             NULL, P);
        fft_small_export_mpn(buf, need, Z, P);
        if (bn > need)
            flint_mpn_zero(buf + need, bn - need);
        sgn = o->t[0].e;
    }

    len = bn;
    while (len > 0 && buf[len - 1] == 0)
        len--;

    fits = fft_small_stage_end(z, zn, buf, len, pool);
    FLINT_ASSERT(fits);
    (void) fits;

    return (len == 0 || !sgn) ? (slong) len : -(slong) len;
}

/* whether [p, p + pn) and [q, q + qn) overlap */
static int
_overlap(nn_srcptr p, slong pn, nn_srcptr q, slong qn)
{
    return (p < q + qn) && (q < p + pn);
}

/* One operand with a single nonzero part: each output part is a single
   product, the two sharing that part. These are the products of
   _fft_small_mpn_mul_shared, whose fused per-product pipeline the
   engine below cannot match: the engine makes the transforms, the
   pointwise products and the inverse transforms separate passes over
   all operands (in the two-prime range this alone costs some 5% over
   the fused pipeline) and forms both products at the plan of the
   longer. Returns 0 (doing nothing) for other shapes, windows not
   starting at limb 0 or too short for the products, and outputs
   overlapping an operand they must not. */
static int
_complex_shared(nn_ptr * z, slong * zlen, slong lo, slong hi,
                const _part_struct * A, const _out_struct * out)
{
    const _term_struct * t0 = out[0].t, * t1 = out[1].t;
    int c, p0, p1, k, first;
    slong n0, n1;

    if (lo != 0 || out[0].nt != 1 || out[1].nt != 1
        || out[0].twice || out[1].twice
        || t0->a == t0->b || t1->a == t1->b)
        return 0;

    /* the common part c and the other part of each product */
    if (t0->a == t1->a || t0->a == t1->b)
        c = t0->a;
    else if (t0->b == t1->a || t0->b == t1->b)
        c = t0->b;
    else
        return 0;
    p0 = (t0->a == c) ? t0->b : t0->a;
    p1 = (t1->a == c) ? t1->b : t1->a;

    n0 = A[p0].n + A[c].n;
    n1 = A[p1].n + A[c].n;
    if (n0 > hi || n1 > hi)
        return 0;

    /* the output written first must not overlap the other product's
       own operand; the common operand is read before either is
       written (see _fft_small_mpn_mul_shared) */
    if (!_overlap(z[0], hi, A[p1].d, A[p1].n))
        first = 0;
    else if (!_overlap(z[1], hi, A[p0].d, A[p0].n))
        first = 1;
    else
        return 0;
    if (_overlap(z[0], hi, z[1], hi))
        return 0;

    if (first == 0)
        _fft_small_mpn_mul_shared(z[0], A[p0].d, A[p0].n,
                                  z[1], A[p1].d, A[p1].n, A[c].d, A[c].n);
    else
        _fft_small_mpn_mul_shared(z[1], A[p1].d, A[p1].n,
                                  z[0], A[p0].d, A[p0].n, A[c].d, A[c].n);

    for (k = 0; k < 2; k++)
    {
        slong n = k ? n1 : n0, len = n;
        flint_mpn_zero(z[k] + n, hi - n);
        while (len > 0 && z[k][len - 1] == 0)
            len--;
        zlen[k] = (len == 0 || !out[k].t[0].e) ? len : -len;
    }

    return 1;
}

static int
_complex_core(nn_ptr * z, slong * zlen, slong lo, slong hi, int mode,
              const _part_struct * A, int squaring)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[6];
    int have[6];
    _out_struct out[2];
    slong an_max, bn_max, zn = hi - lo;
    ulong terms, nslots;
    int any_signed, k, j, nops;

    FLINT_ASSERT(0 <= lo && lo < hi);

    _complex_terms(out, A, squaring);

    if (mode == MODE_SIGNED && !squaring
        && _complex_shared(z, zlen, lo, hi, A, out))
        return 1;

    if (out[0].nt == 0 && out[1].nt == 0)
    {
        for (k = 0; k < 2; k++)
        {
            flint_mpn_zero(z[k], zn);
            if (mode == MODE_SIGNED)
                zlen[k] = 0;
        }
        return 1;
    }

    /* sizes, capacity and signedness */
    an_max = FLINT_MAX(A[0].n, A[1].n);
    bn_max = squaring ? an_max : FLINT_MAX(A[2].n, A[3].n);
    terms = 1;
    any_signed = 0;
    for (k = 0; k < 2; k++)
    {
        ulong tk = out[k].nt + out[k].twice;
        terms = FLINT_MAX(terms, tk);
        out[k].is_signed = (out[k].nt == 2) || (lo > 0);
        if (out[k].nt > 0 && out[k].is_signed)
            any_signed = 1;
    }

    /* the sign-magnitude exports read the whole reconstruction, so
       they never use the cyclic transforms */
    if (!fft_small_plan_init_mpn_window(P, R, an_max, bn_max, terms,
                any_signed, (mode == MODE_WINDOW) ? lo : 0, hi))
        return 0;

    nslots = fft_small_plan_mpn_nslots(P, an_max, bn_max);

    /* storage: the transformed parts, and a slot for an output part when
       no operand of its first product is dead after it */
    nops = 0;
    for (j = 0; j < 6; j++)
        have[j] = (j < 4) && (A[j].n > 0) && !(squaring && j >= 2);
    for (j = 0; j < 4; j++)
        nops += have[j];
    for (k = 0; k < 2; k++)
    {
        out[k].dest = -1;
        if (out[k].nt == 0)
            continue;
        if (!_used_after(out, k, 0, out[k].t[0].a))
            out[k].dest = out[k].t[0].a;
        else if (!_used_after(out, k, 0, out[k].t[0].b))
            out[k].dest = out[k].t[0].b;
        else
        {
            out[k].dest = 4 + k;
            have[4 + k] = 1;
            nops++;
        }
    }

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), nops, 0,
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (j = 0; j < 6; j++)
        if (have[j])
            fft_small_op_init_borrowed(X + j, P, fft_small_op_pool_take(pool));

    for (j = 0; j < 4; j++)
        if (have[j])
            fft_small_fft_mpn(X + j, A[j].d, A[j].n, P);

    for (k = 0; k < 2; k++)
        if (out[k].nt > 0)
            _form_output(X + out[k].dest, out + k, X, P);

    /* the operands are dead: their storage serves the staging */
    for (j = 0; j < 6; j++)
    {
        if (have[j] && (out[0].nt == 0 || j != out[0].dest)
                    && (out[1].nt == 0 || j != out[1].dest))
        {
            fft_small_op_pool_give(pool, X[j].data);
            have[j] = 0;
        }
    }

    for (k = 0; k < 2; k++)
    {
        if (out[k].nt == 0)
        {
            flint_mpn_zero(z[k], zn);
            if (mode == MODE_SIGNED)
                zlen[k] = 0;
            continue;
        }

        if (mode == MODE_WINDOW)
            _export_window(z[k], lo, hi, out + k, X + out[k].dest, terms,
                           nslots, P);
        else
        {
            /* the slots of this output's own products (no aliasing:
               the sign-magnitude modes use the linear plan) */
            ulong ns = 0;
            int j2;
            for (j2 = 0; j2 < out[k].nt; j2++)
                ns = FLINT_MAX(ns, fft_small_plan_mpn_nslots(P,
                        A[out[k].t[j2].a].n, A[out[k].t[j2].b].n));
            zlen[k] = _export_signed(z[k], lo, hi, out + k, X + out[k].dest,
                                     terms, ns, P, pool);
        }

        fft_small_op_pool_give(pool, X[out[k].dest].data);
    }

    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);
    return 1;
}

static int
_mul_entry(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len,
           slong lo, slong hi, int mode,
           nn_srcptr ar, slong arn, int ar_sgn,
           nn_srcptr ai, slong ain, int ai_sgn,
           nn_srcptr br, slong brn, int br_sgn,
           nn_srcptr bi, slong bin, int bi_sgn, int squaring)
{
    _part_struct A[4];
    nn_ptr z[2];
    slong zlen[2];
    int ok;

    _part_set(A + 0, ar, arn, ar_sgn, hi);
    _part_set(A + 1, ai, ain, ai_sgn, hi);
    if (squaring)
    {
        A[2] = A[0];
        A[3] = A[1];
    }
    else
    {
        _part_set(A + 2, br, brn, br_sgn, hi);
        _part_set(A + 3, bi, bin, bi_sgn, hi);
    }

    z[0] = zr;
    z[1] = zi;
    ok = _complex_core(z, zlen, lo, hi, mode, A, squaring);

    if (ok && mode == MODE_SIGNED)
    {
        *zr_len = zlen[0];
        *zi_len = zlen[1];
    }

    return ok;
}

int
_fft_small_mpn_mulmid_complex(nn_ptr zr, nn_ptr zi, slong lo, slong hi,
    nn_srcptr ar, slong arn, int ar_sgn,
    nn_srcptr ai, slong ain, int ai_sgn,
    nn_srcptr br, slong brn, int br_sgn,
    nn_srcptr bi, slong bin, int bi_sgn)
{
    return _mul_entry(zr, NULL, zi, NULL, lo, hi, MODE_WINDOW,
                      ar, arn, ar_sgn, ai, ain, ai_sgn,
                      br, brn, br_sgn, bi, bin, bi_sgn, 0);
}

int
_fft_small_mpn_sqrmid_complex(nn_ptr zr, nn_ptr zi, slong lo, slong hi,
    nn_srcptr ar, slong arn, int ar_sgn,
    nn_srcptr ai, slong ain, int ai_sgn)
{
    return _mul_entry(zr, NULL, zi, NULL, lo, hi, MODE_WINDOW,
                      ar, arn, ar_sgn, ai, ain, ai_sgn,
                      NULL, 0, 0, NULL, 0, 0, 1);
}

int
_fft_small_mpn_mul_complex(nn_ptr zr, slong * zr_len,
    nn_ptr zi, slong * zi_len, slong lo, slong hi,
    nn_srcptr ar, slong arn, int ar_sgn,
    nn_srcptr ai, slong ain, int ai_sgn,
    nn_srcptr br, slong brn, int br_sgn,
    nn_srcptr bi, slong bin, int bi_sgn)
{
    return _mul_entry(zr, zr_len, zi, zi_len, lo, hi, MODE_SIGNED,
                      ar, arn, ar_sgn, ai, ain, ai_sgn,
                      br, brn, br_sgn, bi, bin, bi_sgn, 0);
}

int
_fft_small_mpn_sqr_complex(nn_ptr zr, slong * zr_len,
    nn_ptr zi, slong * zi_len, slong lo, slong hi,
    nn_srcptr ar, slong arn, int ar_sgn,
    nn_srcptr ai, slong ain, int ai_sgn)
{
    return _mul_entry(zr, zr_len, zi, zi_len, lo, hi, MODE_SIGNED,
                      ar, arn, ar_sgn, ai, ain, ai_sgn,
                      NULL, 0, 0, NULL, 0, 0, 1);
}
