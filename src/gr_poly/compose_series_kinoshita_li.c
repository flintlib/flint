/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version. See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* This file was developed with the assistance of Claude Sonnet 4.6. */

/*
    Power series composition f(g(x)) mod x^n via the quasilinear algorithm of

        Yasunori Kinoshita and Baitian Li,
        "Power Series Composition in Near-Linear Time",
        arXiv:2404.05177 (2024).

    Setup
    -----
    Given f (length len1) and g (with g(0) = 0), set:

        P(y)   = f.reverse(n-1)   in R[y]      (P[i] = f[n-1-i])
        Q(x,y) = 1 - y * g(x)    in R[x][y]

    Then  f(g(x)) mod x^n  =  [y^{n-1}]( P(y)/Q(x,y) mod x^n ).

    Algorithm overview
    ------------------
    The algorithm is a two-phase iterative computation.

    At each level k = 0, 1, ..., L-1 (where L = ceil(log2(n))), the
    "current Q" has x-degree < xn_k and actual y-degree ydeg_k-1.
    The Graeffe step replaces Q by V where V(x^2,y) = Q(x,y)*Q(-x,y),
    halving the x-degree and approximately doubling the y-degree.

    DOWNWARD PASS (Graeffe chain):
        Compute Q_0 = initial Q, Q_1 = Graeffe(Q_0), Q_2 = Graeffe(Q_1), ...
        until Q_L has x-degree 1 (base case).  Store all Q_k in a single
        pre-allocated flat array "Qchain".

    BASE CASE:
        Solve the univariate problem: W = P(y) / Q_L(0,y) mod y^n.
        (ylo_arr[L] = 0 always, so the full length-n result is needed.)

    UPWARD PASS (post-multiply chain):
        For k = L-1 down to 0:
            W <- [y-slice of Q_k(-x,y) * W(x^2,y)]
        where the y-slice is chosen so that the final W is [y^{n-1}] of P/Q_0,
        which equals f(g(x)) mod x^n.

    Memory analysis
    ---------------
    The recursive formulation pins O(N) data per live stack frame, and with
    O(log N) frames simultaneously live, uses O(N log N) memory total.

    This iterative formulation separates storage from processing:

      Qchain: holds Q_0, Q_1, ..., Q_L end-to-end.
        - Q_k has x-stride xn_k and y-bound ydeg_k, size = xn_k * ydeg_k.
        - Each entry exceeds n: xn_k * ydeg_k ≈ (n/2^k) * (2^{k+1}-1) > n.
        - With L = ceil(log2(n)) entries, total Qchain size is O(n log n).
          (Measured: ~11n for n=1000, ~18n for n=10000.)

      Working memory (W buffers, KS temporaries): O(n) at every instant.
        - W has x-degree < xn_k and y-width yn_W_k ≈ 2^k, so size ≈ n.
        - KS temporaries scale with the current product length ≈ 4n.

    Total memory is O(n log n), dominated by the Qchain.  The working
    memory outside the Qchain is O(n) at every instant, but the Qchain
    itself must be fully allocated upfront and cannot be reduced below
    O(n log n) without recomputing Graeffe steps on the way back up.

    Bivariate storage layout
    ------------------------
    Q and W are stored as flat arrays throughout the algorithm.  The coefficient
    of x^i y^j in a bivariate with x-degree < xn and y-degree-bound ydeg sits at
    index i + j*xn  (x-contiguous).  This storage format is independent of the
    KS packing used during multiplication.

    Kronecker substitution (KS)
    --------------------------------------
    To multiply two bivariates F (x-degree < a, y-degree < p) and
    G (x-degree < b, y-degree < q), we use the substitution y -> t, x -> t^s
    with stride s = p+q-1.  The coefficient of x^i y^j maps to flat index i*s + j.
    The univariate product has length s*(a+b-1), and _gr_poly_mullow is called
    with an output length trimmed to only the portion that is actually read back.

    In benchmarks, y-major substitution consistently performed better than
    the transposed (x-major) substitution.

    For a bivariate with x-degree < a and y-degree < b packed with stride s,
    the exact (tight) KS array length is (a-1)*s + b: the highest-occupied
    flat index is (a-1)*s + (b-1), so the array length is (a-1)*s + b.

    Downward pass: s = 2*ydeg - 1.  Both inputs have length (xn-1)*s + ydeg.
      Output (rks_len): we unpack even-x columns 0,2,...,2*(xnh-1), each of
      y-height s.  Highest read: 2*(xnh-1)*s + (ydA-1) < s*(2*xnh-1).  -> s*(2*xnh-1).

    Upward pass: s = ydeg + yn_W - 1.  Qneg has length (xn-1)*s + ydeg.
      W(x^2) packs x^{2i} y^j at 2*i*s+j; highest: 2*(xnh-1)*s+(yn_W-1)
      -> W_ks_len = 2*(xnh-1)*s + yn_W.
      Output (R_rks_len): we read R at i*s+jsrc for i<xn, jsrc<yn_W;
      highest: (xn-1)*s + yn_W - 1 -> R_rks_len = (xn-1)*s + yn_W.
      Note: R_rks_len <= W_ks_len + qneg_ks_len - 1, and both inputs are
      shorter than the output (inputs <= R_rks_len because yn_W >= ydeg holds
      at every upward level, provable by induction).
*/

int
_gr_poly_compose_series_kinoshita_li(gr_ptr res, gr_srcptr poly1, slong len1,
                           gr_srcptr poly2, slong len2, slong n,
                           gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong k, L;

    FLINT_ASSERT(len1 >= 1);
    FLINT_ASSERT(len2 >= 1);
    FLINT_ASSERT(n >= 1);

    /* ---- Trivial cases ---- */
    if (len1 == 1 || n == 1)
    {
        status |= gr_set(res, poly1, ctx);
        status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), n - 1, ctx);
        return status;
    }

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    /* ---- Compute level count L = ceil(log2(n)) ---- */
    L = 0;
    { slong tmp = n; while (tmp > 1) { tmp = (tmp + 1) / 2; L++; } }

    /* ---- Precompute level parameters ---- */
    /*
        For level k:
            xn[k]      = x-degree bound of Q at this level
            ylo[k]     = lower y-index of the slice we want
            ydeg[k]    = y-degree bound (= ydegQ_actual) of Q at this level
            qoff[k]    = offset into Qchain where Q_k is stored

        xn[0] = n, ylo[0] = n-1, ydeg[0] = 2 (Q = 1 - y*g has y-degree 1).
        xn[k+1] = ceil(xn[k]/2).
        ylo[k+1] = max(0, ylo[k] - (ydeg[k] - 1)).
        ydeg[k+1] = min(2*ydeg[k] - 1, n).  (y-degree of Graeffe output)

        Level L is the base case (xn[L] = 1).
    */
    slong *xn_arr   = flint_malloc((L + 1) * sizeof(slong));
    slong *ylo_arr  = flint_malloc((L + 1) * sizeof(slong));
    slong *ydeg_arr = flint_malloc((L + 1) * sizeof(slong));
    slong *qoff_arr = flint_malloc((L + 1) * sizeof(slong));

    xn_arr[0]   = n;
    ylo_arr[0]  = n - 1;
    ydeg_arr[0] = 2;

    slong qchain_total = 0;
    for (k = 0; k < L; k++)
    {
        qoff_arr[k]      = qchain_total;
        qchain_total    += xn_arr[k] * ydeg_arr[k];
        xn_arr[k + 1]   = (xn_arr[k] + 1) / 2;
        ylo_arr[k + 1]  = (ylo_arr[k] >= ydeg_arr[k] - 1)
                              ? (ylo_arr[k] - (ydeg_arr[k] - 1)) : 0;
        ydeg_arr[k + 1] = FLINT_MIN(2 * ydeg_arr[k] - 1, n);
    }
    qoff_arr[L]   = qchain_total;
    qchain_total += xn_arr[L] * ydeg_arr[L];

    /* ---- Compute scratch_total: maximum KS working memory across all levels ---- */
    /*
     * Both passes use three temporary buffers laid out consecutively in a
     * single pre-allocated scratch vector, avoiding GR_TMP_INIT_VEC and
     * GR_TMP_CLEAR_VEC calls inside the loops.
     *
     * Downward level k: Qks (qks_len) | Qneg_ks (qks_len) | Rks (rks_len).
     * Upward level k:   Qneg_ks (qneg_ks_len) | W_ks (W_ks_len) | R_rks (R_rks_len).
     *
     * scratch_total = max over all levels of the per-level total size.
     */
    slong scratch_total = 0;
    {
        slong yn_W_pre = n;
        for (k = 0; k < L; k++)
        {
            slong xn   = xn_arr[k];
            slong ydeg = ydeg_arr[k];
            slong xnh  = xn_arr[k + 1];
            slong s       = 2 * ydeg - 1;
            slong qks_len = (xn - 1) * s + ydeg;
            slong rks_len = s * (2 * xnh - 1);
            scratch_total = FLINT_MAX(scratch_total, 2 * qks_len + rks_len);
        }
        for (k = L - 1; k >= 0; k--)
        {
            slong xn   = xn_arr[k];
            slong xnh  = xn_arr[k + 1];
            slong ydeg = ydeg_arr[k];
            slong ylo  = ylo_arr[k];
            slong ylor = ylo_arr[k + 1];
            slong yn_r      = n - ylor;
            slong yslice_lo = ylo - ylor;
            slong s           = ydeg + yn_W_pre - 1;
            slong qneg_ks_len = (xn - 1) * s + ydeg;
            slong W_ks_len    = 2 * (xnh - 1) * s + yn_W_pre;
            slong R_rks_len   = (xn - 1) * s + yn_r - yslice_lo;
            scratch_total = FLINT_MAX(scratch_total, qneg_ks_len + W_ks_len + R_rks_len);
            yn_W_pre = n - ylo;
        }
    }

    /* ---- Allocate Qchain and scratch ---- */
    gr_ptr Qchain;
    GR_TMP_INIT_VEC(Qchain, qchain_total, ctx);

    gr_ptr scratch;
    GR_TMP_INIT_VEC(scratch, scratch_total, ctx);

    /* ---- Initialise Q_0 = 1 - y*g(x) ---- */
    {
        gr_ptr Q0 = GR_ENTRY(Qchain, qoff_arr[0], sz);
        status |= gr_one(GR_ENTRY(Q0, 0, sz), ctx);
        slong glen = FLINT_MIN(len2, n);
        status |= _gr_vec_neg(GR_ENTRY(Q0, n, sz), poly2, glen, ctx);
    }

    /* ================================================================
       DOWNWARD PASS: compute Q_1, Q_2, ..., Q_L via Graeffe steps.

       y-major KS with stride s = 2*ydeg - 1.
       x^i y^j -> flat index i*s + j.
       ================================================================ */
    for (k = 0; k < L; k++)
    {
        slong xn   = xn_arr[k];
        slong ydeg = ydeg_arr[k];
        slong xnh  = xn_arr[k + 1];
        slong ydA  = ydeg_arr[k + 1];

        slong s       = 2 * ydeg - 1;
        slong qks_len = (xn - 1) * s + ydeg;   /* tight: highest pos (xn-1)*s+(ydeg-1) */
        slong rks_len = s * (2 * xnh - 1);      /* even-x cols 0..2*(xnh-1), each of height s */

        gr_srcptr Qk  = GR_ENTRY(Qchain, qoff_arr[k],     sz);
        gr_ptr    Qk1 = GR_ENTRY(Qchain, qoff_arr[k + 1], sz);

        gr_ptr Qks     = GR_ENTRY(scratch, 0,           sz);
        gr_ptr Qneg_ks = GR_ENTRY(scratch, qks_len,     sz);
        gr_ptr Rks     = GR_ENTRY(scratch, 2 * qks_len, sz);

        /* Qks and Qneg_ks must be zeroed before packing: the KS encoding
         * occupies only qks_len of the s*xn positions; gap slots between
         * columns must be zero for the multiply to be correct. */
        status |= _gr_vec_zero(Qks,     qks_len, ctx);
        status |= _gr_vec_zero(Qneg_ks, qks_len, ctx);

        /* Pack Q and Q(-x,y): x^i y^j -> flat index i*s + j */
        for (slong j = 0; j < ydeg; j++)
            for (slong i = 0; i < xn; i++)
            {
                gr_srcptr src = GR_ENTRY(Qk, j * xn + i, sz);
                status |= gr_set(GR_ENTRY(Qks,     i * s + j, sz), src, ctx);
                if (i & 1)
                    status |= gr_neg(GR_ENTRY(Qneg_ks, i * s + j, sz), src, ctx);
                else
                    status |= gr_set(GR_ENTRY(Qneg_ks, i * s + j, sz), src, ctx);
            }

        /* Rks is fully written by mullow; no zeroing needed. */
        status |= _gr_poly_mullow(Rks, Qks, qks_len, Qneg_ks, qks_len, rks_len, ctx);

        /* Unpack even-x columns -> Q_{k+1}.
         * [x^{2i} y^j] lives at flat index 2*i*s + j in the y-major layout. */
        for (slong j = 0; j < ydA; j++)
            for (slong i = 0; i < xnh; i++)
                status |= gr_set(GR_ENTRY(Qk1, j * xnh + i, sz),
                                 GR_ENTRY(Rks, 2 * i * s + j, sz), ctx);
    }

    /* ================================================================
       BASE CASE: solve P(y) / Q_L(0,y) mod y^n.
       Q_L has x-stride = xn_arr[L] = 1, so Q_L[j*1+0] = Q_L[j].
       ylo_arr[L] = 0 always, so the full length-n result fills W.
       ================================================================ */

    /*
     * W is used in-place across all upward steps: at each level the old
     * contents are fully consumed (packed into W_ks) before the new
     * contents are written (unpacked from R_rks), so no temporary is needed.
     * The required capacity is the maximum of xn_k * yn_W_k over all levels,
     * which never exceeds the next power of two >= n.
     */
    slong W_alloc = 1;
    while (W_alloc < n) W_alloc *= 2;

    gr_ptr W;
    GR_TMP_INIT_VEC(W, W_alloc, ctx);

    slong yn_W = n;   /* ylo_arr[L] == 0 always */

    {
        /* Build P(y) = f.reverse(n-1) into scratch, which is at least n elements. */
        gr_ptr P = scratch;
        slong glen = FLINT_MIN(len1, n);
        status |= _gr_poly_reverse(P, poly1, glen, glen, ctx);

        gr_srcptr QL = GR_ENTRY(Qchain, qoff_arr[L], sz);
        slong jmax = FLINT_MIN(ydeg_arr[L], glen);

        status |= _gr_poly_div_series(GR_ENTRY(W, n - glen, sz), P, glen, QL, jmax, glen, ctx);
        /* Low (n - glen) coefficients of W are a priori zero. */
    }

    /* ================================================================
       UPWARD PASS: for k = L-1 down to 0, compute the post-multiply.

       At level k, W has:
           x-stride = xn_arr[k+1]
           y-width  = yn_W  (= n at the start, shrinking toward 1)

       Compute r = Q_k(-x,y) * W(x^2,y) using y-major KS with
       stride s = ydeg + yn_W - 1, then extract the y-slice
       [ylo_arr[k]-ylo_arr[k+1], n-ylo_arr[k+1]) of r into W (in place).

       In-place safety: W is fully read into W_ks before R_rks is unpacked
       back into W, so there is no read-after-write conflict.
       ================================================================ */
    for (k = L - 1; k >= 0; k--)
    {
        slong xn   = xn_arr[k];
        slong xnh  = xn_arr[k + 1];
        slong ydeg = ydeg_arr[k];
        slong ylo  = ylo_arr[k];
        slong ylor = ylo_arr[k + 1];

        slong yn_r = n - ylor;   /* loop bound for jsrc; equals yn_W */

        slong s           = ydeg + yn_W - 1;
        slong qneg_ks_len = (xn - 1) * s + ydeg;       /* tight: highest pos (xn-1)*s+(ydeg-1) */
        slong W_ks_len    = 2 * (xnh - 1) * s + yn_W;  /* tight: highest pos 2*(xnh-1)*s+(yn_W-1) */

        /*
         * We only need product coefficients at positions i*s+jsrc for
         * jsrc in [yslice_lo, yn_r) and i in [0, xn).  The lowest such
         * position is yslice_lo (i=0, jsrc=yslice_lo) and the highest is
         * (xn-1)*s + yn_r - 1.  Using _gr_poly_mulmid with nlo=yslice_lo
         * avoids computing the yslice_lo coefficients at the low end that
         * are never read.  R_rks_len = (xn-1)*s + yn_r - yslice_lo.
         */
        slong yslice_lo   = ylo - ylor;
        slong R_rks_len   = (xn - 1) * s + yn_r - yslice_lo;

        gr_srcptr Qk = GR_ENTRY(Qchain, qoff_arr[k], sz);

        gr_ptr Qneg_ks = GR_ENTRY(scratch, 0,                      sz);
        gr_ptr W_ks    = GR_ENTRY(scratch, qneg_ks_len,            sz);
        gr_ptr R_rks   = GR_ENTRY(scratch, qneg_ks_len + W_ks_len, sz);

        /* Qneg_ks and W_ks must be zeroed before packing: KS gap slots
         * between columns must be zero for the multiply to be correct.
         * R_rks is fully written by mulmid; no zeroing needed. */
        status |= _gr_vec_zero(Qneg_ks, qneg_ks_len, ctx);
        status |= _gr_vec_zero(W_ks,    W_ks_len,    ctx);

        /* Pack Qneg_k: x^i y^j -> flat index i*s + j, with sign flip for odd i */
        for (slong j = 0; j < ydeg; j++)
            for (slong i = 0; i < xn; i++)
            {
                gr_srcptr src = GR_ENTRY(Qk, j * xn + i, sz);
                gr_ptr    dst = GR_ENTRY(Qneg_ks, i * s + j, sz);
                if (i & 1)
                    status |= gr_neg(dst, src, ctx);
                else
                    status |= gr_set(dst, src, ctx);
            }

        /* Pack W(x^2): x^{2i} y^j -> flat index 2*i*s + j  (optimization 3) */
        for (slong j = 0; j < yn_W; j++)
            for (slong i = 0; i < xnh; i++)
                status |= gr_set(GR_ENTRY(W_ks, 2 * i * s + j,   sz),
                                 GR_ENTRY(W,    j * xnh + i,      sz), ctx);

        status |= _gr_poly_mulmid(R_rks, Qneg_ks, qneg_ks_len,
                                         W_ks,    W_ks_len,
                                         yslice_lo, yslice_lo + R_rks_len, ctx);

        /* Unpack y-slice [yslice_lo, yn_r) -> W (in place).
         * R_rks[k] holds product coefficient yslice_lo+k.
         * R[jsrc][i] was at i*s+jsrc; it is now at i*s+jsrc-yslice_lo in R_rks. */
        slong new_yn_W = n - ylo;
        {
            slong jout = 0;
            for (slong jsrc = yslice_lo; jsrc < yn_r; jsrc++, jout++)
                for (slong i = 0; i < xn; i++)
                    status |= gr_set(GR_ENTRY(W, jout * xn + i, sz),
                                     GR_ENTRY(R_rks, i * s + jsrc - yslice_lo, sz), ctx);
        }

        yn_W = new_yn_W;
    }

    /* W[0..n) holds the result f(g(x)) mod x^n. */
    status |= _gr_vec_set(res, W, n, ctx);

    GR_TMP_CLEAR_VEC(W, W_alloc, ctx);
    GR_TMP_CLEAR_VEC(scratch, scratch_total, ctx);
    GR_TMP_CLEAR_VEC(Qchain, qchain_total, ctx);

    flint_free(xn_arr);
    flint_free(ylo_arr);
    flint_free(ydeg_arr);
    flint_free(qoff_arr);

    return status;
}

int
gr_poly_compose_series_kinoshita_li(gr_poly_t res,
                    const gr_poly_t poly1,
                    const gr_poly_t poly2,
                    slong n,
                    gr_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;
    int status;

    if (len2 != 0)
    {
        truth_t is_zero = gr_is_zero(poly2->coeffs, ctx);
        if (is_zero == T_FALSE)   return GR_DOMAIN;
        if (is_zero == T_UNKNOWN) return GR_UNABLE;
    }

    if (len1 == 0 || n == 0)
        return gr_poly_zero(res, ctx);

    if (len2 == 0 || len1 == 1)
        return gr_poly_set_scalar(res, poly1->coeffs, ctx);

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    gr_poly_fit_length(res, lenr, ctx);
    status = _gr_poly_compose_series_kinoshita_li(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, ctx);
    _gr_poly_set_length(res, lenr, ctx);
    _gr_poly_normalise(res, ctx);

    return status;
}

