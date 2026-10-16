/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

/*  Adapted from _gr_poly_compose_series_kinoshita_li (see comments in that
    function for more details about the algorithm) with the only difference that
    we track a common denominator for each polynomial.

    This file was developed with the assistance of Claude Sonnet 4.6.*/

void
_fmpq_poly_compose_series_kinoshita_li(fmpz *res, fmpz_t rden,
                                        const fmpz *poly1, const fmpz_t den1, slong len1,
                                        const fmpz *poly2, const fmpz_t den2, slong len2,
                                        slong n)
{
    slong k, L;

    FLINT_ASSERT(len1 >= 1);
    FLINT_ASSERT(len2 >= 1);
    FLINT_ASSERT(n >= 1);

    /* ---- Trivial cases ---- */
    if (len1 == 1 || n == 1)
    {
        _fmpz_vec_set(res, poly1, 1);
        _fmpz_vec_zero(res + 1, n - 1);
        fmpz_set(rden, den1);
        _fmpq_poly_canonicalise(res, rden, 1);
        return;
    }

    /* ---- Level count L = ceil(log2(n)) ---- */
    L = 0;
    { slong tmp = n; while (tmp > 1) { tmp = (tmp + 1) / 2; L++; } }

    /* ---- Level parameters ---- */
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

    /* ---- Qchain numerators and per-level denominators ---- */
    fmpz *Qchain_num = _fmpz_vec_init(qchain_total);
    fmpz *Qchain_den = _fmpz_vec_init(L + 1);   /* one denominator per level */

    /* ---- Initialise Q_0 = 1 - y*g(x) ---- */
    {
        fmpz *Q0_num = Qchain_num + qoff_arr[0];
        fmpz_set(Q0_num + 0, den2);
        slong glen = FLINT_MIN(len2, n);
        for (slong i = 0; i < glen; i++)
            fmpz_neg(Q0_num + n + i, poly2 + i);
        fmpz_set(Qchain_den + 0, den2);
        _fmpq_poly_canonicalise(Q0_num, Qchain_den + 0, n * 2);
    }

    /* ================================================================
       DOWNWARD PASS: Graeffe steps k = 0, ..., L-1.
       y-major KS with stride s = 2*ydeg - 1.
       x^i y^j maps to flat index i*s + j.
       ================================================================ */
    for (k = 0; k < L; k++)
    {
        slong xn   = xn_arr[k];
        slong ydeg = ydeg_arr[k];
        slong xnh  = xn_arr[k + 1];
        slong ydA  = ydeg_arr[k + 1];

        slong s       = 2 * ydeg - 1;
        slong qks_len = (xn - 1) * s + ydeg;
        slong rks_len = s * (2 * xnh - 1);

        fmpz *Qk_num  = Qchain_num + qoff_arr[k];
        fmpz *Qk1_num = Qchain_num + qoff_arr[k + 1];

        fmpz *Qks_num     = _fmpz_vec_init(qks_len);
        fmpz *Qneg_ks_num = _fmpz_vec_init(qks_len);
        fmpz *Rks_num     = _fmpz_vec_init(rks_len);
        fmpz_t Qks_den, Rks_den;
        fmpz_init_set(Qks_den, Qchain_den + k);
        fmpz_init(Rks_den);

        /* Pack Q and Q(-x,y): x^i y^j -> index i*s + j.
         * Qneg is the same but odd-x entries are negated. */
        for (slong j = 0; j < ydeg; j++)
            for (slong i = 0; i < xn; i++)
            {
                const fmpz *src = Qk_num + j * xn + i;
                if (i & 1)
                {
                    fmpz_set(Qks_num     + i * s + j,  src);
                    fmpz_neg(Qneg_ks_num + i * s + j,  src);
                }
                else
                {
                    fmpz_set(Qks_num     + i * s + j,  src);
                    fmpz_set(Qneg_ks_num + i * s + j,  src);
                }
            }

        /* Product denominator = Qks_den^2 */
        fmpz_mul(Rks_den, Qks_den, Qks_den);

        _fmpz_poly_mullow(Rks_num, Qks_num, qks_len, Qneg_ks_num, qks_len, rks_len);

        _fmpz_vec_clear(Qks_num,     qks_len);
        _fmpz_vec_clear(Qneg_ks_num, qks_len);

        /* Unpack even-x columns -> Q_{k+1}: x^{2i} y^j at index 2*i*s + j.
         * Copy only the xnh*ydA entries that are used, then canonicalise. */
        for (slong j = 0; j < ydA; j++)
            for (slong i = 0; i < xnh; i++)
                fmpz_set(Qk1_num + j * xnh + i, Rks_num + 2 * i * s + j);

        _fmpz_vec_clear(Rks_num, rks_len);

        fmpz_set(Qchain_den + (k + 1), Rks_den);
        _fmpq_poly_canonicalise(Qk1_num, Qchain_den + (k + 1), xnh * ydA);

        fmpz_clear(Qks_den);
        fmpz_clear(Rks_den);
    }

    /* ================================================================
       BASE CASE: W = P / Q_L(0,y) mod y^n.
       Q_L has x-stride 1, so Q_L_num[j] is the j-th y-coefficient.
       ylo_arr[L] = 0 always, so the full length-n result is needed.
       ================================================================ */
    slong W_alloc = 1;
    while (W_alloc < n) W_alloc *= 2;

    fmpz *W_num = _fmpz_vec_init(W_alloc);
    fmpz_t W_den;
    fmpz_init(W_den);
    slong yn_W = n;

    {
        /* ---- Build P_num = f.reverse(n-1), P_den = den1 ---- */
        slong glen = FLINT_MIN(n, len1);
        fmpz *P_num = _fmpz_vec_init(glen);

        _fmpz_poly_reverse(P_num, poly1, glen, glen);

        slong jmax = FLINT_MIN(ydeg_arr[L], glen);
        _fmpq_poly_div_series(W_num + n - glen, W_den,
                              P_num, den1, glen,
                              Qchain_num + qoff_arr[L], Qchain_den + L, jmax,
                              glen);
        /* Low (n - glen) coefficients of W are a priori zero. */

        _fmpz_vec_clear(P_num, glen);
    }

    /* ================================================================
       UPWARD PASS: k = L-1 down to 0.

       At level k, W has x-stride xn_arr[k+1] and y-width yn_W.
       We compute r = Q_k(-x,y) * W(x^2,y), then extract the y-slice
       [yslice_lo, yn_r) into W in place.
       ================================================================ */
    for (k = L - 1; k >= 0; k--)
    {
        slong xn   = xn_arr[k];
        slong xnh  = xn_arr[k + 1];
        slong ydeg = ydeg_arr[k];
        slong ylo  = ylo_arr[k];
        slong ylor = ylo_arr[k + 1];
        slong yn_r = n - ylor;

        slong s           = ydeg + yn_W - 1;
        slong qneg_ks_len = (xn - 1) * s + ydeg;
        slong W_ks_len    = 2 * (xnh - 1) * s + yn_W;
        slong yslice_lo   = ylo - ylor;

        /*
         * We only need product coefficients [yslice_lo, yslice_lo+R_rks_len).
         * _fmpz_poly_mulmid(res, f, flen, g, glen, lo, hi) computes slice [lo, hi)
         * of the full product into res[0..hi-lo-1], so res[k] = product coefficient lo+k.
         */
        slong R_rks_len = (xn - 1) * s + yn_r - yslice_lo;

        fmpz *Qk_num      = Qchain_num + qoff_arr[k];
        fmpz *Qneg_ks_num = _fmpz_vec_init(qneg_ks_len);
        fmpz *W_ks_num    = _fmpz_vec_init(W_ks_len);
        fmpz *R_rks_num   = _fmpz_vec_init(R_rks_len);
        fmpz_t Qneg_den, W_ks_den, R_rks_den;
        fmpz_init_set(Qneg_den, Qchain_den + k);
        fmpz_init_set(W_ks_den, W_den);
        fmpz_init(R_rks_den);

        /* Pack Qneg_k: x^i y^j -> i*s + j, odd-x entries negated. */
        for (slong j = 0; j < ydeg; j++)
            for (slong i = 0; i < xn; i++)
            {
                const fmpz *src = Qk_num + j * xn + i;
                if (i & 1)
                    fmpz_neg(Qneg_ks_num + i * s + j, src);
                else
                    fmpz_set(Qneg_ks_num + i * s + j, src);
            }

        /* Pack W(x^2): x^{2i} y^j -> 2*i*s + j. */
        for (slong j = 0; j < yn_W; j++)
            for (slong i = 0; i < xnh; i++)
                fmpz_set(W_ks_num + 2 * i * s + j, W_num + j * xnh + i);

        fmpz_mul(R_rks_den, Qneg_den, W_ks_den);

        _fmpz_poly_mulmid(R_rks_num, Qneg_ks_num, qneg_ks_len,
                                     W_ks_num,    W_ks_len,
                                     yslice_lo, yslice_lo + R_rks_len);

        _fmpz_vec_clear(Qneg_ks_num, qneg_ks_len);
        _fmpz_vec_clear(W_ks_num,    W_ks_len);

        /*
         * Unpack y-slice [yslice_lo, yn_r) into W in place.
         * R_rks_num[k] holds product coefficient yslice_lo + k.
         * The entry R[jsrc][i] is product coefficient i*s + jsrc,
         * so it sits at R_rks_num[i*s + jsrc - yslice_lo].
         */
        slong new_yn_W = n - ylo;
        {
            slong jout = 0;
            for (slong jsrc = yslice_lo; jsrc < yn_r; jsrc++, jout++)
                for (slong i = 0; i < xn; i++)
                    fmpz_set(W_num + jout * xn + i,
                             R_rks_num + i * s + jsrc - yslice_lo);
        }

        _fmpz_vec_clear(R_rks_num, R_rks_len);

        fmpz_set(W_den, R_rks_den);
        _fmpq_poly_canonicalise(W_num, W_den, new_yn_W * xn);

        fmpz_clear(Qneg_den);
        fmpz_clear(W_ks_den);
        fmpz_clear(R_rks_den);

        yn_W = new_yn_W;
    }

    /* W_num[0..n) / W_den is the result f(g(x)) mod x^n. */
    _fmpz_vec_set(res, W_num, n);
    fmpz_set(rden, W_den);
    /* Already canonical from the last upward step. */

    _fmpz_vec_clear(W_num, W_alloc);
    fmpz_clear(W_den);
    _fmpz_vec_clear(Qchain_num, qchain_total);
    _fmpz_vec_clear(Qchain_den, L + 1);

    flint_free(xn_arr);
    flint_free(ylo_arr);
    flint_free(ydeg_arr);
    flint_free(qoff_arr);
}

void
fmpq_poly_compose_series_kinoshita_li(fmpq_poly_t res,
                                       const fmpq_poly_t poly1,
                                       const fmpq_poly_t poly2,
                                       slong n)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;

    if (len1 == 0 || n == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (len2 != 0 && !fmpz_is_zero(poly2->coeffs))
    {
        flint_throw(FLINT_ERROR, "(fmpq_poly_compose_series_kinoshita_li): "
                "Inner polynomial must have zero constant term.\n");
    }

    if (len2 == 0 || len1 == 1)
    {
        fmpq_poly_set(res, poly1);
        fmpq_poly_truncate(res, 1);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    fmpq_poly_fit_length(res, lenr);
    _fmpq_poly_compose_series_kinoshita_li(res->coeffs, res->den,
                       poly1->coeffs, poly1->den, len1,
                       poly2->coeffs, poly2->den, len2, lenr);
    _fmpq_poly_set_length(res, lenr);
    _fmpq_poly_normalise(res);
}

