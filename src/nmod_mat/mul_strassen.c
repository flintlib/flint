/*
    Copyright (C) 2008, Martin Albrecht
    Copyright (C) 2008, 2009 William Hart.
    Copyright (C) 2010, Fredrik Johansson
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_mat/impl.h"

/*
    Odd dimensions: virtual padding.

    Each dimension d is split as h = ceil(d/2) (top / left blocks) plus
    h - (d mod 2) (bottom / right blocks), and the short blocks are seen as
    padded with a zero row or column, so that the Strassen-Winograd formulas
    apply unchanged to the padded matrices. Nothing is copied: the padding
    only shows in the following places.

    - The temporaries X1, X2 have the full padded size, and the sums of
      blocks that form them are "padded" additions: an ordinary addition on
      the common part, the extra row or column of the larger operand copied
      (or negated), zeros beyond.
    - A product whose left or right operand is a short block of A or B uses
      the real inner dimension or the real rows / columns.
    - The blocks C12, C21, C22 of C lack the padding row or column, while the
      schedule below stores intermediate products in them; following it, no
      missing row or column is ever needed for a real entry of the result:
        . odd n: no result of the left block column (C11, C21) reads C12 or
          C22 (final C11 = P5 + P7, C21 = P1 - P3 - P5 - P6), so the missing
          last column of what they hold is never used;
        . odd m: the last rows of C21 and C22 hold the last rows of P1, then
          of P1 - P3 - P5 (in C11), and P6 (in C21), which only feed the
          missing last rows of the final C21 and C22; the one exception is
          P2 = (A22 - A21)(B22 - B21), held in C22 and added to C12, whose
          last row is zero since both A21 and A22 are padded with a zero
          row;
        . odd k: only the temporaries and the inner dimension of
          P5 = A12 B21 are concerned.
    So odd dimensions cost nothing beyond the padded additions, which
    touch the same entries as the ordinary ones.
*/

/*
    X = P + Q (sub = 0) or X = P - Q (sub = 1), where P and Q have at most
    the dimensions of X and stand for X-sized matrices padded with zeros.
    P may be X itself.
*/
static void
_padded_add(nmod_mat_t X, const nmod_mat_t P, const nmod_mat_t Q, int sub)
{
    slong i, pr, qr, lo, hi;
    nmod_t mod = X->mod;

    for (i = 0; i < X->r; i++)
    {
        nn_ptr x = nmod_mat_entry_ptr(X, i, 0);

        pr = (i < P->r) ? P->c : 0;
        qr = (i < Q->r) ? Q->c : 0;
        lo = FLINT_MIN(pr, qr);
        hi = FLINT_MAX(pr, qr);

        if (sub)
            _nmod_vec_sub(x, nmod_mat_entry_ptr(P, i, 0),
                             nmod_mat_entry_ptr(Q, i, 0), lo, mod);
        else
            _nmod_vec_add(x, nmod_mat_entry_ptr(P, i, 0),
                             nmod_mat_entry_ptr(Q, i, 0), lo, mod);

        if (pr > qr)
        {
            if (P != X)
                _nmod_vec_set(x + lo, nmod_mat_entry_ptr(P, i, lo), hi - lo);
        }
        else if (qr > pr)
        {
            if (sub)
                _nmod_vec_neg(x + lo, nmod_mat_entry_ptr(Q, i, lo), hi - lo,
                              mod);
            else
                _nmod_vec_set(x + lo, nmod_mat_entry_ptr(Q, i, lo), hi - lo);
        }

        if (hi < X->c)
            _nmod_vec_zero(x + hi, X->c - hi);
    }
}

static void _nmod_mat_mul_strassen_rec(nmod_mat_t C, const nmod_mat_t A,
                                       const nmod_mat_t B, slong cutoff);

/*
    The products of one level: with cutoff > 0, those with all dimensions
    at least cutoff use Strassen again (for the tests); otherwise they go
    back to nmod_mat_mul, which chooses.
*/
static void
_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, slong cutoff)
{
    if (cutoff > 0 && A->r >= cutoff && A->c >= cutoff && B->c >= cutoff)
        _nmod_mat_mul_strassen_rec(C, A, B, cutoff);
    else
        nmod_mat_mul(C, A, B);
}

/* The implemented sequence is not Strassen's nor Winograd's, but the sequence
   proposed by Bodrato, which is equivalent to Winograd's, and can be easily
   adapted to compute the square of a matrix. */

static void
_nmod_mat_mul_strassen_rec(nmod_mat_t C, const nmod_mat_t A,
                           const nmod_mat_t B, slong cutoff)
{
    slong a, b, c;
    slong ha, hb, hc;       /* top rows / left columns (padded sizes) */
    slong la, lc;           /* bottom rows / right columns (real sizes) */

    nmod_mat_t A11, A12, A21, A22;
    nmod_mat_t B11, B12, B21, B22;
    nmod_mat_t C11, C12, C21, C22;
    nmod_mat_t X1, X2, W1, W2, W3;

    a = A->r;
    b = A->c;
    c = B->c;

    if (a <= 4 || b <= 4 || c <= 4)
    {
        nmod_mat_mul(C, A, B);
        return;
    }

    ha = (a + 1) / 2;  la = a - ha;
    hb = (b + 1) / 2;
    hc = (c + 1) / 2;  lc = c - hc;

    nmod_mat_window_init(A11, A, 0, 0, ha, hb);
    nmod_mat_window_init(A12, A, 0, hb, ha, b);
    nmod_mat_window_init(A21, A, ha, 0, a, hb);
    nmod_mat_window_init(A22, A, ha, hb, a, b);

    nmod_mat_window_init(B11, B, 0, 0, hb, hc);
    nmod_mat_window_init(B12, B, 0, hc, hb, c);
    nmod_mat_window_init(B21, B, hb, 0, b, hc);
    nmod_mat_window_init(B22, B, hb, hc, b, c);

    nmod_mat_window_init(C11, C, 0, 0, ha, hc);
    nmod_mat_window_init(C12, C, 0, hc, ha, c);
    nmod_mat_window_init(C21, C, ha, 0, a, hc);
    nmod_mat_window_init(C22, C, ha, hc, a, c);

    nmod_mat_init(X1, ha, FLINT_MAX(hb, hc), A->mod.n);
    nmod_mat_init(X2, hb, hc, A->mod.n);

    X1->c = hb;

    /*
        See Jean-Guillaume Dumas, Clement Pernet, Wei Zhou; "Memory
        efficient scheduling of Strassen-Winograd's matrix multiplication
        algorithm"; https://arxiv.org/pdf/0707.2347v3 for reference on the
        used operation scheduling.
    */

    /* P1 = (A12 + A22)(B12 + B22) -> C21 (its last row is not needed) */
    _padded_add(X1, A22, A12, 0);
    _padded_add(X2, B22, B12, 0);
    nmod_mat_window_init(W1, X1, 0, 0, la, hb);
    _mul(C21, W1, X2, cutoff);
    nmod_mat_window_clear(W1);

    /* P2 = (A22 - A21)(B22 - B21) -> C22, whose last row is zero and last
       column not needed */
    _padded_add(X1, A22, A21, 1);
    _padded_add(X2, B22, B21, 1);
    nmod_mat_window_init(W1, X1, 0, 0, la, hb);
    nmod_mat_window_init(W2, X2, 0, 0, hb, lc);
    _mul(C22, W1, W2, cutoff);
    nmod_mat_window_clear(W1);
    nmod_mat_window_clear(W2);

    /* P3 = (A12 + A22 - A21)(B12 + B22 - B21) -> C11 */
    _padded_add(X1, X1, A12, 0);
    _padded_add(X2, X2, B12, 0);
    _mul(C11, X1, X2, cutoff);

    /* P4 = (A12 + A22 - A21 - A11) B12 -> C12 */
    nmod_mat_sub(X1, X1, A11);
    _mul(C12, X1, B12, cutoff);

    /* P5 = A12 B21 -> X1 (inner dimension lb) */
    X1->c = hc;
    _mul(X1, A12, B21, cutoff);

    /* C11 = P3 + P5 */
    nmod_mat_add(C11, C11, X1);

    /* C12 = P4 + P2: the last row of P2 (if a is odd) is zero */
    nmod_mat_window_init(W1, C12, 0, 0, la, lc);
    nmod_mat_add(W1, W1, C22);
    nmod_mat_window_clear(W1);

    /* C12 = P3 + P5 - P4 - P2 (final) */
    nmod_mat_window_init(W1, C11, 0, 0, ha, lc);
    nmod_mat_sub(C12, W1, C12);
    nmod_mat_window_clear(W1);

    /* C11 = P1 - P3 - P5, only needed in its first la rows */
    nmod_mat_window_init(W1, C11, 0, 0, la, hc);
    nmod_mat_sub(W1, C21, W1);
    nmod_mat_window_clear(W1);

    /* P6 = A21 (B12 + B22 - B21 - B11) -> C21 */
    nmod_mat_sub(X2, X2, B11);
    _mul(C21, A21, X2, cutoff);

    nmod_mat_clear(X2);

    /* C21 = P1 - P3 - P5 - P6 and C22 = P1 + P2 - P3 - P5 (final) */
    nmod_mat_window_init(W1, C11, 0, 0, la, hc);
    nmod_mat_sub(C21, W1, C21);
    nmod_mat_window_init(W3, C11, 0, 0, la, lc);
    nmod_mat_add(C22, C22, W3);
    nmod_mat_window_clear(W1);
    nmod_mat_window_clear(W3);

    /* C11 = P5 + P7, P7 = A11 B11 (final) */
    _mul(C11, A11, B11, cutoff);
    nmod_mat_add(C11, X1, C11);

    X1->c = FLINT_MAX(hb, hc);
    nmod_mat_clear(X1);

    nmod_mat_window_clear(A11);
    nmod_mat_window_clear(A12);
    nmod_mat_window_clear(A21);
    nmod_mat_window_clear(A22);

    nmod_mat_window_clear(B11);
    nmod_mat_window_clear(B12);
    nmod_mat_window_clear(B21);
    nmod_mat_window_clear(B22);

    nmod_mat_window_clear(C11);
    nmod_mat_window_clear(C12);
    nmod_mat_window_clear(C21);
    nmod_mat_window_clear(C22);
}

void
_nmod_mat_mul_strassen_cutoff(nmod_mat_t C, const nmod_mat_t A,
                              const nmod_mat_t B, slong cutoff)
{
    _nmod_mat_mul_strassen_rec(C, A, B, cutoff);
}

void
nmod_mat_mul_strassen(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    _nmod_mat_mul_strassen_rec(C, A, B, 0);
}
