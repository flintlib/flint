/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"

/*
    Simultaneously updates two strided vectors:

        (u[i], v[i]) <- (M[0]*u[i] + M[1]*v[i], M[2]*u[i] + M[3]*v[i])

    for i = 0, ..., len - 1, where u[i] and v[i] denote the elements at
    offset i*stride.  With M a Givens rotation this is BLAS's rot; the
    general 2x2 form costs nothing extra and also serves Schur
    reordering, bulge chasing and non-unitary transforms.

    Aliasing.  Each new pair is computed from the old pair, so u[i] and
    v[i] are both inputs and outputs of the same step.  This generic
    implementation never relies on _gr_vec_dot tolerating output/input
    overlap: both dot products go to temporaries which are then swapped
    into place.  That costs one extra swap per element over a kernel
    that can write the second result directly.

    Overriding.  A ring may override this method purely for speed.  It
    is safe to write results directly into u[i], v[i] only if the
    ring's dot product (or whatever the override uses) consumes all
    its inputs before writing its output for all inputs (acf_approx_dot
    does not: its simple version, for short vectors or exponents that
    are not small, writes as it goes).  A machine-float backend would
    implement this as a plain SIMD loop over i.
*/
int
_gr_vec_apply_2x2_generic(gr_ptr u, gr_ptr v, slong stride, slong len, gr_srcptr M, gr_ctx_t ctx)
{
    slong i;
    slong sz = ctx->sizeof_elem;
    gr_ptr t1, t2, w;
    gr_ptr ui, vi;
    gr_method_vec_dot_op dot;
    gr_method_void_unary_op set_shallow;
    gr_method_swap_op swap;
    int status = GR_SUCCESS;

    if (len <= 0)
        return GR_SUCCESS;

    /* Hoist the method lookups out of the loop; this is a hot inner
       kernel and the dispatch is otherwise repeated per element. */
    dot = GR_VEC_DOT_OP(ctx, VEC_DOT);
    set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);
    swap = GR_SWAP_OP(ctx, SWAP);

    GR_TMP_INIT2(t1, t2, ctx);
    w = GR_TMP_ALLOC(2 * sz);

    ui = u;
    vi = v;

    for (i = 0; i < len; i++)
    {
        /* w is a shallow length-2 vector holding (u[i], v[i]); this is
           needed because u[i] and v[i] are in general not adjacent. */
        set_shallow(w, ui, ctx);
        set_shallow(GR_ENTRY(w, 1, sz), vi, ctx);

        status |= dot(t1, NULL, 0, M, w, 2, ctx);
        status |= dot(t2, NULL, 0, GR_ENTRY(M, 2, sz), w, 2, ctx);

        swap(ui, t1, ctx);
        swap(vi, t2, ctx);

        ui = GR_ENTRY(ui, stride, sz);
        vi = GR_ENTRY(vi, stride, sz);
    }

    GR_TMP_FREE(w, 2 * sz);
    GR_TMP_CLEAR2(t1, t2, ctx);

    return status;
}
