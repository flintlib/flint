/*
    Copyright (C) 2025 Albin Ahlbäck
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_MAT_IMPL_H
#define ARB_MAT_IMPL_H

#include "arb_types.h"
#include "arb_mat.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

slong _arb_mat_exp_choose_N(const mag_t norm, slong prec);

/*
    Floating-point types for the approximate (midpoint) functions of
    arb_mat and acb_mat, which are the generic algorithms of gr_mat run
    over a floating-point ring.

    The midpoints are scaled by powers of two and converted to dfloat
    (N doubles, 53N bits) when the working precision plus some padding
    fits, the dynamic range of the entries is small enough and the
    matrix is large enough for dfloat to beat nfloat; otherwise to
    nfloat when the precision allows. If a type fails (an exponent out
    of range, a nonfinite result, a zero pivot, no convergence), the
    next one is tried, and finally arf (acf for complex matrices), which
    has no exponent limits and is not scaled.

    A type is identified by kind: 1 to 4 for dfloat with that many
    doubles, ARB_MAT_APPROX_NFLOAT for nfloat, ARB_MAT_APPROX_ARF for
    arf or acf, and ARB_MAT_APPROX_NONE after arf.
*/

#define ARB_MAT_APPROX_NONE 0
#define ARB_MAT_APPROX_NFLOAT -1
#define ARB_MAT_APPROX_ARF -2

/* operations, for the choice of backend */
#define ARB_MAT_APPROX_OP_LU 0      /* LU, solving, inverse */
#define ARB_MAT_APPROX_OP_TRI 1     /* triangular solving */
#define ARB_MAT_APPROX_OP_EIG 2     /* eigenvalues */

/* Precision padding: dN is used at prec bits when prec + pad <= 53N, since
   the effective accuracy of dN in a long computation is a few bits lower
   than that of arf at 53N bits. */
#define ARB_MAT_APPROX_DFLOAT_PAD(N) (2 * (N) - 2)

/* Largest dynamic range (difference of the largest and smallest exponent
   of the nonzero entries, after scaling) accepted for dfloat: products
   of two entries then keep all N limbs above the subnormal range, with
   some margin. */
#define ARB_MAT_APPROX_DFLOAT_RANGE(N) ((1022 - 53 * (N)) / 2 - 32)

/* Exponents larger than this in absolute value send the input to arf
   directly. */
#define ARB_MAT_APPROX_MAX_EXP (COEFF_MAX / 2)

/* The smallest dimension for which dfloat or nfloat is used. */
#define ARB_MAT_APPROX_MIN_N 4

/* Triangular solving with fewer right-hand sides is done with arf: the
   conversions cost as much as the solving. */
#define ARB_MAT_APPROX_TRI_MIN_RHS 4

int _arb_mat_approx_backend(slong n, slong prec, slong range, int cplx, int op);
int _arb_mat_approx_backend_next(int kind, slong prec);
int _arb_mat_approx_ctx_init(gr_ctx_t ctx, int kind, slong prec, int cplx);

/* Update the exponent range [*emin, *emax] (initially WORD_MAX, WORD_MIN)
   with the nonzero x; return 0 if x is not finite or its exponent is
   too large for the backends. */
int _arb_mat_approx_exp_update(slong * emin, slong * emax, const arf_t x);

/* A dynamic range from the exponent range: zero for no nonzero entries. */
slong _arb_mat_approx_range(slong emin, slong emax);

/* Conversions of one real component (N doubles, an nfloat or an arf)
   with scaling by 2^scale (always 0 for arf); get returns 0 if the value
   is not finite (never for arf). */
int _arb_mat_approx_set_arf(gr_ptr res, const arf_t x, slong scale, int kind, gr_ctx_t ctx);
int _arb_mat_approx_get_arf(arf_t res, gr_srcptr x, slong scale, int kind, gr_ctx_t ctx);
int _arb_mat_approx_is_finite(gr_srcptr x, int kind, gr_ctx_t ctx);

/* The imaginary component of a complex backend element. */
gr_ptr _arb_mat_approx_imag_ptr(gr_ptr x, int kind, gr_ctx_t ctx);

/* Parts of a matrix taking part in a triangular operation. */
#define ARB_MAT_APPROX_ALL 0
#define ARB_MAT_APPROX_LOWER 1          /* i >= j */
#define ARB_MAT_APPROX_STRICT_LOWER 2   /* i > j */
#define ARB_MAT_APPROX_UPPER 3          /* i <= j */
#define ARB_MAT_APPROX_STRICT_UPPER 4   /* i < j */

FLINT_FORCE_INLINE int
_arb_mat_approx_in_part(slong i, slong j, int part)
{
    switch (part)
    {
        case ARB_MAT_APPROX_LOWER: return i >= j;
        case ARB_MAT_APPROX_STRICT_LOWER: return i > j;
        case ARB_MAT_APPROX_UPPER: return i <= j;
        case ARB_MAT_APPROX_STRICT_UPPER: return i < j;
        default: return 1;
    }
}

#ifdef __cplusplus
}
#endif

#endif
