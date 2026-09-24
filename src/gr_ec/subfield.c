/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Counting a curve over F_{p^n} that is really a curve over F_p.

    A curve whose a-invariants all lie in the prime field is the base
    change of a curve over F_p, and the two have the same Frobenius: if
    alpha and beta are the roots of X^2 - t X + p for the curve over F_p,
    then over F_{p^n} the trace is alpha^n + beta^n. Those satisfy the
    linear recurrence

        t_0 = 2,  t_1 = t,  t_k = t t_{k-1} - p t_{k-2},

    so the whole computation over the extension is a point count over F_p
    and n - 1 multiplications of integers the size of p.

    The saving is not small. Counting over F_{p^n} directly costs
    baby-step giant-step at O(q^(1/4)) group operations in F_q, or Schoof
    at a polynomial in log q; counting over F_p costs the same thing in a
    field n times smaller, and then the lift is free. Measured against
    PARI, which does this too, the difference between having it and not
    having it is three to four orders of magnitude by n = 5.

    Only descent to the prime field is attempted. A curve defined over an
    intermediate F_{p^d} with 1 < d < n would want the same treatment, but
    mapping its coefficients down needs an embedding of F_{p^d} in F_{p^n}
    that gr does not currently hand out, whereas the prime field needs
    nothing but gr_get_fmpz.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_ec.h"
#include "impl.h"

/*
    t_m = alpha^m + beta^m for the roots of X^2 - t X + p, by the
    recurrence above.
*/
static void
_trace_of_extension(fmpz_t res, const fmpz_t t, const fmpz_t p, slong m)
{
    fmpz_t a, b, c;
    slong k;

    if (m == 0)
    {
        fmpz_set_ui(res, 2);
        return;
    }

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);

    fmpz_set_ui(a, 2);          /* t_{k-2} */
    fmpz_set(b, t);             /* t_{k-1} */

    for (k = 2; k <= m; k++)
    {
        fmpz_mul(c, t, b);
        fmpz_submul(c, p, a);
        fmpz_swap(a, b);
        fmpz_swap(b, c);
    }

    fmpz_set(res, b);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
}

int
gr_ec_ctx_cardinality_subfield(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ctx_t R1;
    gr_ec_ctx_t E1;
    fmpz_t p, q, t, N1;
    fmpz a[5];
    gr_ptr c;
    slong deg, i;
    int status = GR_SUCCESS, have_ring = 0, have_curve = 0;

    if (gr_ctx_fq_degree(&deg, R) != GR_SUCCESS || deg <= 1)
        return GR_UNABLE;

    fmpz_init(p);
    fmpz_init(q);
    fmpz_init(t);
    fmpz_init(N1);

    for (i = 0; i < 5; i++)
        fmpz_init(a + i);

    if (gr_ctx_fq_prime(p, R) != GR_SUCCESS
            || gr_ctx_cardinality_fmpz(q, R) != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /*
        Every a-invariant has to be an element of the prime field. This is
        the whole test: gr_get_fmpz on a finite field element succeeds
        exactly when it is one, so nothing here needs a Frobenius.
    */
    for (i = 0; i < 5; i++)
        if (gr_get_fmpz(a + i, GR_EC_COEFF(ctx, i), R) != GR_SUCCESS)
        {
            status = GR_UNABLE;
            goto cleanup;
        }

    /* the same curve over F_p */
    if (fmpz_abs_fits_ui(p))
    {
        if (gr_ctx_init_nmod(R1, fmpz_get_ui(p)) != GR_SUCCESS)
        {
            status = GR_UNABLE;
            goto cleanup;
        }
    }
    else
        gr_ctx_init_fmpz_mod(R1, p);

    have_ring = 1;

    /* p came from the base ring, which is a field, so this is a fact */
    GR_IGNORE(gr_ctx_set_is_field(R1, T_TRUE));

    {
        gr_ptr v;
        GR_TMP_INIT_VEC(v, 5, R1);

        for (i = 0; i < 5 && status == GR_SUCCESS; i++)
        {
            c = GR_ENTRY(v, i, R1->sizeof_elem);
            status = gr_set_fmpz(c, a + i, R1);
        }

        if (status == GR_SUCCESS)
            status = gr_ec_ctx_init(E1, R1,
                    GR_ENTRY(v, 0, R1->sizeof_elem),
                    GR_ENTRY(v, 1, R1->sizeof_elem),
                    GR_ENTRY(v, 2, R1->sizeof_elem),
                    GR_ENTRY(v, 3, R1->sizeof_elem),
                    GR_ENTRY(v, 4, R1->sizeof_elem));

        if (status == GR_SUCCESS)
            have_curve = 1;

        GR_TMP_CLEAR_VEC(v, 5, R1);
    }

    if (status != GR_SUCCESS)
        goto cleanup;

    /*
        The count over F_p goes through the ordinary dispatch, so the
        complex multiplication shortcut and everything else applies there.
    */
    status = gr_ec_ctx_cardinality(N1, E1);

    if (status != GR_SUCCESS)
        goto cleanup;

    /* t = p + 1 - #E(F_p), then lift */
    fmpz_add_ui(t, p, 1);
    fmpz_sub(t, t, N1);

    _trace_of_extension(t, t, p, deg);

    fmpz_add_ui(res, q, 1);
    fmpz_sub(res, res, t);

cleanup:
    if (have_curve)
        gr_ec_ctx_clear(E1);

    if (have_ring)
        gr_ctx_clear(R1);

    for (i = 0; i < 5; i++)
        fmpz_clear(a + i);

    fmpz_clear(p);
    fmpz_clear(q);
    fmpz_clear(t);
    fmpz_clear(N1);

    return status;
}
