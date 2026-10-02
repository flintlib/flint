/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_ec.h"
#include "impl.h"

/*
    Arithmetic on the l-torsion of y^2 = f(x) = x^3 + a4 x + a6, carried
    out in F_q[x]/(m) for a modulus m that is the l-division polynomial or
    a factor of it.

    A point is written (u, v y) with u, v in F_q[x]/(m), which is closed
    under the group law because y^2 = f reduces every even power of y. The
    ring is not a field when m is reducible, so an inversion can fail; when
    it does, the failed extended gcd has handed over a proper factor of m,
    and _gr_ec_tors_inv records it in the context for the caller to carry
    on with. Schoof's algorithm uses this with m = psi_l, and SEA will use
    it with m the kernel polynomial of an l-isogeny, which is what this
    file is kept separate for.
*/

void
_gr_ec_tors_init(gr_ec_tors_struct * T, gr_ec_tors_ctx_struct * C)
{
    gr_poly_init(T->u, C->R);
    gr_poly_init(T->v, C->R);
    T->is_inf = 1;
}

void
_gr_ec_tors_clear(gr_ec_tors_struct * T, gr_ec_tors_ctx_struct * C)
{
    gr_poly_clear(T->u, C->R);
    gr_poly_clear(T->v, C->R);
}

int
_gr_ec_tors_set(gr_ec_tors_struct * D, const gr_ec_tors_struct * S, gr_ec_tors_ctx_struct * C)
{
    int status = GR_SUCCESS;
    status |= gr_poly_set(D->u, S->u, C->R);
    status |= gr_poly_set(D->v, S->v, C->R);
    D->is_inf = S->is_inf;
    return status;
}

int
_gr_ec_tors_mulmod(gr_poly_t res, const gr_poly_t a, const gr_poly_t b, gr_ec_tors_ctx_struct * C)
{
    return gr_poly_preinv_mulmod(res, a, b, C->P, C->R);
}

/* 1/a in R, or GR_UNABLE when psi_l turns out to be reducible here */
int
_gr_ec_tors_inv(gr_poly_t res, const gr_poly_t a, gr_ec_tors_ctx_struct * C)
{
    gr_poly_t g, s, t;
    int status = GR_SUCCESS;

    gr_poly_init(g, C->R);
    gr_poly_init(s, C->R);
    gr_poly_init(t, C->R);

    status |= gr_poly_xgcd(g, s, t, a, C->psi, C->R);

    if (status == GR_SUCCESS)
    {
        if (gr_poly_length(g, C->R) != 1)
        {
            /*
                The xgcd has handed us a nontrivial factor of the modulus.
                phi^2 - t phi + q vanishes on the whole l-torsion, so it
                still vanishes modulo this factor, and the caller can start
                again with the smaller modulus. (Recognising that factor as
                a kernel polynomial is exactly what turns Schoof into SEA.)
            */
            if (gr_poly_length(g, C->R) < gr_poly_length(C->psi, C->R)
                    && !C->have_factor)
            {
                if (gr_poly_set(C->factor, g, C->R) == GR_SUCCESS)
                    C->have_factor = 1;
            }

            status = GR_UNABLE;
        }
        else
            status = gr_poly_preinv_rem(res, s, C->P, C->R);
    }

    gr_poly_clear(g, C->R);
    gr_poly_clear(s, C->R);
    gr_poly_clear(t, C->R);

    return status;
}

truth_t
_gr_ec_tors_equal(const gr_ec_tors_struct * A, const gr_ec_tors_struct * B, gr_ec_tors_ctx_struct * C)
{
    if (A->is_inf || B->is_inf)
        return (A->is_inf && B->is_inf) ? T_TRUE : T_FALSE;

    if (gr_poly_equal(A->u, B->u, C->R) != T_TRUE)
        return T_FALSE;

    return gr_poly_equal(A->v, B->v, C->R);
}

int
_gr_ec_tors_neg(gr_ec_tors_struct * D, const gr_ec_tors_struct * S, gr_ec_tors_ctx_struct * C)
{
    int status = GR_SUCCESS;
    status |= gr_poly_set(D->u, S->u, C->R);
    status |= gr_poly_neg(D->v, S->v, C->R);
    D->is_inf = S->is_inf;
    return status;
}

/*
    The chord and tangent law written for (u, v y) with y^2 = f:

      distinct: w = (v2 - v1)/(u2 - u1),  u3 = w^2 f - u1 - u2,
                v3 = w (u1 - u3) - v1
      doubling: w = (3 u^2 + a4) / (2 v f), u3 = w^2 f - 2 u,
                v3 = w (u - u3) - v
*/
int
_gr_ec_tors_add(gr_ec_tors_struct * D, const gr_ec_tors_struct * A, const gr_ec_tors_struct * B,
        gr_ec_tors_ctx_struct * C)
{
    gr_poly_t w, num, den, u3, v3, tmp;
    int status = GR_SUCCESS;

    if (A->is_inf)
        return _gr_ec_tors_set(D, B, C);

    if (B->is_inf)
        return _gr_ec_tors_set(D, A, C);

    gr_poly_init(w, C->R);
    gr_poly_init(num, C->R);
    gr_poly_init(den, C->R);
    gr_poly_init(u3, C->R);
    gr_poly_init(v3, C->R);
    gr_poly_init(tmp, C->R);

    if (gr_poly_equal(A->u, B->u, C->R) == T_TRUE)
    {
        if (gr_poly_equal(A->v, B->v, C->R) != T_TRUE)
        {
            /* B = -A */
            D->is_inf = 1;
            goto cleanup;
        }

        /* doubling */
        status |= _gr_ec_tors_mulmod(tmp, A->u, A->u, C);
        status |= gr_poly_mul_ui(num, tmp, 3, C->R);
        status |= gr_poly_add(num, num, C->a4, C->R);

        status |= _gr_ec_tors_mulmod(den, A->v, C->f, C);
        status |= gr_poly_mul_ui(den, den, 2, C->R);
    }
    else
    {
        status |= gr_poly_sub(num, B->v, A->v, C->R);
        status |= gr_poly_sub(den, B->u, A->u, C->R);
    }

    status |= _gr_ec_tors_inv(tmp, den, C);

    if (status != GR_SUCCESS)
        goto cleanup;

    status |= _gr_ec_tors_mulmod(w, num, tmp, C);

    /* u3 = w^2 f - u1 - u2 */
    status |= _gr_ec_tors_mulmod(u3, w, w, C);
    status |= _gr_ec_tors_mulmod(u3, u3, C->f, C);
    status |= gr_poly_sub(u3, u3, A->u, C->R);
    status |= gr_poly_sub(u3, u3, B->u, C->R);

    /* v3 = w (u1 - u3) - v1 */
    status |= gr_poly_sub(tmp, A->u, u3, C->R);
    status |= _gr_ec_tors_mulmod(v3, w, tmp, C);
    status |= gr_poly_sub(v3, v3, A->v, C->R);

    if (status == GR_SUCCESS)
    {
        status |= gr_poly_set(D->u, u3, C->R);
        status |= gr_poly_set(D->v, v3, C->R);
        D->is_inf = 0;
    }

cleanup:
    gr_poly_clear(w, C->R);
    gr_poly_clear(num, C->R);
    gr_poly_clear(den, C->R);
    gr_poly_clear(u3, C->R);
    gr_poly_clear(v3, C->R);
    gr_poly_clear(tmp, C->R);

    return status;
}

/* k P by a left-to-right binary ladder; k is small (below l) */
int
_gr_ec_tors_mul_ui(gr_ec_tors_struct * D, const gr_ec_tors_struct * S, ulong k, gr_ec_tors_ctx_struct * C)
{
    gr_ec_tors_struct acc;
    int status = GR_SUCCESS;
    slong i;

    _gr_ec_tors_init(&acc, C);

    if (k != 0)
    {
        for (i = FLINT_BIT_COUNT(k) - 1; i >= 0 && status == GR_SUCCESS; i--)
        {
            status |= _gr_ec_tors_add(&acc, &acc, &acc, C);

            if ((k >> i) & 1)
                status |= _gr_ec_tors_add(&acc, &acc, S, C);
        }
    }

    if (status == GR_SUCCESS)
        status = _gr_ec_tors_set(D, &acc, C);

    _gr_ec_tors_clear(&acc, C);

    return status;
}

/*
    The k in [0, l) with target = k base, if there is one.

    Sets *found to 1 and *k to it, or *found to 0 if no k in the range
    works; returns the status of the arithmetic, which is GR_UNABLE when an
    inversion ran into a factor of the modulus (recorded in the context, as
    for every operation here).

    Trying each k in turn would cost a torsion addition per candidate, and
    a torsion addition costs an inversion in F_q[x]/(m) -- an extended gcd
    on polynomials of degree up to (l^2-1)/2, some thirty times a
    multiplication there and by a wide margin the most expensive thing in
    this arithmetic. Baby-step giant-step needs O(sqrt(l)) additions
    instead: write k = i m + j, tabulate j base for j < m, and walk target
    down by m base. Comparing against the table is polynomial equality,
    linear and nearly free, so its m^2 comparisons cost nothing next to the
    additions they replace.

    Schoof asks this for t mod l, with base = phi(P); an Elkies prime asks
    it for the eigenvalue of Frobenius on the kernel, with the same base.
*/
int
_gr_ec_tors_discrete_log(ulong * k, int * found,
        const gr_ec_tors_struct * target, const gr_ec_tors_struct * base,
        ulong l, gr_ec_tors_ctx_struct * C)
{
    slong m = (slong) n_sqrt(l) + 1;
    gr_ec_tors_struct * baby;
    gr_ec_tors_struct step, cur;
    slong i, j;
    int status = GR_SUCCESS;

    *found = 0;

    baby = flint_malloc(m * sizeof(gr_ec_tors_struct));

    for (j = 0; j < m; j++)
        _gr_ec_tors_init(baby + j, C);

    _gr_ec_tors_init(&step, C);
    _gr_ec_tors_init(&cur, C);

    /* baby[j] = j base, starting from the point at infinity */
    for (j = 1; j < m && status == GR_SUCCESS; j++)
        status |= _gr_ec_tors_add(baby + j, baby + j - 1, base, C);

    /* the giant step is -m base */
    status |= _gr_ec_tors_add(&step, baby + m - 1, base, C);
    status |= _gr_ec_tors_neg(&step, &step, C);

    status |= _gr_ec_tors_set(&cur, target, C);

    for (i = 0; i * m < (slong) l && !*found && status == GR_SUCCESS; i++)
    {
        for (j = 0; j < m; j++)
            if (i * m + j < (slong) l
                    && _gr_ec_tors_equal(&cur, baby + j, C) == T_TRUE)
            {
                *k = (ulong) (i * m + j);
                *found = 1;
                break;
            }

        if (!*found)
            status |= _gr_ec_tors_add(&cur, &cur, &step, C);
    }

    _gr_ec_tors_clear(&cur, C);
    _gr_ec_tors_clear(&step, C);

    for (j = 0; j < m; j++)
        _gr_ec_tors_clear(baby + j, C);

    flint_free(baby);

    return status;
}

/*
    F_q[x]/(m): the modulus preconditioned, f = x^3 + a4 x + a6 reduced
    modulo it, and a4 as a constant, which is what the group law reads.
*/
int
_gr_ec_tors_ctx_init(gr_ec_tors_ctx_struct * C, const gr_poly_t m,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    C->R = R;
    C->have_factor = 0;
    gr_poly_init(C->psi, R);
    gr_poly_init(C->factor, R);
    gr_poly_init(C->f, R);
    gr_poly_init(C->a4, R);
    gr_poly_preinv_init(C->P, R);

    status |= gr_poly_set(C->psi, m, R);

    if (status != GR_SUCCESS || gr_poly_length(C->psi, R) < 2)
        return GR_UNABLE;

    status |= gr_poly_preinv_set(C->P, C->psi, R);

    status |= gr_poly_zero(C->f, R);
    status |= gr_poly_set_coeff_si(C->f, 3, 1, R);
    status |= gr_poly_set_coeff_scalar(C->f, 1, GR_EC_A4(ctx), R);
    status |= gr_poly_set_coeff_scalar(C->f, 0, GR_EC_A6(ctx), R);
    status |= gr_poly_preinv_rem(C->f, C->f, C->P, R);

    status |= gr_poly_zero(C->a4, R);
    status |= gr_poly_set_coeff_scalar(C->a4, 0, GR_EC_A4(ctx), R);

    return status;
}

void
_gr_ec_tors_ctx_clear(gr_ec_tors_ctx_struct * C)
{
    gr_poly_clear(C->factor, C->R);
    gr_poly_preinv_clear(C->P, C->R);
    gr_poly_clear(C->a4, C->R);
    gr_poly_clear(C->f, C->R);
    gr_poly_clear(C->psi, C->R);
}

/*
    The generic point P = (x, y) and its Frobenius image
    phi(P) = (x^q, f^((q-1)/2) y), both modulo m.

    The reduction of x matters: once a descent has taken the modulus down
    to a linear factor, x is no longer reduced, and two points with the
    same x-coordinate would otherwise compare unequal and send the group
    law down the "distinct" branch, where it would invert zero.
*/
int
_gr_ec_tors_frobenius(gr_ec_tors_struct * P, gr_ec_tors_struct * phiP,
        const fmpz_t q, gr_ec_tors_ctx_struct * C)
{
    gr_ctx_struct * R = C->R;
    fmpz_t e;
    int status = GR_SUCCESS;

    fmpz_init(e);

    status |= gr_poly_zero(P->u, R);
    status |= gr_poly_set_coeff_si(P->u, 1, 1, R);
    status |= gr_poly_preinv_rem(P->u, P->u, C->P, R);
    status |= gr_poly_one(P->v, R);
    status |= gr_poly_preinv_rem(P->v, P->v, C->P, R);
    P->is_inf = 0;

    status |= gr_poly_preinv_powmod_x_fmpz(phiP->u, q, C->P, R);

    fmpz_sub_ui(e, q, 1);
    fmpz_fdiv_q_ui(e, e, 2);
    status |= gr_poly_preinv_powmod_fmpz_binexp(phiP->v, C->f, e, C->P, R);
    phiP->is_inf = 0;

    fmpz_clear(e);

    return status;
}

/*
    Run step modulo m, and whenever an inversion runs into a proper factor
    of the modulus, start again modulo that factor.

    This is sound for every step that only asks a question about points of
    the subgroup m describes: a relation like phi^2 - t phi + q = 0, or
    phi(P) = lambda P, holds on every point of the subgroup, so it holds
    modulo every factor of m, and the smaller ring gives the same answer.
*/
int
_gr_ec_tors_solve(ulong * res, const gr_poly_t m, ulong l, const fmpz_t q,
        _gr_ec_tors_step_t step, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ec_tors_ctx_struct C;
    gr_poly_t mod;
    slong depth;
    int status = GR_SUCCESS;

    gr_poly_init(mod, R);
    status = gr_poly_set(mod, m, R);

    for (depth = 0; status == GR_SUCCESS && depth < 16; depth++)
    {
        int again;

        status = _gr_ec_tors_ctx_init(&C, mod, ctx);

        if (status == GR_SUCCESS)
            status = step(res, l, q, &C, ctx);

        again = (status != GR_SUCCESS && C.have_factor
                && gr_poly_set(mod, C.factor, R) == GR_SUCCESS);

        _gr_ec_tors_ctx_clear(&C);

        if (!again)
            break;

        status = GR_SUCCESS;
    }

    gr_poly_clear(mod, R);

    return status;
}
