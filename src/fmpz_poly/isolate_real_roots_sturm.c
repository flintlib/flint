/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly/impl.h"
#include "arb.h"

/*
    Real root isolation by bisection using a Sturm sequence.

    The Sturm sequence p_0 = pol, p_1 = pol', ..., p_{n-1} is computed as a
    primitive PRS. We store p_0, p_1 and the pseudo-quotients: with
    lc(p_k)^d p_{k-1} = Q_k p_k + R and R = G_k p_{k+1}, we have

        p_{k+1} = (L_k p_{k-1} - Q_k p_k) / G_k,    L_k = lc(p_k)^d,

    so that the whole sequence can be evaluated at a point using O(n)
    operations. This is done in ball arithmetic, increasing the precision
    and finally falling back on exact arithmetic when some sign cannot be
    determined.

    This is much faster than VCA when the remainders stay small (for
    example, classical orthogonal polynomials, or sparse polynomials with
    close roots), but useless for generic input where the remainders grow;
    we therefore give up (returning 0) as soon as some remainder is too
    large, and the caller falls back on VCA.
*/

typedef struct
{
    slong n;                /* number of polynomials p_0, ..., p_{n-1} */
    slong nQ;               /* Q[1], ..., Q[nQ] have been initialised */
    fmpz_poly_t p0;
    fmpz_poly_t p1;
    fmpz_poly_struct * Q;
    fmpz * L;
    fmpz * G;
    int * sgn_lead;         /* sign of the leading coefficient of p_j */
    int * sgn_const;        /* sign of p_j(0) */
    slong * deg;            /* degree of p_j */
    slong p0_bits;
    /* exact ball copies of the data */
    arb_ptr p0a;
    arb_ptr p1a;
    arb_ptr * Qa;
    arb_ptr La;
    arb_ptr Ga;
}
sturm_seq_struct;

typedef sturm_seq_struct sturm_seq_t[1];

static void
sturm_seq_record(sturm_seq_t S, slong j, const fmpz_poly_t p)
{
    S->sgn_lead[j] = fmpz_sgn(p->coeffs + p->length - 1);
    S->sgn_const[j] = fmpz_sgn(p->coeffs);
    S->deg[j] = p->length - 1;
}

static void
sturm_seq_clear(sturm_seq_t S)
{
    slong j, len = S->p0->length;

    for (j = 1; j <= S->nQ; j++)
        fmpz_poly_clear(S->Q + j);

    if (S->Qa != NULL)
    {
        _arb_vec_clear(S->p0a, S->p0->length);
        _arb_vec_clear(S->p1a, S->p1->length);
        for (j = 1; j <= S->nQ; j++)
            _arb_vec_clear(S->Qa[j], S->Q[j].length);
        _arb_vec_clear(S->La, len);
        _arb_vec_clear(S->Ga, len);
        flint_free(S->Qa);
    }

    fmpz_poly_clear(S->p0);
    fmpz_poly_clear(S->p1);
    flint_free(S->Q);
    _fmpz_vec_clear(S->L, len);
    _fmpz_vec_clear(S->G, len);
    flint_free(S->sgn_lead);
    flint_free(S->sgn_const);
    flint_free(S->deg);
}

/* Assumes len >= 2 and pol[0] != 0. Returns 0 if some remainder
   exceeds max_size, or if pol is not squarefree. */
static int
sturm_seq_init(sturm_seq_t S, const fmpz * pol, slong len, slong max_size)
{
    fmpz_poly_t A, B, R;
    fmpz_t c;
    slong j, k, d, bits;
    int ok = 1;

    fmpz_poly_init2(S->p0, len);
    fmpz_poly_init(S->p1);
    S->Q = flint_malloc(len * sizeof(fmpz_poly_struct));
    S->L = _fmpz_vec_init(len);
    S->G = _fmpz_vec_init(len);
    S->sgn_lead = flint_malloc(len * sizeof(int));
    S->sgn_const = flint_malloc(len * sizeof(int));
    S->deg = flint_malloc(len * sizeof(slong));
    S->nQ = 0;
    S->Qa = NULL;

    /* Replacing pol by -pol does not change the sign variations, and the
       derivative of a polynomial with positive leading coefficient has
       positive leading coefficient, so the normalisation done by
       primitive_part is harmless. */
    _fmpz_poly_primitive_part(S->p0->coeffs, pol, len);
    _fmpz_poly_set_length(S->p0, len);
    fmpz_poly_derivative(S->p1, S->p0);
    fmpz_poly_primitive_part(S->p1, S->p1);

    sturm_seq_record(S, 0, S->p0);
    sturm_seq_record(S, 1, S->p1);
    S->n = 2;
    S->p0_bits = FLINT_ABS(_fmpz_vec_max_bits(S->p0->coeffs, len));

    fmpz_init(c);
    fmpz_poly_init(A);
    fmpz_poly_init(B);
    fmpz_poly_init(R);
    fmpz_poly_set(A, S->p0);
    fmpz_poly_set(B, S->p1);

    while (B->length > 1)
    {
        /* A = p_{k-1}, B = p_k */
        k = S->n - 1;
        d = A->length - B->length + 1;

        fmpz_poly_init(S->Q + k);
        S->nQ = k;
        fmpz_poly_pseudo_divrem_cohen(S->Q + k, R, A, B);

        if (fmpz_poly_is_zero(R))
        {
            /* not squarefree */
            ok = 0;
            break;
        }

        /* Cheap check whether we certainly exceed the bound: the content
           divides the gcd of any two coefficients. */
        bits = FLINT_ABS(_fmpz_vec_max_bits(R->coeffs, R->length));
        if ((double) bits * R->length > (double) max_size)
        {
            for (j = 0; fmpz_is_zero(R->coeffs + j); j++)
                ;
            fmpz_gcd(c, R->coeffs + j, R->coeffs + R->length - 1);
            if ((double) (bits - (slong) fmpz_bits(c)) * R->length > (double) max_size)
            {
                ok = 0;
                break;
            }
        }

        fmpz_pow_ui(S->L + k, B->coeffs + B->length - 1, d);

        /* p_{k+1} = R / c, a positive multiple of -rem(A, B) */
        _fmpz_poly_content(c, R->coeffs, R->length);
        if ((d % 2 == 0) || (fmpz_sgn(B->coeffs + B->length - 1) > 0))
            fmpz_neg(c, c);
        fmpz_set(S->G + k, c);
        _fmpz_vec_scalar_divexact_fmpz(R->coeffs, R->coeffs, R->length, c);

        sturm_seq_record(S, k + 1, R);
        S->n++;

        bits = FLINT_ABS(_fmpz_vec_max_bits(R->coeffs, R->length));
        if ((double) bits * R->length > (double) max_size)
        {
            ok = 0;
            break;
        }

        fmpz_poly_swap(A, B);
        fmpz_poly_swap(B, R);
    }

    fmpz_clear(c);
    fmpz_poly_clear(A);
    fmpz_poly_clear(B);
    fmpz_poly_clear(R);

    if (ok)
    {
        S->p0a = _arb_vec_init(S->p0->length);
        S->p1a = _arb_vec_init(S->p1->length);
        for (j = 0; j < S->p0->length; j++)
            arb_set_fmpz(S->p0a + j, S->p0->coeffs + j);
        for (j = 0; j < S->p1->length; j++)
            arb_set_fmpz(S->p1a + j, S->p1->coeffs + j);

        S->Qa = flint_malloc(len * sizeof(arb_ptr));
        S->La = _arb_vec_init(len);
        S->Ga = _arb_vec_init(len);
        for (k = 1; k <= S->nQ; k++)
        {
            S->Qa[k] = _arb_vec_init(S->Q[k].length);
            for (j = 0; j < S->Q[k].length; j++)
                arb_set_fmpz(S->Qa[k] + j, S->Q[k].coeffs + j);
            arb_set_fmpz(S->La + k, S->L + k);
            arb_set_fmpz(S->Ga + k, S->G + k);
        }
    }

    return ok;
}

/* Horner's rule, skipping runs of zero coefficients. */
static void
_sparse_arb_poly_evaluate(arb_t y, arb_srcptr A, slong len, const arb_t x, slong prec)
{
    slong i, j;
    arb_t t;

    if (len == 0)
    {
        arb_zero(y);
        return;
    }

    arb_init(t);
    arb_set_round(y, A + len - 1, prec);

    for (i = len - 1; i > 0; i = j)
    {
        for (j = i - 1; j > 0 && arb_is_zero(A + j); j--)
            ;

        if (i - j == 1)
        {
            arb_mul(y, y, x, prec);
        }
        else
        {
            arb_pow_ui(t, x, i - j, prec);
            arb_mul(y, y, t, prec);
        }

        arb_add(y, y, A + j, prec);
    }

    arb_clear(t);
}

#define ADD_SIGN(s) \
    do { \
        if ((s) != 0) \
        { \
            if (last != 0 && (s) != last) \
                V++; \
            last = (s); \
        } \
    } while (0)

/* Number of sign variations at x, or -1 if some sign could not be
   determined. Sets *root to whether p_0(x) = 0. If p0x is not NULL,
   it is the exact value of p_0(x). */
static slong
sturm_var_arb(const sturm_seq_t S, const arb_t x, arb_srcptr p0x, slong prec, int * root)
{
    arb_t y0, y1, y2, t;
    slong j, V = 0;
    int s, last = 0;

    arb_init(y0);
    arb_init(y1);
    arb_init(y2);
    arb_init(t);

    if (p0x != NULL)
        arb_set(y0, p0x);
    else
        _sparse_arb_poly_evaluate(y0, S->p0a, S->p0->length, x, prec);
    _sparse_arb_poly_evaluate(y1, S->p1a, S->p1->length, x, prec);

    for (j = 0; j < S->n; j++)
    {
        if (j >= 2)
        {
            /* y2 = (L y0 - Q(x) y1) / G */
            _sparse_arb_poly_evaluate(t, S->Qa[j - 1], S->Q[j - 1].length, x, prec);
            arb_mul(t, t, y1, prec);
            arb_mul(y2, S->La + j - 1, y0, prec);
            arb_sub(y2, y2, t, prec);
            arb_div(y2, y2, S->Ga + j - 1, prec);
            arb_swap(y0, y1);
            arb_swap(y1, y2);
        }

        {
            arb_srcptr y = (j == 0) ? y0 : y1;

            if (arb_is_positive(y))
                s = 1;
            else if (arb_is_negative(y))
                s = -1;
            else if (j == 0 && arb_is_zero(y))
                s = 0;
            else
            {
                V = -1;
                break;
            }
        }

        if (j == 0)
            *root = (s == 0);

        ADD_SIGN(s);
    }

    arb_clear(y0);
    arb_clear(y1);
    arb_clear(y2);
    arb_clear(t);

    return V;
}

static slong
sturm_var_exact(const sturm_seq_t S, const fmpq_t x, int * root)
{
    fmpq_t y0, y1, y2, t;
    slong j, V = 0;
    int s, last = 0;

    fmpq_init(y0);
    fmpq_init(y1);
    fmpq_init(y2);
    fmpq_init(t);

    fmpz_poly_evaluate_fmpq(y0, S->p0, x);
    fmpz_poly_evaluate_fmpq(y1, S->p1, x);

    for (j = 0; j < S->n; j++)
    {
        if (j >= 2)
        {
            fmpz_poly_evaluate_fmpq(t, S->Q + j - 1, x);
            fmpq_mul(t, t, y1);
            fmpq_mul_fmpz(y2, y0, S->L + j - 1);
            fmpq_sub(y2, y2, t);
            fmpq_div_fmpz(y2, y2, S->G + j - 1);
            fmpq_swap(y0, y1);
            fmpq_swap(y1, y2);
        }

        s = fmpq_sgn((j == 0) ? y0 : y1);

        if (j == 0)
            *root = (s == 0);

        ADD_SIGN(s);
    }

    fmpq_clear(y0);
    fmpq_clear(y1);
    fmpq_clear(y2);
    fmpq_clear(t);

    return V;
}

/* Number of sign variations at x = m 2^e. */
static slong
sturm_var(const sturm_seq_t S, const fmpz_t m, slong e, int * root)
{
    arb_t x;
    slong V, prec, prec_start;

    arb_init(x);
    arb_set_fmpz(x, m);
    arb_mul_2exp_si(x, x, e);

    prec_start = S->p0_bits + 2 * fmpz_bits(m) + 64;

    V = -1;
    for (prec = prec_start; prec <= 16 * prec_start && V == -1; prec *= 2)
        V = sturm_var_arb(S, x, NULL, prec, root);

    if (V == -1)
    {
        /* Evaluate p_0 exactly (it may be zero, e.g. if the midpoint is an
           exact root), and retry the others. */
        fmpq_t q;
        arb_t y;

        fmpq_init(q);
        arb_init(y);
        fmpz_set(fmpq_numref(q), m);
        if (e >= 0)
            fmpq_mul_2exp(q, q, e);
        else
            fmpq_div_2exp(q, q, -e);
        fmpz_poly_evaluate_fmpq(q, S->p0, q);
        /* the value is a dyadic number */
        arb_set_fmpz(y, fmpq_numref(q));
        arb_mul_2exp_si(y, y, -(slong) (fmpz_bits(fmpq_denref(q)) - 1));

        for (prec = prec_start; prec <= 4 * prec_start && V == -1; prec *= 2)
            V = sturm_var_arb(S, x, y, prec, root);

        fmpq_clear(q);
        arb_clear(y);
    }

    if (V == -1)
    {
        fmpq_t q;
        fmpq_init(q);
        fmpz_set(fmpq_numref(q), m);
        if (e >= 0)
            fmpq_mul_2exp(q, q, e);
        else
            fmpq_div_2exp(q, q, -e);
        V = sturm_var_exact(S, q, root);
        fmpq_clear(q);
    }

    arb_clear(x);
    return V;
}

static slong
sturm_var_zero(const sturm_seq_t S)
{
    slong j, V = 0;
    int last = 0;

    for (j = 0; j < S->n; j++)
        ADD_SIGN(S->sgn_const[j]);

    return V;
}

static slong
sturm_var_inf(const sturm_seq_t S, int negative)
{
    slong j, V = 0;
    int s, last = 0;

    for (j = 0; j < S->n; j++)
    {
        s = S->sgn_lead[j];
        if (negative && (S->deg[j] % 2 == 1))
            s = -s;
        ADD_SIGN(s);
    }

    return V;
}

typedef struct
{
    fmpz_t c;
    slong k;
    slong Va;
    slong Vb;
    int exact;      /* 1: the exact root c 2^k; 0: the interval (c 2^k, (c+1) 2^k) */
}
bisect_node_struct;

/* Isolates the roots in (c 2^k, (c+1) 2^k) where the endpoints are not
   roots and have sign variations Va and Vb, appending the output in
   increasing order. */
static void
sturm_bisect(fmpq * exact_roots, slong * n_exact, fmpz * c_array,
    slong * k_array, slong * n_interval, const sturm_seq_t S,
    const fmpz_t c, slong k, slong Va, slong Vb)
{
    bisect_node_struct * stack;
    slong alloc, top, i, Vm;
    fmpz_t m;
    int root;

    alloc = 16;
    stack = flint_malloc(alloc * sizeof(bisect_node_struct));
    for (i = 0; i < alloc; i++)
        fmpz_init(stack[i].c);
    fmpz_init(m);

    fmpz_set(stack[0].c, c);
    stack[0].k = k;
    stack[0].Va = Va;
    stack[0].Vb = Vb;
    stack[0].exact = 0;
    top = 1;

    while (top > 0)
    {
        bisect_node_struct * node;

        top--;
        node = stack + top;

        if (node->exact)
        {
            if (exact_roots != NULL)
            {
                fmpz_set(fmpq_numref(exact_roots + *n_exact), node->c);
                fmpz_one(fmpq_denref(exact_roots + *n_exact));
                if (node->k >= 0)
                    fmpq_mul_2exp(exact_roots + *n_exact, exact_roots + *n_exact, node->k);
                else
                    fmpq_div_2exp(exact_roots + *n_exact, exact_roots + *n_exact, -node->k);
            }
            (*n_exact)++;
            continue;
        }

        if (node->Va - node->Vb == 0)
            continue;

        if (node->Va - node->Vb == 1)
        {
            if (c_array != NULL && k_array != NULL)
            {
                fmpz_set(c_array + *n_interval, node->c);
                k_array[*n_interval] = node->k;
            }
            (*n_interval)++;
            continue;
        }

        /* midpoint m 2^(k-1) with m = 2c + 1 */
        fmpz_mul_2exp(m, node->c, 1);
        fmpz_add_ui(m, m, 1);
        k = node->k - 1;
        Va = node->Va;
        Vb = node->Vb;
        Vm = sturm_var(S, m, k, &root);

        /* push right half, the exact root (if any) and the left half */
        if (top + 3 > alloc)
        {
            stack = flint_realloc(stack, 2 * alloc * sizeof(bisect_node_struct));
            for (i = alloc; i < 2 * alloc; i++)
                fmpz_init(stack[i].c);
            alloc *= 2;
        }

        /* note: node may be invalidated by realloc; use m, k, Va, Vb */
        fmpz_set(stack[top].c, m);
        stack[top].k = k;
        stack[top].Va = Vm;
        stack[top].Vb = Vb;
        stack[top].exact = 0;
        top++;

        if (root)
        {
            fmpz_set(stack[top].c, m);
            stack[top].k = k;
            stack[top].exact = 1;
            top++;
        }

        /* the Sturm count V(a) - V(b) counts the roots in (a, b] */
        fmpz_sub_ui(stack[top].c, m, 1);
        stack[top].k = k;
        stack[top].Va = Va;
        stack[top].Vb = Vm + root;
        stack[top].exact = 0;
        top++;
    }

    for (i = 0; i < alloc; i++)
        fmpz_clear(stack[i].c);
    flint_free(stack);
    fmpz_clear(m);
}

int
_fmpz_poly_isolate_real_roots_sturm(fmpq * exact_roots, slong * n_exact,
    fmpz * c_array, slong * k_array, slong * n_interval,
    const fmpz * pol, slong len, int positive_only, slong max_size)
{
    sturm_seq_t S;
    fmpz * tmp;
    slong i, n_zeros, kp, kn, V0;
    fmpz_t c;

    for (n_zeros = 0; n_zeros < len && fmpz_is_zero(pol + n_zeros); n_zeros++)
        ;
    pol += n_zeros;
    len -= n_zeros;

    if (len == 0)
        flint_throw(FLINT_ERROR, "(%s): zero polynomial\n", __func__);

    if (len == 1)
    {
        /* only roots at zero */
        *n_exact = *n_interval = 0;
        if (!positive_only)
        {
            if (exact_roots != NULL)
                for (i = 0; i < n_zeros; i++)
                    fmpq_zero(exact_roots + i);
            *n_exact = n_zeros;
        }
        return 1;
    }

    if (!sturm_seq_init(S, pol, len, max_size))
    {
        sturm_seq_clear(S);
        return 0;
    }

    *n_exact = *n_interval = 0;

    /* strict power-of-two bounds for the positive and negative roots */
    tmp = _fmpz_vec_init(len);
    _fmpz_vec_set(tmp, pol, len);
    kp = _fmpz_poly_scale_positive_roots_0_1(tmp, len);
    _fmpz_vec_set(tmp, pol, len);
    for (i = 1; i < len; i += 2)
        fmpz_neg(tmp + i, tmp + i);
    kn = _fmpz_poly_scale_positive_roots_0_1(tmp, len);
    _fmpz_vec_clear(tmp, len);

    V0 = sturm_var_zero(S);
    fmpz_init(c);

    if (!positive_only)
    {
        if (kn != WORD_MIN)
        {
            fmpz_set_si(c, -1);
            sturm_bisect(exact_roots, n_exact, c_array, k_array, n_interval,
                S, c, kn, sturm_var_inf(S, 1), V0);
        }

        if (exact_roots != NULL)
            for (i = 0; i < n_zeros; i++)
                fmpq_zero(exact_roots + *n_exact + i);
        *n_exact += n_zeros;
    }

    if (kp != WORD_MIN)
    {
        fmpz_zero(c);
        sturm_bisect(exact_roots, n_exact, c_array, k_array, n_interval,
            S, c, kp, V0, sturm_var_inf(S, 0));
    }

    fmpz_clear(c);
    sturm_seq_clear(S);
    return 1;
}
