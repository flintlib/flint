/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "thread_support.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "factor_impl.h"

FLINT_DLL slong gr_poly_factor_threaded_cutoff = GR_POLY_FACTOR_THREADED_CUTOFF;

/* Should we compute iterated Frobenius powers x^(q^i) mod v using modular
   composition with a precomputed matrix (cost: ~ n^2 + 2 sqrt(n) mulmods
   per step) rather than by repeated exponentiation (cost: ~ log2(q) mulmods
   per step)? Threshold taken from nmod_poly_factor_distinct_deg. */
static int
_use_composition(const fmpz_t q, slong lenv)
{
    return fmpz_bits(q) > ((n_sqrt(lenv - 1) + 1) * 3) / 4;
}

int
_gr_poly_iterated_frobenius_preinv(gr_poly_struct * h, slong n,
    const gr_poly_t frob, const gr_poly_preinv_t P,
    const fmpz_t q, gr_ctx_t ctx)
{
    slong i;
    int status = GR_SUCCESS;

    if (n <= 0)
        return GR_SUCCESS;

    status |= gr_poly_gen(h + 0, ctx);
    status |= gr_poly_preinv_rem(h + 0, h + 0, P, ctx);

    if (n == 1)
        return status;

    status |= gr_poly_set(h + 1, frob, ctx);

    if (n == 2 || P->lenf <= 2)
    {
        for (i = 2; i < n; i++)
            status |= gr_poly_set(h + i, h + 1, ctx);
        return status;
    }

    if (_use_composition(q, P->lenf))
    {
        /* Doubling: given h[1..2^i], compute h[2^i+1..2^(i+1)] as the
           compositions h[1..2^i] o h[2^i] using a single vector
           composition (one matrix product). */
        slong done = 1;   /* h[1..done] computed */

        while (done < n - 1 && status == GR_SUCCESS)
        {
            slong num = FLINT_MIN(done, n - 1 - done);
            status |= gr_poly_preinv_compose_mod_brent_kung_vec(h + done + 1, h + 1, num, num, h + done, P, ctx);
            done += num;
        }
    }
    else
    {
        for (i = 2; i < n && status == GR_SUCCESS; i++)
            status |= gr_poly_preinv_powmod_fmpz_sliding(h + i, h + i - 1, q, 0, P, ctx);
    }

    return status;
}

int
gr_poly_iterated_frobenius_preinv(gr_poly_struct * h, slong n,
    const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    gr_poly_t frob;
    fmpz_t q;
    int status = GR_SUCCESS;

    if (n <= 0)
        return GR_SUCCESS;

    fmpz_init(q);
    gr_poly_init(frob, ctx);

    status |= _gr_poly_factor_ff_info(q, NULL, NULL, ctx);

    if (status == GR_SUCCESS)
    {
        if (n >= 2)
            status |= gr_poly_preinv_powmod_x_fmpz(frob, q, P, ctx);

        status |= _gr_poly_iterated_frobenius_preinv(h, n, frob, P, q, ctx);
    }

    fmpz_clear(q);
    gr_poly_clear(frob, ctx);
    return status;
}

/*
    Distinct degree factorization of a monic squarefree polynomial v
    using the baby-step giant-step algorithm of Kaltofen and Shoup.
    Requires the precomputed modulus P for v and frob = x^q mod v.
*/
int
_gr_poly_factor_distinct_deg_with_frob(gr_poly_vec_t fac, fmpz_vec_t degs,
    const gr_poly_t poly, const gr_poly_preinv_t polyP, const gr_poly_t frob,
    const fmpz_t q, gr_ctx_t ctx)
{
    gr_poly_t f, g, v, tmp;
    gr_poly_preinv_t P;
    int have_P = 0;
    const gr_poly_preinv_struct * PP = polyP;
    gr_poly_struct * h, * H, * I;
    slong i, j, l, m, n, d;
    gr_mat_t HH, HHH;
    double beta;
    int status = GR_SUCCESS;
    int have_HH = 0;

    n = poly->length - 1;

    if (n == 1)
    {
        status |= gr_poly_vec_append(fac, poly, ctx);
        fmpz_vec_append_ui(degs, 1);
        return status;
    }

    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);

    gr_poly_init(f, ctx);
    gr_poly_init(g, ctx);
    gr_poly_init(v, ctx);
    gr_poly_init(tmp, ctx);

    h = flint_malloc((2 * m + l + 1) * sizeof(gr_poly_struct));
    H = h + (l + 1);
    I = H + m;

    for (i = 0; i < 2 * m + l + 1; i++)
        gr_poly_init(h + i, ctx);

    status |= gr_poly_set(v, poly, ctx);

    /* compute baby steps: h[i] = x^{q^i} mod v */
    status |= _gr_poly_iterated_frobenius_preinv(h, l + 1, frob, polyP, q, ctx);
    GR_POLY_FACTOR_CHECK_STATUS()

    /* compute coarse distinct-degree factorisation */
    status |= gr_poly_set(H + 0, h + l, ctx);
    gr_mat_init(HH, n_sqrt(v->length - 1) + 1, v->length - 1, ctx);
    have_HH = 1;
    status |= gr_poly_preinv_precompute_matrix(HH, H + 0, PP, ctx);
    GR_POLY_FACTOR_CHECK_STATUS()

    d = 1;
    for (j = 0; j < m; j++)
    {
        /* compute giant steps: H[j] = x^{q^(lj)} mod v */
        if (j > 0)
        {
            if (I[j - 1].length > 1)
            {
                /* v has been replaced by a proper divisor; reduce the matrix
                   of powers of H[0] accordingly */
                status |= _gr_poly_reduce_matrix_mod_poly(HHH, HH, v, ctx);
                gr_mat_swap(HH, HHH, ctx);
                gr_mat_clear(HHH, ctx);
                status |= gr_poly_rem(tmp, H + j - 1, v, ctx);
                status |= gr_poly_preinv_compose_mod_brent_kung_precomp(H + j, tmp, HH, PP, ctx);
            }
            else
            {
                status |= gr_poly_preinv_compose_mod_brent_kung_precomp(H + j, H + j - 1, HH, PP, ctx);
            }

            GR_POLY_FACTOR_CHECK_STATUS()
        }

        /* compute interval polynomials */
        status |= gr_poly_one(I + j, ctx);
        for (i = l - 1; (i >= 0) && (2 * d <= v->length - 1); i--, d++)
        {
            status |= gr_poly_rem(tmp, h + i, v, ctx);
            status |= gr_poly_sub(tmp, H + j, tmp, ctx);
            status |= gr_poly_preinv_mulmod(I + j, tmp, I + j, PP, ctx);
        }

        /* compute F_j = f^{[j*l+1]} * ... * f^{[j*l+l]}; stored in I_j
           (the product may vanish mod v when all remaining factors of v
           have degrees in this interval) */
        if (I[j].length == 0)
            status |= gr_poly_set(I + j, v, ctx);
        else
        {
            status |= gr_poly_gcd(I + j, v, I + j, ctx);
        }
        GR_POLY_FACTOR_CHECK(I + j)

        if (I[j].length > 1)
        {
            status |= gr_poly_divexact(v, v, I + j, ctx);
            if (!have_P)
            {
                gr_poly_preinv_init(P, ctx);
                have_P = 1;
            }
            status |= gr_poly_preinv_set(P, v, ctx);
            PP = P;
            GR_POLY_FACTOR_CHECK(v)
        }

        if (v->length - 1 < 2 * d)
            break;
    }

    if (v->length > 1)
    {
        fmpz_vec_append_ui(degs, v->length - 1);
        status |= gr_poly_vec_append(fac, v, ctx);
    }

    /* compute fine distinct-degree factorisation */
    for (j = 0; j < m; j++)
    {
        if (I[j].length - 1 > (j + 1) * l || j == 0)
        {
            status |= gr_poly_set(g, I + j, ctx);

            for (i = l - 1; i >= 0 && (g->length > 1); i--)
            {
                /* compute f^{[l*(j+1)-i]} */
                status |= gr_poly_sub(tmp, H + j, h + i, ctx);
                status |= gr_poly_gcd(f, g, tmp, ctx);
                GR_POLY_FACTOR_CHECK(f)

                if (f->length > 1)
                {
                    status |= gr_poly_divexact(g, g, f, ctx);
                    GR_POLY_FACTOR_CHECK(g)

                    fmpz_vec_append_ui(degs, l * (j + 1) - i);
                    status |= gr_poly_vec_append(fac, f, ctx);
                }
            }
        }
        else if (I[j].length > 1)
        {
            fmpz_vec_append_ui(degs, I[j].length - 1);
            status |= gr_poly_vec_append(fac, I + j, ctx);
        }
    }

cleanup:
    gr_poly_clear(f, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_clear(v, ctx);
    gr_poly_clear(tmp, ctx);

    if (have_P)
        gr_poly_preinv_clear(P, ctx);

    if (have_HH)
        gr_mat_clear(HH, ctx);

    for (i = 0; i < 2 * m + l + 1; i++)
        gr_poly_clear(h + i, ctx);
    flint_free(h);

    return status;
}

/*
    Multithreaded distinct degree factorization.

    The serial algorithm computes the giant steps H[j] = x^(q^(l(j+1))) mod v
    one after the other and reduces v as soon as an interval polynomial
    splits off factors. Here the giant steps are computed in blocks of B
    independent compositions: with the matrix of powers of
    x^(q^(lB)) = H[B-1] precomputed once (and reduced whenever v shrinks),
    H[j + B] = H[j] o H[B-1] for the B values j of a block are independent,
    as are the B interval polynomials. Within a block, v is only reduced
    after the gcd of v with the product of the block's interval
    polynomials, so slightly more work is done than in the serial version
    when factors of small degree are present; the fine distinct degree
    factorization is also done in parallel over the intervals.

    Matrix products are left to gr_mat_mul (which may use threads
    internally); the parallel sections only contain polynomial
    arithmetic, using flint_parallel_do.
*/

typedef struct
{
    gr_poly_struct * H;             /* giant steps */
    gr_poly_struct * I;             /* interval polynomials */
    const gr_poly_struct * h;       /* baby steps */
    const gr_mat_struct * HH;       /* matrix of powers of H[B-1] mod v */
    const gr_poly_struct * v;
    const gr_poly_preinv_struct * PP;
    slong l, B, j0, d0, nmax;
    gr_ctx_struct * ctx;
    int * status;
}
ddf_block_args_t;

/* H[j0 + i] = H[j0 + i - B] o H[B-1] mod v */
static void
_ddf_giant_worker(slong i, void * arg)
{
    ddf_block_args_t * a = arg;
    gr_poly_t tmp;
    int status = GR_SUCCESS;
    slong j = a->j0 + i;

    gr_poly_init(tmp, a->ctx);
    status |= gr_poly_preinv_rem(tmp, a->H + j - a->B, a->PP, a->ctx);
    status |= gr_poly_preinv_compose_mod_brent_kung_precomp(a->H + j, tmp, a->HH, a->PP, a->ctx);
    gr_poly_clear(tmp, a->ctx);
    a->status[i] = status;
}

/* I[j0 + i] = prod_k (H[j0 + i] - h[k]) over the degrees d covered by
   the interval, 2d <= deg v */
static void
_ddf_interval_worker(slong i, void * arg)
{
    ddf_block_args_t * a = arg;
    gr_poly_t tmp, hk;
    int status = GR_SUCCESS;
    slong j = a->j0 + i, k, d;

    gr_poly_init(tmp, a->ctx);
    gr_poly_init(hk, a->ctx);
    status |= gr_poly_one(a->I + j, a->ctx);

    for (k = a->l - 1, d = a->d0 + i * a->l; k >= 0 && 2 * d <= a->nmax; k--, d++)
    {
        status |= gr_poly_preinv_rem(hk, a->h + k, a->PP, a->ctx);
        status |= gr_poly_sub(tmp, a->H + j, hk, a->ctx);
        status |= gr_poly_preinv_mulmod(a->I + j, tmp, a->I + j, a->PP, a->ctx);
    }

    gr_poly_clear(tmp, a->ctx);
    gr_poly_clear(hk, a->ctx);
    a->status[i] = status;
}

typedef struct
{
    const gr_poly_struct * H;
    const gr_poly_struct * I;
    const gr_poly_struct * h;
    gr_poly_vec_struct * out;       /* per interval: factors found */
    slong * out_deg;                /* per interval: their degrees, l slots each */
    slong l;
    gr_ctx_struct * ctx;
    int * status;
}
ddf_fine_args_t;

static void
_ddf_fine_worker(slong j, void * arg)
{
    ddf_fine_args_t * a = arg;
    gr_poly_t g, f, tmp;
    slong i, num = 0;
    int status = GR_SUCCESS;

    gr_poly_init(g, a->ctx);
    gr_poly_init(f, a->ctx);
    gr_poly_init(tmp, a->ctx);
    gr_poly_vec_set_length(a->out + j, 0, a->ctx);

    if (a->I[j].length - 1 > (j + 1) * a->l || j == 0)
    {
        status |= gr_poly_set(g, a->I + j, a->ctx);

        for (i = a->l - 1; i >= 0 && g->length > 1 && status == GR_SUCCESS; i--)
        {
            /* factors of degree l (j + 1) - i */
            status |= gr_poly_sub(tmp, a->H + j, a->h + i, a->ctx);
            status |= gr_poly_gcd(f, g, tmp, a->ctx);

            if (f->length > 1)
            {
                status |= gr_poly_divexact(g, g, f, a->ctx);
                gr_poly_vec_append_swap(a->out + j, f, a->ctx);
                a->out_deg[j * a->l + num] = a->l * (j + 1) - i;
                num++;
            }
            else if (f->length < 1)
            {
                status |= GR_UNABLE;
            }
        }
    }
    else if (a->I[j].length > 1)
    {
        status |= gr_poly_vec_append(a->out + j, a->I + j, a->ctx);
        a->out_deg[j * a->l] = a->I[j].length - 1;
    }

    gr_poly_clear(g, a->ctx);
    gr_poly_clear(f, a->ctx);
    gr_poly_clear(tmp, a->ctx);
    a->status[j] = status;
}

int
_gr_poly_factor_distinct_deg_with_frob_threaded(gr_poly_vec_t fac, fmpz_vec_t degs,
    const gr_poly_t poly, const gr_poly_preinv_t polyP, const gr_poly_t frob,
    const fmpz_t q, gr_ctx_t ctx)
{
    gr_poly_t v, tmp, G;
    gr_poly_preinv_t P;
    const gr_poly_preinv_struct * PP = polyP;
    gr_poly_struct * h, * H, * I;
    slong i, j, l, m, n, d, B, j0, nthreads;
    gr_mat_t HH, HHH;
    double beta;
    int status = GR_SUCCESS;
    int have_HH = 0, have_P = 0;
    int * wstatus;

    n = poly->length - 1;

    if (n == 1)
    {
        status |= gr_poly_vec_append(fac, poly, ctx);
        fmpz_vec_append_ui(degs, 1);
        return status;
    }

    nthreads = flint_get_num_available_threads();
    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);
    B = FLINT_MIN(m, FLINT_MAX(2, nthreads));

    gr_poly_init(v, ctx);
    gr_poly_init(tmp, ctx);
    gr_poly_init(G, ctx);
    h = flint_malloc((2 * m + l + 1) * sizeof(gr_poly_struct));
    H = h + (l + 1);
    I = H + m;
    for (i = 0; i < 2 * m + l + 1; i++)
        gr_poly_init(h + i, ctx);
    wstatus = flint_calloc(FLINT_MAX(B, m), sizeof(int));

    status |= gr_poly_set(v, poly, ctx);

    /* baby steps: h[i] = x^{q^i} mod v (the vector composition used
       inside is parallel) */
    status |= _gr_poly_iterated_frobenius_preinv(h, l + 1, frob, polyP, q, ctx);
    GR_POLY_FACTOR_CHECK_STATUS()

    /* first block of giant steps: H[0] = h[l], H[j] = H[j-1] o H[0] */
    status |= gr_poly_set(H + 0, h + l, ctx);
    gr_mat_init(HH, n_sqrt(n) + 1, n, ctx);
    have_HH = 1;
    status |= gr_poly_preinv_precompute_matrix(HH, H + 0, PP, ctx);
    for (j = 1; j < B && status == GR_SUCCESS; j++)
        status |= gr_poly_preinv_compose_mod_brent_kung_precomp(H + j, H + j - 1, HH, PP, ctx);
    GR_POLY_FACTOR_CHECK_STATUS()

    if (m > B)
    {
        /* matrix of powers of x^(q^(lB)) = H[B-1] for the subsequent
           blocks: H[j + B] = H[j] o H[B-1] */
        status |= gr_poly_preinv_precompute_matrix(HH, H + B - 1, PP, ctx);
        GR_POLY_FACTOR_CHECK_STATUS()
    }

    d = 1;
    for (j0 = 0; j0 < m && status == GR_SUCCESS; j0 += B)
    {
        slong nb = FLINT_MIN(B, m - j0);
        ddf_block_args_t args;

        args.H = H; args.I = I; args.h = h; args.HH = HH; args.v = v; args.PP = PP;
        args.l = l; args.B = B; args.j0 = j0; args.d0 = d; args.nmax = v->length - 1;
        args.ctx = ctx; args.status = wstatus;

        if (j0 > 0)
        {
            /* giant steps of this block, in parallel */
            if (nb > 1 && gr_ctx_is_threadsafe(ctx) == T_TRUE)
                flint_parallel_do(_ddf_giant_worker, &args, nb, -1, FLINT_PARALLEL_UNIFORM);
            else
                for (i = 0; i < nb; i++)
                    _ddf_giant_worker(i, &args);
            for (i = 0; i < nb; i++)
                status |= wstatus[i];
            GR_POLY_FACTOR_CHECK_STATUS()
        }

        /* interval polynomials, in parallel */
        if (nb > 1 && gr_ctx_is_threadsafe(ctx) == T_TRUE)
            flint_parallel_do(_ddf_interval_worker, &args, nb, -1, FLINT_PARALLEL_UNIFORM);
        else
            for (i = 0; i < nb; i++)
                _ddf_interval_worker(i, &args);
        for (i = 0; i < nb; i++)
            status |= wstatus[i];
        GR_POLY_FACTOR_CHECK_STATUS()

        /* gcd of v with the product of the block's interval polynomials,
           then refine per interval */
        status |= gr_poly_one(G, ctx);
        for (i = 0; i < nb; i++)
            status |= gr_poly_preinv_mulmod(G, G, I + j0 + i, PP, ctx);
        /* the product may vanish mod v (all remaining factors covered) */
        if (G->length == 0)
            status |= gr_poly_set(G, v, ctx);
        else
            status |= gr_poly_gcd(G, v, G, ctx);
        GR_POLY_FACTOR_CHECK(G)

        if (G->length > 1)
        {
            for (i = 0; i < nb; i++)
            {
                if (I[j0 + i].length == 0)
                    status |= gr_poly_set(I + j0 + i, G, ctx);
                else
                    status |= gr_poly_gcd(I + j0 + i, I + j0 + i, G, ctx);
                if (I[j0 + i].length > 1)
                    status |= gr_poly_divexact(G, G, I + j0 + i, ctx);
            }
            GR_POLY_FACTOR_CHECK_STATUS()

            /* v /= product of the intervals' factors */
            for (i = 0; i < nb; i++)
                if (I[j0 + i].length > 1)
                    status |= gr_poly_divexact(v, v, I + j0 + i, ctx);
            GR_POLY_FACTOR_CHECK(v)

            if (v->length > 1)
            {
                if (!have_P)
                {
                    gr_poly_preinv_init(P, ctx);
                    have_P = 1;
                }
                status |= gr_poly_preinv_set(P, v, ctx);
                PP = P;

                status |= _gr_poly_reduce_matrix_mod_poly(HHH, HH, v, ctx);
                gr_mat_swap(HH, HHH, ctx);
                gr_mat_clear(HHH, ctx);
            }
        }
        else
        {
            for (i = 0; i < nb; i++)
                status |= gr_poly_one(I + j0 + i, ctx);
        }

        /* degrees covered so far; once more than half the degree of the
           remaining v is covered, v is irreducible (or 1) */
        d += nb * l;
        if (v->length - 1 < 2 * d)
            break;
    }

    if (status == GR_SUCCESS && v->length > 1)
    {
        fmpz_vec_append_ui(degs, v->length - 1);
        status |= gr_poly_vec_append(fac, v, ctx);
    }

    /* fine distinct degree factorization, in parallel over the intervals */
    if (status == GR_SUCCESS)
    {
        ddf_fine_args_t fargs;
        gr_poly_vec_struct * out = flint_malloc(m * sizeof(gr_poly_vec_struct));
        slong * out_deg = flint_malloc(m * l * sizeof(slong));

        for (j = 0; j < m; j++)
            gr_poly_vec_init(out + j, 0, ctx);

        fargs.H = H; fargs.I = I; fargs.h = h; fargs.out = out; fargs.out_deg = out_deg;
        fargs.l = l; fargs.ctx = ctx; fargs.status = wstatus;

        if (m > 1 && gr_ctx_is_threadsafe(ctx) == T_TRUE && nthreads > 1)
            flint_parallel_do(_ddf_fine_worker, &fargs, m, -1, FLINT_PARALLEL_UNIFORM);
        else
            for (j = 0; j < m; j++)
                _ddf_fine_worker(j, &fargs);

        for (j = 0; j < m; j++)
        {
            status |= wstatus[j];
            for (i = 0; i < out[j].length; i++)
            {
                fmpz_vec_append_ui(degs, out_deg[j * l + i]);
                gr_poly_vec_append_swap(fac, out[j].entries + i, ctx);
            }
        }

        for (j = 0; j < m; j++)
            gr_poly_vec_clear(out + j, ctx);
        flint_free(out);
        flint_free(out_deg);
    }

cleanup:
    gr_poly_clear(v, ctx);
    gr_poly_clear(tmp, ctx);
    gr_poly_clear(G, ctx);
    if (have_P)
        gr_poly_preinv_clear(P, ctx);
    if (have_HH)
        gr_mat_clear(HH, ctx);
    for (i = 0; i < 2 * m + l + 1; i++)
        gr_poly_clear(h + i, ctx);
    flint_free(h);
    flint_free(wstatus);

    return status;
}

int
gr_poly_factor_distinct_deg(gr_poly_vec_t fac, fmpz_vec_t degs,
    const gr_poly_t poly, gr_ctx_t ctx)
{
    gr_poly_t v, frob;
    gr_poly_preinv_t P;
    fmpz_t q;
    int status = GR_SUCCESS;

    gr_poly_vec_set_length(fac, 0, ctx);
    fmpz_vec_set_length(degs, 0);

    if (poly->length <= 1)
        return GR_DOMAIN;

    fmpz_init(q);
    gr_poly_init(v, ctx);
    gr_poly_preinv_init(P, ctx);
    gr_poly_init(frob, ctx);

    status |= _gr_poly_factor_ff_info(q, NULL, NULL, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;

    status |= gr_poly_make_monic(v, poly, ctx);
    GR_POLY_FACTOR_CHECK(v)

    if (v->length == 2)
    {
        status |= gr_poly_vec_append(fac, v, ctx);
        fmpz_vec_append_ui(degs, 1);
        goto cleanup;
    }

    status |= gr_poly_preinv_set(P, v, ctx);
    status |= gr_poly_preinv_powmod_x_fmpz(frob, q, P, ctx);
    GR_POLY_FACTOR_CHECK_STATUS()

    if (v->length > gr_poly_factor_threaded_cutoff && flint_get_num_available_threads() > 1 && gr_ctx_is_threadsafe(ctx) == T_TRUE)
        status |= _gr_poly_factor_distinct_deg_with_frob_threaded(fac, degs, v, P, frob, q, ctx);
    else
        status |= _gr_poly_factor_distinct_deg_with_frob(fac, degs, v, P, frob, q, ctx);

cleanup:
    fmpz_clear(q);
    gr_poly_clear(v, ctx);
    gr_poly_preinv_clear(P, ctx);
    gr_poly_clear(frob, ctx);

    return status;
}
