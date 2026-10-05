/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_poly_mat.h"
#include "fmpz_mat.h"
#include "mpoly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_q.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

/*
    Dense arithmetic for towers of number fields whose moduli are monic
    univariate integer polynomials, each in its own variable (a single
    algebraic generator, as in Q(zeta_N) or Q(sqrt(-163)), or a tensor
    product of such fields: roots of unity of coprime orders, radicals of
    integers), with no transcendental generators. Elements with integer
    denominators are multiplied as dense univariate integer polynomials
    through Kronecker substitution (the variable of step k becoming
    x^(S_k), with S_1 = 1, S_(k+1) = S_k (2 d_k - 1), so that the product
    of two reduced elements has no carries between the digits), then
    reduced digit by digit by the moduli; in a single generator, inverses
    are computed by the extended Euclidean algorithm over Q. This is the
    representation of number field elements as nf_elem in Calcium, and
    replaces the sparse multivariate products (heap multiplication) and
    the reductions of the general flat arithmetic.
*/

/* (GR_TOWER_OPT_DENSE_LIMIT: the size limit of the dense products, entries of the product array) */

int
_gr_tower_flat_dense_applicable(gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong k, L = 1;

    if (T->length == 0 || T->num_trans != 0 || T->num_gens != T->length ||
        !GR_TOWER_BASE_IS_CONSTS(T) || !GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL))
        return 0;

    gr_tower_flat_ensure(F);

    for (k = 1; k <= T->length; k++)
    {
        slong d = gr_tower_step_degree(T, k);
        (void) _gr_tower_flat_ideal_elem(F, k);
        if (F->ideal_univar[k - 1] == NULL)
            return 0;
        L *= (2 * d - 1);
        if (L > GR_TOWER_OPTION(F->T, GR_TOWER_OPT_DENSE_LIMIT))
            return 0;
    }

    return 1;
}

/* fmpz vectors in temporary storage (TMP_ALLOC: on the stack when
   small): zero entries on init; clearing frees the mpz entries only */
static fmpz *
_tmp_vec_init(slong n, void * mem)
{
    fmpz * v = (fmpz *) mem;
    slong i;
    for (i = 0; i < n; i++)
        fmpz_init(v + i);
    return v;
}

static void
_tmp_vec_clear(fmpz * v, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        fmpz_clear(v + i);
}

/* whether some algebraic generator of F is a root of unity (a product by
   a polynomial with several terms may then be a monomial: zeta_5^4 =
   -1 - zeta_5 - zeta_5^2 - zeta_5^3) */
static int
_any_root_of_unity(gr_tower_flat_t F)
{
    slong k;
    for (k = 1; k <= F->T->length; k++)
        if (GR_TOWER_STEP(F->T, k - 1)->def_kind == GR_TOWER_ROOT_OF_UNITY)
            return 1;
    return 0;
}

/* packs the reduced polynomial P into the dense vector A (entries at
   sum e_k S_k); returns 0 if P is not reduced or involves another
   variable */
static int _pack_generic(fmpz * A, const fmpz_mpoly_t P, const slong * var, const slong * deg, const slong * S, slong m, gr_tower_flat_t F, slong modn);

/* (exponents of at most one word per field: the digits are read from
   the packed exponents directly, and the monomial rebuilt from them must
   be the original one, i.e. no other variable occurs) */
static int
_pack_ex(fmpz * A, const fmpz_mpoly_t P, const slong * var, const slong * deg, const slong * S, slong m, gr_tower_flat_t F, slong modn)
{
    flint_bitcnt_t bits = P->bits;
    const mpoly_ctx_struct * minfo = F->mctx->minfo;
    slong N, t, k;
    slong * off, * shift;
    ulong * one, * tmp, mask;
    int ok = 1;
    TMP_INIT;

    if (bits > FLINT_BITS)
        return _pack_generic(A, P, var, deg, S, m, F, modn);

    TMP_START;
    N = mpoly_words_per_exp_sp(bits, minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    off = TMP_ALLOC(sizeof(slong) * 2 * m);
    shift = off + m;
    one = TMP_ALLOC(sizeof(ulong) * N * (m + 1));
    tmp = one + N * m;
    for (k = 0; k < m; k++)
        mpoly_gen_monomial_offset_shift_sp(one + N * k, off + k, shift + k, var[k], bits, minfo);

    if (minfo->ord == ORD_LEX)
    {
        /* (lex: the fields of the other variables must be zero) */
        slong w;
        mpoly_monomial_zero(tmp, N);
        for (k = 0; k < m; k++)
            tmp[off[k]] |= mask << shift[k];
        for (t = 0; t < P->length && ok; t++)
        {
            const ulong * exp = P->exps + N * t;
            slong idx = 0;
            for (w = 0; w < N; w++)
                if (exp[w] & ~tmp[w])
                    ok = 0;
            for (k = 0; k < m && ok; k++)
            {
                ulong e = (exp[off[k]] >> shift[k]) & mask;
                if (e >= (ulong) deg[k])
                    ok = 0;
                idx += (slong) e * S[k];
            }
            if (ok)
                fmpz_set(A + ((modn != 0) ? idx % modn : idx), P->coeffs + t);
        }
    }
    else
    for (t = 0; t < P->length && ok; t++)
    {
        const ulong * exp = P->exps + N * t;
        slong idx = 0;
        mpoly_monomial_zero(tmp, N);
        for (k = 0; k < m; k++)
        {
            ulong e = (exp[off[k]] >> shift[k]) & mask;
            if (e >= (ulong) deg[k])
            {
                ok = 0;
                break;
            }
            idx += (slong) e * S[k];
            mpoly_monomial_madd(tmp, tmp, e, one + N * k, N);
        }
        if (ok && !mpoly_monomial_equal(tmp, exp, N))
            ok = 0;
        if (ok)
            fmpz_set(A + ((modn != 0) ? idx % modn : idx), P->coeffs + t);
    }

    TMP_END;
    return ok;
}

static int
_pack_generic(fmpz * A, const fmpz_mpoly_t P, const slong * var, const slong * deg, const slong * S, slong m, gr_tower_flat_t F, slong modn)
{
    slong t, k, nvars = F->mctx->minfo->nvars;
    ulong * e = flint_malloc(sizeof(ulong) * nvars);
    int * used = flint_calloc(nvars, sizeof(int));
    int ok = 1;

    for (k = 0; k < m; k++)
        used[var[k]] = 1;

    for (t = 0; t < P->length && ok; t++)
    {
        slong i, idx = 0;
        fmpz_mpoly_get_term_exp_ui(e, P, t, F->mctx);
        for (i = 0; i < nvars && ok; i++)
            if (e[i] != 0 && !used[i])
                ok = 0;
        for (k = 0; k < m && ok; k++)
        {
            if (e[var[k]] >= (ulong) deg[k])
                ok = 0;
            else
                idx += (slong) e[var[k]] * S[k];
        }
        if (ok)
            fmpz_set(A + ((modn != 0) ? idx % modn : idx), P->coeffs + t);
    }

    flint_free(e);
    flint_free(used);
    return ok;
}

static int
_pack(fmpz * A, const fmpz_mpoly_t P, const slong * var, const slong * deg, const slong * S, slong m, gr_tower_flat_t F)
{
    return _pack_ex(A, P, var, deg, S, m, F, 0);
}

#define CYCLIC_LIMIT 100000
#define DENSE_REM_CUTOFF 64
#define DENSE_REM_NNZ 16

static slong
_nnz(const fmpz_poly_struct * M)
{
    slong i, n = 0;
    for (i = 0; i < M->length; i++)
        n += !fmpz_is_zero(M->coeffs + i);
    return n;
}

/* whether M = 1 + x + ... + x^d (d >= 1) */
static int
_is_all_ones(const fmpz_poly_struct * M)
{
    slong i;
    if (M->length < 2)
        return 0;
    for (i = 0; i < M->length; i++)
        if (!fmpz_is_one(M->coeffs + i))
            return 0;
    return 1;
}

/* reduces the dense product array A (digits 0 .. 2 d_k - 2) by the monic
   moduli, then writes it as a polynomial */
static void
_fold_array(fmpz * A, slong L, const slong * deg, const slong * S, slong m, gr_tower_flat_t F)
{
    slong k;
    /* variable by variable: for each block of the higher digits, the
       digit e_k from the top (2 d_k - 2 in the Kronecker layout, whose
       extents S_(k+1) / S_k are 2 d_k - 1) down to d_k, the lower digits
       forming a contiguous vector of length S_k:
       A[e_k] -> A[e_k - d_k + j] - M_j A[e_k] */
    for (k = 0; k < m; k++)
    {
        const fmpz_poly_struct * M = F->ideal_univar[k];
        slong d = deg[k], Sk = S[k], blk = S[k + 1], hi, ek, j;
        slong top = blk / Sk - 1;

        if (_is_all_ones(M))
        {
            /* the cyclotomic polynomial 1 + x + ... + x^d of a prime
               q = d + 1: x^q = 1 (one addition per coefficient), then
               x^d = -(1 + ... + x^(d-1)) */
            for (hi = 0; hi < L; hi += blk)
            {
                for (ek = top; ek >= d + 1; ek--)
                {
                    fmpz * src = A + hi + ek * Sk;
                    _fmpz_vec_add(A + hi + (ek - d - 1) * Sk, A + hi + (ek - d - 1) * Sk, src, Sk);
                    _fmpz_vec_zero(src, Sk);
                }
                {
                    fmpz * src = A + hi + d * Sk;
                    if (!_fmpz_vec_is_zero(src, Sk))
                    {
                        for (j = 0; j < d; j++)
                            _fmpz_vec_sub(A + hi + j * Sk, A + hi + j * Sk, src, Sk);
                        _fmpz_vec_zero(src, Sk);
                    }
                }
            }
            continue;
        }

        for (hi = 0; hi < L; hi += blk)
        {
            for (ek = top; ek >= d; ek--)
            {
                fmpz * src = A + hi + ek * Sk;
                if (_fmpz_vec_is_zero(src, Sk))
                    continue;
                for (j = 0; j < d; j++)
                    if (!fmpz_is_zero(M->coeffs + j))
                        _fmpz_vec_scalar_submul_fmpz(A + hi + (ek - d + j) * Sk, src, Sk, M->coeffs + j);
                _fmpz_vec_zero(src, Sk);
            }
        }
    }

}

/* the terms of the reduced array A (digits 0 .. d_k - 1; the entries
   are moved out, and A is left zero) as the polynomial R: the packed
   exponents are sums of multiples of the generators' monomials */
static void
_unpack(fmpz_mpoly_t R, fmpz * A, const slong * var, const slong * deg, const slong * S, slong m, gr_tower_flat_t F)
{
    const mpoly_ctx_struct * minfo = F->mctx->minfo;
    flint_bitcnt_t bits;
    slong k, N, len, idx;
    slong * digit, * didx;
    ulong * one, * cmpmask, * delta, * cur;
    int sorted = 1;
    ulong tot = 0;
    TMP_INIT;

    TMP_START;
    /* (room for the total degree, for the degree orderings) */
    for (k = 0; k < m; k++)
        tot += deg[k] - 1;
    bits = mpoly_fix_bits(FLINT_MAX(MPOLY_MIN_BITS, 1 + FLINT_BIT_COUNT(tot)), minfo);
    N = mpoly_words_per_exp(bits, minfo);

    /* the number of terms */
    len = 0;
    {
        slong top = 0;
        for (k = 0; k < m; k++)
            top += (deg[k] - 1) * S[k];
        for (idx = 0; idx <= top; idx++)
            len += !fmpz_is_zero(A + idx);
    }

    fmpz_mpoly_fit_length_reset_bits(R, len, bits, F->mctx);
    /* one[k]: the monomial of the variable of step k; delta[k]: the change
       of the monomial (modulo 2^(N FLINT_BITS)) and didx[k] of the index
       when digit k is decremented and the lower digits wrap around */
    one = TMP_ALLOC(sizeof(ulong) * N * (2 * m + 2));
    delta = one + N * m;
    cur = delta + N * m;
    cmpmask = cur + N;
    didx = TMP_ALLOC(sizeof(slong) * 2 * m);
    digit = didx + m;
    for (k = 0; k < m; k++)
        mpoly_gen_monomial_sp(one + N * k, var[k], bits, minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, minfo);
    for (k = 0; k < m; k++)
    {
        slong j;
        mpoly_monomial_zero(delta + N * k, N);
        didx[k] = -S[k];
        for (j = 0; j < k; j++)
        {
            mpoly_monomial_madd(delta + N * k, delta + N * k, deg[j] - 1, one + N * j, N);
            didx[k] += (deg[j] - 1) * S[j];
        }
        /* (subtracting one[k] as the two's complement) */
        mpn_sub_n(delta + N * k, delta + N * k, one + N * k, N);
    }

    /* the terms, by decreasing index (a mixed-radix counter, digit 0
       fastest) */
    {
        slong t = 0;
        mpoly_monomial_zero(cur, N);
        idx = 0;
        for (k = 0; k < m; k++)
        {
            digit[k] = deg[k] - 1;
            mpoly_monomial_madd(cur, cur, deg[k] - 1, one + N * k, N);
            idx += (deg[k] - 1) * S[k];
        }
        for (;;)
        {
            if (!fmpz_is_zero(A + idx))
            {
                ulong * exp = R->exps + N * t;
                mpoly_monomial_set(exp, cur, N);
                fmpz_swap(R->coeffs + t, A + idx);
                fmpz_zero(A + idx);
                if (t > 0 && sorted && !mpoly_monomial_gt(exp - N, exp, N, cmpmask))
                    sorted = 0;
                t++;
            }
            for (k = 0; k < m; k++)
            {
                if (digit[k] > 0)
                {
                    digit[k]--;
                    break;
                }
                digit[k] = deg[k] - 1;
            }
            if (k == m)
                break;
            if (N == 1)
                cur[0] += delta[k];
            else
                mpn_add_n(cur, cur, delta + N * k, N);
            idx += didx[k];
        }
        _fmpz_mpoly_set_length(R, t, F->mctx);
    }

    if (!sorted)
        fmpz_mpoly_sort_terms(R, F->mctx);
    TMP_END;
}

static void
_reduce_unpack(fmpz_mpoly_t R, fmpz * A, slong L, const slong * var, const slong * deg, const slong * S, slong m, gr_tower_flat_t F)
{
    _fold_array(A, L, deg, S, m, F);
    _unpack(R, A, var, deg, S, m, F);
}

/* x = num / den with den an integer: cancels the content */
static void
_canonicalise_int_den(fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    fmpz_t g, dc;
    fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(x);
    fmpz_mpoly_struct * den = fmpz_mpoly_q_denref(x);

    if (fmpz_mpoly_is_zero(num, F->mctx))
    {
        fmpz_mpoly_one(den, F->mctx);
        return;
    }

    if (fmpz_mpoly_is_one(den, F->mctx))
        return;

    fmpz_init(g);
    fmpz_init(dc);
    fmpz_mpoly_get_fmpz(dc, den, F->mctx);
    {
        /* gcd of the denominator and the coefficients, stopping at 1 */
        slong i;
        fmpz_abs(g, dc);
        for (i = 0; i < num->length && !fmpz_is_one(g); i++)
            fmpz_gcd(g, g, num->coeffs + i);
    }
    if (fmpz_sgn(dc) < 0)
        fmpz_neg(g, g);
    if (!fmpz_is_one(g))
    {
        fmpz_mpoly_scalar_divexact_fmpz(num, num, g, F->mctx);
        fmpz_divexact(dc, dc, g);
        fmpz_mpoly_set_fmpz(den, dc, F->mctx);
    }
    fmpz_clear(g);
    fmpz_clear(dc);
}

typedef struct
{
    ulong version;      /* ideal version of F */
    ulong layout;       /* layout version of F */
    slong m;
    slong W;            /* the number of words of the matrices V_k and V_k^(-1) */
    slong R;            /* roots per prime: d_1 + ... + d_m */
    ulong step;         /* the candidates are 1 mod step */
    ulong next;         /* the next candidate (decreasing) */
    slong num;
    slong alloc;
    ulong * primes;
    ulong * roots;      /* for each prime, the roots of the moduli modulo it */
    ulong * cyc;        /* cyc[k] = n if the modulus of step k is the cyclotomic polynomial of order n, else 0 */
    int exhausted;
}
modp_struct;

void
_gr_tower_flat_modp_clear(gr_tower_flat_t F)
{
    modp_struct * M = (modp_struct *) F->modp;
    if (M != NULL)
    {
        flint_free(M->primes);
        flint_free(M->roots);
        flint_free(M->cyc);
        flint_free(M);
        F->modp = NULL;
    }
}

static modp_struct *
_modp_get(gr_tower_flat_t F)
{
    slong m = F->T->length;
    modp_struct * M = (modp_struct *) F->modp;
    slong k;
    ulong top;

    if (M != NULL && (M->version != F->ideal_version || M->layout != F->layout_version || M->m != m))
    {
        _gr_tower_flat_modp_clear(F);
        M = NULL;
    }

    if (M == NULL)
    {
        M = flint_malloc(sizeof(modp_struct));
        M->version = F->ideal_version;
        M->layout = F->layout_version;
        M->m = m;
        M->W = 0;
        M->R = 0;
        /* (a heuristic for the search: n-th roots of unity and radicals
           split only if p = 1 mod n; the roots are checked anyway) */
        M->step = 2;
        M->cyc = flint_calloc(m, sizeof(ulong));
        for (k = 0; k < m; k++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_STEP(F->T, k);
            const fmpz_poly_struct * Mk = F->ideal_univar[k];
            slong d = gr_tower_step_degree(F->T, k + 1);
            ulong n = 1;

            M->W += 2 * d * d;
            M->R += d;
            if (g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param > 1)
            {
                n = g->def_param;
                if (n_euler_phi(n) == (ulong) d)
                {
                    fmpz_poly_t cp;
                    fmpz_poly_init(cp);
                    fmpz_poly_cyclotomic(cp, n);
                    if (fmpz_poly_equal(cp, Mk))
                        M->cyc[k] = n;
                    fmpz_poly_clear(cp);
                }
            }
            else if (Mk->length == d + 1 && _fmpz_vec_is_zero(Mk->coeffs + 1, d - 1))
            {
                /* x^d - c: the d-th roots of unity must be in F_p; for
                   d = 2, p = 1 mod 8|c| makes c a square (quadratic
                   reciprocity) */
                n = d;
                if (d == 2 && fmpz_bits(Mk->coeffs) <= FLINT_BITS / 4)
                    n = 8 * FLINT_ABS(fmpz_get_si(Mk->coeffs));
            }
            if (n > 1)
            {
                ulong t;
                t = M->step / n_gcd(M->step, n);
                if (n < (UWORD(1) << (FLINT_BITS / 4)) && t * (double) n < (double) (UWORD(1) << (FLINT_BITS / 2)))
                    M->step = t * n;
            }
        }
        /* (primes below 2^(FLINT_BITS - 4): dot products of moderate
           length accumulate in two words) */
        top = UWORD(1) << (FLINT_BITS - 4);
        M->next = ((top - 1) / M->step) * M->step + 1;
        M->num = 0;
        M->alloc = 0;
        M->primes = NULL;
        M->roots = NULL;
        M->exhausted = 0;
        F->modp = M;
    }

    return M;
}

/*
    Dense coefficient vectors of numerators in the variable v alone
    (lexicographic contexts, one-word exponents): the entries of A must be
    zero on input; returns the length (0 for zero), or -1 if x has a term
    of degree >= len or the layout does not apply.
*/
static slong
_univar_unpack(fmpz * A, slong len, const fmpz_mpoly_t x, slong v, const fmpz_mpoly_ctx_t mctx)
{
    slong N, off, shift, t, n = 0;
    ulong mask, e;

    if (x->length == 0)
        return 0;
    if (x->bits > FLINT_BITS || mctx->minfo->ord != ORD_LEX)
        return -1;

    N = mpoly_words_per_exp_sp(x->bits, mctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, v, x->bits, mctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - x->bits);

    for (t = 0; t < x->length; t++)
    {
        e = (x->exps[N * t + off] >> shift) & mask;
        if (e >= (ulong) len)
            return -1;
        fmpz_set(A + e, x->coeffs + t);
        if ((slong) e + 1 > n)
            n = e + 1;
    }
    return n;
}

/* res = sum A[e] x_v^e (the entries of A are moved out: zero on return) */
static void
_univar_pack(fmpz_mpoly_t res, fmpz * A, slong len, slong v, const fmpz_mpoly_ctx_t mctx)
{
    slong N, off, shift, e, k, nnz = 0;
    flint_bitcnt_t bits;

    for (e = 0; e < len; e++)
        nnz += !fmpz_is_zero(A + e);

    bits = mpoly_fix_bits(FLINT_MAX(MPOLY_MIN_BITS, 1 + FLINT_BIT_COUNT(len)), mctx->minfo);
    fmpz_mpoly_fit_length_reset_bits(res, nnz, bits, mctx);
    N = mpoly_words_per_exp_sp(bits, mctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, v, bits, mctx->minfo);

    for (e = len - 1, k = 0; e >= 0; e--)
    {
        if (fmpz_is_zero(A + e))
            continue;
        mpoly_monomial_zero(res->exps + N * k, N);
        res->exps[N * k + off] = ((ulong) e) << shift;
        fmpz_swap(res->coeffs + k, A + e);
        k++;
    }
    _fmpz_mpoly_set_length(res, k, mctx);
}

#define UNIVAR_STACK 64

/* a single generator: dense arithmetic on coefficient vectors (on the
   stack for small degrees), the product folded by the (sparse or dense)
   monic modulus */
static int
_mul_dense_univar(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_tower_flat_t F)
{
    slong v = GR_TOWER_FLAT_VAR(F, 1), d = gr_tower_step_degree(F->T, 1), i, j, n, lenA, lenB, lenC;
    const fmpz_poly_struct * M = F->ideal_univar[0];
    fmpz Abuf[UNIVAR_STACK], Bbuf[UNIVAR_STACK], Cbuf[2 * UNIVAR_STACK];
    fmpz * A, * B, * C;
    fmpz_t dx, dy;
    int heap = (d > UNIVAR_STACK), ok;

    if (heap)
    {
        A = _fmpz_vec_init(d);
        B = _fmpz_vec_init(d);
        C = _fmpz_vec_init(2 * d);
    }
    else
    {
        A = Abuf;
        B = Bbuf;
        C = Cbuf;
        for (i = 0; i < d; i++)
        {
            fmpz_init(A + i);
            fmpz_init(B + i);
        }
        for (i = 0; i < 2 * d; i++)
            fmpz_init(C + i);
    }

    /* (the numerators involve the variable of the generator only: with a
       single algebraic generator and no transcendental one, the other
       variables of the context are unused) */
    lenA = _univar_unpack(A, d, fmpz_mpoly_q_numref(x), v, F->mctx);
    lenB = (lenA > 0) ? _univar_unpack(B, d, fmpz_mpoly_q_numref(y), v, F->mctx) : -1;
    ok = (lenA > 0 && lenB > 0);

    if (ok)
    {
        fmpz_init(dx);
        fmpz_init(dy);
        fmpz_mpoly_get_fmpz(dx, fmpz_mpoly_q_denref(x), F->mctx);
        fmpz_mpoly_get_fmpz(dy, fmpz_mpoly_q_denref(y), F->mctx);
        fmpz_mul(dx, dx, dy);

        lenC = lenA + lenB - 1;
        if (fmpz_mpoly_q_numref(y)->length == 1 || fmpz_mpoly_q_numref(x)->length == 1)
        {
            /* a monomial factor (a power of a root of unity): a shift */
            const fmpz * big = (fmpz_mpoly_q_numref(y)->length == 1) ? A : B;
            const fmpz * mon = (big == A) ? B : A;
            slong lenbig = (big == A) ? lenA : lenB;
            slong e = ((big == A) ? lenB : lenA) - 1;
            _fmpz_vec_scalar_mul_fmpz(C + e, big, lenbig, mon + e);
        }
        else if (lenA >= lenB)
            _fmpz_poly_mul(C, A, lenA, B, lenB);
        else
            _fmpz_poly_mul(C, B, lenB, A, lenA);

        n = lenC;
        if (n > d && d >= DENSE_REM_CUTOFF && _nnz(M) > DENSE_REM_NNZ)
        {
            /* a large dense modulus: fast division (the sequential
               folding, like the precomputed powers of nf_elem, is
               quadratic; measured faster below degree 64 for random
               dense moduli with small coefficients) */
            /* (R of length n, as _fmpz_poly_rem requires) */
            fmpz * R = _fmpz_vec_init(n);
            _fmpz_poly_rem(R, C, n, M->coeffs, M->length);
            _fmpz_vec_swap(C, R, d);
            _fmpz_vec_zero(C + d, n - d);
            _fmpz_vec_clear(R, n);
            n = d;
        }
        else if (_is_all_ones(M) && n > d)
        {
            /* x^(d+1) = 1, then x^d = -(1 + ... + x^(d-1)) */
            for (i = n - 1; i >= d + 1; i--)
            {
                fmpz_add(C + i - d - 1, C + i - d - 1, C + i);
                fmpz_zero(C + i);
            }
            if (!fmpz_is_zero(C + d))
            {
                for (i = 0; i < d; i++)
                    fmpz_sub(C + i, C + i, C + d);
                fmpz_zero(C + d);
            }
            n = d;
        }
        for (i = n - 1; i >= d; i--)
        {
            if (fmpz_is_zero(C + i))
                continue;
            for (j = 0; j < d; j++)
                if (!fmpz_is_zero(M->coeffs + j))
                    fmpz_submul(C + i - d + j, C + i, M->coeffs + j);
            fmpz_zero(C + i);
        }

        _univar_pack(fmpz_mpoly_q_numref(res), C, FLINT_MIN(n, d), v, F->mctx);
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), dx, F->mctx);
        _canonicalise_int_den(res, F);
        fmpz_clear(dx);
        fmpz_clear(dy);
    }

    if (heap)
    {
        _fmpz_vec_clear(A, d);
        _fmpz_vec_clear(B, d);
        _fmpz_vec_clear(C, 2 * d);
    }
    else
    {
        for (i = 0; i < d; i++)
        {
            fmpz_clear(A + i);
            fmpz_clear(B + i);
        }
        for (i = 0; i < 2 * d; i++)
            fmpz_clear(C + i);
    }
    return ok;
}

/*
    Whether the polynomial Y (reduced) is c times a product of powers
    a_k^(e_k) of the generators with 0 <= e_k < q_k for roots of unity of
    order q_k (0 <= e_k < d_k otherwise): a monomial, or a tensor product
    in which the power of a root of unity of odd prime power order
    q = p^f with e >= d = q - q/p appears reduced as
    -a^(e - d) (1 + a^(q/p) + ... + a^((p-2) q/p)), as zeta_n^j does in the
    canonical decomposition (zeta_105 = zeta_3^2 zeta_5^2 zeta_7^3 has
    eight terms). Sets ey and c.
*/
#define ROOT_MONOMIAL_MAX_TERMS 256

static int
_root_monomial(slong * ey, fmpz_t c, const fmpz_mpoly_t Y, const slong * var, const slong * deg, slong m, gr_tower_flat_t F)
{
    slong t = Y->length, nvars = F->mctx->minfo->nvars, i, j, k, prod = 1, nexp = 0;
    ulong * e, * E;
    int ok = 1;

    if (t == 0 || t > ROOT_MONOMIAL_MAX_TERMS)
        return 0;
    for (i = 1; i < t; i++)
        if (!fmpz_equal(Y->coeffs + i, Y->coeffs))
            return 0;

    e = flint_malloc(sizeof(ulong) * nvars);
    E = flint_malloc(sizeof(ulong) * t * m);
    for (i = 0; i < t && ok; i++)
    {
        fmpz_mpoly_get_term_exp_ui(e, Y, i, F->mctx);
        for (j = 0; j < nvars && ok; j++)
        {
            for (k = 0; k < m; k++)
                if (var[k] == j)
                    break;
            if (k == m)
            {
                if (e[j] != 0)
                    ok = 0;
            }
            else
            {
                E[i * m + k] = e[j];
                if (e[j] >= (ulong) deg[k])
                    ok = 0;
            }
        }
    }
    flint_free(e);

    if (ok && t == 1)
    {
        for (k = 0; k < m; k++)
            ey[k] = E[k];
        fmpz_set(c, Y->coeffs);
        flint_free(E);
        return 1;
    }

    if (ok)
    {
        modp_struct * M = _modp_get(F);

        for (k = 0; k < m && ok; k++)
        {
            ulong lo = E[k], q = M->cyc[k], p, s, v;
            slong nk = 0, ndist;

            for (i = 1; i < t; i++)
                lo = FLINT_MIN(lo, E[i * m + k]);

            /* the distinct values along the axis */
            ndist = 0;
            {
                ulong hi = lo;
                for (i = 0; i < t; i++)
                    hi = FLINT_MAX(hi, E[i * m + k]);
                if (hi == lo)
                    ndist = 1;
            }

            if (ndist == 1)
            {
                ey[k] = lo;
                continue;
            }

            /* an expanded power: q = p^f odd */
            if (q == 0 || q % 2 == 0)
            {
                ok = 0;
                break;
            }
            {
                n_factor_t fac;
                n_factor_init(&fac);
                n_factor(&fac, q, 1);
                if (fac.num != 1)
                {
                    ok = 0;
                    break;
                }
                p = fac.p[0];
            }
            s = q / p;
            for (i = 0; i < t && ok; i++)
            {
                v = E[i * m + k] - lo;
                if (v % s != 0 || v / s > p - 2)
                    ok = 0;
            }
            /* all p - 1 values occur: counted through the product below */
            nk = p - 1;
            prod *= nk;
            nexp++;
            ey[k] = lo + deg[k];
            if (ey[k] >= (slong) q)
                ok = 0;
        }

        /* (distinct terms within the product set, as many as its size:
           the whole tensor product) */
        if (ok && prod != t)
            ok = 0;

        if (ok)
        {
            fmpz_set(c, Y->coeffs);
            if (nexp % 2)
                fmpz_neg(c, c);
        }
    }

    flint_free(E);
    return ok;
}

/*
    x times the monomial y = c * prod a_k^(e_k) (several generators): on
    the compact array of x (D entries, no Kronecker room), each fibre
    along axis k is multiplied by a_k^(e_k) and reduced by M_k (for
    1 + x + ... + x^d, as x^(d+1) = 1 then one subtraction). O(D) per
    axis for sparse moduli, against O(L) = O(D prod (2 - 1/d_k)) for the
    shift in the Kronecker array followed by the folds.
*/
static int
_mul_monomial_compact(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong m = T->length, k, j, D = 1, maxd = 1;
    slong * var, * deg, * C, * ey;
    fmpz * A, * t;
    fmpz_t c, dx;
    int ok;
    TMP_INIT;

    /* (most products are not by monomials: a cheap exclusion first) */
    if (fmpz_mpoly_q_numref(y)->length > 1 && !_any_root_of_unity(F))
        return 0;

    TMP_START;
    var = TMP_ALLOC(sizeof(slong) * 4 * m);
    deg = var + m;
    C = var + 2 * m;
    ey = var + 3 * m;
    for (k = 0; k < m; k++)
    {
        var[k] = GR_TOWER_FLAT_VAR(F, k + 1);
        deg[k] = gr_tower_step_degree(T, k + 1);
        C[k] = D;
        D *= deg[k];
        maxd = FLINT_MAX(maxd, deg[k]);
    }

    /* y = c prod a_k^(e_k) */
    fmpz_init(c);
    ok = _root_monomial(ey, c, fmpz_mpoly_q_numref(y), var, deg, m, F);
    if (!ok)
    {
        fmpz_clear(c);
        TMP_END;
        return 0;
    }

    A = _fmpz_vec_init(D);
    t = _fmpz_vec_init(3 * maxd);
    ok = ok && _pack(A, fmpz_mpoly_q_numref(x), var, deg, C, m, F);

    for (k = 0; k < m && ok; k++)
    {
        const fmpz_poly_struct * M = F->ideal_univar[k];
        slong d = deg[k], a = ey[k], blk = C[k] * d, base, off, i;
        int ones = _is_all_ones(M);

        if (a == 0)
            continue;

        for (base = 0; base < D; base += blk)
        {
            for (off = 0; off < C[k]; off++)
            {
                fmpz * f = A + base + off;

                /* t = fibre * x^a */
                for (i = 0; i < d; i++)
                    fmpz_swap(t + a + i, f + i * C[k]);

                if (ones)
                {
                    for (i = d + a - 1; i >= d + 1; i--)
                    {
                        fmpz_add(t + i - d - 1, t + i - d - 1, t + i);
                        fmpz_zero(t + i);
                    }
                    if (!fmpz_is_zero(t + d))
                    {
                        for (i = 0; i < d; i++)
                            fmpz_sub(t + i, t + i, t + d);
                        fmpz_zero(t + d);
                    }
                }
                else
                {
                    for (i = d + a - 1; i >= d; i--)
                    {
                        if (fmpz_is_zero(t + i))
                            continue;
                        for (j = 0; j < d; j++)
                            if (!fmpz_is_zero(M->coeffs + j))
                                fmpz_submul(t + i - d + j, t + i, M->coeffs + j);
                        fmpz_zero(t + i);
                    }
                }

                for (i = 0; i < d; i++)
                    fmpz_swap(t + i, f + i * C[k]);
            }
        }
    }

    if (ok)
    {
        if (!fmpz_is_one(c))
            _fmpz_vec_scalar_mul_fmpz(A, A, D, c);
        fmpz_init(dx);
        fmpz_mpoly_get_fmpz(dx, fmpz_mpoly_q_denref(x), F->mctx);
        fmpz_mpoly_get_fmpz(c, fmpz_mpoly_q_denref(y), F->mctx);
        fmpz_mul(dx, dx, c);
        _unpack(fmpz_mpoly_q_numref(res), A, var, deg, C, m, F);
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), dx, F->mctx);
        _canonicalise_int_den(res, F);
        fmpz_clear(dx);
    }

    _fmpz_vec_clear(A, D);
    _fmpz_vec_clear(t, 3 * maxd);
    fmpz_clear(c);
    TMP_END;
    return ok;
}

/*
    Products in a compositum of cyclotomic fields Q(zeta_q1, ..., zeta_qm)
    with pairwise coprime orders q_k (the canonical decomposition of
    Q(zeta_n), n = q_1 ... q_m): with x_k -> x^(c_k), c_k the CRT
    idempotents (c_k = 1 mod q_k, 0 mod q_j), the tensor product of the
    Z[x_k] / (x_k^q_k - 1) is Z[x] / (x^n - 1), so the product is one
    cyclic convolution of length n (Good-Thomas), mapped back digit by
    digit (e_k = e mod q_k) and reduced by the Phi_(q_k). The Kronecker
    array has prod (2 d_k - 1) entries instead: 4389 against 1155 for
    Q(zeta_1155).
*/

static int
_mul_cyclic(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_tower_flat_t F)
{
    modp_struct * M = _modp_get(F);
    slong m = F->T->length, k, j, n = 1, lenA, lenB, e, idx;
    slong * var, * deg, * c, * Q, * r;
    fmpz * A, * B, * C;
    fmpz_t dx, dy;

    for (k = 0; k < m; k++)
    {
        if (M->cyc[k] == 0)
            return 0;
        for (j = 0; j < k; j++)
            if (n_gcd(M->cyc[j], M->cyc[k]) != 1)
                return 0;
        n *= M->cyc[k];
        if (n > CYCLIC_LIMIT)
            return 0;
    }

    /* (only if shorter than the Kronecker array) */
    {
        slong L = 1;
        for (k = 0; k < m; k++)
            L *= 2 * gr_tower_step_degree(F->T, k + 1) - 1;
        if (2 * n - 1 >= L)
            return 0;
    }

    var = flint_malloc(sizeof(slong) * 5 * m + sizeof(slong));
    deg = var + m;
    c = var + 2 * m;
    r = var + 3 * m;
    Q = var + 4 * m;
    Q[0] = 1;
    for (k = 0; k < m; k++)
    {
        ulong q = M->cyc[k], u = n / q;
        var[k] = GR_TOWER_FLAT_VAR(F, k + 1);
        deg[k] = gr_tower_step_degree(F->T, k + 1);
        c[k] = u * n_invmod(u % q, q);
        Q[k + 1] = Q[k] * q;
    }

    A = _fmpz_vec_init(n);
    B = _fmpz_vec_init(n);
    if (!_pack_ex(A, fmpz_mpoly_q_numref(x), var, deg, c, m, F, n) ||
        !_pack_ex(B, fmpz_mpoly_q_numref(y), var, deg, c, m, F, n))
    {
        _fmpz_vec_clear(A, n);
        _fmpz_vec_clear(B, n);
        flint_free(var);
        return 0;
    }

    lenA = n;
    lenB = n;
    while (lenA > 0 && fmpz_is_zero(A + lenA - 1)) lenA--;
    while (lenB > 0 && fmpz_is_zero(B + lenB - 1)) lenB--;
    C = _fmpz_vec_init(2 * n);
    if (lenA >= lenB)
        _fmpz_poly_mul(C, A, lenA, B, lenB);
    else
        _fmpz_poly_mul(C, B, lenB, A, lenA);
    for (e = 0; e + n < lenA + lenB - 1; e++)
    {
        fmpz_add(C + e, C + e, C + e + n);
        fmpz_zero(C + e + n);
    }

    /* back to the digits: A[sum (e mod q_k) Q_k] = C[e] */
    for (k = 0; k < m; k++)
        r[k] = 0;
    idx = 0;
    for (e = 0; e < n; e++)
    {
        fmpz_swap(A + idx, C + e);
        for (k = 0; k < m; k++)
        {
            r[k]++;
            idx += Q[k];
            if (r[k] == (slong) M->cyc[k])
            {
                r[k] = 0;
                idx -= Q[k + 1];
            }
        }
    }

    _fmpz_vec_zero(C, 2 * n);   /* (holds the old entries of A) */
    _fold_array(A, n, deg, Q, m, F);

    fmpz_init(dx);
    fmpz_init(dy);
    fmpz_mpoly_get_fmpz(dx, fmpz_mpoly_q_denref(x), F->mctx);
    fmpz_mpoly_get_fmpz(dy, fmpz_mpoly_q_denref(y), F->mctx);
    fmpz_mul(dx, dx, dy);
    _unpack(fmpz_mpoly_q_numref(res), A, var, deg, Q, m, F);
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), dx, F->mctx);
    _canonicalise_int_den(res, F);
    fmpz_clear(dx);
    fmpz_clear(dy);

    _fmpz_vec_clear(A, n);
    _fmpz_vec_clear(B, n);
    _fmpz_vec_clear(C, 2 * n);
    flint_free(var);
    return 1;
}

/*
    Whether the product of x and y should rather be computed sparsely: with
    several generators, the Kronecker array has L = prod (2 d_k - 1)
    entries whatever the number of terms, while the sparse product has
    len(x) len(y) terms before the reduction (8-sqrt tower, 12-term
    elements: 4.8 ms dense against ~0.05 ms sparse).
*/
int
_gr_tower_flat_mul_prefers_sparse(const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_tower_flat_t F)
{
    slong k, L = 1, m = F->T->length;
    double t;

    /* (a monomial factor is a shift of the packed array) */
    if (m <= 1 || fmpz_mpoly_q_numref(x)->length <= 1 || fmpz_mpoly_q_numref(y)->length <= 1)
        return 0;
    for (k = 1; k <= m; k++)
        L *= 2 * gr_tower_step_degree(F->T, k) - 1;
    t = (double) fmpz_mpoly_q_numref(x)->length * (double) fmpz_mpoly_q_numref(y)->length;
    return 8 * t < (double) L;
}

int
_gr_tower_flat_mul_dense(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong m = T->length, k, L, Lin;
    slong * var, * deg, * S;
    fmpz * A, * B;
    fmpz_t dx, dy;
    int ok;
    TMP_INIT;

    if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), F->mctx) ||
        !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(y), F->mctx))
        return 0;

    if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), F->mctx) || fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(y), F->mctx))
    {
        fmpz_mpoly_q_zero(res, F->mctx);
        return 1;
    }

    if (m == 1)
        return _mul_dense_univar(res, x, y, F);

    /* (a power of a root of unity, or a monomial: fibre by fibre) */
    if (fmpz_mpoly_q_numref(y)->length <= fmpz_mpoly_q_numref(x)->length)
    {
        if (_mul_monomial_compact(res, x, y, F) || _mul_monomial_compact(res, y, x, F))
            return 1;
    }
    else
    {
        if (_mul_monomial_compact(res, y, x, F) || _mul_monomial_compact(res, x, y, F))
            return 1;
    }

    if (_mul_cyclic(res, x, y, F))
        return 1;

    TMP_START;
    var = TMP_ALLOC(sizeof(slong) * m);
    deg = TMP_ALLOC(sizeof(slong) * m);
    S = TMP_ALLOC(sizeof(slong) * (m + 1));
    S[0] = 1;
    Lin = 1;
    for (k = 0; k < m; k++)
    {
        var[k] = GR_TOWER_FLAT_VAR(F, k + 1);
        deg[k] = gr_tower_step_degree(T, k + 1);
        S[k + 1] = S[k] * (2 * deg[k] - 1);
        Lin += (deg[k] - 1) * S[k];
    }
    L = S[m];

    A = _tmp_vec_init(Lin, TMP_ALLOC(sizeof(fmpz) * Lin));
    B = _tmp_vec_init(Lin, TMP_ALLOC(sizeof(fmpz) * Lin));
    ok = _pack(A, fmpz_mpoly_q_numref(x), var, deg, S, m, F) &&
         _pack(B, fmpz_mpoly_q_numref(y), var, deg, S, m, F);

    if (ok)
    {
        fmpz * C;
        slong lenA = Lin, lenB = Lin, lenC;

        while (lenA > 0 && fmpz_is_zero(A + lenA - 1)) lenA--;
        while (lenB > 0 && fmpz_is_zero(B + lenB - 1)) lenB--;
        lenC = lenA + lenB - 1;

        C = _tmp_vec_init(FLINT_MAX(L, lenC), TMP_ALLOC(sizeof(fmpz) * FLINT_MAX(L, lenC)));
        if (fmpz_mpoly_q_numref(y)->length == 1 || fmpz_mpoly_q_numref(x)->length == 1)
        {
            /* a monomial factor: a shift of the packed array */
            const fmpz * big = (fmpz_mpoly_q_numref(y)->length == 1) ? A : B;
            const fmpz * mon = (big == A) ? B : A;
            slong lenbig = (big == A) ? lenA : lenB;
            slong e = ((big == A) ? lenB : lenA) - 1;
            _fmpz_vec_scalar_mul_fmpz(C + e, big, lenbig, mon + e);
        }
        else if (lenA >= lenB)
            _fmpz_poly_mul(C, A, lenA, B, lenB);
        else
            _fmpz_poly_mul(C, B, lenB, A, lenA);

        fmpz_init(dx);
        fmpz_init(dy);
        fmpz_mpoly_get_fmpz(dx, fmpz_mpoly_q_denref(x), F->mctx);
        fmpz_mpoly_get_fmpz(dy, fmpz_mpoly_q_denref(y), F->mctx);
        fmpz_mul(dx, dx, dy);

        _reduce_unpack(fmpz_mpoly_q_numref(res), C, L, var, deg, S, m, F);
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), dx, F->mctx);
        _canonicalise_int_den(res, F);

        fmpz_clear(dx);
        fmpz_clear(dy);
        _tmp_vec_clear(C, FLINT_MAX(L, lenC));
    }

    _tmp_vec_clear(A, Lin);
    _tmp_vec_clear(B, Lin);
    TMP_END;
    return ok;
}

/*
    Inverse: in a single generator, the extended Euclidean algorithm over
    Q with the modulus (a closed formula for a quadratic modulus); with
    several generators, a linear system with the multiplication matrix. Returns 1 on success (0 if not applicable, or if
    the element is a zero divisor).
*/
static int _inv_dense_multi(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t F);

static int _inv_dense_impl(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t F);
static int _inv_dense_modular(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t F);

/* fields of at least this degree cache their last inverse */
#define INV_CACHE_DEGREE 8

/*
    (Fraction-free elimination divides many elements by the same pivot:
    the last inverse is kept, keyed by the element and the versions of
    the context, the ideal and the variable layout.)
*/
int
_gr_tower_flat_inv_dense(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    slong D;
    fmpz_mpoly_q_t xc;
    int ok;

    if (F->inv_cache != NULL && F->inv_cache_mctx == F->mctx &&
        F->inv_cache_version[0] == F->ideal_version &&
        F->inv_cache_version[1] == F->layout_version &&
        fmpz_mpoly_q_equal(F->inv_cache, x, F->mctx))
    {
        fmpz_mpoly_q_set(res, F->inv_cache + 1, F->mctx);
        return 1;
    }

    D = gr_tower_degree(F->T);

    if (D < INV_CACHE_DEGREE)
        return _inv_dense_impl(res, x, F);

    fmpz_mpoly_q_init(xc, F->mctx);
    fmpz_mpoly_q_set(xc, x, F->mctx);
    ok = _inv_dense_impl(res, x, F);
    if (ok)
    {
        if (F->inv_cache == NULL)
            F->inv_cache = flint_malloc(2 * sizeof(fmpz_mpoly_q_struct));
        else
        {
            fmpz_mpoly_q_clear(F->inv_cache, F->inv_cache_mctx);
            fmpz_mpoly_q_clear(F->inv_cache + 1, F->inv_cache_mctx);
        }
        F->inv_cache_mctx = F->mctx;
        F->inv_cache_version[0] = F->ideal_version;
        F->inv_cache_version[1] = F->layout_version;
        fmpz_mpoly_q_init(F->inv_cache, F->mctx);
        fmpz_mpoly_q_init(F->inv_cache + 1, F->mctx);
        fmpz_mpoly_q_swap(F->inv_cache, xc, F->mctx);
        fmpz_mpoly_q_set(F->inv_cache + 1, res, F->mctx);
    }
    fmpz_mpoly_q_clear(xc, F->mctx);
    return ok;
}

static int
_inv_dense_impl(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong v;
    fmpz_poly_t P;
    fmpq_poly_t A, M, G, S, U;
    fmpz_t d;
    int ok;

    if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), F->mctx) ||
        fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), F->mctx))
        return 0;

    /* several generators: modulo primes at which the moduli split, else
       the linear system with the multiplication matrix */
    if (T->length > 1)
    {
        slong alg = GR_TOWER_OPTION(T, GR_TOWER_OPT_INV_DENSE_ALG);
        if (alg != 2 && _inv_dense_modular(res, x, F))
            return 1;
        if (alg == 1)
            return 0;
        return _inv_dense_multi(res, x, F);
    }

    if (GR_TOWER_OPTION(T, GR_TOWER_OPT_INV_DENSE_ALG) == 1)
        return _inv_dense_modular(res, x, F);

    v = GR_TOWER_FLAT_VAR(F, 1);
    {
        /* (the numerator must involve the variable of the generator only) */
        int * used = flint_malloc(sizeof(int) * F->mctx->minfo->nvars);
        slong i;
        ok = 1;
        fmpz_mpoly_used_vars(used, fmpz_mpoly_q_numref(x), F->mctx);
        for (i = 0; i < F->mctx->minfo->nvars; i++)
            if (used[i] && i != v)
                ok = 0;
        flint_free(used);
        if (!ok)
            return 0;
    }

    fmpz_poly_init(P);
    fmpq_poly_init(A);
    fmpq_poly_init(M);
    fmpq_poly_init(G);
    fmpq_poly_init(S);
    fmpq_poly_init(U);
    fmpz_init(d);

    fmpz_mpoly_get_fmpz_poly(P, fmpz_mpoly_q_numref(x), v, F->mctx);
    fmpz_mpoly_get_fmpz(d, fmpz_mpoly_q_denref(x), F->mctx);

    /* a quadratic modulus t^2 + p t + q: (a + b t)^(-1) =
       (a - b p - b t) / (a^2 - a b p + b^2 q) */
    if (fmpz_poly_degree(F->ideal_univar[0]) == 2 && fmpz_poly_length(P) <= 2)
    {
        const fmpz * M2 = F->ideal_univar[0]->coeffs;
        fmpz_t a, b, n, t;
        fmpz_poly_t R;
        fmpz_init(a); fmpz_init(b); fmpz_init(n); fmpz_init(t);
        fmpz_poly_init(R);
        fmpz_poly_get_coeff_fmpz(a, P, 0);
        fmpz_poly_get_coeff_fmpz(b, P, 1);
        /* n = a (a - b p) + b^2 q */
        fmpz_mul(t, b, M2 + 1);
        fmpz_sub(t, a, t);                 /* a - b p */
        fmpz_mul(n, a, t);
        fmpz_mul(a, b, b);
        fmpz_addmul(n, a, M2 + 0);
        ok = !fmpz_is_zero(n);
        if (ok)
        {
            /* 1/x = d (a - b p - b t) / n */
            fmpz_poly_set_coeff_fmpz(R, 0, t);
            fmpz_neg(b, b);
            fmpz_poly_set_coeff_fmpz(R, 1, b);
            fmpz_poly_scalar_mul_fmpz(R, R, d);
            if (fmpz_sgn(n) < 0)
            {
                fmpz_neg(n, n);
                fmpz_poly_neg(R, R);
            }
            /* (common content) */
            fmpz_poly_content(t, R);
            fmpz_gcd(t, t, n);
            if (!fmpz_is_one(t))
            {
                fmpz_poly_scalar_divexact_fmpz(R, R, t);
                fmpz_divexact(n, n, t);
            }
            fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_q_numref(res), R, v, F->mctx);
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), n, F->mctx);
        }
        fmpz_clear(a); fmpz_clear(b); fmpz_clear(n); fmpz_clear(t);
        fmpz_poly_clear(R);
        fmpz_poly_clear(P);
        fmpq_poly_clear(A);
        fmpq_poly_clear(M);
        fmpq_poly_clear(G);
        fmpq_poly_clear(S);
        fmpq_poly_clear(U);
        fmpz_clear(d);
        return ok;
    }

    fmpq_poly_set_fmpz_poly(A, P);
    fmpq_poly_set_fmpz_poly(M, F->ideal_univar[0]);

    ok = (fmpq_poly_degree(A) < fmpq_poly_degree(M));
    if (ok)
    {
        fmpq_poly_xgcd(G, S, U, A, M);
        ok = fmpq_poly_is_one(G);
    }

    if (ok)
    {
        /* 1/x = d S */
        fmpq_poly_scalar_mul_fmpz(S, S, d);
        fmpq_poly_get_numerator(P, S);
        fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_q_numref(res), P, v, F->mctx);
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), fmpq_poly_denref(S), F->mctx);
    }

    fmpz_poly_clear(P);
    fmpq_poly_clear(A);
    fmpq_poly_clear(M);
    fmpq_poly_clear(G);
    fmpq_poly_clear(S);
    fmpq_poly_clear(U);
    fmpz_clear(d);
    return ok;
}

/* ------------------------------------------------------------------------ */
/* polynomials and matrices                                                   */
/* ------------------------------------------------------------------------ */

typedef struct
{
    slong m;
    slong * var;
    slong * deg;
    slong * S;      /* S[k]: stride of the variable of step k + 1; S[m] = L */
    slong L;        /* extent of a product of two reduced elements */
    slong Lin;      /* extent of a reduced element */
    slong cyc_n;    /* n > 0: the cyclic (Good-Thomas) layout of length n (Lin = n, L = 2 n - 1) */
    slong * c;      /* (cyclic layout) the CRT idempotents: strides modulo n */
    slong * Q;      /* (cyclic layout) the strides of the digits after the product */
    ulong * q;      /* (cyclic layout) the orders */
}
dense_layout_struct;

static void
_layout_init(dense_layout_struct * D, gr_tower_flat_t F)
{
    slong k, m = F->T->length;
    D->m = m;
    D->var = flint_malloc(sizeof(slong) * m);
    D->deg = flint_malloc(sizeof(slong) * m);
    D->S = flint_malloc(sizeof(slong) * (m + 1));
    D->S[0] = 1;
    D->Lin = 1;
    for (k = 0; k < m; k++)
    {
        D->var[k] = GR_TOWER_FLAT_VAR(F, k + 1);
        D->deg[k] = gr_tower_step_degree(F->T, k + 1);
        D->S[k + 1] = D->S[k] * (2 * D->deg[k] - 1);
        D->Lin += (D->deg[k] - 1) * D->S[k];
    }
    D->L = D->S[m];
    D->cyc_n = 0;
    D->c = NULL;
    D->Q = NULL;
    D->q = NULL;
}

/* the layout for products: the cyclic one for composita of cyclotomic
   fields of coprime orders when it is shorter (see _mul_cyclic) */
static void
_layout_init_products(dense_layout_struct * D, gr_tower_flat_t F)
{
    modp_struct * M;
    slong k, j, n = 1;

    _layout_init(D, F);
    if (D->m < 2)
        return;

    M = _modp_get(F);
    for (k = 0; k < D->m; k++)
    {
        if (M->cyc[k] == 0)
            return;
        for (j = 0; j < k; j++)
            if (n_gcd(M->cyc[j], M->cyc[k]) != 1)
                return;
        n *= M->cyc[k];
        if (n > CYCLIC_LIMIT)
            return;
    }
    if (2 * n - 1 >= D->L)
        return;

    D->cyc_n = n;
    D->c = flint_malloc(sizeof(slong) * D->m);
    D->Q = flint_malloc(sizeof(slong) * (D->m + 1));
    D->q = flint_malloc(sizeof(ulong) * D->m);
    D->Q[0] = 1;
    for (k = 0; k < D->m; k++)
    {
        ulong qk = M->cyc[k], u = n / qk;
        D->q[k] = qk;
        D->c[k] = u * n_invmod(u % qk, qk);
        D->Q[k + 1] = D->Q[k] * qk;
    }
    D->Lin = n;
    D->L = 2 * n - 1;
}

static void
_layout_clear(dense_layout_struct * D)
{
    flint_free(D->var);
    flint_free(D->deg);
    flint_free(D->S);
    flint_free(D->c);
    flint_free(D->Q);
    flint_free(D->q);
}

/* the lcm of the (integer) denominators of the n elements x[i * stride]
   into D; returns 0 if some denominator is not an integer */
static int
_lcm_dens(fmpz_t D, fmpz_mpoly_q_struct * const * x, slong n, gr_tower_flat_t F)
{
    slong i;
    fmpz_t d;
    fmpz_init(d);
    fmpz_one(D);
    for (i = 0; i < n; i++)
    {
        if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x[i]), F->mctx))
        {
            fmpz_clear(d);
            return 0;
        }
        fmpz_mpoly_get_fmpz(d, fmpz_mpoly_q_denref(x[i]), F->mctx);
        fmpz_lcm(D, D, d);
    }
    fmpz_clear(d);
    return 1;
}

/* packs x * D (D a multiple of the denominator of x) at A */
static int
_pack_scaled(fmpz * A, const fmpz_mpoly_q_t x, const fmpz_t D, const dense_layout_struct * DL, gr_tower_flat_t F)
{
    fmpz_t s;
    slong j;
    int ok;

    if (DL->cyc_n)
        ok = _pack_ex(A, fmpz_mpoly_q_numref(x), DL->var, DL->deg, DL->c, DL->m, F, DL->cyc_n);
    else
        ok = _pack(A, fmpz_mpoly_q_numref(x), DL->var, DL->deg, DL->S, DL->m, F);
    if (ok)
    {
        fmpz_init(s);
        fmpz_mpoly_get_fmpz(s, fmpz_mpoly_q_denref(x), F->mctx);
        fmpz_divexact(s, D, s);
        if (!fmpz_is_one(s))
            for (j = 0; j < DL->Lin; j++)
                fmpz_mul(A + j, A + j, s);
        fmpz_clear(s);
    }
    return ok;
}

/* res = the reduced element packed at C (extent L, destroyed) / den */
static void
_unpack_div(fmpz_mpoly_q_t res, fmpz * C, const fmpz_t den, const dense_layout_struct * DL, gr_tower_flat_t F)
{
    if (DL->cyc_n)
    {
        /* modulo x^n - 1, then back to the digits e mod q_k */
        slong n = DL->cyc_n, e, k, idx = 0;
        slong * r = flint_calloc(DL->m, sizeof(slong));
        fmpz * T = _fmpz_vec_init(n);
        for (e = 0; e < n - 1; e++)
            fmpz_add(C + e, C + e, C + e + n);
        for (e = 0; e < n; e++)
        {
            fmpz_swap(T + idx, C + e);
            for (k = 0; k < DL->m; k++)
            {
                r[k]++;
                idx += DL->Q[k];
                if (r[k] == (slong) DL->q[k])
                {
                    r[k] = 0;
                    idx -= DL->Q[k + 1];
                }
            }
        }
        _fold_array(T, n, DL->deg, DL->Q, DL->m, F);
        _unpack(fmpz_mpoly_q_numref(res), T, DL->var, DL->deg, DL->Q, DL->m, F);
        _fmpz_vec_clear(T, n);
        flint_free(r);
    }
    else
        _reduce_unpack(fmpz_mpoly_q_numref(res), C, DL->L, DL->var, DL->deg, DL->S, DL->m, F);
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), den, F->mctx);
    _canonicalise_int_den(res, F);
}

/*
    res[0 .. n-1] = the low n coefficients of the product of the
    polynomials with coefficients x[0 .. len1-1] and y[0 .. len2-1] (all
    in the current context of F): Kronecker substitution of the
    polynomial variable into the dense packing (stride L), a single
    integer polynomial multiplication, and reductions of the n product
    coefficients. res must not alias the inputs. Returns 0 if not
    applicable (nothing written).
*/
int
_gr_tower_flat_poly_mullow_dense(fmpz_mpoly_q_struct * const * res, fmpz_mpoly_q_struct * const * x, slong len1,
    fmpz_mpoly_q_struct * const * y, slong len2, slong n, gr_tower_flat_t F)
{
    dense_layout_struct DL;
    fmpz_t D1, D2, D;
    fmpz * A, * B, * C;
    slong i, L, lenA, lenB, lenC;
    int ok;

    if (len1 < 1 || len2 < 1 || n < 1 || !_gr_tower_flat_dense_applicable(F))
        return 0;

    _layout_init_products(&DL, F);
    L = DL.L;
    if ((len1 + len2) * (double) L > 1e8)
    {
        _layout_clear(&DL);
        return 0;
    }

    fmpz_init(D1);
    fmpz_init(D2);
    fmpz_init(D);

    ok = _lcm_dens(D1, x, len1, F) && _lcm_dens(D2, y, len2, F);

    A = _fmpz_vec_init(len1 * L);
    B = _fmpz_vec_init(len2 * L);
    for (i = 0; i < len1 && ok; i++)
        ok = _pack_scaled(A + i * L, x[i], D1, &DL, F);
    for (i = 0; i < len2 && ok; i++)
        ok = _pack_scaled(B + i * L, y[i], D2, &DL, F);

    if (ok)
    {
        slong nout = FLINT_MIN(n, len1 + len2 - 1);
        lenA = len1 * L;
        lenB = len2 * L;
        while (lenA > 0 && fmpz_is_zero(A + lenA - 1)) lenA--;
        while (lenB > 0 && fmpz_is_zero(B + lenB - 1)) lenB--;

        lenC = FLINT_MAX(n, len1 + len2) * L;
        C = _fmpz_vec_init(lenC);
        if (lenA > 0 && lenB > 0)
        {
            slong lp = FLINT_MIN(lenA + lenB - 1, nout * L);
            if (lenA >= lenB)
                _fmpz_poly_mullow(C, A, lenA, B, lenB, lp);
            else
                _fmpz_poly_mullow(C, B, lenB, A, lenA, lp);
        }

        fmpz_mul(D, D1, D2);
        for (i = 0; i < n; i++)
        {
            if (i < nout)
                _unpack_div(res[i], C + i * L, D, &DL, F);
            else
                fmpz_mpoly_q_zero(res[i], F->mctx);
        }
        _fmpz_vec_clear(C, lenC);
    }

    _fmpz_vec_clear(A, len1 * L);
    _fmpz_vec_clear(B, len2 * L);
    fmpz_clear(D1);
    fmpz_clear(D2);
    fmpz_clear(D);
    _layout_clear(&DL);
    return ok;
}

/*
    C = A B for matrices given as arrays of pointers (row-major) to
    elements in the current context of F: common denominators by rows of
    A and columns of B, the entries packed densely as integer polynomials,
    one fmpz_poly_mat multiplication, and reductions of the entries.
    C must not alias the inputs. Returns 0 if not applicable.
*/
int
_gr_tower_flat_mat_mul_dense(fmpz_mpoly_q_struct * const * C, fmpz_mpoly_q_struct * const * A,
    fmpz_mpoly_q_struct * const * B, slong r, slong s, slong c, gr_tower_flat_t F)
{
    dense_layout_struct DL;
    fmpz_poly_mat_t PA, PB, PC;
    fmpz * rowden, * colden;
    fmpz * tmp;
    fmpz_t den;
    slong i, j, L;
    int ok = 1;

    if (r < 1 || s < 1 || c < 1 || !_gr_tower_flat_dense_applicable(F))
        return 0;

    _layout_init_products(&DL, F);
    L = DL.L;

    rowden = _fmpz_vec_init(r);
    colden = _fmpz_vec_init(c);
    fmpz_init(den);

    for (i = 0; i < r && ok; i++)
    {
        fmpz_mpoly_q_struct ** row = flint_malloc(sizeof(fmpz_mpoly_q_struct *) * s);
        for (j = 0; j < s; j++)
            row[j] = A[i * s + j];
        ok = _lcm_dens(rowden + i, (fmpz_mpoly_q_struct * const *) row, s, F);
        flint_free(row);
    }
    for (j = 0; j < c && ok; j++)
    {
        fmpz_mpoly_q_struct ** col = flint_malloc(sizeof(fmpz_mpoly_q_struct *) * s);
        for (i = 0; i < s; i++)
            col[i] = B[i * c + j];
        ok = _lcm_dens(colden + j, (fmpz_mpoly_q_struct * const *) col, s, F);
        flint_free(col);
    }

    fmpz_poly_mat_init(PA, r, s);
    fmpz_poly_mat_init(PB, s, c);
    fmpz_poly_mat_init(PC, r, c);
    tmp = _fmpz_vec_init(L);

    for (i = 0; i < r && ok; i++)
        for (j = 0; j < s && ok; j++)
        {
            fmpz_poly_struct * p = fmpz_poly_mat_entry(PA, i, j);
            _fmpz_vec_zero(tmp, DL.Lin);
            ok = _pack_scaled(tmp, A[i * s + j], rowden + i, &DL, F);
            if (ok)
            {
                fmpz_poly_fit_length(p, DL.Lin);
                _fmpz_vec_set(p->coeffs, tmp, DL.Lin);
                _fmpz_poly_set_length(p, DL.Lin);
                _fmpz_poly_normalise(p);
            }
        }
    for (i = 0; i < s && ok; i++)
        for (j = 0; j < c && ok; j++)
        {
            fmpz_poly_struct * p = fmpz_poly_mat_entry(PB, i, j);
            _fmpz_vec_zero(tmp, DL.Lin);
            ok = _pack_scaled(tmp, B[i * c + j], colden + j, &DL, F);
            if (ok)
            {
                fmpz_poly_fit_length(p, DL.Lin);
                _fmpz_vec_set(p->coeffs, tmp, DL.Lin);
                _fmpz_poly_set_length(p, DL.Lin);
                _fmpz_poly_normalise(p);
            }
        }

    if (ok)
    {
        fmpz_poly_mat_mul(PC, PA, PB);

        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
            {
                fmpz_poly_struct * p = fmpz_poly_mat_entry(PC, i, j);
                _fmpz_vec_zero(tmp, L);
                _fmpz_vec_set(tmp, p->coeffs, FLINT_MIN(p->length, L));
                fmpz_mul(den, rowden + i, colden + j);
                _unpack_div(C[i * c + j], tmp, den, &DL, F);
            }
    }

    _fmpz_vec_clear(tmp, L);
    fmpz_poly_mat_clear(PA);
    fmpz_poly_mat_clear(PB);
    fmpz_poly_mat_clear(PC);
    _fmpz_vec_clear(rowden, r);
    _fmpz_vec_clear(colden, c);
    fmpz_clear(den);
    _layout_clear(&DL);
    return ok;
}

/* ------------------------------------------------------------------------ */
/* inverses modulo primes at which the moduli split                           */
/* ------------------------------------------------------------------------ */

/*
    With all moduli split into distinct linear factors modulo p, the
    algebra modulo p is F_p^D (evaluation at the tuples of roots: one
    Vandermonde matrix per generator, applied along its axis of the
    compact array of coefficients), in which inversion is pointwise. For
    an integer numerator X, the adjugate A = N(X) / X (N the norm, the
    determinant of the multiplication matrix, so that A has integer
    coefficients) and N(X) are computed modulo primes of this kind and
    combined by the Chinese remainder theorem until the product of the
    primes exceeds a bound for the coefficients of X A - N (which is
    then zero, no exact check being needed). The cost per prime is
    O(D (d_1 + ... + d_m)), against O(D^3) for a linear system with the
    multiplication matrix. The primes found for a tower (with the
    Vandermonde matrices and their inverses) are kept in F->modp.
*/

/* the roots of the monic f of degree d if it splits into distinct
   linear factors */
static int
_nmod_poly_split_roots(ulong * r, const nmod_poly_t f, slong d)
{
    nmod_t mod = f->mod;
    ulong p = mod.n;
    slong i;
    int ok;

    if (d == 1)
    {
        r[0] = nmod_neg(f->coeffs[0], mod);
        return 1;
    }

    if (d == 2 && p > 2)
    {
        ulong b = f->coeffs[1], c = f->coeffs[0], disc, sq, h;
        disc = nmod_sub(nmod_mul(b, b, mod), nmod_mul(4 % p, c, mod), mod);
        if (disc == 0)
            return 0;
        sq = n_sqrtmod(disc, p);
        if (sq == 0)
            return 0;
        h = nmod_inv(2, mod);
        r[0] = nmod_mul(nmod_sub(sq, b, mod), h, mod);
        r[1] = nmod_mul(nmod_sub(nmod_neg(sq, mod), b, mod), h, mod);
        return 1;
    }

    /* x^p = x modulo f */
    {
        nmod_poly_t xp, finv;
        nmod_poly_init_mod(xp, mod);
        nmod_poly_init_mod(finv, mod);
        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series(finv, finv, f->length);
        nmod_poly_powmod_x_ui_preinv(xp, p, f, finv);
        ok = (xp->length == 2 && xp->coeffs[0] == 0 && xp->coeffs[1] == 1);
        nmod_poly_clear(xp);
        nmod_poly_clear(finv);
    }

    if (ok)
    {
        nmod_poly_factor_t fac;
        nmod_poly_factor_init(fac);
        nmod_poly_roots(fac, f, 0);
        ok = (fac->num == d);
        for (i = 0; i < d && ok; i++)
            r[i] = nmod_neg(fac->p[i].coeffs[0], mod);
        nmod_poly_factor_clear(fac);
    }

    return ok;
}

/* the primitive n-th roots of unity modulo p = 1 mod n (p prime) */
static void
_cyclotomic_roots(ulong * r, ulong n, nmod_t mod)
{
    ulong p = mod.n, a, z, t;
    n_factor_t fac;
    slong i, j;
    int prim;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    for (a = 2; ; a++)
    {
        z = nmod_pow_ui(a, (p - 1) / n, mod);
        prim = 1;
        for (i = 0; i < fac.num && prim; i++)
            if (nmod_pow_ui(z, n / fac.p[i], mod) == 1)
                prim = 0;
        if (prim)
            break;
    }

    t = 1;
    for (i = 1, j = 0; i < (slong) n; i++)
    {
        t = nmod_mul(t, z, mod);
        if (n_gcd(i, n) == 1)
            r[j++] = t;
    }
}

/* appends a prime at which all moduli split into distinct linear factors
   (with the roots); returns 0 after max_tries candidates without one */
static int
_modp_extend(modp_struct * M, const dense_layout_struct * DL, gr_tower_flat_t F, slong max_tries)
{
    slong tries, k, i;
    ulong p;
    nmod_poly_t f;
    ulong * roots;
    int ok;

    if (M->exhausted)
        return 0;

    roots = flint_malloc(sizeof(ulong) * M->R);

    for (tries = 0; tries < max_tries; tries++)
    {
        slong off = 0;

        p = M->next;
        if (p <= M->step || p < (UWORD(1) << (FLINT_BITS - 6)))
        {
            M->exhausted = 1;
            break;
        }
        M->next -= M->step;

        if (!n_is_prime(p))
            continue;

        ok = 1;
        for (k = 0; k < DL->m && ok; k++)
        {
            slong d = DL->deg[k];

            nmod_poly_init(f, p);
            fmpz_poly_get_nmod_poly(f, F->ideal_univar[k]);
            if (nmod_poly_degree(f) != d)
                ok = 0;
            else if (M->cyc[k] != 0 && p % M->cyc[k] == 1)
                _cyclotomic_roots(roots + off, M->cyc[k], f->mod);
            else if (!_nmod_poly_split_roots(roots + off, f, d))
                ok = 0;
            nmod_poly_clear(f);
            off += d;
        }

        if (ok)
        {
            if (M->num == M->alloc)
            {
                M->alloc = FLINT_MAX(8, 2 * M->alloc);
                M->primes = flint_realloc(M->primes, sizeof(ulong) * M->alloc);
                M->roots = flint_realloc(M->roots, sizeof(ulong) * M->alloc * M->R);
            }
            M->primes[M->num] = p;
            for (i = 0; i < M->R; i++)
                M->roots[M->num * M->R + i] = roots[i];
            M->num++;
            flint_free(roots);
            return 1;
        }
    }

    /* (no more searching for this tower) */
    M->exhausted = 1;
    flint_free(roots);
    return 0;
}

/* the matrices V_k (V[i][e] = r_i^e) and V_k^(-1) ([x^e] f / ((x - r_i)
   f'(r_i)) in column i) for the roots modulo p; returns 0 if two roots
   coincide */
static int
_modp_mats(ulong * mats, const ulong * roots, const dense_layout_struct * DL, gr_tower_flat_t F, nmod_t mod)
{
    slong k, i, j, e, off = 0, roff = 0, maxd = 1;
    ulong p = mod.n;
    nmod_poly_t f;
    ulong * q;
    int ok = 1;

    for (k = 0; k < DL->m; k++)
        maxd = FLINT_MAX(maxd, DL->deg[k]);
    q = flint_malloc(sizeof(ulong) * maxd);

    for (k = 0; k < DL->m && ok; k++)
    {
        slong d = DL->deg[k];
        ulong * V = mats + off, * Vi = mats + off + d * d;

        nmod_poly_init_mod(f, mod);
        fmpz_poly_get_nmod_poly(f, F->ideal_univar[k]);

        for (i = 0; i < d && ok; i++)
        {
            ulong r = roots[roff + i], rpre = n_mulmod_precomp_shoup(r, p), t = 1, fr, w, wpre;

            for (j = 0; j < d; j++)
            {
                V[i * d + j] = t;
                t = n_mulmod_shoup(r, t, rpre, p);
            }

            /* the quotient f / (x - r) by synthetic division, and its value
               at r, f'(r) */
            q[d - 1] = 1;
            for (e = d - 1; e > 0; e--)
                q[e - 1] = nmod_add(f->coeffs[e], n_mulmod_shoup(r, q[e], rpre, p), mod);
            fr = 0;
            for (e = d - 1; e >= 0; e--)
                fr = nmod_add(n_mulmod_shoup(r, fr, rpre, p), q[e], mod);
            if (fr == 0)
            {
                ok = 0;
                break;
            }
            w = nmod_inv(fr, mod);
            wpre = n_mulmod_precomp_shoup(w, p);
            for (e = 0; e < d; e++)
                Vi[e * d + i] = n_mulmod_shoup(w, q[e], wpre, p);
        }

        nmod_poly_clear(f);
        off += 2 * d * d;
        roff += d;
    }

    flint_free(q);
    return ok;
}

/* v <- the matrix A (d x d) applied along the axis of stride C of the
   compact array v of length D */
static void
_axis_apply(ulong * v, slong D, slong C, slong d, const ulong * A, nmod_t mod, ulong * t)
{
    slong blk = C * d, base, off, i, e;
    dot_params_t params = _nmod_vec_dot_params(d, mod);

    for (base = 0; base < D; base += blk)
    {
        for (off = 0; off < C; off++)
        {
            ulong * src = v + base + off;
            for (e = 0; e < d; e++)
                t[e] = src[e * C];
            for (i = 0; i < d; i++)
                src[i * C] = _nmod_vec_dot(A + i * d, t, d, mod, params);
        }
    }
}

/* v <- the adjugate N(X) / X modulo p, *N <- N(X) mod p, from v = X mod p;
   returns 0 if X is not a unit modulo p */
static int
_modp_adj(ulong * v, ulong * N, const ulong * mats, const dense_layout_struct * DL, const slong * C, slong D, nmod_t mod, ulong * t)
{
    slong k, j, off;
    ulong * pre;
    ulong acc, inv;

    /* evaluation */
    for (k = 0, off = 0; k < DL->m; k++)
    {
        slong d = DL->deg[k];
        _axis_apply(v, D, C[k], d, mats + off, mod, t);
        off += 2 * d * d;
    }

    /* the norm and the pointwise inverses (one inversion) */
    pre = flint_malloc(sizeof(ulong) * D);
    acc = 1;
    for (j = 0; j < D; j++)
    {
        if (v[j] == 0)
        {
            flint_free(pre);
            return 0;
        }
        pre[j] = acc;
        acc = nmod_mul(acc, v[j], mod);
    }
    *N = acc;
    inv = nmod_inv(acc, mod);
    for (j = D - 1; j >= 0; j--)
    {
        ulong w = nmod_mul(inv, pre[j], mod);   /* 1 / v[j] */
        inv = nmod_mul(inv, v[j], mod);
        v[j] = nmod_mul(w, acc, mod);            /* N / v[j] */
    }
    flint_free(pre);

    /* interpolation */
    for (k = 0, off = 0; k < DL->m; k++)
    {
        slong d = DL->deg[k];
        _axis_apply(v, D, C[k], d, mats + off + d * d, mod, t);
        off += 2 * d * d;
    }

    return 1;
}

/*
    A bound (in bits) for the growth of the coefficients by the folds:
    the reduction of a product array with entries at most c has entries
    at most c R_1 ... R_m, where R_k is the largest over i of
    sum_{e < 2 d_k - 1} |[x^i] (x^e mod M_k)|.
*/
static slong
_fold_growth_bits(const dense_layout_struct * DL, gr_tower_flat_t F)
{
    slong k, e, i, bits = 0;

    for (k = 0; k < DL->m; k++)
    {
        const fmpz_poly_struct * Mk = F->ideal_univar[k];
        slong d = DL->deg[k];
        fmpz * r = _fmpz_vec_init(d + 1);
        fmpz * sum = _fmpz_vec_init(d);
        fmpz_t t, R;

        fmpz_init(t);
        fmpz_init(R);
        fmpz_one(r);    /* x^0 */
        for (e = 0; e < 2 * d - 1; e++)
        {
            for (i = 0; i < d; i++)
            {
                fmpz_abs(t, r + i);
                fmpz_add(sum + i, sum + i, t);
            }
            /* r <- x r mod M_k (monic) */
            for (i = d; i > 0; i--)
                fmpz_swap(r + i, r + i - 1);
            fmpz_zero(r);
            if (!fmpz_is_zero(r + d))
            {
                for (i = 0; i < d; i++)
                    fmpz_submul(r + i, r + d, Mk->coeffs + i);
                fmpz_zero(r + d);
            }
        }
        for (i = 0; i < d; i++)
            if (fmpz_cmp(sum + i, R) > 0)
                fmpz_set(R, sum + i);
        bits += fmpz_bits(R);

        _fmpz_vec_clear(r, d + 1);
        _fmpz_vec_clear(sum, d);
        fmpz_clear(t);
        fmpz_clear(R);
    }

    return bits;
}

/* the number of candidates examined for a new prime before giving up
   (for good: a failed search marks the tower); in proportion to the
   degree D, as the alternative costs O(D^3) */
#define MODP_MAX_TRIES(D) FLINT_MIN(100000, 1000 + 100 * (D))

static int
_inv_dense_modular(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    dense_layout_struct DL;
    modp_struct * M;
    slong D = 1, k, j, i;
    slong * C;
    fmpz * X, * A;
    fmpz_t N, P, t;
    ulong * v, * tmp, * mats, Np;
    slong extra_bits = 0, abits;
    int ok, done = 0;

    _layout_init(&DL, F);
    C = flint_malloc(sizeof(slong) * DL.m);
    for (k = 0; k < DL.m; k++)
    {
        C[k] = D;
        D *= DL.deg[k];
    }

    X = _fmpz_vec_init(D);
    A = _fmpz_vec_init(D);
    fmpz_init(N);
    fmpz_init(P);
    fmpz_init(t);
    v = flint_malloc(sizeof(ulong) * D);
    tmp = flint_malloc(sizeof(ulong) * 2 * (1 + FLINT_MAX(1, D)));

    ok = _pack(X, fmpz_mpoly_q_numref(x), DL.var, DL.deg, C, DL.m, F);
    M = _modp_get(F);
    mats = flint_malloc(sizeof(ulong) * M->W);
    fmpz_one(P);

    /* X A - N is 0 modulo every prime used; it is 0 once P exceeds twice
       a bound for its coefficients: ||X||_1 ||A||_oo G + |N| */
    if (ok)
    {
        fmpz_t s1;
        fmpz_init(s1);
        for (j = 0; j < D; j++)
        {
            fmpz_abs(t, X + j);
            fmpz_add(s1, s1, t);
        }
        extra_bits = fmpz_bits(s1) + _fold_growth_bits(&DL, F) + 2;
        fmpz_clear(s1);
    }

    /* residues for batches of primes (doubling), combined by the
       Chinese remainder theorem after each batch */
    {
        ulong * R = NULL, * pr = NULL, * g;
        slong cnt = 0, target = 4;
        /* primes at which the norm of X vanishes: for a unit, at most
           as many as the bits of a Hadamard bound for the norm (a zero
           divisor, for which every prime fails, is thus detected) */
        slong fails = 0, max_fails = 16 + D * (extra_bits + FLINT_BIT_COUNT(D));

        g = flint_malloc(sizeof(ulong) * 1);
        i = 0;
        while (ok && !done)
        {
            fmpz_comb_t comb;
            fmpz_comb_temp_t ct;

            R = flint_realloc(R, sizeof(ulong) * target * (D + 1));
            pr = flint_realloc(pr, sizeof(ulong) * target);
            g = flint_realloc(g, sizeof(ulong) * target);

            while (cnt < target)
            {
                ulong p;
                nmod_t mod;

                if (i == M->num && !_modp_extend(M, &DL, F, MODP_MAX_TRIES(D)))
                {
                    ok = 0;
                    break;
                }
                p = M->primes[i];
                nmod_init(&mod, p);
                _fmpz_vec_get_nmod_vec(v, X, D, mod);
                if (_modp_mats(mats, M->roots + i * M->R, &DL, F, mod) &&
                    _modp_adj(v, &Np, mats, &DL, C, D, mod, tmp))
                {
                    for (j = 0; j < D; j++)
                        R[cnt * (D + 1) + j] = v[j];
                    R[cnt * (D + 1) + D] = Np;
                    pr[cnt] = p;
                    cnt++;
                }
                else if (++fails > max_fails)
                {
                    ok = 0;
                    break;
                }
                i++;
            }

            if (!ok)
                break;

            fmpz_comb_init(comb, pr, cnt);
            fmpz_comb_temp_init(ct, comb);
            abits = 0;
            for (j = 0; j <= D; j++)
            {
                slong u;
                for (u = 0; u < cnt; u++)
                    g[u] = R[u * (D + 1) + j];
                fmpz_multi_CRT_ui((j < D) ? A + j : N, g, comb, ct, 1);
                if (j < D)
                    abits = FLINT_MAX(abits, (slong) fmpz_bits(A + j));
            }
            fmpz_comb_temp_clear(ct);
            fmpz_comb_clear(comb);

            fmpz_one(P);
            for (j = 0; j < cnt; j++)
                fmpz_mul_ui(P, P, pr[j]);

            /* (P exceeds the bound: proved) */
            if ((slong) fmpz_bits(P) > FLINT_MAX(abits + extra_bits, (slong) fmpz_bits(N) + 1) + 1)
                done = 1;
            else
                target *= 2;
        }

        flint_free(R);
        flint_free(pr);
        flint_free(g);
    }

    if (ok && done)
    {
        /* 1 / x = dx / X = dx A / N */
        fmpz * B = _fmpz_vec_init(DL.L);
        for (j = 0; j < D; j++)
        {
            slong u = j, idx = 0;
            for (k = 0; k < DL.m; k++)
            {
                idx += (u % DL.deg[k]) * DL.S[k];
                u /= DL.deg[k];
            }
            fmpz_swap(B + idx, A + j);
        }
        fmpz_mpoly_get_fmpz(t, fmpz_mpoly_q_denref(x), F->mctx);
        if (!fmpz_is_one(t))
            _fmpz_vec_scalar_mul_fmpz(B, B, DL.L, t);
        _reduce_unpack(fmpz_mpoly_q_numref(res), B, DL.L, DL.var, DL.deg, DL.S, DL.m, F);
        if (fmpz_sgn(N) < 0)
        {
            fmpz_neg(N, N);
            fmpz_mpoly_neg(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), F->mctx);
        }
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), N, F->mctx);
        _canonicalise_int_den(res, F);
        _fmpz_vec_clear(B, DL.L);
    }

    _fmpz_vec_clear(X, D);
    _fmpz_vec_clear(A, D);
    fmpz_clear(N);
    fmpz_clear(P);
    fmpz_clear(t);
    flint_free(v);
    flint_free(tmp);
    flint_free(mats);
    flint_free(C);
    _layout_clear(&DL);
    return ok && done;
}

/* ------------------------------------------------------------------------ */
/* inverses in products of fields (several generators)                       */
/* ------------------------------------------------------------------------ */


/*
    The inverse of x in a tower of several generators (the dense case):
    the multiplication matrix of x on the monomial basis (columns: x
    times the basis monomials, as shifts of the packed array followed by
    the reductions) and one integer linear system. Returns 0 if not
    applicable or if x is a zero divisor.
*/
static int
_inv_dense_multi(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    dense_layout_struct DL;
    slong D = 1, i, j, k, L, nvars = F->mctx->minfo->nvars;
    fmpz * A, * C;
    fmpz_mat_t M, B, X;
    fmpz_t den, dx;
    slong * digit, * idx;
    ulong * e;
    int ok;

    D = gr_tower_degree(F->T);
    if (D > GR_TOWER_OPTION(F->T, GR_TOWER_OPT_INV_DENSE_DEGREE_LIMIT))
        return 0;

    _layout_init(&DL, F);
    L = DL.L;
    A = _fmpz_vec_init(L);
    C = _fmpz_vec_init(L);
    digit = flint_calloc(DL.m, sizeof(slong));
    idx = flint_malloc(sizeof(slong) * D);

    ok = _pack(A, fmpz_mpoly_q_numref(x), DL.var, DL.deg, DL.S, DL.m, F);

    /* the packed indices of the basis monomials (mixed radix, digit 0
       fastest) */
    for (j = 0; j < D; j++)
    {
        slong t = j;
        idx[j] = 0;
        for (k = 0; k < DL.m; k++)
        {
            idx[j] += (t % DL.deg[k]) * DL.S[k];
            t /= DL.deg[k];
        }
    }

    fmpz_mat_init(M, D, D);
    fmpz_mat_init(B, D, 1);
    fmpz_mat_init(X, D, 1);
    fmpz_init(den);
    fmpz_init(dx);

    for (j = 0; j < D && ok; j++)
    {
        /* x times the basis monomial j: a shift */
        _fmpz_vec_zero(C, L);
        for (i = 0; i < D; i++)
            fmpz_set(C + idx[i] + idx[j], A + idx[i]);
        _fold_array(C, L, DL.deg, DL.S, DL.m, F);
        for (i = 0; i < D; i++)
            fmpz_set(fmpz_mat_entry(M, i, j), C + idx[i]);
    }

    if (ok)
    {
        fmpz_mpoly_get_fmpz(dx, fmpz_mpoly_q_denref(x), F->mctx);
        fmpz_set(fmpz_mat_entry(B, 0, 0), dx);   /* the basis monomial 1 */
        ok = fmpz_mat_solve(X, den, M, B);
    }

    if (ok)
    {
        fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(res);
        e = flint_calloc(nvars, sizeof(ulong));
        fmpz_mpoly_zero(num, F->mctx);
        for (i = 0; i < D; i++)
        {
            slong t = i;
            if (fmpz_is_zero(fmpz_mat_entry(X, i, 0)))
                continue;
            for (k = 0; k < DL.m; k++)
            {
                e[DL.var[k]] = t % DL.deg[k];
                t /= DL.deg[k];
            }
            fmpz_mpoly_push_term_fmpz_ui(num, fmpz_mat_entry(X, i, 0), e, F->mctx);
        }
        fmpz_mpoly_sort_terms(num, F->mctx);
        if (fmpz_sgn(den) < 0)
        {
            fmpz_neg(den, den);
            fmpz_mpoly_neg(num, num, F->mctx);
        }
        fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), den, F->mctx);
        _canonicalise_int_den(res, F);
        flint_free(e);
    }

    fmpz_mat_clear(M);
    fmpz_mat_clear(B);
    fmpz_mat_clear(X);
    fmpz_clear(den);
    fmpz_clear(dx);
    _fmpz_vec_clear(A, L);
    _fmpz_vec_clear(C, L);
    flint_free(digit);
    flint_free(idx);
    _layout_clear(&DL);
    return ok;
}
