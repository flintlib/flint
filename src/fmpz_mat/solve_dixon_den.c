/*
    Copyright (C) 2011, 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "perm.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

/* Tuning parameters for the rational linear solvers. The cost of all
   operations is estimated in units of one word (limb) multiply-add;
   an a x b limb integer multiplication is modelled as
   min(mul_quad a b, mul_C (a + b)) such units. Not public API. */
typedef struct
{
    double mul_quad;      /* schoolbook multiplication constant */
    double mul_C;         /* quasi-linear multiplication constant */
    double nmod_op;       /* one modular multiply-add, in units */
    double comb_op;       /* per-prime per-entry cost of multi_mod / multi_CRT, in units */
    double switch_factor; /* lifting work required before a digit-size change,
                             as a multiple of the estimated cost of the change */
    double probe_frac;    /* max fraction of the lifting work spent on probes */
    slong blind_steps;    /* attempt without probe information for this many
                             initial steps */
} fmpz_mat_solve_tuning_struct;

static const fmpz_mat_solve_tuning_struct flint_fmpz_mat_solve_tuning = {
    0.8,    /* mul_quad */
    64.0,   /* mul_C */
    2.0,    /* nmod_op */
    2.0,    /* comb_op */
    2.0,    /* switch_factor */
    0.1,    /* probe_frac */
    8,      /* blind_steps */
};

#define FMPZ_MAT_SOLVE_MUL_COST(ka, kz) \
    FLINT_MIN(flint_fmpz_mat_solve_tuning.mul_quad * (double) (ka) * (double) (kz), \
              flint_fmpz_mat_solve_tuning.mul_C * ((double) (ka) + (double) (kz)))

/* Print compact per-step diagnostics (parameters, lifting steps, probes,
   reconstruction attempts, digit-size changes). */
#define FMPZ_MAT_SOLVE_DIXON_VERBOSE 0

#if FMPZ_MAT_SOLVE_DIXON_VERBOSE
#include <stdio.h>
#define DIXON_VERB(...) flint_printf(__VA_ARGS__)
#else
#define DIXON_VERB(...) do { } while (0)
#endif
#include "thread_support.h"


/*
    Attempt to reconstruct Z / den from x mod m, where the entries of Z
    are bounded by N in absolute value and den (a common denominator)
    is at most D. den_init is a known divisor of the denominator (e.g.
    from a probe). Entries are first tried as t = den * x mod m in
    symmetric representation; if |t| <= N we accept t directly (this
    is unique when 2 N D' < m for the remaining denominator bound D').
    Otherwise we run a bounded rational reconstruction with numerator
    bound N and denominator bound D / den, and on success rescale the
    entries processed so far.

    Returns 1 on success, 0 if some entry could not be reconstructed.
*/
int
_fmpz_mat_reconstruct_matwise(fmpz_mat_t Z, fmpz_t den,
        const fmpz_mat_t x, const fmpz_t m,
        const fmpz_t N, const fmpz_t D, const fmpz_t den_init)
{
    fmpz_t t, u, d2, Dq, hm;
    fmpz_preinvn_t minv;
    slong i, j, k, l, r, c;
    int success = 1, use_preinv;

    r = x->r;
    c = x->c;

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(d2);
    fmpz_init(Dq);
    fmpz_init(hm);

    fmpz_set(den, den_init);
    fmpz_tdiv_q_2exp(hm, m, 1);   /* m/2, for symmetric reduction */

    /* For a large fixed modulus a precomputed inverse makes each reduction
       den x mod m (a dividend of about twice the size of m) noticeably
       cheaper (measured ~1.2x at 10^5 bits, ~1.6x at 10^6 bits); below
       that the plain division is as fast or faster. */
    use_preinv = fmpz_bits(m) >= 16384;
    if (use_preinv)
        fmpz_preinvn_init(minv, m);

    for (i = 0; i < r && success; i++)
    {
        for (j = 0; j < c; j++)
        {
            if (fmpz_is_one(den))
                fmpz_set(t, fmpz_mat_entry(x, i, j));
            else
                fmpz_mul(t, den, fmpz_mat_entry(x, i, j));
            if (use_preinv)
                fmpz_fdiv_r_preinvn(t, t, m, minv);
            else
                fmpz_mod(t, t, m);
            if (fmpz_cmp(t, hm) > 0)
                fmpz_sub(t, t, m);

            if (fmpz_cmpabs(t, N) <= 0)
            {
                fmpz_set(fmpz_mat_entry(Z, i, j), t);
                continue;
            }

            /* remaining denominator bound */
            fmpz_fdiv_q(Dq, D, den);
            if (fmpz_cmp_ui(Dq, 2) < 0)
            {
                success = 0;
                break;
            }
            if (fmpz_sgn(t) < 0)
                fmpz_add(t, t, m);
            success = _fmpq_reconstruct_fmpz_2(u, d2, t, m, N, Dq);
            if (!success)
                break;
            fmpz_set(fmpz_mat_entry(Z, i, j), u);
            /* rescale previously reconstructed entries by d2 */
            for (k = 0; k < i; k++)
                for (l = 0; l < c; l++)
                    fmpz_mul(fmpz_mat_entry(Z, k, l), fmpz_mat_entry(Z, k, l), d2);
            for (l = 0; l < j; l++)
                fmpz_mul(fmpz_mat_entry(Z, i, l), fmpz_mat_entry(Z, i, l), d2);
            fmpz_mul(den, den, d2);
        }
    }

    if (use_preinv)
        fmpz_preinvn_clear(minv);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(d2);
    fmpz_clear(Dq);
    fmpz_clear(hm);
    return success;
}

/*
    R = A reduced into the symmetric residue system mod the modulus of ctx.
    fmpz_mod_mat_set_fmpz_mat uses the precomputed inverse in ctx, which is
    noticeably faster than plain division for large fixed moduli.
*/
static void
_dixon_mat_smod(fmpz_mat_t R, const fmpz_mat_t A, const fmpz_mod_ctx_t ctx, const fmpz_t halfq)
{
    slong i, j;
    fmpz_mod_mat_set_fmpz_mat(R, A, ctx);
    for (i = 0; i < R->r; i++)
        for (j = 0; j < R->c; j++)
            if (fmpz_cmp(fmpz_mat_entry(R, i, j), halfq) > 0)
                fmpz_sub(fmpz_mat_entry(R, i, j), fmpz_mat_entry(R, i, j), fmpz_mod_ctx_modulus(ctx));
}

/* Cost of verifying A Z = den B explicitly, per entry of Z, relative to
   multiplying A by a word-size matrix: the cheaper of schoolbook
   multiplication (limbs(Z) word operations per word of A) and quasi-linear
   multiplication. */
static double
_fmpz_mat_solve_verify_cost(slong bits_A, slong bits_Z)
{
    slong KA = FLINT_MAX(1, bits_A / FLINT_BITS + 1);
    slong KZ = FLINT_MAX(1, bits_Z / FLINT_BITS + 1);
    return FMPZ_MAT_SOLVE_MUL_COST(KA, KZ) / (double) KA;
}

/*
    Estimated cost, in word-operation units, of an n x n by n x cols
    matrix product with entries of ka and kb limbs: the minimum of
    entrywise integer multiplication and the multimodular algorithm
    (products modulo ~ 64 (ka + kb) / 58 word-size primes, plus the
    reduction of the inputs and the reconstruction of the output).
*/
static double
_fmpz_mat_solve_matmul_cost(slong n, slong cols, slong ka, slong kb)
{
    double np = (FLINT_BITS * (double) (ka + kb)) / (NMOD_MAT_OPTIMAL_MODULUS_BITS - 1) + 1.0;
    double direct = (double) n * FMPZ_MAT_SOLVE_MUL_COST(ka, kb);
    double mm = (double) n * np * flint_fmpz_mat_solve_tuning.nmod_op
        + 3.0 * np * flint_fmpz_mat_solve_tuning.comb_op;
    return (double) n * cols * FLINT_MIN(direct, mm);
}

/*
    Cancel any common factor of den and the entries of Z. The reconstruction
    keys the entries to the denominator it is given, which on the probe path
    is the projection denominator; if that is a proper multiple of the true
    common denominator (possible when a probe result is spurious, which
    carries a 2^-30 margin but is not impossible), the reconstruction and
    the certificate below are both still valid, but for a non-minimal
    denominator. Cancelling makes den the exact common denominator of the
    solution unconditionally, so that den always divides det(A). The
    common case exits after one or two entries with gcd 1.
*/
static void
_fmpz_mat_solve_canonicalise(fmpz_mat_t Z, fmpz_t den)
{
    fmpz_t g;
    slong i;

    if (fmpz_is_one(den))
        return;

    fmpz_init_set(g, den);

    for (i = 0; i < Z->r && !fmpz_is_one(g); i++)
        _fmpz_vec_content_chained(g, fmpz_mat_row(Z, i), Z->c, g);

    if (!fmpz_is_one(g))
    {
        fmpz_divexact(den, den, g);
        fmpz_mat_scalar_divexact_fmpz(Z, Z, g);
    }

    fmpz_clear(g);
}

/*
    Cheap test, from the probe alone, of whether a full reconstruction
    attempt can succeed at the current modulus: the numerator of the probe
    (a random combination of all entries with weights < 2^4) gives an
    estimate of the size of the numerators of the solution; the attempt is
    worthwhile if either the size certificate can hold or an explicit
    verification would be cheaper than lifting to the certificate. Called
    before the (expensive, for wide right hand sides) reconstruction of
    the full solution modulo m.
*/
int
_fmpz_mat_solve_attempt_worthwhile(const fmpz_t m, const fmpz_t snum,
        const fmpz_t sden, const fmpz_t Amax, const fmpz_t Bmax,
        slong n, slong cols, slong step_bits, double step_cost)
{
    slong bits_Z, mbits = fmpz_bits(m), cert_bits, steps_needed;

    /* numerators after clearing the probe denominator: the probe numerator
       is a sum of n cols terms with weights < 2^4 */
    bits_Z = fmpz_bits(snum) - 4 - FLINT_BIT_COUNT(n * cols);
    if (bits_Z < 0)
        bits_Z = 0;

    /* certificate: n ||A|| 2^bits_Z + den ||B|| < m */
    cert_bits = FLINT_MAX(fmpz_bits(Amax) + bits_Z, fmpz_bits(sden) + fmpz_bits(Bmax)) + 2;
    if (cert_bits <= mbits)
        return 1;

    steps_needed = (cert_bits - mbits) / step_bits + 1;
    return _fmpz_mat_solve_verify_cost(fmpz_bits(Amax), bits_Z + fmpz_bits(sden))
        < (double) steps_needed * step_cost;
}

/*
    Try to obtain a certified solution Z / den of A Z = den B from the
    solution x mod m (m being either a prime power from p-adic lifting or a
    product of primes), given a candidate common denominator sden (from a
    probe) and the global bounds N, D. Amax = n ||A||, Bmax = ||B||.

    Since A Z = den B mod m always holds for a reconstructed Z / den, the
    solution is certified when n ||A|| ||Z|| + den ||B|| < m. If that is
    not attainable at this precision we still reconstruct (with the
    uniqueness bound only) and verify explicitly when that is estimated to
    be cheaper than extending the modulus by the missing number of bits;
    step_bits is the number of bits added per step by the caller and
    step_cost the cost of a step in units of "multiply A by a word-size
    matrix" (may be a small fraction, e.g. for one modular solve);
    den_margin_bits the number of extra bits allowed in the denominator
    beyond sden.

    Returns 1 if Z / den is a certified solution, 0 otherwise.
*/
int
_fmpz_mat_solve_reconstruct_attempt(fmpz_mat_t Z, fmpz_t den,
        const fmpz_mat_t x, const fmpz_t m, const fmpz_t N, const fmpz_t D,
        const fmpz_t sden, const fmpz_t Amax, const fmpz_t Bmax,
        const fmpz_mat_t A, const fmpz_mat_t B, slong step_bits, double step_cost,
        slong den_margin_bits)
{
    fmpz_t Ncert, Dcert, t;
    int r = 0;
    slong n = A->r, cols = B->c;

    fmpz_init(Ncert);
    fmpz_init(Dcert);
    fmpz_init(t);

    /* Dcert = min(D, 2^margin sden): margin for entries whose denominator
       is a proper multiple of the probe denominator (the caller grows the
       margin after failed attempts). With no information about the
       denominator (sden = 1 and margin 0), split the modulus evenly. */
    if (fmpz_is_one(sden) && den_margin_bits == 0)
    {
        fmpz_one(Dcert);
        fmpz_mul_2exp(Dcert, Dcert, (fmpz_bits(m) - 1) / 2);
    }
    else
        fmpz_mul_2exp(Dcert, sden, den_margin_bits);
    if (fmpz_cmp(Dcert, D) > 0)
        fmpz_set(Dcert, D);

    /* uniqueness bound: 2 Ncert Dcert < m */
    fmpz_mul_2exp(t, Dcert, 1);
    fmpz_add_ui(t, t, 1);
    fmpz_fdiv_q(Ncert, m, t);
    if (fmpz_cmp(Ncert, N) > 0)
        fmpz_set(Ncert, N);
    /* certificate bound, if attainable at this precision */
    fmpz_mul(t, Dcert, Bmax);
    fmpz_sub(t, m, t);
    if (fmpz_sgn(t) > 0)
    {
        fmpz_fdiv_q(t, t, Amax);
        fmpz_sub_ui(t, t, 1);
        if (fmpz_cmp(Ncert, t) > 0)
            fmpz_set(Ncert, t);
    }

    if (fmpz_sgn(Ncert) > 0 && fmpz_cmp(Dcert, sden) >= 0 &&
        _fmpz_mat_reconstruct_matwise(Z, den, x, m, Ncert, Dcert, sden))
    {
        slong bits_Z;

        _fmpz_mat_solve_canonicalise(Z, den);
        bits_Z = FLINT_ABS(fmpz_mat_max_bits(Z));

        /* certificate with the actual sizes */
        fmpz_one(t);
        fmpz_mul_2exp(t, t, bits_Z);
        fmpz_mul(t, t, Amax);
        fmpz_addmul(t, den, Bmax);
        if (fmpz_cmp(t, m) < 0)
        {
            r = 1;
        }
        else
        {
            /* explicit check when cheaper than extending the modulus */
            slong steps_needed = (fmpz_bits(t) - fmpz_bits(m)) / step_bits + 1;
            if (_fmpz_mat_solve_verify_cost(fmpz_bits(Amax), bits_Z + fmpz_bits(den))
                < (double) steps_needed * step_cost)
            {
                fmpz_mat_t AZ, dB;
                fmpz_mat_init(AZ, n, cols);
                fmpz_mat_init(dB, n, cols);
                /* for small Z the classical product (an mpn_mul_1-type
                   product per entry) beats the multimodular algorithm,
                   which would reduce A modulo bits(A)/59 primes */
                slong KA = FLINT_MAX(1, fmpz_bits(Amax) / FLINT_BITS + 1);
                slong KZ = (bits_Z + fmpz_bits(den)) / FLINT_BITS + 1;
                if (KZ <= 4 && KA >= 8)
                    fmpz_mat_mul_classical(AZ, A, Z);
                else
                    fmpz_mat_mul(AZ, A, Z);
                fmpz_mat_scalar_mul_fmpz(dB, B, den);
                r = fmpz_mat_equal(AZ, dB);
                fmpz_mat_clear(AZ);
                fmpz_mat_clear(dB);
            }
        }
    }

    fmpz_clear(Ncert);
    fmpz_clear(Dcert);
    fmpz_clear(t);
    return r;
}

/*
    Multimodular matrix-times-(narrow matrix) product with a fixed left
    matrix whose residues are precomputed. Used for the p^k lifting where
    the left matrices (A^-1 mod q and A) have huge entries and the right
    hand sides are narrow: the cost per product is n^2 cols num_primes
    word operations plus O(n cols) reductions and CRTs, instead of
    n^2 cols products of multi-limb integers.
*/
typedef struct
{
    slong np;
    ulong * primes;
    nmod_mat_struct * M;
    slong n, cols;
    nmod_mat_struct * in_mod;
    nmod_mat_struct * out_mod;
} dixon_mm_struct;

static void
dixon_mm_init(dixon_mm_struct * mm, const fmpz_mat_t X, slong bits, slong cols)
{
    slong i, k, n = X->r;
    fmpz_t prod;

    mm->n = n;
    mm->cols = cols;
    mm->np = bits / (NMOD_MAT_OPTIMAL_MODULUS_BITS - 1) + 2;
    mm->primes = flint_malloc(sizeof(ulong) * mm->np);
    fmpz_init(prod);
    fmpz_one(prod);
    mm->primes[0] = n_nextprime(UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS, 0);
    fmpz_mul_ui(prod, prod, mm->primes[0]);
    for (i = 1; i < mm->np; i++)
    {
        mm->primes[i] = n_nextprime(mm->primes[i - 1], 0);
        fmpz_mul_ui(prod, prod, mm->primes[i]);
    }
    /* trim primes beyond what is needed */
    while (mm->np > 1)
    {
        fmpz_divexact_ui(prod, prod, mm->primes[mm->np - 1]);
        if ((slong) fmpz_bits(prod) <= bits + 1)
        {
            fmpz_mul_ui(prod, prod, mm->primes[mm->np - 1]);
            break;
        }
        mm->np--;
    }
    fmpz_clear(prod);

    mm->M = flint_malloc(sizeof(nmod_mat_struct) * mm->np);
    mm->in_mod = flint_malloc(sizeof(nmod_mat_struct) * mm->np);
    mm->out_mod = flint_malloc(sizeof(nmod_mat_struct) * mm->np);
    for (k = 0; k < mm->np; k++)
    {
        nmod_mat_init(mm->M + k, n, X->c, mm->primes[k]);
        nmod_mat_init(mm->in_mod + k, X->c, cols, mm->primes[k]);
        nmod_mat_init(mm->out_mod + k, n, cols, mm->primes[k]);
    }
    /* reduction of all entries modulo all primes at once */
    fmpz_mat_multi_mod_ui((nmod_mat_t *) mm->M, mm->np, X);
}

static void
dixon_mm_clear(dixon_mm_struct * mm)
{
    slong k;
    for (k = 0; k < mm->np; k++)
    {
        nmod_mat_clear(mm->M + k);
        nmod_mat_clear(mm->in_mod + k);
        nmod_mat_clear(mm->out_mod + k);
    }
    flint_free(mm->M);
    flint_free(mm->in_mod);
    flint_free(mm->out_mod);
    flint_free(mm->primes);
}

static void
dixon_mm_mul_worker(slong k, void * arg)
{
    dixon_mm_struct * mm = (dixon_mm_struct *) arg;
    nmod_mat_mul(mm->out_mod + k, mm->M + k, mm->in_mod + k);
}

/* out = X * in exactly (requires |X in| below half the product of primes).
   The reduction and reconstruction are threaded internally; the products
   for the different primes are independent and run in parallel here. */
static void
dixon_mm_mul(fmpz_mat_t out, dixon_mm_struct * mm, const fmpz_mat_t in)
{
    fmpz_mat_multi_mod_ui((nmod_mat_t *) mm->in_mod, mm->np, in);
    /* thread_limit 1 makes flint_parallel_do a serial loop */
    flint_parallel_do(dixon_mm_mul_worker, mm, mm->np,
        FLINT_MIN(FLINT_MAX(1, (slong) ((double) mm->np * mm->n * mm->n * mm->cols / FMPZ_MAT_SOLVE_MIN_WORK_PER_THREAD)), 1024),
        FLINT_PARALLEL_UNIFORM);
    fmpz_mat_multi_CRT_ui(out, (nmod_mat_t *) mm->out_mod, mm->np, 1);
}

/* residues of A y for the word-size phase */
typedef struct
{
    nmod_mat_t * A_mod;
    nmod_mat_t * Ay_mod;
    const nmod_mat_struct * y_mod;
} dixon_ay_args;

static void
dixon_ay_worker(slong j, void * arg)
{
    dixon_ay_args * a = (dixon_ay_args *) arg;
    nmod_mat_t y;
    /* a shallow copy with the modulus of this prime */
    *y = *a->y_mod;
    nmod_mat_set_mod(y, a->A_mod[j]->mod.n);
    nmod_mat_mul(a->Ay_mod[j], a->A_mod[j], y);
}

static int _fmpz_mat_solve_dixon_use_mm(slong n, slong bits);

/*
    Cheap probe reconstruction of the rational s mod m. First try the
    candidate denominators 1 and den_prev (the denominator found by the
    previous probe): if c s mod m in symmetric representation is at least
    2^30 times smaller than m, accept num = c s mod m, den = c (the
    probability of a spurious acceptance is ~2^-30 per probe, and any
    wrong probe result is caught by the certificate). This costs O(bits)
    and, for solutions with small denominators, terminates at about
    bits(num) + bits(den) rather than the 2 (bits(num) + bits(den))
    needed by the balanced reconstruction, which is used as fallback.
*/
int
_fmpz_mat_solve_probe_reconstruct(fmpz_t num, fmpz_t den, const fmpz_t s,
                         const fmpz_t m, const fmpz_t den_prev, int have_prev)
{
    fmpz_t t, Nb, Db;
    slong margin = 30, mbits = fmpz_bits(m);
    int c;

    fmpz_init(t);
    fmpz_init(Nb);
    fmpz_init(Db);
    for (c = 0; c < 2; c++)
    {
        if (c == 0)
            fmpz_set(t, s);
        else
        {
            if (!have_prev || fmpz_is_one(den_prev))
                break;
            fmpz_mul(t, s, den_prev);
        }
        fmpz_mod(t, t, m);
        if ((slong) fmpz_bits(t) + margin > mbits)
        {
            fmpz_sub(t, t, m);   /* symmetric representative */
            if ((slong) fmpz_bits(t) + margin > mbits)
                continue;
        }
        fmpz_swap(num, t);
        if (c == 0)
            fmpz_one(den);
        else
            fmpz_set(den, den_prev);
        fmpz_clear(t);
        fmpz_clear(Nb);
        fmpz_clear(Db);
        return 1;
    }
    /* Fallback: the true value num / den appears among the convergents of
       the continued fraction of s / m as soon as m > 2^30 num den,
       whatever the relative sizes of num and den (it is the convergent
       preceding a partial quotient of ~ m / (2 num den), and a partial
       quotient of >= 30 bits identifies it with false-positive
       probability ~2^-30 per probe; wrong probe results are caught by the
       certificate). This detects solutions with unbalanced numerator and
       denominator at about half the precision that any fixed split of the
       modulus would need. The quotients are computed with the subquadratic
       fmpq_get_cfrac, and the convergent itself is then recovered by a
       rational reconstruction with the split read off from the quotient
       sizes. */
    {
        fmpq_t xq, remq;
        fmpz * cf;
        slong nterms, k, jj, J, qual, prefix, prefJ, split;

        c = 0;
        fmpz_mod(t, s, m);
        fmpq_init(xq);
        fmpq_init(remq);
        fmpz_set(fmpq_numref(xq), t);
        fmpz_set(fmpq_denref(xq), m);
        fmpq_canonicalise(xq);
        nterms = fmpq_cfrac_bound(xq);
        cf = _fmpz_vec_init(nterms);
        k = fmpq_get_cfrac(cf, remq, xq, nterms);

        J = -1;
        qual = 30;
        prefix = 0;
        prefJ = 0;
        for (jj = 1; jj < k; jj++)
        {
            if ((slong) fmpz_bits(cf + jj) >= qual)
            {
                qual = fmpz_bits(cf + jj);
                J = jj;
                prefJ = prefix;
            }
            prefix += fmpz_bits(cf + jj);
        }

        if (J >= 0)
        {
            /* denominator of the convergent before position J has about
               prefJ bits */
            split = prefJ + FLINT_BIT_COUNT(J + 1) + 8;
            if (split + 32 <= (slong) fmpz_bits(m))
            {
                fmpz_one(Db);
                fmpz_mul_2exp(Db, Db, split);
                fmpz_one(Nb);
                fmpz_mul_2exp(Nb, Nb, fmpz_bits(m) - split - 2);
                c = _fmpq_reconstruct_fmpz_2(num, den, t, m, Nb, Db);
                /* keep only results identified with the 2^30 margin */
                if (c && fmpz_bits(num) + fmpz_bits(den) + 30 > fmpz_bits(m))
                    c = 0;
            }
        }

        _fmpz_vec_clear(cf, nterms);
        fmpq_clear(xq);
        fmpq_clear(remq);
    }
    fmpz_clear(Nb);
    fmpz_clear(Db);
    fmpz_clear(t);
    return c;
}

/* x = sum over filled levels, high level = low digits */
static void
_dixon_materialise(fmpz_mat_t x, fmpz_mat_struct * lev, const int * filled,
                   slong nlev, const fmpz * pow2, const fmpz_mat_t x0, const fmpz_t pshift)
{
    slong j;
    int first = 1;
    for (j = nlev - 1; j >= 0; j--)
    {
        if (!filled[j])
            continue;
        if (first)
        {
            fmpz_mat_set(x, lev + j);
            first = 0;
        }
        else
        {
            /* x holds digits below; lev[j] goes on top: x += lev[j] * p^(#digits in x) */
            /* the number of digits in x is the sum of 2^k over filled k > j;
               shift lev[j] by that amount: multiply by product of pow2[k] */
            slong k;
            fmpz_mat_t t;
            fmpz_mat_init_set(t, lev + j);
            for (k = j + 1; k < nlev; k++)
                if (filled[k])
                    fmpz_mat_scalar_mul_fmpz(t, t, pow2 + k);
            fmpz_mat_add(x, x, t);
            fmpz_mat_clear(t);
        }
    }
    if (first)
        fmpz_mat_zero(x);
    /* digits lifted before a change of modulus */
    if (!fmpz_is_one(pshift))
    {
        fmpz_mat_scalar_mul_fmpz(x, x, pshift);
        fmpz_mat_add(x, x, x0);
    }
}

/* Returns 0 if max_steps > 0 lifting steps were exhausted before a
   solution was found (Z, den are then undefined), else 1. */
static int
_fmpz_mat_solve_dixon_den_lu(fmpz_mat_t Z, fmpz_t den,
                    const fmpz_mat_t A, const fmpz_mat_t B,
                    const nmod_mat_t LU_in, const slong * P_in, ulong p,
                    const fmpz_t N, const fmpz_t D, slong max_steps, slong kexp)
{
    const nmod_mat_struct * LU = LU_in;
    const slong * P = P_in;
    nmod_mat_t LU2;
    slong * P2;
    int have_lu2 = 0;
    int found = 1;
    fmpz_t q;            /* lifting modulus q = p^kexp */
    fmpz_mat_t Ainv_q;   /* A^-1 mod q (only when kexp > 1) */
    dixon_mm_struct mm_inv, mm_A;
    int use_mm = 0, amod_alloc = 0, use_classical;
    slong kexp_target, kexp_first = kexp;
    fmpz_mat_t x0;
    fmpz_t pshift;
    /* Word-size steps before switching to p^k. The Newton inverse and
       residue tables cost a few n^3 multimodular products, i.e. roughly
       3n/cols word-size steps in operation count, but word-size steps
       use cheap mul_1 products so in practice the crossover is later;
       n + 2n/cols (empirical) bounds the overhead on medium-size
       solutions while still switching early enough on long lifts. */
    /* number of initial steps during which no digit-size change is
       considered (the residual has not reached its steady state yet) */
    slong warmup, bits_A;
    double newton_est;   /* estimated cost of the switch, in seconds */
    fmpz_mod_ctx_t qctx;
    fmpz_t halfq;
    int qctx_init = 0;
    fmpz_t bound, ppow, prod, Amax, Bmax, s, snum, sden,
           snum_prev, sden_prev, t, s_blk, ppow_blk, ppow_base;
    /* attempts reconstruct into Zt, swapped into Z only on success, so
       that Z may alias A or B (as fmpz_mat_inv does through
       fmpz_mat_solve) */
    fmpz_mat_t x, d, Ay, Zt;
    /* binary-counter accumulation of p-adic digit blocks: level j holds
       a block of 2^j digits (or is empty); pow2[j] = p^(2^j) */
#define DIXON_MAX_LEVELS 64
    fmpz_mat_struct lev[DIXON_MAX_LEVELS];
    fmpz pow2[DIXON_MAX_LEVELS];
    int filled[DIXON_MAX_LEVELS];
    slong nlev = 0, npow = 0;
    fmpz_mat_t cur;
    int x_valid = 0;
    ulong * crt_primes;
    ulong * u;
    slong den_margin = 16;
    nmod_mat_t * A_mod;
    nmod_mat_t d_mod, y_mod, Ainv;
    nmod_mat_t * Ay_mod;
    int have_inv = 0, use_fmpz_mul;
    fmpz_mat_t y;
    slong i, j, k, n, cols, num_primes, next_full, last_probe, probe_gap, ieff;
    /* estimated work done since the start / since the last digit-size
       change, in word-operation units (see fmpz_mat_solve_tuning_struct) */
    double work_switch, step_work;
    int probe_ok = 0, probe_stable;
    flint_rand_t state;

    n = A->r;
    cols = B->c;
    warmup = 8;
    bits_A = FLINT_ABS(fmpz_mat_max_bits(A));
    newton_est = 0.0;

    P2 = _perm_init(n);
    fmpz_init(halfq);
    fmpz_init(bound);
    fmpz_init(ppow);
    fmpz_init(q);
    fmpz_init(pshift);
    fmpz_one(pshift);
    fmpz_mat_init(x0, B->r, B->c);
    kexp = 1;
    fmpz_set_ui(q, p);
    /* with residue tables the digits may be doubled up to q ~ 2 bits(A),
       beyond which the cost per bit no longer improves; without them
       (fmpz_mat_mul products) the first digit size q ~ bits(A) is already
       in the flat region */
    /* the digit sizes were chosen in units of full-size (~2^59) digits;
       convert to units of the actual prime */
    if (kexp_first > 1)
        kexp_first = FLINT_MAX(2, (kexp_first * NMOD_MAT_OPTIMAL_MODULUS_BITS) / FLINT_BIT_COUNT(p));
    kexp_target = kexp_first;
    if (kexp_first > 1 && _fmpz_mat_solve_dixon_use_mm(n, FLINT_ABS(fmpz_mat_max_bits(A))))
        kexp_target = FLINT_MAX(kexp_first, 2 * FLINT_ABS(fmpz_mat_max_bits(A)) / ((slong) FLINT_BIT_COUNT(p) - 1));
    fmpz_init(prod);
    fmpz_init(Amax);
    fmpz_init(Bmax);
    fmpz_init(s);
    fmpz_init(snum);
    fmpz_init(sden);
    fmpz_init(snum_prev);
    fmpz_init(sden_prev);
    fmpz_init(t);
    fmpz_init(s_blk);
    fmpz_init(ppow_blk);
    fmpz_init(ppow_base);
    fmpz_one(ppow_blk);
    fmpz_one(ppow_base);

    fmpz_mat_init(x, n, cols);
    fmpz_mat_init(Zt, n, cols);
    fmpz_mat_init(cur, n, cols);
    fmpz_mat_init(Ay, n, cols);
    fmpz_mat_init(y, n, cols);
    fmpz_mat_init_set(d, B);

    /* Modulus needed for guaranteed reconstruction: m > 2 N D. */
    fmpz_mul(bound, N, D);
    fmpz_mul_2exp(bound, bound, 1);

    /* Max-norms for the size certificate. */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (fmpz_cmpabs(Amax, fmpz_mat_entry(A, i, j)) < 0)
                fmpz_abs(Amax, fmpz_mat_entry(A, i, j));
    for (i = 0; i < n; i++)
        for (j = 0; j < cols; j++)
            if (fmpz_cmpabs(Bmax, fmpz_mat_entry(B, i, j)) < 0)
                fmpz_abs(Bmax, fmpz_mat_entry(B, i, j));
    fmpz_mul_ui(Amax, Amax, n);   /* n * ||A|| */

    /* Random probe vector with small entries; s = u^T x[:,0] mod p^i is
       maintained incrementally. */
    flint_rand_init(state);
    /* random projection over all entries of the solution, so that the
       probe denominator is the common denominator (a single column may
       have a proper divisor of it, e.g. for structured inverses) */
    u = (ulong *) flint_malloc(sizeof(ulong) * n * cols);
    for (i = 0; i < n * cols; i++)
        u[i] = n_randint(state, UWORD(1) << 4) + 1;

    crt_primes = fmpz_mat_dixon_get_crt_primes(&num_primes, A, UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS);
    DIXON_VERB("dixon: n %wd cols %wd bits(A) %wd bits(B) %wd p %wu kexp_first %wd kexp_target %wd num_primes %wd bound %wd bits\n",
        n, cols, bits_A, FLINT_ABS(fmpz_mat_max_bits(B)), p, kexp_first, kexp_target, num_primes, fmpz_bits(bound));

    /* The multimodular product only pays off for large n with entries of
       moderate size (a few hundred to ~1000 bits), where nmod_mat_mul is
       efficient. For huge entries the direct product is a word-by-limb
       multiplication per entry and beats num_primes modular products plus
       CRT by a wide margin; for wide RHS with few primes fmpz_mat_mul is
       also better. */
    use_fmpz_mul = (num_primes >= 8) || (n < 32) || (cols >= 8 && num_primes <= 3);
    /* for word-size y and entries beyond a few hundred bits, the classical
       product (an mpn_mul_1 per entry) beats fmpz_mat_mul's multimodular
       algorithm, which re-reduces A on every call */
    use_classical = (num_primes >= 8);

    A_mod = (nmod_mat_t *) flint_malloc(sizeof(nmod_mat_t) * num_primes);
    if (!use_fmpz_mul)
    {
        amod_alloc = 1;
        for (j = 0; j < num_primes; j++)
        {
            nmod_mat_init(A_mod[j], n, n, crt_primes[j]);
            fmpz_mat_get_nmod_mat(A_mod[j], A);
        }
    }

    Ay_mod = (nmod_mat_t *) flint_malloc(sizeof(nmod_mat_t) * num_primes);
    for (j = 0; j < num_primes; j++)
        nmod_mat_init(Ay_mod[j], n, cols, crt_primes[j]);
    nmod_mat_init(d_mod, n, cols, p);
    nmod_mat_init(y_mod, n, cols, p);

    fmpz_one(ppow);
    i = 1; /* working with p^i */
    next_full = 1;
    last_probe = 0;
    probe_gap = 1;
    work_switch = 0.0;
    step_work = 1.0;
    probe_stable = 0;

    while (1)
    {
        /* Digit size policy. Switch from word-size digits to p^k digits,
           and later double k, only when (a) the residual d has shrunk to
           the size it keeps for the rest of the lifting (~ n ||A||); before
           that, e.g. when ||B|| >> ||A||, each word-size step is exactly
           what is needed and the residue tables would be sized for the
           huge initial d; and (b) the lifting since the last change has
           already cost twice the estimated cost of the change (Newton
           products of n x n matrices with q-size entries plus residue
           tables), which bounds the overhead to a factor 1.5 for
           solutions of any size while delaying long lifts only by a small
           fraction of their cost. Doubling k step by step rather than
           jumping to the target avoids paying for a large Newton inverse
           that a medium-size solution never uses. */
        if (kexp < kexp_target && i > warmup
            && FLINT_ABS(fmpz_mat_max_bits(d)) <= bits_A + (slong) FLINT_BIT_COUNT(n) + 2 * FLINT_BITS)
        {
            /* first jump to the digit size where p^k steps are clearly
               cheaper per bit than word-size steps (with residue tables the
               per-step overhead is set by bits(A), so small k does not pay),
               then double */
            slong kexp_new = (kexp == 1) ? kexp_first : FLINT_MIN(kexp_target, 2 * kexp);

            /* Estimated cost of the change, in word-operation units:
               the Newton iteration costs ~2 products of n x n matrices
               with entries of the new digit size (plus a geometric tail),
               and with residue tables the reduction of A and A^-1 mod q
               modulo ~ (bits(A) + bits(q)) / 59 primes each. */
            {
                slong Kq = (kexp_new * FLINT_BIT_COUNT(p)) / FLINT_BITS + 1;
                double prods = (kexp == 1) ? 3.0 : 2.0;
                newton_est = prods * _fmpz_mat_solve_matmul_cost(n, n, Kq, Kq);
                if (_fmpz_mat_solve_dixon_use_mm(n, bits_A))
                {
                    slong npt = (bits_A + kexp_new * FLINT_BIT_COUNT(p)) / (NMOD_MAT_OPTIMAL_MODULUS_BITS - 1) + 2;
                    newton_est += 2.0 * (double) n * n * npt * flint_fmpz_mat_solve_tuning.comb_op;
                }
            }
            /* A change must be amortised by the lifting it accelerates.
               When the probe has succeeded, it predicts how much lifting
               remains (up to the certificate), and the change is taken
               only if that remaining work exceeds switch_factor times the
               cost of the change. Before the probe succeeds, the solution
               is at least as large as the current modulus, so the work
               done so far is a lower bound on the work remaining and is
               used in its place. */
            {
                double gate_work;
                if (probe_ok)
                {
                    slong rem_bits = bits_A + FLINT_MAX(0, (slong) fmpz_bits(snum_prev) - 4 - (slong) FLINT_BIT_COUNT(n * cols))
                        + (slong) fmpz_bits(sden_prev) + FLINT_BIT_COUNT(n) + FLINT_BITS
                        - (slong) fmpz_bits(ppow);
                    gate_work = (double) FLINT_MAX(0, rem_bits) / fmpz_bits(q) * step_work;
                }
                else
                    gate_work = work_switch;
                if (gate_work >= flint_fmpz_mat_solve_tuning.switch_factor * newton_est)

            {
                fmpz_t qq;
                fmpz_mat_t E;
                slong m = kexp, m2;

                /* fold the digits lifted so far into x0 */
                {
                    fmpz_mat_t tmp;
                    fmpz_mat_init(tmp, n, cols);
                    _dixon_materialise(tmp, lev, filled, nlev, pow2, x0, pshift);
                    fmpz_mat_swap(x0, tmp);
                    fmpz_mat_clear(tmp);
                }
                fmpz_set(pshift, ppow);
                for (j = 0; j < nlev; j++) { fmpz_mat_clear(lev + j); fmpz_clear(pow2 + j); }
                nlev = 0; npow = 0; x_valid = 0;
                use_fmpz_mul = 1;
                probe_gap = 1;
                next_full = i;

                if (kexp == 1)
                {
                    /* A^-1 mod p from the LU */
                    nmod_mat_t Ainv1, PB;
                    nmod_mat_init(Ainv1, n, n, p);
                    nmod_mat_init(PB, n, n, p);
                    for (k = 0; k < n; k++)
                        nmod_mat_entry(PB, k, P[k]) = UWORD(1);
                    nmod_mat_solve_tril(Ainv1, LU, PB, 1);
                    nmod_mat_solve_triu(Ainv1, LU, Ainv1, 0);
                    nmod_mat_clear(PB);
                    fmpz_mat_init(Ainv_q, n, n);
                    fmpz_mat_set_nmod_mat(Ainv_q, Ainv1);   /* symmetric residues */
                    nmod_mat_clear(Ainv1);
                }
                else if (use_mm)
                {
                    dixon_mm_clear(&mm_inv);
                    dixon_mm_clear(&mm_A);
                    use_mm = 0;
                }

                /* Newton iteration X <- X + X (I - A X) from A^-1 mod p^kexp
                   to A^-1 mod p^kexp_new */
                fmpz_mat_init(E, n, n);
                fmpz_init(qq);
                while (m < kexp_new)
                {
                    fmpz_mod_ctx_t qqctx;
                    fmpz_t halfqq;
                    m2 = FLINT_MIN(2 * m, kexp_new);
                    fmpz_ui_pow_ui(qq, p, m2);
                    fmpz_mod_ctx_init(qqctx, qq);
                    fmpz_init(halfqq);
                    fmpz_tdiv_q_2exp(halfqq, qq, 1);
                    fmpz_mat_mul(E, A, Ainv_q);
                    fmpz_mat_neg(E, E);
                    for (k = 0; k < n; k++)
                        fmpz_add_ui(fmpz_mat_entry(E, k, k), fmpz_mat_entry(E, k, k), 1);
                    _dixon_mat_smod(E, E, qqctx, halfqq);
                    fmpz_mat_mul(E, Ainv_q, E);
                    fmpz_mat_add(Ainv_q, Ainv_q, E);
                    _dixon_mat_smod(Ainv_q, Ainv_q, qqctx, halfqq);
                    fmpz_mod_ctx_clear(qqctx);
                    fmpz_clear(halfqq);
                    m = m2;
                }
                fmpz_mat_clear(E);
                fmpz_clear(qq);

                kexp = kexp_new;
                fmpz_ui_pow_ui(q, p, kexp);
                if (qctx_init)
                    fmpz_mod_ctx_clear(qctx);
                fmpz_mod_ctx_init(qctx, q);
                fmpz_tdiv_q_2exp(halfq, q, 1);
                qctx_init = 1;

                /* Precompute residues of A^-1 mod q and A for multimodular
                   matrix-vector products. |d| stays below max(current |d|,
                   n ||A||) + 1, |y| < q/2. */
                {
                    slong bits_d = FLINT_MAX(FLINT_ABS(fmpz_mat_max_bits(d)), (slong) fmpz_bits(Amax)) + 2;
                    slong bits_inv = fmpz_bits(q) + bits_d + FLINT_BIT_COUNT(n) + 2;
                    slong bits_Aq = fmpz_bits(Amax) + fmpz_bits(q) + 2;
                    use_mm = _fmpz_mat_solve_dixon_use_mm(n, bits_A);
                    if (use_mm)
                    {
                        dixon_mm_init(&mm_inv, Ainv_q, bits_inv, cols);
                        dixon_mm_init(&mm_A, A, bits_Aq, cols);
                    }
                }

                work_switch = 0.0;
                DIXON_VERB("  digit size now p^%wd (%wd bits) at i %wd, use_mm %d (np %wd), est %.3g units\n",
                    kexp, fmpz_bits(q), i, use_mm, use_mm ? mm_inv.np : 0, newton_est);
            }
            }
        }

        /* Once the lifting has run for a while, forming the explicit
           inverse (2/3 n^3) pays for itself against the extra cost of
           triangular solves (~0.4 n^2 cols per step). */
        if (kexp > 1)
        {
            /* y = A^-1 d mod q, in symmetric representation */
            if (use_mm)
                dixon_mm_mul(y, &mm_inv, d);
            else
                fmpz_mat_mul(y, Ainv_q, d);
            _dixon_mat_smod(y, y, qctx, halfq);
        }
        else if (!have_inv && i * cols >= n / 2)
        {
            nmod_mat_t PB;
            /* The lifting started with a small prime (cheap LU, and the
               triangular solves cost the same per bit as with a full-size
               prime). For a long lift with the explicit inverse a full-size
               prime is ~20% cheaper per bit, so re-factor A modulo a
               ~2^59 prime here (cost ~ 1/3 of the inverse we form anyway)
               and continue the lifting with the new modulus. */
            if (FLINT_BIT_COUNT(p) < NMOD_MAT_OPTIMAL_MODULUS_BITS && !flint_fmpz_mat_force_small_primes)
            {
                ulong p2;
                nmod_mat_init(LU2, n, n, 2);
                /* A is nonsingular here (we are lifting a solution), so no
                   singularity certification is wanted */
                p2 = _fmpz_mat_find_good_prime_and_lu2(LU2, P2, A, NULL,
                        NMOD_MAT_OPTIMAL_MODULUS_BITS, 1);
                if (p2 != 0)
                {
                    fmpz_mat_t tmp;
                    fmpz_mat_init(tmp, n, cols);
                    _dixon_materialise(tmp, lev, filled, nlev, pow2, x0, pshift);
                    fmpz_mat_swap(x0, tmp);
                    fmpz_mat_clear(tmp);
                    fmpz_set(pshift, ppow);
                    for (j = 0; j < nlev; j++) { fmpz_mat_clear(lev + j); fmpz_clear(pow2 + j); }
                    nlev = 0; npow = 0; x_valid = 0;
                    p = p2;
                    fmpz_set_ui(q, p);
                    LU = LU2;
                    P = P2;
                    nmod_mat_set_mod(d_mod, p);
                    nmod_mat_set_mod(y_mod, p);
                    have_lu2 = 1;
                }
                else
                    nmod_mat_clear(LU2);
            }
            nmod_mat_init(Ainv, n, n, p);
            nmod_mat_init(PB, n, n, p);
            for (k = 0; k < n; k++)
                nmod_mat_entry(PB, k, P[k]) = UWORD(1);
            nmod_mat_solve_tril(Ainv, LU, PB, 1);
            nmod_mat_solve_triu(Ainv, LU, Ainv, 0);
            nmod_mat_clear(PB);
            have_inv = 1;
        }

        /* y = A^(-1) * d  (mod p) */
        if (kexp > 1)
        {
        }
        else if (have_inv)
        {
            fmpz_mat_get_nmod_mat(d_mod, d);
            nmod_mat_mul(y_mod, Ainv, d_mod);
        }
        else
        {
            /* via LU: y = U^-1 L^-1 (P d) */
            for (k = 0; k < n; k++)
                for (j = 0; j < cols; j++)
                    nmod_mat_entry(d_mod, k, j) = fmpz_fdiv_ui(fmpz_mat_entry(d, P[k], j), p);
            nmod_mat_solve_tril(y_mod, LU, d_mod, 1);
            nmod_mat_solve_triu(y_mod, LU, y_mod, 0);
        }
        if (kexp == 1)
            fmpz_mat_set_nmod_mat_unsigned(y, y_mod);

        /* Accumulate digit y into the binary-counter structure. Doing
           x += y * p^i directly would cost O(i) per entry per step, i.e.
           quadratic in the number of steps overall. */
        fmpz_mat_set(cur, y);
        x_valid = 0;
        for (j = 0; ; j++)
        {
            if (j >= nlev)
            {
                fmpz_mat_init(lev + j, n, cols);
                filled[j] = 0;
                nlev++;
                fmpz_init(pow2 + j);
                if (j == 0)
                    fmpz_set(pow2 + j, q);
                else
                    fmpz_mul(pow2 + j, pow2 + j - 1, pow2 + j - 1);
                npow++;
            }
            if (!filled[j])
            {
                fmpz_mat_swap(lev + j, cur);
                filled[j] = 1;
                break;
            }
            /* merge: level j holds the lower 2^j digits */
            fmpz_mat_scalar_mul_fmpz(cur, cur, pow2 + j);
            fmpz_mat_add(cur, cur, lev + j);
            filled[j] = 0;
        }

        /* probe: s += p^i * (u . y), over all entries of y; u < 2^4 and
           the digits y are < 2^64 */
        fmpz_zero(t);
        for (k = 0; k < n; k++)
            for (j = 0; j < cols; j++)
                fmpz_addmul_ui(t, fmpz_mat_entry(y, k, j), u[k * cols + j]);
        /* s += t * p^i, accumulated blockwise (s_blk * p^i_blk) to avoid
           multiplying by the full-length p^i every step */
        fmpz_addmul(s_blk, t, ppow_blk);
        fmpz_mul(ppow_blk, ppow_blk, q);

        /* Estimated work of the step just performed, in word-operation
           units, for the probe and digit-size schedules. */
        {
            slong KA = bits_A / FLINT_BITS + 1;
            slong Kq = fmpz_bits(q) / FLINT_BITS + 1;
            double tun_nmod = flint_fmpz_mat_solve_tuning.nmod_op;
            if (kexp == 1)
            {
                /* triangular solves or gemv, plus A y */
                step_work = 2.0 * (double) n * n * cols * tun_nmod;
                if (use_mm || !use_fmpz_mul)
                    step_work += (double) n * n * cols * num_primes * tun_nmod;
                else
                    step_work += (double) n * n * cols * FLINT_MIN((double) KA, FMPZ_MAT_SOLVE_MUL_COST(KA, 1));
            }
            else if (use_mm)
            {
                slong np = mm_inv.np;
                step_work = 2.0 * (double) n * n * cols * np * tun_nmod
                    + 4.0 * (double) n * cols * np * flint_fmpz_mat_solve_tuning.comb_op;
            }
            else
                step_work = 2.0 * _fmpz_mat_solve_matmul_cost(n, cols, Kq, Kq);
            work_switch += step_work;
        }

        /* ppow = p^(i+1) */
        fmpz_mul(ppow, ppow, q);
        /* number of steps of the current digit size that would give the
           current precision; schedules below are expressed in these units
           so that they stay meaningful after a change of digit size */
        ieff = fmpz_bits(ppow) / fmpz_bits(q);

        if (fmpz_cmp(ppow, bound) > 0)
        {
            /* guaranteed: reconstruct with the global bounds */
            DIXON_VERB("  full bound reached at i %wd (%wd bits)\n", i, fmpz_bits(ppow));
            fmpz_one(t);
            if (!x_valid) { _dixon_materialise(x, lev, filled, nlev, pow2, x0, pshift); x_valid = 1; }
            /* The reconstruction is unique for m > 2 N D and the true
               solution satisfies the bounds, so this cannot fail; report
               failure rather than a wrong answer if it ever does. */
            if (!_fmpz_mat_reconstruct_matwise(Zt, den, x, ppow, N, D, t))
            {
                found = 0;
                goto dixon_done;
            }
            fmpz_mat_swap(Z, Zt);
            goto dixon_done;
        }

        /* Cheap probe reconstruction (balanced bounds). The hgcd costs
           O(M(i log p)), so probe at most ~16 times per doubling of i. */
        if (i - last_probe >= probe_gap)
        {
            last_probe = i;
            fmpz_addmul(s, s_blk, ppow_base);
            fmpz_zero(s_blk);
            fmpz_one(ppow_blk);
            fmpz_set(ppow_base, ppow);
            probe_ok = _fmpz_mat_solve_probe_reconstruct(snum, sden, s, ppow, sden_prev, probe_ok);
            /* Keep the total probing cost to a small fraction of the
               lifting: a probe costs a few multiplications at the size of
               the modulus. Probing every i/16 steps also bounds the number
               of wasted steps. */
            {
                slong Km = i * FLINT_BIT_COUNT(p) / FLINT_BITS + 1;
                /* a probe is dominated by the continued fraction of s/m,
                   which costs about 5 log2(limbs) full-size multiplications
                   (measured) */
                double probe_cost = 5.0 * FLINT_BIT_COUNT(Km) * FMPZ_MAT_SOLVE_MUL_COST(Km, Km);
                slong gap2 = (slong) (probe_cost / (flint_fmpz_mat_solve_tuning.probe_frac * step_work + 1.0));
                probe_gap = FLINT_MAX(ieff / 16, 1);
                probe_gap = FLINT_MAX(probe_gap, FLINT_MIN(gap2, ieff + 1));
            }
            DIXON_VERB("  probe i %wd (%wd bits): ok %d num %wd bits den %wd bits\n",
                i, fmpz_bits(ppow), probe_ok, probe_ok ? fmpz_bits(snum) : -1, probe_ok ? fmpz_bits(sden) : -1);
            /* every probe acceptance path (candidate denominator or
               continued-fraction convergent) carries a 2^-30 margin, so a
               single successful probe is confident enough to attempt on;
               agreement with the previous probe only resets the back-off */
            if (probe_ok)
                probe_stable = 1;
            else
            {
                probe_stable = 0;
                next_full = i;   /* reset backoff when the probe changes */
            }
            fmpz_swap(snum, snum_prev);
            fmpz_swap(sden, sden_prev);
        }

        /* An attempt costs O(n cols i) fmpz operations, a step costs
           O(n^2 cols) or more, so while i is small an attempt is cheap
           relative to a step: try on any successful probe rather than
           waiting for two agreeing probes. */
        /* While the modulus is small an attempt is cheap: try even without
           a successful probe (the probe, a weighted sum, needs a larger
           modulus than the entries), reconstructing with balanced bounds
           and relying on the certificate / verification. */
        if (((probe_ok && (probe_stable || ieff <= 8)) || ieff <= 8) && i >= next_full)
        {
            /* at small moduli a spurious probe result is likely (about half
               of all residues have a small-fraction reconstruction), so
               ignore the probe there as well */
            int blind = !probe_ok || ieff <= 8;
            if (!blind && !_fmpz_mat_solve_attempt_worthwhile(ppow, snum_prev, sden_prev, Amax, Bmax, n, cols, fmpz_bits(q) - 1, 1.0))
                goto dixon_no_attempt;
            if (!x_valid) { _dixon_materialise(x, lev, filled, nlev, pow2, x0, pshift); x_valid = 1; }
            if (blind)
                fmpz_one(t);
            DIXON_VERB("  attempt i %wd (%wd bits): blind %d sden %wd bits margin %wd\n",
                i, fmpz_bits(ppow), blind, blind ? 0 : fmpz_bits(sden_prev), blind ? 0 : den_margin);
            if (_fmpz_mat_solve_reconstruct_attempt(Zt, den, x, ppow, N, D,
                    blind ? t : sden_prev, Amax, Bmax, A, B, fmpz_bits(q) - 1, 1, blind ? 0 : den_margin))
            {
                fmpz_mat_swap(Z, Zt);
                DIXON_VERB("  success: den %wd bits, X %wd bits\n", fmpz_bits(den), FLINT_ABS(fmpz_mat_max_bits(Z)));
                goto dixon_done;
            }
            DIXON_VERB("  attempt failed\n");
            /* failed with a stable probe: allow a larger denominator next time */
            if (probe_stable && den_margin < (slong) fmpz_bits(D))
                den_margin *= 2;
            dixon_no_attempt: ;
            /* failed: back off */
            next_full = i + FLINT_MAX(1, ieff / 4);
        }

        i++;
        if (max_steps > 0 && i * kexp > max_steps)
        {
            found = 0;
            goto dixon_done;
        }

        /* d = (d - Ay) / p */
        if (use_mm)
            dixon_mm_mul(Ay, &mm_A, y);
        else if (kexp == 1 && use_classical)
            fmpz_mat_mul_classical(Ay, A, y);
        else if (use_fmpz_mul || kexp > 1)
            fmpz_mat_mul(Ay, A, y);
        else
        {
            /* residues of A y for each prime, then one multimodular
               reconstruction (one multimodular
               reconstruction over all primes) */
            {
                dixon_ay_args ay_args;
                ay_args.A_mod = A_mod;
                ay_args.Ay_mod = Ay_mod;
                ay_args.y_mod = y_mod;
                flint_parallel_do(dixon_ay_worker, &ay_args, num_primes,
                    FLINT_MIN(FLINT_MAX(1, (slong) ((double) num_primes * n * n * cols / FMPZ_MAT_SOLVE_MIN_WORK_PER_THREAD)), 1024),
                    FLINT_PARALLEL_UNIFORM);
            }
            fmpz_mat_multi_CRT_ui(Ay, Ay_mod, num_primes, 1);
        }

        fmpz_mat_sub(d, d, Ay);
        if (kexp == 1)
            fmpz_mat_scalar_divexact_ui(d, d, p);
        else
            fmpz_mat_scalar_divexact_fmpz(d, d, q);
    }

dixon_done:

    flint_rand_clear(state);
    flint_free(u);

    if (have_inv)
        nmod_mat_clear(Ainv);
    if (have_lu2)
        nmod_mat_clear(LU2);
    _perm_clear(P2);
    nmod_mat_clear(y_mod);
    nmod_mat_clear(d_mod);
    for (j = 0; j < num_primes; j++)
        nmod_mat_clear(Ay_mod[j]);
    flint_free(Ay_mod);

    if (amod_alloc)
        for (j = 0; j < num_primes; j++)
            nmod_mat_clear(A_mod[j]);
    flint_free(A_mod);
    flint_free(crt_primes);

    if (qctx_init)
        fmpz_mod_ctx_clear(qctx);
    fmpz_clear(halfq);
    fmpz_clear(bound);
    fmpz_clear(ppow);
    fmpz_clear(q);
    fmpz_clear(pshift);
    fmpz_mat_clear(x0);
    if (kexp > 1)
        fmpz_mat_clear(Ainv_q);
    if (use_mm)
    {
        dixon_mm_clear(&mm_inv);
        dixon_mm_clear(&mm_A);
    }
    fmpz_clear(prod);
    fmpz_clear(Amax);
    fmpz_clear(Bmax);
    fmpz_clear(s);
    fmpz_clear(snum);
    fmpz_clear(sden);
    fmpz_clear(snum_prev);
    fmpz_clear(sden_prev);
    fmpz_clear(t);
    fmpz_clear(s_blk);
    fmpz_clear(ppow_blk);
    fmpz_clear(ppow_base);

    for (j = 0; j < nlev; j++)
        fmpz_mat_clear(lev + j);
    for (j = 0; j < npow; j++)
        fmpz_clear(pow2 + j);
    fmpz_mat_clear(cur);
    fmpz_mat_clear(d);
    fmpz_mat_clear(x);
    fmpz_mat_clear(Zt);
    fmpz_mat_clear(Ay);
    fmpz_mat_clear(y);
    return found;
}

#undef DIXON_MAX_LEVELS

/* Multimodular matvecs with precomputed residues cost O(n^2 num_primes)
   per product plus O(n M(bits) log) for reducing and reconstructing the
   vector entries; the latter dominates for small n and huge entries. */
static int
_fmpz_mat_solve_dixon_use_mm(slong n, slong bits)
{
    return (n >= 24 && bits <= 20000) || n >= 100;
}

/* Number of word-size digits per lifting step. */
static slong
_fmpz_mat_solve_dixon_choose_k(const fmpz_mat_t A, const fmpz_mat_t B)
{
    (void) B;
    slong bits = FLINT_ABS(fmpz_mat_max_bits(A));
    /* Empirical: with multimodular matvecs (large n) the optimum is
       q ~ bits(A)/2 for moderate sizes, decreasing towards bits(A)/8 for
       very large entries since Newton (~3 full-precision matrix products)
       grows faster than the lifting; without them (small n, where
       fmpz_mat_mul is used) larger digits pay off once the integer
       multiplication is subquadratic, q ~ bits(A). */
    if (bits < 150)
        return 1;
    if (_fmpz_mat_solve_dixon_use_mm(A->r, bits))
        return FLINT_MAX(1, bits / (50 + bits / 150));
    return FLINT_MAX(1, bits / (NMOD_MAT_OPTIMAL_MODULUS_BITS - 1));
}

/* max_steps > 0 limits the number of lifting steps; returns -1 if the
   limit was reached without finding the solution. */
static int
_fmpz_mat_solve_dixon_den_limited(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B, slong max_steps)
{
    int found = 1;
    nmod_mat_t LU;
    slong * P;
    fmpz_t N, D;
    ulong p;

    if (!fmpz_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_solve_dixon_den). Non-square system matrix.\n");
    }

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
    {
        fmpz_one(den);
        return 1;
    }

    fmpz_init(N);
    fmpz_init(D);

    fmpz_mat_solve_bound(N, D, A, B);

    nmod_mat_init(LU, A->r, A->r, 2);
    P = _perm_init(A->r);
    /* start with a ~2^30 prime: the LU and the triangular solves are about
       2x cheaper than with a ~2^59 prime, at the same cost per lifted bit;
       a full-size prime is substituted later if the lifting turns out to
       be long (see the explicit-inverse step in the core) */
    p = _fmpz_mat_find_good_prime_and_lu2(LU, P, A, D, 30, 0);

    if (p != 0)
    {
        slong kexp = _fmpz_mat_solve_dixon_choose_k(A, B);
        found = _fmpz_mat_solve_dixon_den_lu(X, den, A, B, LU, P, p, N, D, max_steps, kexp);
    }

    nmod_mat_clear(LU);
    _perm_clear(P);
    fmpz_clear(N);
    fmpz_clear(D);

    if (p == 0)
        return 0;
    return found ? 1 : -1;
}

/*
    Solve A X = den B for nonsingular A, given a prime p for which A is
    nonsingular modulo p and the corresponding LU factorisation (as
    produced by nmod_mat_lu). Since no prime search is involved, this
    entry point neither needs nor performs any singularity certification.
    Returns 1 on success and 0 if the lifting was exhausted.
*/
int
_fmpz_mat_solve_dixon_den_given_lu(fmpz_mat_t X, fmpz_t den,
                    const fmpz_mat_t A, const fmpz_mat_t B,
                    const nmod_mat_t LU, const slong * P, ulong p)
{
    fmpz_t N, D;
    int found;

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
    {
        fmpz_one(den);
        return 1;
    }

    fmpz_init(N);
    fmpz_init(D);
    fmpz_mat_solve_bound(N, D, A, B);

    found = _fmpz_mat_solve_dixon_den_lu(X, den, A, B, LU, P, p, N, D,
                WORD_MAX, _fmpz_mat_solve_dixon_choose_k(A, B));

    fmpz_clear(N);
    fmpz_clear(D);

    return found;
}

int
fmpz_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_mat_solve_dixon_den_limited(X, den, A, B, 0);
}

