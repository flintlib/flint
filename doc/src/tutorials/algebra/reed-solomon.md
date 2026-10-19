---
title: Decoding Reed–Solomon codes
summary: A complete Reed–Solomon encoder and decoder over a prime field, from the generator polynomial to syndromes, the Berlekamp–Massey algorithm, root finding and Forney's formula, scaled up to codes of length 32768 with fast evaluation.
order: 3
---

Reed–Solomon codes [@RS1960] protect data against corruption by adding
redundant symbols computed from polynomial arithmetic over a finite field.
They are used in CDs, DVDs and QR codes over $\mathbb{F}_{2^8}$, and the
PDF417 barcode format uses them over the prime field $\mathbb{F}_{929}$,
which makes it a natural setting for [](nmod_poly). This tutorial writes a
complete encoder and decoder. A codeword of length $N$ carries $K$ message
symbols and $N - K = 2T$ check symbols, and the decoder corrects any $T$
or fewer corrupted symbols, without knowing where they are. Along the way
we use polynomial division, multipoint evaluation, the
[Berlekamp–Massey algorithm](nmod_poly#berlekamp-massey), root finding and
power series arithmetic.

Throughout, $p$ is a prime and $\alpha$ an element of order at least $N$ in
$\mathbb{F}_p^\times$; we start with $p = 929$, $\alpha = 3$ (a
primitive root), $N = 40$, $K = 24$ and $T = 8$.

## Step 1: encoding {#encoding}

A word $c = (c_0, \ldots, c_{N-1})$ is identified with the polynomial
$c(x) = \sum c_i x^i$. The code consists of the polynomials of degree less
than $N$ that vanish at $\alpha, \alpha^2, \ldots, \alpha^{2T}$, that is,
the multiples of the *generator polynomial*
$$
g(x) = (x - \alpha)(x - \alpha^2) \cdots (x - \alpha^{2T}).
$$
Two different codewords differ in at least $2T + 1$ positions: their
difference is a nonzero multiple of $g$ of degree less than $N$, and a
polynomial with at most $2T$ nonzero coefficients cannot vanish at $2T$
distinct powers of $\alpha$ (the corresponding Vandermonde matrix is
invertible). So $T$ errors can be corrected.

In *systematic* encoding the message $m(x)$ of length $K$ appears
unchanged in the high coefficients: $c = m x^{2T} - (m x^{2T} \bmod g)$,
which is divisible by $g$. FLINT gives us $g$ with
`nmod_poly_product_roots_nmod_vec` and the remainder with `nmod_poly_rem`:

```c example name=tutorial_rs_encode
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

#define P 929       /* the field GF(929) of PDF417 barcodes */
#define ALPHA 3     /* a primitive root modulo 929 */
#define N 40        /* codeword length */
#define K 24        /* message length */
#define T 8         /* errors that can be corrected: N - K = 2T */

int main(void)
{
    nmod_t mod;
    nmod_poly_t g, m, c;
    ulong pts[2 * T];
    slong i;

    nmod_init(&mod, P);
    nmod_poly_init_mod(g, mod);
    nmod_poly_init_mod(m, mod);
    nmod_poly_init_mod(c, mod);

    /* generator polynomial g = (x - a)(x - a^2)...(x - a^2T) */
    pts[0] = ALPHA;
    for (i = 1; i < 2 * T; i++)
        pts[i] = nmod_mul(pts[i - 1], ALPHA, mod);
    nmod_poly_product_roots_nmod_vec(g, pts, 2 * T);

    /* message: the symbols 1, 2, ..., K */
    for (i = 0; i < K; i++)
        nmod_poly_set_coeff_ui(m, i, i + 1);

    /* systematic encoding: c = m x^2T - (m x^2T mod g) */
    nmod_poly_shift_left(c, m, 2 * T);
    nmod_poly_rem(g, c, g);
    nmod_poly_sub(c, c, g);

    flint_printf("check symbols:");
    for (i = 0; i < 2 * T; i++)
        flint_printf(" %wu", nmod_poly_get_coeff_ui(c, i));
    flint_printf("\n");

    {
        ulong s[2 * T];
        nmod_poly_evaluate_nmod_vec(s, c, pts, 2 * T);
        flint_printf("c(a^i) = 0 for i = 1..2T: %d\n", _nmod_vec_is_zero(s, 2 * T));
    }

    nmod_poly_clear(g);
    nmod_poly_clear(m);
    nmod_poly_clear(c);
    return 0;
}
```
```output
check symbols: 382 684 808 907 581 514 884 383 661 205 350 74 800 634 107 338
c(a^i) = 0 for i = 1..2T: 1
```

The check that $c$ vanishes at $\alpha, \ldots, \alpha^{2T}$ uses
`nmod_poly_evaluate_nmod_vec`, which for short polynomials evaluates point
by point and for long ones uses a subproduct tree
([](nmod_poly/multipoint)).

## Step 2: syndromes and the error locator {#locator}

Suppose the received word is $r = c + e$, where the error $e$ is nonzero in
$\nu \le T$ positions $e_1, \ldots, e_\nu$ with values $Y_1, \ldots,
Y_\nu$. Write $X_k = \alpha^{e_k}$ for the *error locators*. Since $c$
vanishes at the powers of $\alpha$, the *syndromes*
$$
s_i = r(\alpha^i) = e(\alpha^i) = \sum_{k=1}^{\nu} Y_k X_k^i, \qquad 1 \le
i \le 2T,
$$
depend only on the error. All of them are zero if and only if no error is
detected.

The sequence $s_i$ is a sum of $\nu$ geometric sequences, so it satisfies
the linear recurrence whose characteristic polynomial is
$$
V(x) = (x - X_1)(x - X_2) \cdots (x - X_\nu),
$$
and no shorter one. Finding the shortest recurrence of a sequence is
exactly what the Berlekamp–Massey algorithm [@Ber1968] [@Mas1969] does,
and $2\nu \le 2T$ terms suffice to determine a recurrence of order $\nu$.
FLINT's `nmod_berlekamp_massey_t` implements it as a partial extended
Euclidean algorithm [@SKHN1975], switching to the
[half-GCD algorithm](algorithms/polynomials/half-gcd) for long sequences,
and returns $V$ (up to a constant factor) with
`nmod_berlekamp_massey_V_poly`. Its documentation explains the second
polynomial $R$: if $\deg R \ge \deg V$, the recurrence does not explain
the whole sequence, which means that more than $T$ errors occurred.

The error locators are the roots of $V$, which are distinct and nonzero.
`nmod_poly_find_distinct_nonzero_roots` finds them by splitting $V$ with
GCDs with $(x + \delta)^{(p-1)/2} - 1$ [@Rab1980b], and fails if $V$ does
not have $\deg V$ distinct roots in $\mathbb{F}_p$ — another sign of too
many errors. Classical decoders search all $N$ powers of $\alpha$ instead
(Chien search [@Chi1964]), which is fine for short codes. The position
$e_k$ is the discrete logarithm of $X_k$ to the base $\alpha$, which we
look up in a table.

## Step 3: the error values {#forney}

With $\Lambda(x) = \prod_k (1 - X_k x)$, the reversal of $V$, and the
syndrome series $S(x) = \sum_{i=0}^{2T-1} s_{i+1} x^i$, we have
$$
S(x) \equiv \sum_k \frac{Y_k X_k}{1 - X_k x} \pmod{x^{2T}},
\qquad
\Omega(x) = S(x) \Lambda(x) \bmod x^{2T} = \sum_k Y_k X_k \prod_{j \ne k}
(1 - X_j x).
$$
Evaluating $\Omega$ and $\Lambda'(x) = -\sum_k X_k \prod_{j \ne k} (1 -
X_j x)$ at $X_k^{-1}$ gives Forney's formula [@For1965]
$$
Y_k = -\frac{\Omega(X_k^{-1})}{\Lambda'(X_k^{-1})}.
$$
In FLINT this is a truncated product (`nmod_poly_mullow`), a derivative
and two evaluations per error. The complete decoder corrects 3 and 6
errors, and detects that 9 errors are too many:

```c example name=tutorial_rs_decode
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

#define P 929
#define ALPHA 3
#define N 40
#define K 24
#define T 8

/* Decode the received word r (length N) in place. Returns the number of
   corrected errors, or -1 if more than T errors were detected. */
slong rs_decode(nmod_poly_t r, nmod_t mod, int verbose)
{
    ulong pts[2 * T], s[2 * T], X[T], inv, num, den, logtab[P];
    nmod_berlekamp_massey_t B;
    nmod_poly_t V, Lambda, S, Omega, dLambda;
    slong i, nu, result;

    /* syndromes s_i = r(a^i), i = 1..2T */
    pts[0] = ALPHA;
    for (i = 1; i < 2 * T; i++)
        pts[i] = nmod_mul(pts[i - 1], ALPHA, mod);
    nmod_poly_evaluate_nmod_vec(s, r, pts, 2 * T);
    if (_nmod_vec_is_zero(s, 2 * T))
        return 0;

    /* error locator: minimal polynomial of the syndrome sequence */
    nmod_berlekamp_massey_init(B, mod.n);
    nmod_berlekamp_massey_add_points(B, s, 2 * T);
    nmod_berlekamp_massey_reduce(B);

    nmod_poly_init_mod(V, mod);
    nmod_poly_init_mod(Lambda, mod);
    nmod_poly_init_mod(S, mod);
    nmod_poly_init_mod(Omega, mod);
    nmod_poly_init_mod(dLambda, mod);
    nmod_poly_make_monic(V, nmod_berlekamp_massey_V_poly(B));
    nu = nmod_poly_degree(V);
    result = -1;

    if (verbose)
    {
        flint_printf("  V = "); nmod_poly_print_pretty(V, "x"); flint_printf("\n");
    }

    /* the roots of V are the error locators X_k = a^(e_k) */
    if (nu > T || nmod_poly_degree(nmod_berlekamp_massey_R_poly(B)) >= nu
               || !nmod_poly_find_distinct_nonzero_roots(X, V))
        goto cleanup;

    /* Forney: Lambda = rev(V), Omega = S Lambda mod x^2T,
       Y_k = -Omega(1/X_k) / Lambda'(1/X_k) */
    nmod_poly_reverse(Lambda, V, nu + 1);
    for (i = 0; i < 2 * T; i++)
        nmod_poly_set_coeff_ui(S, i, s[i]);
    nmod_poly_mullow(Omega, S, Lambda, 2 * T);
    nmod_poly_derivative(dLambda, Lambda);

    for (i = 0; i < P; i++)
        logtab[i] = N;                      /* positions >= N are invalid */
    for (i = 0, num = 1; i < N; i++, num = nmod_mul(num, ALPHA, mod))
        logtab[num] = i;

    for (i = 0; i < nu; i++)
    {
        slong e = logtab[X[i]];
        if (e >= N)
            goto cleanup;                   /* locator outside the word */
        inv = nmod_inv(X[i], mod);
        num = nmod_poly_evaluate_nmod(Omega, inv);
        den = nmod_poly_evaluate_nmod(dLambda, inv);
        num = nmod_neg(nmod_div(num, den, mod), mod);
        if (verbose)
            flint_printf("  error %wu at position %wd\n", num, e);
        nmod_poly_set_coeff_ui(r, e,
            nmod_sub(nmod_poly_get_coeff_ui(r, e), num, mod));
    }
    result = nu;

cleanup:
    nmod_berlekamp_massey_clear(B);
    nmod_poly_clear(V);
    nmod_poly_clear(Lambda);
    nmod_poly_clear(S);
    nmod_poly_clear(Omega);
    nmod_poly_clear(dLambda);
    return result;
}

int main(void)
{
    nmod_t mod;
    nmod_poly_t g, c, r;
    ulong pts[2 * T];
    slong i, nerr;
    slong pos[9] = { 3, 7, 11, 17, 22, 25, 30, 38, 39 };
    ulong val[9] = { 100, 1, 928, 5, 77, 300, 2, 40, 9 };

    nmod_init(&mod, P);
    nmod_poly_init_mod(g, mod);
    nmod_poly_init_mod(c, mod);
    nmod_poly_init_mod(r, mod);

    pts[0] = ALPHA;
    for (i = 1; i < 2 * T; i++)
        pts[i] = nmod_mul(pts[i - 1], ALPHA, mod);
    nmod_poly_product_roots_nmod_vec(g, pts, 2 * T);

    for (i = 0; i < K; i++)
        nmod_poly_set_coeff_ui(c, i, i + 1);
    nmod_poly_shift_left(c, c, 2 * T);
    nmod_poly_rem(r, c, g);
    nmod_poly_sub(c, c, r);

    for (nerr = 3; nerr <= 9; nerr += 3)
    {
        nmod_poly_set(r, c);
        for (i = 0; i < nerr; i++)
            nmod_poly_set_coeff_ui(r, pos[i],
                nmod_add(nmod_poly_get_coeff_ui(r, pos[i]), val[i], mod));
        flint_printf("%wd errors:\n", nerr);
        i = rs_decode(r, mod, nerr == 3);
        flint_printf("  decoder returns %wd, codeword restored: %d\n", i, nmod_poly_equal(r, c));
    }

    nmod_poly_clear(g);
    nmod_poly_clear(c);
    nmod_poly_clear(r);
    return 0;
}
```
```output
3 errors:
  V = x^3+865*x^2+618*x+68
  error 1 at position 7
  error 100 at position 3
  error 928 at position 11
  decoder returns 3, codeword restored: 1
6 errors:
  decoder returns 6, codeword restored: 1
9 errors:
  decoder returns -1, codeword restored: 0
```

The monic $V$ printed for three errors is $(x - 3^7)(x - 3^3)(x - 3^{11}) =
(x - 329)(x - 27)(x - 637)$. With more than $T$ errors, detection is
likely but not certain: the received word can also be within distance $T$
of another codeword, and is then "corrected" to the wrong one. The
`goto cleanup` checks catch the other cases: a locator polynomial of
degree above $T$, one that does not explain all syndromes, one without
enough roots, or a root that is not a position of the word.

## Step 4: long codes {#long}

Over $\mathbb{F}_{65537}$, with $\alpha = 9 = 3^2$ of order $2^{15}$, we
can build codes of length $N = 32768$. With $2T = 4096$ check symbols the
decoder handles 2048 errors, and every step must now be quasi-linear:

- the generator polynomial of degree 4096 and the encoding remainder use
  fast multiplication and [Newton division](nmod_poly/divrem);
- the syndromes are the values of $r$ at $1, \alpha, \ldots, \alpha^{2T}$,
  a geometric progression with ratio $\alpha = 3^2$, so
  `nmod_poly_evaluate_geometric_nmod_vec_fast` computes all of them with
  one middle product, by Bluestein's trick [@Blu1970], instead of a
  subproduct tree (see [](nmod_poly/multipoint#geometric));
- the Berlekamp–Massey algorithm runs on 4096 terms with half-GCD steps;
- `nmod_poly_find_distinct_nonzero_roots` splits the locator of degree
  2048, where a Chien search would evaluate it at 32768 points;
- Forney's formula is evaluated at all inverse locators at once with
  `nmod_poly_evaluate_nmod_vec`, after inverting all locators with a single
  modular inversion (`_nmod_vec_invert`, see
  [](algorithms/arithmetic/batch-inversion)).

```c example name=tutorial_rs_long
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

#define P 65537     /* Fermat prime 2^16 + 1 */
#define R 3         /* primitive root; the code uses alpha = R^2 */
#define N 32768     /* order of alpha = 9 */
#define T 2048      /* N - K = 2T */

int main(void)
{
    flint_rand_t state;
    nmod_t mod;
    nmod_poly_t g, c, r, V, Lambda, S, Omega, dLambda;
    nmod_berlekamp_massey_t B;
    nn_ptr pts, s, X, Xinv, num, den;
    slong * logtab;
    slong i, nu;
    ulong alpha, t;

    flint_rand_init(state);
    nmod_init(&mod, P);
    nmod_poly_init_mod(g, mod);
    nmod_poly_init_mod(c, mod);
    nmod_poly_init_mod(r, mod);
    nmod_poly_init_mod(V, mod);
    nmod_poly_init_mod(Lambda, mod);
    nmod_poly_init_mod(S, mod);
    nmod_poly_init_mod(Omega, mod);
    nmod_poly_init_mod(dLambda, mod);
    pts = _nmod_vec_init(2 * T + 1);
    s = _nmod_vec_init(2 * T + 1);
    X = _nmod_vec_init(2 * T);
    Xinv = _nmod_vec_init(2 * T);
    num = _nmod_vec_init(2 * T);
    den = _nmod_vec_init(2 * T);
    logtab = flint_malloc(P * sizeof(slong));

    /* encode a random message */
    alpha = nmod_mul(R, R, mod);
    pts[0] = alpha;
    for (i = 1; i < 2 * T; i++)
        pts[i] = nmod_mul(pts[i - 1], alpha, mod);
    nmod_poly_product_roots_nmod_vec(g, pts, 2 * T);
    nmod_poly_rand(c, state, N - 2 * T);
    nmod_poly_shift_left(c, c, 2 * T);
    nmod_poly_rem(r, c, g);
    nmod_poly_sub(c, c, r);

    /* T errors at distinct random positions */
    nmod_poly_set(r, c);
    for (i = 0; i < T; i++)
    {
        slong e;
        do e = n_randint(state, N);
        while (nmod_poly_get_coeff_ui(r, e) != nmod_poly_get_coeff_ui(c, e));
        t = 1 + n_randint(state, P - 1);
        nmod_poly_set_coeff_ui(r, e, nmod_add(nmod_poly_get_coeff_ui(r, e), t, mod));
    }

    /* syndromes r(alpha^i), i = 0..2T, at a geometric progression */
    nmod_poly_evaluate_geometric_nmod_vec_fast(s, r, R, 2 * T + 1);

    /* error locator from s_1..s_2T */
    nmod_berlekamp_massey_init(B, P);
    nmod_berlekamp_massey_add_points(B, s + 1, 2 * T);
    nmod_berlekamp_massey_reduce(B);
    nmod_poly_set(V, nmod_berlekamp_massey_V_poly(B));
    nu = nmod_poly_degree(V);
    flint_printf("degree of the error locator: %wd\n", nu);

    if (!nmod_poly_find_distinct_nonzero_roots(X, V))
        flint_printf("decoding failure\n");
    else
    {
        /* Forney's formula at all error locators at once */
        nmod_poly_reverse(Lambda, V, nu + 1);
        nmod_poly_fit_length(S, 2 * T);
        _nmod_vec_set(S->coeffs, s + 1, 2 * T);
        _nmod_poly_set_length(S, 2 * T);
        _nmod_poly_normalise(S);
        nmod_poly_mullow(Omega, S, Lambda, 2 * T);
        nmod_poly_derivative(dLambda, Lambda);
        _nmod_vec_invert(Xinv, X, nu, mod);
        nmod_poly_evaluate_nmod_vec(num, Omega, Xinv, nu);
        nmod_poly_evaluate_nmod_vec(den, dLambda, Xinv, nu);

        for (i = 0, t = 1; i < N; i++, t = nmod_mul(t, alpha, mod))
            logtab[t] = i;
        for (i = 0; i < nu; i++)
        {
            slong e = logtab[X[i]];
            ulong y = nmod_neg(nmod_div(num[i], den[i], mod), mod);
            nmod_poly_set_coeff_ui(r, e, nmod_sub(nmod_poly_get_coeff_ui(r, e), y, mod));
        }
        flint_printf("codeword restored: %d\n", nmod_poly_equal(r, c));
    }

    nmod_berlekamp_massey_clear(B);
    nmod_poly_clear(g); nmod_poly_clear(c); nmod_poly_clear(r);
    nmod_poly_clear(V); nmod_poly_clear(Lambda); nmod_poly_clear(S);
    nmod_poly_clear(Omega); nmod_poly_clear(dLambda);
    _nmod_vec_clear(pts); _nmod_vec_clear(s); _nmod_vec_clear(X);
    _nmod_vec_clear(Xinv); _nmod_vec_clear(num); _nmod_vec_clear(den);
    flint_free(logtab);
    flint_rand_clear(state);
    return 0;
}
```
```output
degree of the error locator: 2048
codeword restored: 1
```

The whole program takes a few tens of milliseconds. A decoder that is
quadratic in $N$ would take seconds.

## Going further {#further}

- The error locator can also be found from the extended Euclidean
  algorithm applied to $x^{2T}$ and $S(x)$, stopped when the remainder
  has degree less than $T$ [@SKHN1975]; then the remainder is $\Omega$ up
  to a constant, and $\Lambda$ is the cofactor of $S$. `nmod_poly_hgcd`
  performs the first half of this computation in one step.
- Gao's decoder [@Gao2003] works directly with the received values in the
  evaluation view of the code, using interpolation
  (`nmod_poly_interpolate_nmod_vec_fast`) and a partial GCD.
- *Erasures*, errors at known positions, are handled by multiplying the
  syndrome series by the polynomial $\prod (1 - \alpha^{e} x)$ of the
  known positions; then up to $2T$ erasures can be corrected.
