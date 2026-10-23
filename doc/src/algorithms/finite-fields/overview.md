---
title: Arithmetic in finite fields
kind: topic
summary: FLINT represents $\mathbb{F}_{p^d}$ by Zech logarithms, by polynomials over a word-size or a multiprecision prime field, or, for $d = 1$, by residues; this page compares the representations, says how `fq_default` chooses among them, and lists the algorithms for reduction, inversion, powering, Frobenius, norm, trace, embeddings and defining polynomials.
order: 0
---

A finite field with $q = p^d$ elements is unique up to isomorphism, but
computing in it requires a concrete model. FLINT, like most systems,
uses the model $\mathbb{F}_q = \mathbb{F}_p[x]/(f)$ for an irreducible
polynomial $f$ of degree $d$, the *defining polynomial*, and stores
elements either as polynomials of degree less than $d$ or, for small
fields, as discrete logarithms with respect to a root of $f$. The
standard references for the mathematics are Lidl and Niederreiter
[@LN1997] and, for the algorithms, von zur Gathen and Gerhard
[@GG2003] and Shoup [@Sho2009].

## The representations {#representations}

| Module | Elements | Characteristic | Product | Sum |
|---|---|---|---|---|
| `fq_zech` | discrete logarithm $k$ of $a = g^k$, a word | $q = p^d < 2^{\bits}$ | addition modulo $q - 1$ | one lookup in a table of $q$ words |
| `fq_nmod` | polynomial of degree $< d$, `nmod_poly` | $p < 2^{\bits}$ | polynomial product and reduction | $d$ additions modulo $p$ |
| `fq` | polynomial of degree $< d$, `fmpz_poly` | any $p$ | polynomial product and reduction | $d$ multiprecision additions |
| `nmod` | residue in $[0, p)$ | $p < 2^{\bits}$, $d = 1$ | modular multiplication | modular addition |
| `fmpz_mod` | residue in $[0, p)$ | any $p$, $d = 1$ | modular multiplication | modular addition |

The three extension field modules `fq_zech`, `fq_nmod` and `fq` share
their interface: their headers and most of their code are generated from
the same templates (`fq_templates.h`, `fq_poly_templates.h` and so on),
and the vector, matrix, polynomial, factorisation and embedding modules
exist for each of them. `fq_default` wraps all five representations
behind one interface and chooses the representation when a field is
created ([](fq_default/getting-started)):

- Zech logarithms if $d > 1$ and $\operatorname{bits}(p) \cdot d \le 16$,
  so that $q < 2^{16}$ and the tables of $2q + p$ words fit in about a
  megabyte, in the cache of current processors;
- `fq_nmod` for other extension fields with $p < 2^{\bits}$;
- `fq` for extension fields with larger $p$;
- `nmod` and `fmpz_mod` for prime fields.

The Zech representation needs a *primitive* defining polynomial, one
whose root generates the cyclic group $\mathbb{F}_q^\times$; when a
non-primitive polynomial is given, `fq_default` falls back to `fq_nmod`.

## Costs {#costs}

In the polynomial representations a product costs a product of two
polynomials of degree less than $d$ over $\mathbb{F}_p$ (or over $\Z$, for
`fq`) and a reduction modulo $f$, that is $O(M(d))$ operations in
$\mathbb{F}_p$, where $M(d)$ is the cost of polynomial multiplication;
for the small $d$ that are common in practice this is $O(d^2)$ word
operations plus the allocation and normalisation overhead of the
polynomial types. An inverse costs an extended Euclidean algorithm in
$\mathbb{F}_p[x]$, $O(d^2)$ operations. In the Zech representation a
product, quotient, inverse or power costs one or two word operations,
and a sum one table lookup, but building the tables costs $q$ products
when the context is created.

The following table shows timings measured with the installed library on
the x86-64 machine used to build this documentation, through the
`fq_default` interface: one multiplication and one addition of field
elements, and the product of two polynomials of length 1000 over the
field.

| Field | `fq_zech` | `fq_nmod` | `fq` |
|---|---|---|---|
| $\mathbb{F}_{3^8}$, mul and add | 8 ns | 330 ns | 740 ns |
| $\mathbb{F}_{3^8}$, polynomial product | 0.6 ms | 0.7 ms | 2.1 ms |
| $\mathbb{F}_{2^{12}}$, mul and add | 10 ns | 470 ns | 1300 ns |
| $\mathbb{F}_{2^{12}}$, polynomial product | 1.2 ms | 1.3 ms | 4.1 ms |
| $\mathbb{F}_{101^3}$, mul and add | 150 ns | 100 ns | 160 ns |
| $\mathbb{F}_{101^3}$, polynomial product | 0.3 ms | 0.2 ms | 0.4 ms |

Zech logarithms are faster by more than an order of magnitude for
element arithmetic in small fields, but lose their advantage when the
tables no longer fit in the cache ($q \approx 10^6$ in the last rows),
and for long polynomials, whose products are computed by converting the
coefficients to the polynomial representation and using fast
multiplication over $\Z/p\Z$ in all three modules ([](fq_poly/mul)).

## Algorithms {#algorithms}

| Task | Method | Functions | Page |
|---|---|---|---|
| Reduction modulo $f$ | elimination with the nonzero terms of $f$ if it has at most 5; Newton division with a precomputed inverse otherwise | `_fq_reduce`, `_fq_nmod_reduce` | [](algorithms/finite-fields/reduction) |
| Sums of products | accumulate unreduced products, reduce once | `_fq_vec_dot`, `fq_mat_mul_classical`, `fq_poly_mul_classical` | [](algorithms/finite-fields/reduction) |
| Inversion | extended Euclidean algorithm in $\mathbb{F}_p[x]$ | `fq_inv`, `fq_nmod_inv` | [](algorithms/polynomials/gcd) |
| Powering | binary exponentiation, exponent reduced modulo $q - 1$ | `fq_pow`, `fq_nmod_pow` | |
| Zech arithmetic | table of $Z(k)$ with $1 + g^k = g^{Z(k)}$ | `fq_zech_add`, `fq_zech_ctx_init_fq_nmod_ctx` | [](algorithms/finite-fields/zech-logarithms) |
| Frobenius map | exponentiation by $p^e$ | `fq_frobenius` | [](algorithms/finite-fields/frobenius-norm-trace) |
| Trace | traces of $1, x, \ldots, x^{d-1}$ by Newton's identities | `fq_trace` | [](algorithms/finite-fields/frobenius-norm-trace) |
| Norm | resultant $\operatorname{Res}(f, a)$ | `fq_norm` | [](algorithms/finite-fields/frobenius-norm-trace) |
| Square roots | Tonelli–Shanks in $\mathbb{F}_q^\times$; Frobenius inverse for $p = 2$ | `fq_sqrt` | [](algorithms/finite-fields/frobenius-norm-trace) |
| Polynomial multiplication | Kronecker substitution, substitution $y = x^{2d-1}$, reordering | `fq_poly_mul` | [](fq_poly/mul), [](algorithms/polynomials/kronecker) |
| Matrix multiplication | Kronecker substitution into `fmpz_mat` | `fq_mat_mul_KS` | [](algorithms/linear-algebra/matrix-multiplication) |
| Factorisation, roots | generic `gr_poly` algorithms | `fq_poly_factor`, `fq_poly_roots` | [](algorithms/polynomials/factoring-finite-fields) |
| Embeddings $\mathbb{F}_{p^m} \to \mathbb{F}_{p^n}$ | root of the defining polynomial, linear algebra and the trace | `fq_embed_gens`, `fq_embed_matrices` | [](algorithms/finite-fields/embeddings) |
| Defining polynomials | Conway polynomials from a table; sparse or minimal-weight irreducible, or random primitive, polynomials | `fq_ctx_init`, `fq_nmod_ctx_init_ui`, `fq_zech_ctx_init_ui` | [](algorithms/finite-fields/defining-polynomials) |

Polynomials, matrices and power series over $\mathbb{F}_q$ mostly use the
generic algorithms of `gr_poly` and `gr_mat`, with the multiplication
of the respective module: division and power series by Newton iteration
([](algorithms/polynomials/newton)), GCDs by the half-GCD algorithm
([](algorithms/polynomials/half-gcd)), modular composition by the
Brent–Kung algorithm ([](algorithms/polynomials/composition)), and
linear algebra by Gaussian elimination
([](algorithms/linear-algebra/gaussian-elimination)).
