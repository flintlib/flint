---
title: Choosing a representation for exact numbers
kind: topic
summary: FLINT represents exact real and complex numbers as rationals, elements of a fixed number field, individual algebraic numbers (qqbar) or elements of formal fields with transcendental generators (ca). This page compares what each representation can decide, what it costs, and which algorithms the main operations use.
order: 0
---

An exact number is one that is represented symbolically, so that the
value is determined without error and can be approximated to any
precision. FLINT has four exact representations of increasing
generality, and they differ less in what they can express than in what
they can decide: whether two numbers are equal, whether a number is
real, positive or rational. The more general the representation, the
more expensive arithmetic becomes and the weaker the guarantees.

## The representations {#representations}

| Representation | Numbers | Element stored as | Equality | Typical cost |
|---|---|---|---|---|
| [](fmpq) | $\mathbb{Q}$ | reduced fraction | canonical | quasi-linear in the bit size |
| [](nf_elem) | a fixed number field $\mathbb{Q}(a)$ | polynomial in $a$ of degree less than $[\mathbb{Q}(a) : \mathbb{Q}]$ | canonical | polynomial arithmetic modulo the defining polynomial |
| [](qqbar) | $\overline{\mathbb{Q}}$ | minimal polynomial and isolating ball | decidable | resultants and factorisation; degrees multiply |
| [](ca) | real and complex numbers built from $\mathbb{Q}$ by arithmetic, roots, $\exp$, $\log$, $\pi$, trigonometric and some special functions | element of a formal field $\mathbb{Q}(a_1, \ldots, a_n)$ reduced modulo proved relations | semi-decidable: answers may be unknown | multivariate rational functions, relation searches, numerical evaluation |

The two exact algebraic representations behave very differently. A
[](nf_elem) is cheap, but all numbers must belong to one field fixed in
advance, which is not always known and whose degree can be large. A
[](qqbar) number carries its own minimal polynomial, so any two algebraic
numbers can be combined; arithmetic is self-contained and canonical, but
an operation may need a resultant and its factorisation, and the degree
of the result can be as large as the degree of the compositum, so long
computations become expensive (see [](algorithms/exact/qqbar-arithmetic)).

A [](ca) number represents its value as an element of a formal field
whose generators are algebraic numbers or transcendental extension
numbers such as $\pi$, $e^{\sqrt 2}$ or $\log 3$. Field arithmetic is
formal and cheap compared with resultants, several algebraic generators
can coexist without computing a primitive element, and transcendental
numbers can be represented at all. The price is that equality is not
canonical: a zero test combines numerical evaluation with exact
reductions and can fail (see [](algorithms/exact/zero-testing)). When
all generators are algebraic, `ca` can fall back on [](qqbar) and
decide equality, within the degree limit `CA_OPT_QQBAR_DEG_LIMIT`.

The [](gr) interface provides contexts for all of them (see
[](gr_domains#algebraic) and [](gr_domains#real-complex)): `qqbar`
through `gr_ctx_init_real_qqbar` and `gr_ctx_init_complex_qqbar`, and
`ca` through `gr_ctx_init_real_ca`, `gr_ctx_init_complex_ca`, their
algebraic variants and `gr_ctx_init_complex_extended_ca`. Generic
polynomial and matrix algorithms then work over exact numbers and return
`GR_UNABLE` when a zero test fails. Symbolic expressions ([](fexpr)) are
not numbers but trees; they are used to construct, print and convert
exact numbers.

## Which algorithm is used {#methods}

| Task | Method | Functions |
|---|---|---|
| Create an algebraic number | factor the defining polynomial, isolate the roots, certify with interval Newton | `qqbar_roots_fmpz_poly`; [](algorithms/exact/qqbar-representation#creating) |
| Refine an enclosure | certified Newton iteration, recomputation of all roots as fallback | `qqbar_get_acb`, `qqbar_cache_enclosure`; [](algorithms/exact/qqbar-representation#refining) |
| Compare algebraic numbers | minimal polynomials, then enclosures | `qqbar_equal`, `qqbar_cmp_re`; [](algorithms/exact/qqbar-representation#predicates) |
| Add or multiply algebraic numbers | closed forms for quadratics, subfield and cyclotomic detection, guessing with proof, composed polynomials and factor selection | `qqbar_add`, `qqbar_mul`, `qqbar_binary_op`; [](algorithms/exact/qqbar-arithmetic) |
| Create a field with given generators | ideal of relations: minimal polynomials, logarithmic and multiplicative relations, Gröbner basis | `ca_field_build_ideal`; [](algorithms/exact/field-towers) |
| Arithmetic in `ca` | formal arithmetic in the common field of the operands, reduction modulo the ideal | `ca_add`, `ca_mul`, `ca_merge_fields`; [](algorithms/exact/field-towers) |
| Decide $x = 0$ | numerical evaluation, algebraic evaluation, normal form rewriting, factoring | `ca_check_is_zero` ([](ca/check-is-zero)); [](algorithms/exact/zero-testing) |
| Simplify, recognise closed forms | canonical number field forms, rewriting rules, LLL with exact verification | `ca_condense_field`, `qqbar_guess`, `qqbar_get_fexpr_formula`; [](algorithms/exact/simplification) |
| Roots of polynomials | squarefree factorisation, rational factors via [](qqbar), formulas up to degree 3 | `ca_poly_roots` ([](ca_poly/roots)) |
| Linear algebra | simplest-first pivoting, fraction-free and division-free methods, rational and number field specialisations | `ca_mat_det`, `ca_mat_rank`, `ca_mat_eigenvalues`; [](algorithms/exact/linear-algebra) |

## Practical advice {#advice}

- Use [](fmpq) or [](nf_elem) when the numbers are known to lie in
  $\mathbb{Q}$ or in one known number field: arithmetic is much faster
  than with `qqbar` or `ca`.
- Use [](qqbar) for algebraic numbers when exact answers to all
  predicates are required and degrees stay moderate. Computations whose
  intermediate degrees explode, such as sums of many square roots, can be
  cheaper in [](ca), which keeps them in a field with several quadratic
  generators.
- Use [](ca) whenever transcendental numbers occur, and check the
  `truth_t` results of predicates: `T_UNKNOWN` is a legitimate answer.
  The options of the context (see [](ca/options)) bound the work spent on
  zero tests and relation searches.
- For numerical values only, [](arb) and [](acb) are much faster;
  `ca_get_acb` and `qqbar_get_acb` convert exact numbers to balls of any
  precision.
