---
title: The design of generic rings
kind: topic
summary: An overview of the generic rings framework: run-time ring objects with method tables, headerless elements, explicit status codes and three-valued predicates, and how these choices make generic code both rigorous and fast.
order: 0
---

A generic algorithm is written once and runs over any ring that
provides the operations it needs. Computer algebra systems achieve this
in different ways: C++ libraries use templates, which compile a separate
copy of each algorithm for each element type; Julia's Nemo and
AbstractAlgebra [@FHHJ2017] rely on just-in-time compilation of generic
code for each type; interpreted systems dispatch on dynamically typed
objects. FLINT's framework [`gr`](gr), introduced in FLINT 3 and
described in [@Joh2025gr], works in plain C with run-time dispatch: a
ring is an object, operations are looked up in a table of function
pointers, and the critical kernels are specialised by hand rather than
by a compiler.

## Rings as objects {#objects}

A ring is a context object `gr_ctx_t` holding the size of an element, a
pointer to a method table indexed by the enumeration `gr_method`, an
implementation code, and a few words of parameters (a modulus, a
precision, a base ring). Elements carry no type information: an element
is `sizeof_elem` bytes of the ring's own element type, passed as a
`void` pointer together with its context. This has three consequences
([](algorithms/generic/method-tables)):

* Element arrays are packed and identical to the arrays of the
  type-specific modules, so a `gr` function and an `nmod_poly` function
  can work on the same coefficient array without conversion. Hundreds
  of type-specific functions are now wrappers around generic
  implementations.
* A call costs one table lookup and one indirect call, a few cycles.
  This is significant only for very cheap elements, and is removed
  entirely by overriding vector operations, at which point generic and
  type-specific code run at the same speed.
* Creating a ring is cheap: most constructors allocate nothing, so an
  algorithm can create a context for each of millions of primes in a
  multimodular computation.

Rings are typed strictly. There is no implicit coercion: $1/2$ in $\Z$ is
an error, not a rational number, and combining elements of different
rings requires an explicit conversion (`gr_set_other`) or a mixed
operation (`gr_mul_other`). Properties such as "commutative",
"integral domain" or "field" are predicates of the context rather than
classes in a hierarchy, and generic algorithms branch on them.

## Specialisation in three levels {#levels}

Operations are organised in three levels: elementwise operations
($z = x + y$), vector operations ($z_i = x_i + y_i$, dot products), and
operations of higher complexity (polynomial and matrix
multiplication, division, GCD, linear algebra). Every method has a
default implemented in terms of the level below or the same level, so a
new ring needs only about fifteen methods to support everything; and a
ring gains speed by overriding the methods of the levels that matter,
most importantly dot products. Parametric types such as `mpn_mod`
(moduli of 2 to 16 limbs) and `nfloat` (floating point of 1 to 66 limbs)
choose a fixed element size at run time and switch to kernels for that
size internally, which gives them shallow elements and fixed-length
arithmetic without compiling a copy of every algorithm per size
([](algorithms/generic/method-tables#parametric)).

## Rigour: enclosures, truth values and status codes {#rigour}

The framework uses the word "ring" in its strict sense: a context that
claims to be a ring must satisfy the axioms exactly. Inexact
representations are allowed only with enclosure semantics: an element
of an inexact ring (a real ball, a power series with an error term)
represents a set containing the exact value, and every operation must
return a set containing every possible exact result
([](algorithms/generic/inexact)). Floating-point contexts, which do not
satisfy this, report that they are not rings.

Because elements may be inexact, and because some questions are
undecidable or too expensive (is this exact real number zero? is this
large modulus prime?), predicates return a `truth_t` with a third value
`T_UNKNOWN`, and generic algorithms must handle it. Polynomials over
inexact rings are only *weakly* normalised: a leading coefficient that
might be zero is kept, so the length bounds the degree, and operations
that need the exact degree fail cleanly
([](algorithms/generic/inexact#weak-normalisation)).

Failure is reported by status codes rather than by exceptions or junk
values. `GR_DOMAIN` says the result does not exist; `GR_UNABLE` says it
was not computed. The distinction lets algorithms choose fallbacks: a
fast numerical method that is unable to decide is followed by a slower
exact one, while a domain error is final
([](algorithms/generic/errors)). Rings that pretend to be fields turn
the discovery of a zero divisor into a recoverable failure, which gives
dynamic evaluation [@DDDD1985] at no extra cost
([](algorithms/generic/errors#dynamic-evaluation)).

## Correctness by testing {#testing}

Every ring is tested against the algebraic laws it claims: associativity,
distributivity, the consistency of the many variants of each operation,
aliasing, and the interaction with vector and matrix operations, using
random elements from deliberately skewed distributions. The test
framework is generic as well, so a new ring is tested by one call to
`gr_test_ring`, and the tests remain valid for inexact rings because a
test fails only on a provable contradiction
([](algorithms/generic/testing)).

## Where the framework is used {#use}

| Module | Generic algorithms |
|---|---|
| [`gr_vec`](gr_vec) | elementwise arithmetic, dot products, sums and products by binary splitting |
| [`gr_poly`](gr_poly) | multiplication, division and Newton iteration, half-GCD, resultants, composition, evaluation and interpolation, root finding and factoring over finite fields |
| [`gr_mat`](gr_mat) | Strassen–Winograd multiplication, LU and fraction-free LU, solving, determinants, characteristic polynomials, eigenvalues, matrix functions |
| [`gr_series`](gr_series), [`gr_mpoly`](gr_mpoly), [`gr_ore_poly`](gr_ore_poly) | power series with error terms, sparse multivariate polynomials, Ore polynomials |
| [`gr_special`](gr_special) | special functions with generic fallbacks (factorials, binomials, Fibonacci, Bernoulli and Stirling numbers, Chebyshev polynomials) |
| [`gr_dft`](gr_dft) | discrete Fourier transforms of any length |

The overview pages [](algorithms/polynomials/overview) and
[](algorithms/linear-algebra/overview) say which generic algorithm each
type-specific function uses.
