---
title: Randomised testing of ring axioms
summary: How ring implementations are tested by checking algebraic identities on random elements from skewed distributions, why such tests stay valid for inexact rings, and how generic algorithms are tested over many rings at once.
complexity: One test of an identity costs a few ring operations; gr_test_ring runs about 90 tests, each a given number of times.
order: 40
---

A ring implementation consists of hundreds of methods, many of them
specialised for speed, and generic algorithms combine them in ways that
no implementer anticipates. FLINT tests rings the way algebra defines
them: by checking the identities a ring must satisfy on random elements,
in the style of property-based testing [@ClaHug2000]. The test framework
is itself generic, so one function, `gr_test_ring`, tests any ring
[@Joh2025gr].

## Identities as tests {#identities}

Each test draws random elements, evaluates both sides of an identity
with the methods of the ring, and compares the results with `gr_equal`.
The identities checked by `gr_test_ring` include:

* the ring axioms: associativity of $+$ and $\cdot$, commutativity of
  $\cdot$ (if the ring claims to be commutative), left and right
  distributivity, $x - y = x + (-y)$, and the properties of $0$ and $1$;
* consistency of variants: `gr_add_ui`, `gr_add_si`, `gr_add_fmpz` and
  `gr_add_fmpq` agree with `gr_add` after conversion, and likewise for
  subtraction, multiplication, `addmul`, `submul`, division, exact
  division and powers, and the `_other` variants agree with conversion
  followed by the plain operation;
* aliasing: each operation gives the same result when the output is one
  of the inputs;
* division: $(x / y) \cdot y = x$ and $(x y) / y = x$ when the quotients
  exist, $x^{-1} x = 1$, $(x^{-1})^{-1} = x$, consistency of
  `gr_is_invertible`, `gr_div_nonunique` and `gr_divides`, and no zero
  divisors in a ring that claims to be an integral domain;
* powers: $x^{a+b} = x^a x^b$, $(cx)^a = c^a x^a$ for an integer $c$,
  $(xy)^a = x^a y^a$ in commutative rings, and the same for `fmpz` and
  general exponents;
* square roots, `gr_cmp` and `gr_cmpabs` in ordered rings, real and
  complex parts, numerators and denominators, canonical associates, and
  in unique factorisation domains GCD, LCM and factorisation;
* conversions and printing: round trips through `ulong`, `slong`,
  `fmpz`, `fmpq`, `double`, strings and symbolic expressions;
* vector operations against the elementwise ones, dot products (also
  strided, if overridden) against loops, and associativity of classical
  matrix multiplication.

Testing the variants against each other is what makes specialisation
safe: a ring that overrides `_gr_vec_dot` or `gr_mul_si` for speed is
checked against its own plain `gr_mul` and `gr_add`.

## Soundness for inexact and partial rings {#soundness}

A test reports a failure only when every operation succeeded and
`gr_equal` returns `T_FALSE`, that is, when the two sides are *proven*
different. An operation that returns `GR_DOMAIN` or `GR_UNABLE`, or a
comparison that returns `T_UNKNOWN`, counts as a pass (with the flag
`GR_TEST_ALWAYS_ABLE`, an unable result counts as a failure instead).

::: proposition {#prop-sound-test}
If a ring satisfies the enclosure semantics of
[](algorithms/generic/inexact#def-enclosure) and its predicates are
sound, then `gr_test_ring` never reports a failure for it. Conversely,
any reported failure exhibits elements for which a method violates
either an identity of the ring or the enclosure semantics.
:::

::: proof
Both sides of each identity are computed by programs that enclose the
exact values of the two sides
([](algorithms/generic/inexact#thm-composition)), and the exact values
are equal because the identity holds in the ring. Two enclosures of the
same value intersect, so a sound `gr_equal` cannot return `T_FALSE`.
:::

This is what lets the same tests run on real and complex balls, power
series with error terms, exact real numbers with undecidable equality,
and rings that implement only part of the interface. It also explains
why floating-point contexts are tested by `gr_test_floating_point` and
the approximate tests (`gr_test_approx_unary_op`, `gr_test_approx_dot`
and others), which compare with a more precise reference ring up to a
relative tolerance: floating-point results are not enclosures, so exact
identities fail.

## Random elements {#random}

Uniformly random elements are poor at finding bugs: a random 64-bit
residue almost never triggers a carry into a third word, and a random
large integer almost never has the bit pattern that exposes a
normalisation error. The `randtest` methods are therefore deliberately
skewed. `n_randtest` and `fmpz_randtest` produce numbers with long runs
of zero and one bits, small and extreme values; ball generators produce
exact and inexact balls with exponents of varying size; polynomial and matrix
generators produce sparse and short objects as well as full ones; and
`_gr_vec_randtest` makes each entry zero with probability $1/2$. The
default `gr_randtest_small` returns small combinations $a g + b$ of the
generator, for tests whose cost grows quickly with the size of the
elements. If a bug manifests on a fraction $\varepsilon$ of the inputs
drawn by a test, then $N$ iterations miss it with probability
$(1 - \varepsilon)^N \le e^{-\varepsilon N}$, so skewing the
distribution to raise $\varepsilon$ is far more effective than more
iterations.

## Testing generic algorithms over many rings {#many-rings}

Generic algorithms are tested the same way, over rings drawn by
`gr_ctx_init_random` and its variants: the tests of `gr_poly` and
`gr_mat` check, for example, that products computed by different
algorithms agree, that $q b + r = a$ after division, and that a GCD
divides both inputs, over integers, residue rings, finite fields, number
fields, balls, exact real numbers and composite rings. The same test
then exercises both the algorithm and the specialised methods of every
ring it runs on.

Two further tools target error handling. The debug wrapper
`gr_ctx_init_debug` can make the wrapped methods randomly return
`GR_UNABLE` and predicates randomly return `T_UNKNOWN`, which tests that
algorithms propagate failures instead of returning wrong results; and
with `GR_DEBUG_WRAP` it stores elements in boxes that record their
context, which detects elements used with the wrong context. Randomised
testing cannot prove an implementation correct; formal verification of
generic algorithms against the ring axioms would be the next step, and
is discussed in [@Joh2025gr].
