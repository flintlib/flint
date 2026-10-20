---
title: Testing, bugs and checking results
summary: How FLINT is tested, why bugs that give wrong results still occur, how to check results independently, debug builds with assertions, and how to report a bug.
order: 2
---

FLINT is extensively tested, but it is not formally verified: neither
the implementations nor, in most cases, the proofs behind them have been
checked by a machine. Bugs that make a function return a mathematically
wrong result are rare, but they are found and fixed from time to time.
Results that matter should therefore be checked independently.

## How FLINT is tested {#test-suite}

The test suite consists of about 4000 test programs, mostly one for
each public function or group of functions, run by `make check` (see [](guide/getting-started/building#testing)).
Most tests are randomised: they generate many inputs with the
`randtest` functions, which favour special cases such as zero and
numbers with long runs of zero and one bits, and check properties
instead of fixed answers. Typical checks compare two
algorithms for the same operation (Karatsuba against classical
multiplication), verify identities ($(a b) c = a (b c)$, $q b + r = a$
with $\deg r < \deg b$), check the results with aliased arguments, and,
for ball arithmetic, check that the result contains the exact value
computed independently. The environment variable
`FLINT_TEST_MULTIPLIER` scales the number of iterations. The test suite
is run in continuous integration on several platforms and
configurations, including builds with assertions enabled.

Random testing finds most bugs, but not all of them. A defect may only
show for inputs that random generation is unlikely to produce: a
particular combination of lengths, a modulus in a narrow range, or a
coincidence that occurs once in thousands of cases. For example,
`nmod_poly_mul` in FLINT 3.6 and earlier could return wrong products
modulo $n \le 23$ for some lengths above about 370. The bug was found by
reading the code, not by the test suite, which did not generate such
inputs.

## Checking results {#checking}

Independent checks catch bugs in FLINT, and also mistakes in using it,
such as a violated precondition or a wrong aliasing assumption. Useful
techniques:

- *Verify the defining property.* A factorisation can be multiplied
  out and its factors tested for primality or irreducibility; a
  solution $x$ of $A x = b$ can be substituted; a GCD must divide both
  inputs; an inverse modulo $n$ multiplied by the original must give 1.
  Verifying is usually much cheaper than computing.
- *Compute in two ways.* Use a different algorithm (most operations
  have explicitly named variants, such as `nmod_poly_mul_classical`), a
  different representation (`fmpz_poly` against `gr_poly` over
  `gr_ctx_init_fmpz`), or another computer algebra system.
- *Check random specialisations.* Two polynomials or matrices that
  should be equal can be compared at random points modulo a random
  prime; by the Schwartz–Zippel lemma [@Sch1980] [@Zip1979] a
  difference is detected with high probability, at a fraction of the
  cost of the full computation.
- *Vary the parameters.* In ball arithmetic, repeat the computation at
  a higher precision: the two balls must overlap, and the second should
  usually be smaller. With multimodular or randomised algorithms,
  change the random seed or the number of threads.
- *Check identities and known values.* Special functions satisfy
  functional equations and have known values at special points;
  number-theoretic functions have known values in tables such as the
  OEIS.

The following program checks a product computed by `nmod_poly_mul` in
two ways: by evaluating both sides at random points, and by recomputing
it with the classical algorithm.

```c example name=guide_correctness_check
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>

int main(void)
{
    flint_rand_t state;
    nmod_poly_t a, b, c, d;
    ulong n = 1000003, x;
    slong i, len = 2000;
    int ok = 1;

    flint_rand_init(state);
    nmod_poly_init(a, n);
    nmod_poly_init(b, n);
    nmod_poly_init(c, n);
    nmod_poly_init(d, n);

    nmod_poly_randtest(a, state, len);
    nmod_poly_randtest(b, state, len + 1);
    nmod_poly_mul(c, a, b);

    /* c(x) = a(x) b(x) at random points */
    for (i = 0; i < 10; i++)
    {
        x = n_randint(state, n);
        if (nmod_poly_evaluate_nmod(c, x) !=
            nmod_mul(nmod_poly_evaluate_nmod(a, x),
                     nmod_poly_evaluate_nmod(b, x), c->mod))
            ok = 0;
    }
    flint_printf("evaluation check: %s\n", ok ? "passed" : "FAILED");

    /* the same product with another algorithm */
    nmod_poly_mul_classical(d, a, b);
    flint_printf("classical algorithm agrees: %s\n",
        nmod_poly_equal(c, d) ? "yes" : "NO");

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    nmod_poly_clear(d);
    flint_rand_clear(state);
    return 0;
}
```

```output
evaluation check: passed
classical algorithm agrees: yes
```

If the computed product $c$ differs from $ab$, the nonzero polynomial
$c - ab$ has degree at most 4000, so it vanishes at a random point
modulo the prime $n = 1000003$ with probability at most $4000/n$, and
ten points miss the error with probability below $10^{-20}$. Modulo a
small prime this check is weak (modulo 23 every polynomial takes the
same values as its remainder modulo $x^{23} - x$), and one should
evaluate over an extension field or compare with another algorithm, as
the second check does.

## Assertions and debug builds {#assertions}

FLINT contains many internal consistency checks and checks of
preconditions written with the macro `FLINT_ASSERT`. They are compiled
only when FLINT is configured with `--enable-assert`, which defines
`FLINT_WANT_ASSERT` in the installed `flint-config.h`; in a normal build
`FLINT_ASSERT` expands to nothing and costs nothing. An assertion
build is several times slower, and is meant for development and for
finding errors in programs that use FLINT: a failed assertion stops the
program (through the C `assert` macro) at the first point where an
invariant is violated, which is often much earlier than the point where
a wrong result or a crash becomes visible. Some preconditions of public
functions, such as compatible matrix dimensions in `fmpz_mat_mul` and
`nmod_mat_mul`, are checked only in this way. In assertion builds the
`TMP_ALLOC` scratch memory is also taken from the heap, so that memory
checkers see it.

The configure option `--enable-debug` (the default) only adds debugging
information (`-g`); it does not enable assertions. To debug a program
that uses FLINT, combine an assertion build with `-O0 -g` in `CFLAGS`
and run the program under Valgrind or with the address sanitizer
(`--enable-asan`); see [](guide/getting-started/building#debugging).

## Reporting bugs {#reporting}

If a result looks wrong, first reduce the problem to a small complete
program, check that the input satisfies the documented preconditions,
and try it with an assertion build. Then report it as described in
[](guide/development/bug-reporting). Reports of wrong results are
especially welcome, including those found by comparing with other
software.
