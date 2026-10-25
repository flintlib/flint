---
title: Rigorous numerics
short: Numerics
summary: Programs that compute with balls: constants, special functions, roots, integrals and zeta zeros, with guaranteed error bounds and adaptive precision.
order: 5
---

These tutorials compute with real and complex numbers in ball arithmetic
([](arb), [](acb) and the modules built on them). A ball
$[m \pm r]$ is a floating-point midpoint $m$ with an error bound $r$, and
every operation returns a ball that is guaranteed to contain the exact
result for all inputs in the input balls. A program built from such
operations therefore produces a mathematically proved enclosure of the
true answer: the digits it prints are correct, not estimates. What ball
arithmetic does not guarantee is that the enclosure is narrow. The
recurring theme of these tutorials is how to get narrow enclosures at a
reasonable cost:

- choosing the working precision, and increasing it when the result is
  not accurate enough (the precision is an input; the accuracy is
  measured on the output);
- adding the errors that ball arithmetic cannot see, such as the
  truncation error of a series or of an iteration, with explicit bounds;
- writing formulas that do not overestimate the error, and recognising
  when a wide ball is caused by the formula rather than by the problem;
- using the library functions that do all of this internally for
  standard problems: constants, special functions, polynomial roots,
  root isolation for real functions, integration, zeros of the zeta
  function.

The first tutorials introduce the basic techniques with constants and
simple iterations; the later ones apply them to roots, integrals and the
Riemann zeta function. The last one, [](tutorials/numerics/padic-interpolation),
transfers the theme to $p$-adic numbers, whose balls $x + O(p^N)$ raise
the same question of how much precision an algorithm can guarantee. Two of them develop a complete rigorous
computation from scratch: [](tutorials/numerics/series) and
[](tutorials/numerics/interval-newton).

Each tutorial also describes the related programs in FLINT's `examples`
directory: what they compute, their options and sample output. They are
built with `make examples` (see
[](tutorials/getting-started/example-programs)); the sample output shown
was produced with `build/examples/<name>`, and the timing and memory
lines it contains vary between machines.

For background, see the user guide chapter [](guide/types/ball-arithmetic) for
the semantics of balls and the conventions of the ball arithmetic
modules, and the overviews of the algorithms for
[ball arithmetic](algorithms/ball-arithmetic/overview) and for
[special functions](algorithms/special-functions/overview).
