# Generic rings in C (experimental)

Generic Flint-style rings in C using void pointers + context objects.

The idea is to support fully generic recursive constructions of rings
(polynomials, fractions, power series, matrices, etc.) over arbitrary
user-defined base rings.

## Principles/goals/benefits

* Small code size, fast compilation.
* Possible to pack data efficiently (down to 1 byte / element).
* Plain C, similar interface to existing Flint code.
* Data layouts backwards compatible with most existing Flint types.
* Support all unusual cases in Flint/Arb/Calcium uniformly (error handling,
  inexact rings, noncomputable rings, context objects), with a uniform
  interface.
* Support fast stack based allocation of temporary variables and arrays
  (with safe limits).
* Fully in-place operations.
* Allow fast (a few cycles) runtime construction of context objects.
  (May be on the order of hundreds of cycles when constructing a new
  method table; this could be optimized.)
* Possibility to use generic default methods or provide optimized versions
  (e.g. for vector operations).

## Disadvantages

* Potentially a few cycles overhead for each dynamic function lookup and
  function call. Precludes intra-method compiler optimizations and
  inlining (however, providing custom vector methods should partially
  compensate for this).
* Less compiler protection against type errors.
* Some function call signatures need to change in ways that make
  the interface less convenient (for uniformity).

