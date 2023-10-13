.. _examples:

**Examples**
===============================================================================

Example programs
-------------------------------------------------------------------------------

FLINT comes with example programs to demonstrate current and future FLINT
features. To build the example programs run::

    make examples

The example programs are built in the ``build/examples`` directory.

For Arb and Calcium there are separate example pages
:ref:`examples-arb` and :ref:`examples-calcium`. Below are some
general examples.

- ``partitions`` Demonstrates the partition counting code, e.g.
  ``build/examples/partitions 1000000000`` will compute the number of
  partitions of ``10^9``.

- ``delta_qexp`` Computes the `n`-th term of the delta function, e.g.
  ``build/examples/delta_qexp 1000000`` will compute the one million-th
  term of the `q`-expansion of delta.

- ``crt`` Demonstrates the integer Chinese Remainder code, e.g.
  ``build/examples/crt 10382788`` will build up the given integer from its
  value mod various primes.

- ``multi_crt`` Demonstrates the fast tree version of the integer Chinese
  Remainder code, e.g. ``build/examples/multi_crt 100493287498239 13`` will
  build up the given integer from its value mod the given number of primes.

- ``stirling_matrix`` Generates Stirling number matrices of the first and
  second kind and computes their product, which should come out as the
  identity matrix. The matrices are printed to standard output. For example
  ``build/examples/stirling_matrix 10`` does this with 10 by 10 matrices.

- ``fmpz_poly_factor_zassenhaus`` Demonstrates the factorisation of a small
  polynomial. A larger polynomial is also provided on disk and a small
  (obvious) change to the example program will read this file instead of
  using the hard coded polynomial.

- ``padic`` Gives examples of the usage of many functions in the padic
  module.

- ``fmpz_poly_q`` Gives a very simple example of the ``fmpz_poly_q`` module.

- ``fmpq_poly`` Gives a very simple example of the ``fmpq_poly`` module.
