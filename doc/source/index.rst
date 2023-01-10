.. generic-rings documentation master file, created by
   sphinx-quickstart on Sun Mar 13 11:37:54 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

generic-rings documentation
=========================================

These modules provide a plain C implementation of generic rings
intended for use with Flint, Antic, Arb and Calcium.
The idea is to support fully generic recursive constructions of rings
(polynomials, fractions, power series, matrices, etc.) over arbitrary
user-defined base rings (including all Flint/Antic/Arb/Calcium types)
with the ability to add specializations (e.g. ``fmpz_mat`` functions
instead of generic functions for matrices over ``fmpz``).
This code is *experimental*.

Contents
-------------------------------------------------------------------------------

.. toctree::
   :maxdepth: 2

   gr.rst
   gr_domains.rst
   gr_special.rst
   gr_mat.rst
   gr_poly.rst
   gr_mpoly.rst
   python_flint.rst

Design
-------------------------------------------------------------------------------

We use element pointers + context objects holding method pointers.
Vectors of elements can be packed densely without indirection.

Principles/goals/benefits
...............................................................................

* Plain C, similar interface to existing Flint code.
* Small code size, fast compilation.
* Possible to pack data efficiently (down to 1 byte / element).
* Data layouts backwards compatible with most existing Flint element types,
  and mostly also with polynomials, matrices, etc.
* Support all unusual cases in Flint/Arb/Calcium uniformly (error handling,
  inexact rings, noncomputable rings, context objects), with a uniform
  interface.
* Support fast stack based allocation of temporary variables and arrays
  (with safe limits).
* Fully in-place operations.
* Very fast (a few cycles) runtime construction of rings / context objects.
  (May be hundreds of cycles when constructing a new
  method table, but this is rarely needed on the fly.)
* Possibility to use generic default methods or provide optimized versions
  (e.g. for vector operations).

Disadvantages
...............................................................................

* Runtime generics in C means very little compiler protection against
  type errors; rigorous testing is crucial.
* Some function signatures become more complicated in order to provide
  uniform error handling.
* A few cycles overhead for each dynamic function lookup and
  function call. Function pointer generics also preclude intra-method
  compiler optimizations and inlining.
  Overloading vector methods should partially compensate for this.

Possible applications
...............................................................................

* At minimum, this could be useful to bootstrap new types: we only
  need to provide some basic methods (``init``, ``clear``, ``set_fmpz``,
  ``randtest``, ``add``, ``mul``, ``equal``, etc.) and the generic
  code provides fallbacks for everything else (which can be overloaded
  later for better performance). Importantly, the generic code also
  provides unit tests including tests for boring and easily borked
  things like type conversions.
* This could simplify interfacing from other libraries and languages;
  it is enough to wrap the generic interface once, without having
  to wrap all the methods for every single Flint type (with all
  their quirks and subtle differences in interface).
* There is already plenty of ad-hoc generic code in Flint
  and its descendants:

  * Function pointer generics are used, for example, in some ``mpoly`` code and in ``fmpz_mod``.

  * Template-based generics for the different ``fq`` implementations.

  * Type union + switch statement generics for ``fq_default`` and ``nf``.

  * Some Arb functions use ugly hacks to distinguish between constants and power series, etc.

  The generics module could potentially replace such code and also
  eliminate a lot of other copy-pasted boilerplate.
  By specializing methods at runtime for rings with different parameters,
  it should even be possible to improve performance in some cases.


.. raw:: latex

    \newpage
