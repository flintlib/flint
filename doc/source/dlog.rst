.. _dlog:

**dlog.h** -- discrete logarithms mod ulong primes
===============================================================================

This module implements discrete logarithms, with the application
to Dirichlet characters in mind.

In particular, this module defines a :type:`dlog_precomp_t` structure
permitting to describe a discrete log problem  in some subgroup
of `(\mathbb Z/p^e \mathbb Z)^\times` for primepower moduli `p^e`,
and store precomputed data for faster computation of several such
discrete logarithms.

When initializing this data, the user provides both a group description and the expected
number of subsequent discrete logarithms calls. The choice of algorithm and
the amount of stored data depend both on the structure of the group and this number.

No particular effort has been made towards single discrete logarithm
computation. Currently only machine size primepower moduli
are supported.

Types, macros and constants
-------------------------------------------------------------------------------

.. macro:: DLOG_NONE

   Return value when the discrete logarithm does not exist

.. type:: dlog_precomp_struct

.. type:: dlog_precomp_t

   Structure for discrete logarithm precomputed data.

   A :type:`dlog_precomp_t` is defined as an array of length one of type
   :type:`dlog_precomp_struct`, permitting a :type:`dlog_precomp_t` to be passed by
   reference.

Single evaluation
-------------------------------------------------------------------------------

.. function:: ulong dlog_once(ulong b, ulong a, const nmod_t mod, ulong n)

   Return `x` such that `b = a^x` in `(\mathbb Z/mod \mathbb Z)^\times`,
   where *a* is known to have order *n*.

Precomputations
-------------------------------------------------------------------------------

.. function:: void dlog_precomp_n_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num)

   Precompute data for *num* discrete logarithms evaluations in the subgroup generated
   by *a* modulo *mod*, where *a* is known to have order *n*.

.. function:: ulong dlog_precomp(const dlog_precomp_t pre, ulong b)

   Return `\log(b)` for the group described in *pre*.

.. function:: void dlog_precomp_clear(dlog_precomp_t pre)

   Clears *t*.

Specialized versions of :func:`dlog_precomp_n_init` are available when specific information
is known about the group:

.. function:: void dlog_precomp_modpe_init(dlog_precomp_t pre, ulong a, ulong p, ulong e, ulong pe, ulong num)

   Assume that *a* generates the group of residues modulo *pe* equal `p^e` for
   prime *p*.

.. function:: void dlog_precomp_p_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong num)

   Assume that *a* has prime order *p*.

.. function:: void dlog_precomp_pe_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong e, ulong pe, ulong num)

   Assume that *a* has primepower order *pe* `p^e`.

.. function:: void dlog_precomp_small_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num)

   Make a complete lookup table of size *n*.
   If *mod* is small, this is done using an element-indexed array (see
   :type:`dlog_table_t`), otherwise with a sorted array allowing binary search.

Vector evaluations
-------------------------------------------------------------------------------

These functions compute all logarithms of successive integers `1\dots n`.

.. function:: void dlog_vec_fill(ulong * v, ulong nv, ulong x)

   Sets values *v[k]* to *x* for all *k* less than *nv*.

.. function:: void dlog_vec_set_not_found(ulong * v, ulong nv, nmod_t mod)

   Sets values *v[k]* to :macro:`DLOG_NONE` for all *k* not coprime to *mod*.

.. function:: void dlog_vec(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

   Sets *v[k]* to `\log(k,a)` times value *va*  for `0\leq k < nv`, where *a*
   has order *na*. *va* should be *1* for usual log computation.

.. function:: void dlog_vec_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

   Same parameters as before, but adds `\log(k,a)\times v_a`
   to *v[k]* and reduce modulo *order* instead of replacing the value. Indices
   *k* such that *v[k]* equals *DLOG_NONE* are ignored.

Depending on the relative size of *nv* and *na*, these two *dlog_vec* functions
call one of the following functions.

.. function:: void dlog_vec_loop(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

.. function:: void dlog_vec_loop_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

   Perform a complete loop of size *na* on powers of *a* to fill the logarithm
   values, discarding powers outside the bounds of *v*. This requires no
   discrete logarithm computation.

.. function:: void dlog_vec_eratos(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

.. function:: void dlog_vec_eratos_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

   Compute discrete logarithms of prime numbers less than *nv* and propagate to composite numbers.

.. function:: void dlog_vec_sieve_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

.. function:: void dlog_vec_sieve(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

   Compute the discrete logarithms of the first few prime numbers, then
   use them as a factor base to obtain the logarithms of larger primes
   by sieving techniques.

   In the the present implementation, the full index-calculus method is not
   implemented.

Internal discrete logarithm strategies
-------------------------------------------------------------------------------

Several discrete logarithms strategies are implemented:

- Complete lookup table for small groups.

- Baby-step giant-step table.

combined with mathematical reductions:

- Pohlig-Hellman decomposition (Chinese remainder decomposition on the
  order of the group and base `p` decomposition for primepower order).

- p-adic log for primepower modulus `p^e`.

The *dlog_precomp* structure makes recursive use of the following
method-specific structures.

Complete table
...............................................................................

.. type:: dlog_table_struct

.. type:: dlog_table_t

   Structure for complete lookup table.

.. function:: ulong dlog_table_init(dlog_table_t t, ulong a, ulong mod)

   Initialize a table of powers of *a* modulo *mod*, storing all elements
   in an array of size *mod*.

.. function:: void dlog_table_clear(dlog_table_t t)

   Clears *t*.

.. function:: ulong dlog_table(dlog_table_t t, ulong b)

   Return `\log(b,a)` using the precomputed data *t*.

Baby-step giant-step table
...............................................................................

.. type:: dlog_bsgs_struct

.. type:: dlog_bsgs_t

   Structure for Baby-Step Giant-Step decomposition.

.. function:: ulong dlog_bsgs_init(dlog_bsgs_t t, ulong a, ulong mod, ulong n, ulong m)

   Initialize *t* and store the first *m* powers of *a* in a sorted array. The
   return value is a rought measure of the cost of each logarithm using this
   table.
   The user should take `m\approx\sqrt{kn}` to compute k logarithms in a group of size n.

.. function:: void dlog_bsgs_clear(dlog_bsgs_t t)

   Clears *t*.

.. function:: ulong dlog_bsgs(dlog_bsgs_t t, ulong b)

   Return `\log(b,a)` using the precomputed data *t*.

Prime-power modulus decomposition
...............................................................................

.. type:: dlog_modpe_struct

.. type:: dlog_modpe_t

   Structure for discrete logarithm modulo primepower `p^e`.

   A :type:`dlog_modpe_t` is defined as an array of length one of type
   :type:`dlog_modpe_struct`, permitting a :type:`dlog_modpe_t` to be passed by
   reference.

.. function:: ulong dlog_modpe_init(dlog_modpe_t t, ulong a, ulong p, ulong e, ulong pe, ulong num)

.. function:: void dlog_modpe_clear(dlog_modpe_t t)

   Clears *t*.

.. function:: ulong dlog_modpe(dlog_modpe_t t, ulong b)

   Return `\log(b,a)` using the precomputed data *t*.

CRT decomposition
...............................................................................

.. type:: dlog_crt_struct

.. type:: dlog_crt_t

   Structure for discrete logarithm for groups of composite order.
   A :type:`dlog_crt_t` is defined as an array of length one of type
   :type:`dlog_crt_struct`, permitting a :type:`dlog_crt_t` to be passed by
   reference.

.. function:: ulong dlog_crt_init(dlog_crt_t t, ulong a, ulong mod, ulong n, ulong num)

   Precompute data for *num* evaluations of discrete logarithms in base *a* modulo *mod*,
   where *a* has composite order *n*, using chinese remainder decomposition.

.. function:: void dlog_crt_clear(dlog_crt_t t)

   Clears *t*.

.. function:: ulong dlog_crt(dlog_crt_t t, ulong b)

   Return `\log(b,a)` using the precomputed data *t*.

padic decomposition
...............................................................................

.. type:: dlog_power_struct

.. type:: dlog_power_t

   Structure for discrete logarithm for groups of primepower order.
   A :type:`dlog_power_t` is defined as an array of length one of type
   :type:`dlog_power_struct`, permitting a :type:`dlog_power_t` to be passed by
   reference.

.. function:: ulong dlog_power_init(dlog_power_t t, ulong a, ulong mod, ulong p, ulong e, ulong num)

   Precompute data for *num* evaluations of discrete logarithms in base *a* modulo *mod*,
   where *a* has prime power order *pe* equals `p^e`, using decomposition in base *p*.

.. function:: void dlog_power_clear(dlog_power_t t)

   Clears *t*.

.. function:: ulong dlog_power(dlog_power_t t, ulong b)

   Return `\log(b,a)` using the precomputed data *t*.

Pollard rho method
...............................................................................

.. type:: dlog_rho_struct

.. type:: dlog_rho_t

   Structure for discrete logarithm using Pollard rho.
   A :type:`dlog_rho_t` is defined as an array of length one of type
   :type:`dlog_rho_struct`, permitting a :type:`dlog_rho_t` to be passed by
   reference.

.. function:: ulong dlog_rho_init(dlog_rho_t t, ulong a, ulong mod, ulong n, ulong num)

   Initialize random walks for evaluations of discrete logarithms in base *a* modulo *mod*,
   where *a* has order *n*.

.. function:: void dlog_rho_clear(dlog_rho_t t)

   Clears *t*.

.. function:: ulong dlog_rho(dlog_rho_t t, ulong b)

   Return `\log(b,a)` by the rho method in the group described by *t*.
