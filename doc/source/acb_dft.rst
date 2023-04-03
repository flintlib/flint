.. _acb-dft:

**acb_dft.h** -- Discrete Fourier transform
===================================================================================

*Warning: the interfaces in this module are experimental and may change
without notice.*

All functions support aliasing.

Let *G* be a finite abelian group, and `\chi` a character of *G*.
For any map `f:G\to\mathbb C`, the discrete fourier transform
`\hat f:\hat G\to \mathbb C` is defined by

.. math::

   \hat f(\chi) = \sum_{x\in G}\overline{\chi(x)}f(x)

Note that by the inversion formula

.. math::

   \widehat{\hat f}(\chi) = \#G\times f(\chi^{-1})

it is straightforward to recover `f` from its DFT `\hat f`.

Main DFT functions
-------------------------------------------------------------------------------

If `G=\mathbb Z/n\mathbb Z`, we compute the DFT according to the usual convention

.. math::

   w_x = \sum_{y\bmod n} v_y e^{-\frac{2i \pi}nxy}

.. function:: void acb_dft(acb_ptr w, acb_srcptr v, slong n, slong prec)

   Set *w* to the DFT of *v* of length *len*, using an automatic choice
   of algorithm.

.. function:: void acb_dft_inverse(acb_ptr w, acb_srcptr v, slong n, slong prec)

   Compute the inverse DFT of *v* into *w*.

If several computations are to be done on the same group, the FFT scheme
should be reused.

.. type:: acb_dft_pre_struct

.. type:: acb_dft_pre_t

    Stores a fast DFT scheme on :math:`\mathbb Z/n\mathbb Z`
    as a recursive decomposition into simpler DFT
    with some tables of roots of unity.

    An *acb_dft_pre_t* is defined as an array of *acb_dft_pre_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void acb_dft_precomp_init(acb_dft_pre_t pre, slong len, slong prec)

   Initializes the fast DFT scheme of length *len*, using an automatic choice of
   algorithms depending on the factorization of *len*.

   The length *len* is stored as *pre->n*.

.. function:: void acb_dft_precomp_clear(acb_dft_pre_t pre)

   Clears *pre*.

.. function:: void acb_dft_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)

   Computes the DFT of the sequence *v* into *w* by applying the precomputed scheme
   *pre*. Both *v* and *w* must have length *pre->n*.

.. function:: void acb_dft_inverse_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)

   Compute the inverse DFT of *v* into *w*.

DFT on products
-------------------------------------------------------------------------------

A finite abelian group is isomorphic to a product of cyclic components

.. math::

   G = \bigoplus_{i=1}^r \mathbb Z/n_i\mathbb Z

Characters are product of component characters and the DFT reads

.. math::

   \hat f(x_1,\dots x_r) = \sum_{y_1\dots y_r} f(y_1,\dots y_r)
   e^{-2i \pi \sum\frac{x_i y_i}{n_i}}

We assume that `f` is given by a vector of length `\prod n_i` corresponding
to a lexicographic ordering of the values `y_1,\dots y_r`, and the computation
returns the same indexing for values of `\hat f`.

.. function:: void acb_dirichlet_dft_prod(acb_ptr w, acb_srcptr v, slong * cyc, slong num, slong prec)

   Computes the DFT on the group product of *num* cyclic components of sizes *cyc*. Assume the entries
   of *v* are indexed according to lexicographic ordering of the cyclic
   components.

.. type:: acb_dft_prod_struct

.. type:: acb_dft_prod_t

    Stores a fast DFT scheme on a product of cyclic groups.

    An *acb_dft_prod_t* is defined as an array of *acb_dft_prod_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void acb_dft_prod_init(acb_dft_prod_t t, slong * cyc, slong num, slong prec)

   Stores in *t* a DFT scheme for the product of *num* cyclic components whose sizes are given in the array *cyc*.

.. function:: void acb_dft_prod_clear(acb_dft_prod_t t)

   Clears *t*.

.. function:: void acb_dirichlet_dft_prod_precomp(acb_ptr w, acb_srcptr v, const acb_dft_prod_t prod, slong prec)

   Sets *w* to the DFT of *v*. Assume the entries are lexicographically
   ordered according to the product of cyclic groups initialized in *t*.

Convolution
-------------------------------------------------------------------------------

For functions `f` and `g` on `G` we consider the convolution

.. math::

   (f \star g)(x) = \sum_{y\in G} f(x-y)g(y)

.. function:: void acb_dft_convol_naive(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)

.. function:: void acb_dft_convol_rad2(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)

.. function:: void acb_dft_convol(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)

   Sets *w* to the convolution of *f* and *g* of length *len*.

   The *naive* version simply uses the definition.

   The *rad2* version embeds the sequence into a power of 2 length and
   uses the formula

   .. math::

      \widehat{f \star g}(\chi) = \hat f(\chi)\hat g(\chi)

   to compute it using three radix 2 FFT.

   The default version uses radix 2 FFT unless *len* is a product of small
   primes where a non padded FFT is faster.

FFT algorithms
-------------------------------------------------------------------------------

Fast Fourier transform techniques allow to compute efficiently
all values `\hat f(\chi)` by reusing common computations.

Specifically, if `H\triangleleft G` is a subgroup of size `M` and index
`[G:H]=m`, then writing `f_x(h)=f(xh)` the translate of `f` by representatives
`x` of `G/H`, one has a decomposition

.. math::

   \hat f(\chi) = \sum_{x\in G/H} \overline{\chi(x)} \hat{f_x}(\chi_{H})

so that the DFT on `G` can be computed using `m` DFT  on `H` (of
appropriate translates of `f`), then `M` DFT on `G/H`, one for
each restriction `\chi_{H}`.

This decomposition can be done recursively.

Naive algorithm
...............................................................................

.. function:: void acb_dft_naive(acb_ptr w, acb_srcptr v, slong n, slong prec)

   Computes the DFT of *v* into *w*, where *v* and *w* have size *n*,
   using the naive `O(n^2)` algorithm.

.. type:: acb_dft_naive_struct

.. type:: acb_dft_naive_t

.. function:: void acb_dft_naive_init(acb_dft_naive_t t, slong len, slong prec)

.. function:: void acb_dft_naive_clear(acb_dft_naive_t t)

   Stores a table of roots of unity in *t*.
   The length *len* is stored as *t->n*.

.. function:: void acb_dft_naive_precomp(acb_ptr w, acb_srcptr v, const acb_dft_naive_t t, slong prec)

   Sets *w* to the DFT of *v* of size *t->n*, using the naive algorithm data *t*.

CRT decomposition
...............................................................................

.. function:: void acb_dft_crt(acb_ptr w, acb_srcptr v, slong n, slong prec)

   Computes the DFT of *v* into *w*, where *v* and *w* have size *len*,
   using CRT to express `\mathbb Z/n\mathbb Z` as a product of cyclic groups.

.. type:: acb_dft_crt_struct

.. type:: acb_dft_crt_t

.. function:: void acb_dft_crt_init(acb_dft_crt_t t, slong len, slong prec)

.. function:: void acb_dft_crt_clear(acb_dft_crt_t t)

   Initialize a CRT decomposition of `\mathbb Z/n\mathbb Z` as a direct product
   of cyclic groups.
   The length *len* is stored as *t->n*.

.. function:: void acb_dft_crt_precomp(acb_ptr w, acb_srcptr v, const acb_dft_crt_t t, slong prec)

   Sets *w* to the DFT of *v* of size *t->n*, using the CRT decomposition scheme *t*.

Cooley-Tukey decomposition
...............................................................................

.. function:: void acb_dft_cyc(acb_ptr w, acb_srcptr v, slong n, slong prec)

   Computes the DFT of *v* into *w*, where *v* and *w* have size *n*,
   using each prime factor of `m` of `n` to decompose with
   the subgroup `H=m\mathbb Z/n\mathbb Z`.

.. type:: acb_dft_cyc_struct

.. type:: acb_dft_cyc_t

.. function:: void acb_dft_cyc_init(acb_dft_cyc_t t, slong len, slong prec)

.. function:: void acb_dft_cyc_clear(acb_dft_cyc_t t)

   Initialize a decomposition of `\mathbb Z/n\mathbb Z` into cyclic subgroups.
   The length *len* is stored as *t->n*.

.. function:: void acb_dft_cyc_precomp(acb_ptr w, acb_srcptr v, const acb_dft_cyc_t t, slong prec)

   Sets *w* to the DFT of *v* of size *t->n*, using the cyclic decomposition scheme *t*.

Radix 2 decomposition
...............................................................................

.. function:: void acb_dft_rad2(acb_ptr w, acb_srcptr v, int e, slong prec)

   Computes the DFT of *v* into *w*, where *v* and *w* have size `2^e`,
   using a radix 2 FFT.

.. function:: void acb_dft_inverse_rad2(acb_ptr w, acb_srcptr v, int e, slong prec)

   Computes the inverse DFT of *v* into *w*, where *v* and *w* have size `2^e`,
   using a radix 2 FFT.

.. type:: acb_dft_rad2_struct

.. type:: acb_dft_rad2_t

.. function:: void acb_dft_rad2_init(acb_dft_rad2_t t, int e, slong prec)

.. function:: void acb_dft_rad2_clear(acb_dft_rad2_t t)

   Initialize and clear a radix 2 FFT of size `2^e`, stored as *t->n*.

.. function:: void acb_dft_rad2_precomp(acb_ptr w, acb_srcptr v, const acb_dft_rad2_t t, slong prec)

   Sets *w* to the DFT of *v* of size *t->n*, using the precomputed radix 2 scheme *t*.

Bluestein transform
...............................................................................

.. function:: void acb_dft_bluestein(acb_ptr w, acb_srcptr v, slong n, slong prec)

   Computes the DFT of *v* into *w*, where *v* and *w* have size *n*,
   by conversion to a radix 2 one using Bluestein's convolution trick.

.. type:: acb_dft_bluestein_struct

.. type:: acb_dft_bluestein_t

   Stores a Bluestein scheme for some length *n* : that is a :type:`acb_dft_rad2_t` of size
   `2^e \geq 2n-1` and a size *n* array of convolution factors.

.. function:: void acb_dft_bluestein_init(acb_dft_bluestein_t t, slong len, slong prec)
 
.. function:: void acb_dft_bluestein_clear(acb_dft_bluestein_t t)

   Initialize and clear a Bluestein scheme to compute DFT of size *len*.

.. function:: void acb_dft_bluestein_precomp(acb_ptr w, acb_srcptr v, const acb_dft_bluestein_t t, slong prec)

   Sets *w* to the DFT of *v* of size *t->n*, using the precomputed Bluestein scheme *t*.

