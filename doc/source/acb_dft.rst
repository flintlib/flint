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

.. function:: void acb_dft(acb_ptr w, acb_srcptr v, slong len, slong prec)

   Set *w* to the DFT of *v* of length *len*, using an automatic choice
   of algorithm.

.. function:: void acb_dft_inverse(acb_ptr w, acb_srcptr v, slong len, slong prec)

   Compute the inverse DFT of *v* into *w*.

If several computations are to be done on the same group, the FFT scheme
should be reused.

The precomputed schemes build their internal tables of roots of unity at the
precision passed to their *init* function. Build a scheme at (at least) the
precision at which you intend to evaluate transforms; *acb_dft* and the
precomputed schemes then return results accurate to about *prec* bits.

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

    These functions are thin wrappers around the :doc:`gr_dft <gr_dft>`
    module, which computes DFTs over generic rings. For complex balls
    the transform is carried out in fixed-point arithmetic with
    rigorous error bounds whenever the input permits, falling back to
    ball arithmetic otherwise, and the precomputation object is the
    :type:`gr_dft_acb_pre_t` plan itself.

Obsolete functions
-------------------------------------------------------------------------------

The remaining functions of this module (product DFTs, convolutions,
and direct access to the naive, CRT, cyclic, radix-2 and Bluestein
algorithms) have been removed. Product DFTs with complex ball input
and output are provided by :func:`gr_dft_acb_prod`; the individual
algorithms, transforms over other rings and further functionality are
available in the :doc:`gr_dft <gr_dft>` module, and cyclic
convolutions are easily expressed through forward and inverse
transforms.

