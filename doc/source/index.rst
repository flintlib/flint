Arb - a C library for arbitrary-precision ball arithmetic
=============================================================

.. only:: latex

    Introduction
    ::::::::::::

    Welcome to Arb's documentation!
    Arb is a C library for rigorous real and complex arithmetic with arbitrary precision.
    Arb tracks numerical errors automatically using
    *ball arithmetic*, a form of interval arithmetic based on a midpoint-radius
    representation.
    On top of this, Arb provides a wide range of mathematical functionality, including polynomials,
    power series, matrices, integration, root-finding, and transcendental functions.
    Arb is designed with efficiency as a primary goal, and is usually competitive with or faster
    than other arbitrary-precision packages.
    The code is thread-safe, portable, and extensively tested.

    Arb is free software distributed under the
    GNU Lesser General Public License (LGPL), version 2.1 or later
    (see :ref:`license`).

    The git repository is https://github.com/fredrik-johansson/arb/

    Arb is developed by `Fredrik Johansson <http://fredrikj.net/>`_
    (fredrik.johansson@gmail.com), with help from many
    contributors (see :ref:`credits`).
    Questions and discussion about Arb are welcome on the
    `flint-devel <https://groups.google.com/d/forum/flint-devel>`_ mailing list.
    There is also an `issue tracker <https://github.com/fredrik-johansson/arb/issues>`_
    for bug reports and feature requests.
    Development progress is sometimes covered on
    `Fredrik's blog <http://fredrikj.net/blog/>`_.

    This documentation is available in HTML format at http://arblib.org and in
    PDF format at http://arblib.org/arb.pdf.
    This edition of the documentation was updated
    |today| and describes Arb |version|.
    Documentation for :ref:`specific release versions <history>`
    is also available in PDF format.

.. only:: html

        .. image:: _static/banner.jpg
            :align: center

    Welcome to Arb's documentation!
    Arb is a C library for rigorous real and complex arithmetic with arbitrary precision.
    Arb tracks numerical errors automatically using
    *ball arithmetic*, a form of interval arithmetic based on a midpoint-radius
    representation.
    On top of this, Arb provides a wide range of mathematical functionality, including polynomials,
    power series, matrices, integration, root-finding, and many transcendental functions.
    Arb is designed with efficiency as a primary goal, and is usually competitive with or faster
    than other arbitrary-precision packages.
    The code is thread-safe, portable, and extensively tested.

    Arb is free software distributed under the
    GNU Lesser General Public License (LGPL), version 2.1 or later
    (see :ref:`license`).

    The git repository is https://github.com/fredrik-johansson/arb/

    Arb is developed by `Fredrik Johansson <http://fredrikj.net/>`_
    (fredrik.johansson@gmail.com), with help from many
    contributors (see :ref:`credits`).
    Questions and discussion about Arb are welcome on the
    `flint-devel <https://groups.google.com/d/forum/flint-devel>`_ mailing list.
    There is also an `issue tracker <https://github.com/fredrik-johansson/arb/issues>`_
    for bug reports and feature requests.
    Development progress is sometimes covered on
    `Fredrik's blog <http://fredrikj.net/blog/>`_.

    This documentation is available in HTML format at http://arblib.org and in
    PDF format at http://arblib.org/arb.pdf.
    This edition of the documentation was updated
    |today| and describes Arb |version|.
    Documentation for :ref:`specific release versions <history>`
    is also available in PDF format.

General information
::::::::::::::::::::

.. toctree::
   :maxdepth: 2

   overview.rst
   setup.rst
   using.rst
   issues.rst
   contributing.rst
   credits.rst

Example programs
::::::::::::::::::::

.. toctree::
   :maxdepth: 2

   examples.rst

Floating-point numbers
::::::::::::::::::::::::::::::::::::

Arb uses two custom floating-point types in its implementation of ball
arithmetic. The radius of a ball is represented using the type *mag_t* which is
unsigned and has a fixed precision. The midpoint is represented using the
type *arf_t* which has arbitrary precision.

.. toctree::
   :maxdepth: 2

   mag.rst
   arf.rst

Real and complex numbers
::::::::::::::::::::::::::::::::::::

Real numbers (*arb_t*) are represented as midpoint-radius intervals,
also known as balls. Complex numbers (*acb_t*) are represented in rectangular
form, with *arb_t* balls for the real and imaginary parts.

.. toctree::
   :maxdepth: 2

   arb.rst
   acb.rst

Polynomials and power series
::::::::::::::::::::::::::::::::::::

These modules implement dense univariate polynomials with real and complex
coefficients. Truncated power series are supported via methods acting
on polynomials, without introducing a separate power series type.

.. toctree::
   :maxdepth: 2

   arb_poly.rst
   acb_poly.rst

.. toctree::
   :maxdepth: 2

   arb_fmpz_poly.rst

Transforms
::::::::::::::::::::::::::::::::::::

.. toctree::
   :maxdepth: 2

   acb_dft.rst

Matrices
::::::::::::::::::::::::::::::::::::

These modules implement dense matrices with real and complex coefficients.
Rudimentary linear algebra is supported.

.. toctree::
   :maxdepth: 2

   arb_mat.rst
   acb_mat.rst

Special functions
::::::::::::::::::::::::::::::::::::

These modules implement mathematical functions with complexity
that goes beyond the basics covered directly in the *arb* and *acb*
modules.

.. toctree::
   :maxdepth: 2

   acb_hypgeom.rst
   arb_hypgeom.rst
   acb_elliptic.rst
   acb_modular.rst
   dirichlet.rst
   acb_dirichlet.rst
   bernoulli.rst
   hypgeom.rst
   partitions.rst

Calculus
::::::::::::::::::::::::::::::::::::

Using ball arithmetic, it is possible to do rigorous root-finding and
integration (among other operations)
with generic functions. This code should be considered experimental.

.. toctree::
   :maxdepth: 2

   arb_calc.rst
   acb_calc.rst

Wrappers
::::::::::::::::::::::::::::::::::::

Floating-point wrappers for Arb functions.

.. toctree::
   :maxdepth: 2

   arb_fpwrap.rst


Extra utility modules
::::::::::::::::::::::::::::::::::::

Mainly for internal use.

.. toctree::
   :maxdepth: 1

   fmpzi.rst
   acf.rst
   double_interval.rst
   fmpz_extras.rst
   bool_mat.rst
   dlog.rst
   fmpr.rst

Supplementary algorithm notes
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Here, we give extra proofs, error bounds, and formulas that would be too
lengthy to reproduce in the documentation for each module.

.. toctree::
   :maxdepth: 1

   formulas.rst
   constants.rst
   gamma.rst
   hurwitz.rst
   polylogarithms.rst
   hypergeometric.rst
   agm.rst

Version history
:::::::::::::::::::::::::::::::::

.. toctree::
   :maxdepth: 1

   history.rst

