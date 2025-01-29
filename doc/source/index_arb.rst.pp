.. _index-arb:

**Real and complex numbers (Arb)** : *detailed table of contents*
=================================================================

.. only:: html

    General information
    ::::::::::::::::::::

    .. toctree::
       :maxdepth: 2

       overview.rst
       using.rst
       issues.rst

    Example programs
    ::::::::::::::::::::

    .. toctree::
       :maxdepth: 2

       examples_arb.rst

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
       acb_theta.rst
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
