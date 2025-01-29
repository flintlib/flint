.. _overview:

Feature overview
===============================================================================

Ball arithmetic, also known as mid-rad interval arithmetic, is an
extension of floating-point arithmetic in which an error bound is
attached to each variable. This allows computing with real and
complex numbers in a mathematically rigorous way.

With plain floating-point arithmetic, the user must do an error analysis
to guarantee that results are correct. Manual error analysis is time-consuming
and bug-prone. Ball arithmetic effectively makes error analysis
automatic.

In traditional (inf-sup) interval arithmetic, both endpoints of an interval
`[a,b]` are full-precision numbers, which makes interval arithmetic
twice as expensive as floating-point arithmetic.
In ball arithmetic, only the midpoint *m* of an interval `[m \pm r]`
is a full-precision number, and a few bits suffice for the radius *r*.
At high precision, ball arithmetic is therefore not more expensive than
plain floating-point arithmetic.

Joris van der Hoeven's paper [Hoe2009]_ is a good introduction to the subject.

Other implementations of ball arithmetic include
`iRRAM <http://irram.uni-trier.de/>`_ and
`Mathemagix <http://www.mathemagix.org/www/mmdoc/doc/html/main/index.en.html>`_.
Arb differs from earlier implementations in technical aspects of the
implementation, which makes certain computations more efficient.
It also provides a more comprehensive low-level interface, giving
the user full access to the internals. Finally, it implements a wider
range of transcendental functions, covering a large portion of the
special functions in standard reference works such as [NIST2012]_.

The ball arithmetic routines in FLINT (formerly the standalone Arb
library) are designed for computer algebra and computational number
theory, but may be useful in any area demanding
reliable or precise numerical computing.
The contents include:

* A module (:ref:`arf <arf>`) for correctly rounded arbitrary-precision
  floating-point arithmetic. Arb's floating-point numbers have a few special
  features, such as arbitrary-size exponents (useful for combinatorics and
  asymptotics) and dynamic allocation (facilitating implementation of hybrid
  integer/floating-point and mixed-precision algorithms).

* A module (:ref:`mag <mag>`) for representing magnitudes (error bounds)
  more efficiently than with an arbitrary-precision floating-point type.

* A module (:ref:`arb <arb>`) for real ball arithmetic, where a ball is
  implemented as an *arf* midpoint and a *mag* radius.

* A module (:ref:`acb <acb>`) for complex numbers in rectangular form,
  represented as pairs of real balls.

* Modules (:ref:`arb_poly <arb-poly>`, :ref:`acb_poly <acb-poly>`)
  for polynomials or power series over the real and complex numbers,
  implemented using balls as coefficients,
  with asymptotically fast polynomial multiplication and
  many other operations.

* Modules (:ref:`arb_mat <arb-mat>`, :ref:`acb_mat <acb-mat>`)
  for matrices over the real and complex numbers,
  implemented using balls as coefficients.
  At the moment, only rudimentary linear algebra operations are provided.

* Functions for high-precision evaluation of various
  mathematical constants and special functions, implemented using
  ball arithmetic with rigorous error bounds.
