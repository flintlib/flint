Calcium
===================================

.. only:: html

    .. image:: _static/ca2.svg
        :align: center

**Calcium** (pronounced "kalkium") is a C library for
exact computation
with real and complex numbers.
It is capable of rigorously deciding the truth of any
constant relation involving algebraic
numbers and many relations involving transcendental numbers, for example

.. math ::

    \frac{\log(\sqrt{2}+\sqrt{3})}{\log(5+2\sqrt{6})} = \frac{1}{2}, \quad
    i^{\,i} = \exp\left(\frac{\pi}{ \left({\left(\sqrt{-2}\right)}^{\sqrt{2}}\right)^{\sqrt{2}}}\right), \quad
    10^{-30} < \frac{640320^3 + 744}{e^{\pi \sqrt{163}}} - 1 < 10^{-29}.

Calcium is free software (LGPL). It depends on
`GMP <https://gmplib.org/>`_, `MPFR <https://mpfr.org/>`_,
`Flint <http://flintlib.org/>`_, `Antic <https://github.com/wbhart/antic/>`_
and `Arb <http://arblib.org/>`_.

* Source code: https://github.com/fredrik-johansson/calcium
* Bug reports and feature requests: https://github.com/fredrik-johansson/calcium/issues
* Mailing list: `flint-devel <https://groups.google.com/d/forum/flint-devel>`_

This documentation is available in HTML format at http://fredrikj.net/calcium/
and in PDF format at http://fredrikj.net/calcium/calcium.pdf.
This edition of the documentation was updated
|today| and describes Calcium |version|.

General information
---------------------------

.. toctree::
   :maxdepth: 2

   introduction.rst
   setup.rst
   examples.rst

General modules
----------------------

.. toctree::
   :maxdepth: 2

   calcium.rst

Numbers
----------------------

.. toctree::
   :maxdepth: 2

   ca.rst
   ca_vec.rst

Matrices and polynomials
------------------------

.. toctree::
   :maxdepth: 2

   ca_poly.rst
   ca_mat.rst

Field and extension number constructions
----------------------------------------

These modules are used internally by the :type:`ca_t` type
to construct towers of algebraic and transcendental number fields.
The user does not normally need to use these modules directly
outside of advanced applications requiring inspection of the
symbolic representations of numbers.

.. toctree::
   :maxdepth: 2

   ca_ext.rst
   ca_field.rst

Basic algebraic structures
--------------------------------------------

The following modules implement useful exact structures
independently of the :type:`ca_t` type. They are used
internally in Calcium, but the interfaces are stable and
use in appropriate external applications is encouraged.

.. toctree::
   :maxdepth: 2

   fmpz_mpoly_q.rst
   qqbar.rst
   utils_flint.rst


Python interface
----------------------

.. toctree::
   :maxdepth: 1

   pycalcium.rst

Credits and references
----------------------

.. toctree::
   :maxdepth: 2

   bibliography.rst

Version history
----------------------

.. toctree::
   :maxdepth: 1

   history.rst


.. only:: html

    Indices and tables
    ----------------------

    * :ref:`genindex`

