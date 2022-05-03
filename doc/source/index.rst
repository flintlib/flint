.. FLINT documentation master file, created by
   sphinx-quickstart on Fri Nov 16 21:59:21 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

FLINT: Fast Library for Number Theory
======================================

Welcome to FLINT's documentation! FLINT is a C library for doing number theory,
maintained by William Hart and Fredrik Johansson.

* Source code on GitHub: https://github.com/wbhart/flint2
* Issue tracker: https://github.com/wbhart/flint2/issues
* Mailing list: https://groups.google.com/group/flint-devel

FLINT is free software distributed under the
GNU Lesser General Public License (LGPL), version 2.1 or later.

Introduction
------------

.. toctree::
   :maxdepth: 1

   introduction.rst
   building.rst
   bug_reporting.rst
   contributors.rst
   examples.rst
   macros.rst
   memory.rst
   portability.rst
   threading.rst

General utilities
-----------------

.. toctree::
   :maxdepth: 1

   flint.rst
   profiler.rst
   thread_pool.rst
   perm.rst
   mpoly.rst

Integers
---------------

.. toctree::
   :maxdepth: 1

   ulong_extras.rst
   fmpz.rst
   fmpz_vec.rst
   fmpz_factor.rst
   fmpz_mat.rst
   fmpz_lll.rst
   fmpz_poly.rst
   fmpz_poly_mat.rst
   fmpz_poly_factor.rst
   fmpz_mpoly.rst
   fmpz_mpoly_factor.rst
   long_extras.rst
   longlong.rst
   mpn_extras.rst
   aprcl.rst
   arith.rst
   fft.rst
   qsieve.rst

Rational numbers
----------------

.. toctree::
   :maxdepth: 1

   fmpq.rst
   fmpq_vec.rst
   fmpq_mat.rst
   fmpq_poly.rst
   fmpq_mpoly_factor.rst
   fmpq_mpoly.rst
   fmpz_poly_q.rst

Integers mod n
---------------

.. toctree::
   :maxdepth: 1

   nmod.rst
   nmod_vec.rst
   nmod_mat.rst
   nmod_poly.rst
   nmod_poly_mat.rst
   nmod_poly_factor.rst
   nmod_mpoly.rst
   nmod_mpoly_factor.rst
   fmpz_mod.rst
   fmpz_mod_vec.rst
   fmpz_mod_mat.rst
   fmpz_mod_poly.rst
   fmpz_mod_poly_factor.rst
   fmpz_mod_mpoly.rst
   fmpz_mod_mpoly_factor.rst

Finite fields
---------------

.. toctree::
   :maxdepth: 1

   fq.rst
   fq_default.rst
   fq_vec.rst
   fq_mat.rst
   fq_default_mat.rst
   fq_poly.rst
   fq_default_poly.rst
   fq_poly_factor.rst
   fq_default_poly_factor.rst
   fq_embed.rst

.. toctree::
   :maxdepth: 1

   fq_nmod.rst
   fq_nmod_vec.rst
   fq_nmod_mat.rst
   fq_nmod_poly.rst
   fq_nmod_poly_factor.rst
   fq_nmod_embed.rst
   fq_nmod_mpoly.rst
   fq_nmod_mpoly_factor.rst

.. toctree::
   :maxdepth: 1

   fq_zech.rst
   fq_zech_vec.rst
   fq_zech_mat.rst
   fq_zech_poly.rst
   fq_zech_poly_factor.rst
   fq_zech_embed.rst

p-adic numbers
---------------

.. toctree::
   :maxdepth: 1

   padic.rst
   padic_poly.rst
   padic_mat.rst
   qadic.rst

Floating-point support code
-----------------------------------

.. toctree::
   :maxdepth: 1

   double_extras.rst
   d_vec.rst
   d_mat.rst
   mpf_vec.rst
   mpf_mat.rst
   mpfr_vec.rst
   mpfr_mat.rst

C++ Interface
-----------------------------------

.. toctree::
   :maxdepth: 1

   flintxx.rst

References
----------------

.. toctree::
   :maxdepth: 1

   references.rst
