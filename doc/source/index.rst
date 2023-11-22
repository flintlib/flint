.. FLINT documentation master file, created by
   sphinx-quickstart on Fri Nov 16 21:59:21 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

FLINT: Fast Library for Number Theory
======================================

Welcome to FLINT's documentation! FLINT is a C library for doing number theory.

* Website: https://flintlib.org
* Source code on GitHub: https://github.com/flintlib/flint
* Issue tracker: https://github.com/flintlib/flint/issues
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
   contributing.rst
   contributors.rst
   examples.rst
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
   mpoly.rst
   machine_vectors.rst

Generic rings
-----------------------------------------------------------------------

.. toctree::
   :maxdepth: 1

   gr.rst
   gr_implementing.rst
   gr_domains.rst
   gr_generic.rst
   gr_special.rst
   gr_vec.rst
   gr_mat.rst
   gr_poly.rst
   gr_mpoly.rst

.. only:: not latex

   .. toctree::
      :maxdepth: 1

      index_generic.rst


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
   fft_small.rst
   qsieve.rst


.. only:: not latex
	  
   .. toctree::
      :maxdepth: 1

      index_integers.rst

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
   fmpz_mpoly_q.rst

.. only:: not latex
	  
   .. toctree::
      :maxdepth: 1
		 
      index_rationals.rst

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

.. only:: not latex

   .. toctree::
      :maxdepth: 1

      index_integers_mod.rst

Groups and other structures
---------------------------

.. toctree::
   :maxdepth: 1

   perm.rst
   qfb.rst
   dirichlet.rst
   dlog.rst
   bool_mat.rst

Number fields and algebraic numbers
-----------------------------------

.. toctree::
   :maxdepth: 1

   nf.rst
   nf_elem.rst
   fmpzi.rst
   qqbar.rst

Real and complex numbers
----------------------------------------

.. toctree::
   :maxdepth: 1

   overview.rst
   using.rst
   issues.rst
   examples_arb.rst
   mag.rst
   arf.rst
   acf.rst
   arb.rst
   acb.rst
   arb_poly.rst
   acb_poly.rst
   arb_fmpz_poly.rst
   acb_dft.rst
   arb_mat.rst
   acb_mat.rst
   acb_hypgeom.rst
   arb_hypgeom.rst
   acb_elliptic.rst
   acb_modular.rst
   acb_theta.rst
   acb_dirichlet.rst
   bernoulli.rst
   hypgeom.rst
   partitions.rst
   arb_calc.rst
   acb_calc.rst
   arb_fpwrap.rst
   fmpz_extras.rst
   formulas.rst
   constants.rst
   gamma.rst
   hurwitz.rst
   polylogarithms.rst
   hypergeometric.rst
   agm.rst

.. only:: not latex

    .. toctree::
       :maxdepth: 1

       index_arb.rst

Exact real and complex numbers
----------------------------------------

.. toctree::
   :maxdepth: 1

   introduction_calcium.rst
   examples_calcium.rst
   calcium.rst
   ca.rst
   ca_vec.rst
   ca_poly.rst
   ca_mat.rst
   ca_ext.rst
   ca_field.rst
   fexpr.rst
   fexpr_builtin.rst

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
   double_interval.rst
   d_vec.rst
   d_mat.rst
   mpf_vec.rst
   mpf_mat.rst
   mpfr_vec.rst
   mpfr_mat.rst

Interfaces
-----------------------------------

.. toctree::
   :maxdepth: 1

   python_flint.rst

References
-----------------------------------

.. toctree::
   :maxdepth: 1

   references.rst

Version history
-----------------------------------

.. toctree::
   :maxdepth: 1

   history.rst
