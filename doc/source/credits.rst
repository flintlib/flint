.. _credits:

Credits and references
===============================================================================

.. _license:

License
-------------------------------------------------------------------------------

Arb is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

Arb is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with Arb (see the LICENSE file in the root of the Arb source
directory).  If not, see http://www.gnu.org/licenses/.

Versions of Arb up to and including 2.8 were distributed under
the GNU General Public License (GPL), not the LGPL. The switch to
the LGPL applies retroactively; i.e. users may redistribute older versions
of Arb under the LGPL if they prefer.

Authors
-------------------------------------------------------------------------------

Fredrik Johansson is the main author. The project was started in 2012
as a numerical extension of FLINT, and the initial design was heavily based
on FLINT 2.0 (with particular credit to Bill Hart and Sebastian Pancratz).

The following authors have developed major new features.

* Pascal Molin - discrete Fourier transform (DFT), Dirichlet characters, Dirichlet L-functions, discrete logarithm computation
* Alex Griffing - sinc function, matrix trace, improved matrix squaring, boolean matrices, improved structured matrix exponentials, Cholesky decomposition, miscellaneous patches
* Marc Mezzarobba - fast evaluation of Legendre polynomials, work on Arb interface in Sage, bug reports, feedback

Several people have contributed patches, bug reports, or substantial feedback.
This list (ordered by time of first contribution) is probably incomplete.

* Bill Hart - build system, Windows 64 support, design of FLINT
* Sebastian Pancratz - divide-and-conquer polynomial composition algorithm (taken from FLINT)
* The MPFR development team - Arb includes two-limb multiplication code taken from MPFR
* Jonathan Bober - original code for Dirichlet characters, C++ compatibility fixes
* Yuri Matiyasevich - feedback about the zeta function and root-finding code
* Abhinav Baid - dot product and norm functions
* Ondřej Čertík - bug reports, feedback
* Andrew Booker - bug reports, feedback
* Francesco Biscani - C++ compatibility fixes, feedback
* Clemens Heuberger - work on Arb interface in Sage, feedback
* Ricky Farr - convenience functions, feedback
* Marcello Seri - fix for static builds on OS X
* Tommy Hofmann - matrix transpose, comparison, other utility methods, Julia interface
* Alexander Kobel - documentation and code cleanup patches
* Hrvoje Abraham - patches for MinGW compatibility
* Julien Puydt - soname versioning support, bug reports, Debian testing
* Jeroen Demeyer - patch for major bug on PPC64
* Isuru Fernando - continuous integration setup, support for cmake and MSVC builds
* François Bissey - build system patches
* Jean-Pierre Flori - code simplifications for Gauss periods, feedback
* arbguest - preconditioned linear algebra algorithms
* Ralf Stephan - return exact real parts in acos and acosh
* Vincent Delecroix - various feedback and patches, work on Sage interface
* D.H.J Polymath - Riemann xi function, Riemann zeta zeros
* Joel Dahne - feedback and improvements for Legendre functions
* Gianfranco Costamagna - bug reports, Debian testing
* Julian Rüth - serialization support
* Michael Orlitzky - build system patches
* David Berghaus - aliased window matrix multiplication
* Albin Ahlbäck - uniformly distributed random numbers
* Daniel Schultz - derivative of Weierstrass elliptic function
* Matthias Gessinger - Graeffe transforms
* David Harvey - modular computation of Bernoulli numbers (code taken from Harvey's bernmm package)
* Erik Postma - improved handling of infinities

Funding
-------------------------------------------------------------------------------

From 2012 to July 2014, Fredrik's work on Arb was supported by
Austrian Science Fund FWF Grant Y464-N18 (Fast Computer Algebra
for Special Functions).
During that period, he was a PhD student (and briefly a postdoc) at
RISC, Johannes Kepler University, Linz, supervised by Manuel Kauers.

From September 2014 to October 2015, Fredrik was a postdoc at
INRIA Bordeaux and Institut de Mathématiques de Bordeaux,
in the LFANT project-team headed by Andreas Enge. During that period,
Fredrik's work on Arb was supported
by ERC Starting Grant ANTICS 278537 (Algorithmic Number Theory in
Computer Science) http://cordis.europa.eu/project/rcn/101288_en.html
Since October 2015, Fredrik is a CR2 researcher in the LFANT team,
funded by INRIA.

Software
-------------------------------------------------------------------------------

The following software has been helpful in the development of Arb.

* GMP (Torbjörn Granlund and others), http://gmplib.org
* MPIR (Brian Gladman, Jason Moxham, William Hart and others), http://mpir.org
* MPFR (Guillaume Hanrot, Vincent Lefèvre, Patrick Pélissier, Philippe Théveny, Paul Zimmermann and others), http://mpfr.org
* FLINT (William Hart, Sebastian Pancratz, Andy Novocin, Fredrik Johansson, David Harvey and others), http://flintlib.org
* Sage (William Stein and others), http://sagemath.org
* Pari/GP (The Pari group), http://pari.math.u-bordeaux.fr/
* SymPy (Ondřej Čertík, Aaron Meurer and others), http://sympy.org
* mpmath (Fredrik Johansson and others), http://mpmath.org
* Mathematica (Wolfram Research), http://www.wolfram.com/mathematica
* HolonomicFunctions (Christoph Koutschan), http://www.risc.jku.at/research/combinat/software/HolonomicFunctions/
* Sphinx (George Brandl and others), http://sphinx.pocoo.org
* CM (Andreas Enge), http://www.multiprecision.org/index.php?prog=cm
* ore_algebra (Manuel Kauers, Maximilian Jaroschek, Fredrik Johansson), http://www.risc.jku.at/research/combinat/software/ore_algebra/

Citing Arb
-------------------------------------------------------------------------------

To cite Arb in a scientific paper, the following reference can be used:

\F. Johansson. "Arb: efficient arbitrary-precision midpoint-radius interval arithmetic", *IEEE Transactions on Computers*, 66(8):1281-1292, 2017. DOI: `10.1109/TC.2017.2690633 <https://doi.org/10.1109/TC.2017.2690633>`_.

In BibTeX format::

  @article{Johansson2017arb,
    author = {F. Johansson},
    title = {Arb: efficient arbitrary-precision midpoint-radius interval arithmetic},
    journal = {IEEE Transactions on Computers},
    year = {2017},
    volume = {66},
    issue = {8},
    pages = {1281--1292},
    doi = {10.1109/TC.2017.2690633},
  }

Alternatively, the Arb manual or website can be cited directly.

The *IEEE Transactions on Computers* paper supersedes the following extended abstract,
which is now outdated:

\F. Johansson. "Arb: a C library for ball arithmetic", *ACM Communications in Computer Algebra*, 47(4):166-169, 2013.

