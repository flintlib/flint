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

Bibliography
-------------------------------------------------------------------------------

(In the PDF edition, this section is empty. See the bibliography listing at the end of the document.)

.. [Ari2011] \J. Arias de Reyna, "High precision computation of Riemann’s zeta function by the Riemann-Siegel formula, I", Mathematics of Computation 80 (2011), 995-1009

.. [Ari2012] \J. Arias de Reyna, "Programs for Riemann's zeta function", (J. A. J. van Vonderen, Ed.) *Leven met getallen : liber amicorum ter gelegenheid van de pensionering van Herman te Riele* CWI (2012) 102-112, https://ir.cwi.nl/pub/19724

.. [Arn2010] \J. Arndt, *Matters Computational*, Springer (2010), http://www.jjj.de/fxt/#fxtbook

.. [BBC1997] \D. H. Bailey, J. M. Borwein and R. E. Crandall, "On the Khintchine constant", Mathematics of Computation 66 (1997) 417-431

.. [Blo2009] \R. Bloemen, "Even faster zeta(2n) calculation!", https://web.archive.org/web/20141101133659/http://xn--2-umb.com/09/11/even-faster-zeta-calculation

.. [BBC2000] \J. Borwein, D. M. Bradley and R. E. Crandall, "Computational strategies for the Riemann zeta function", Journal of Computational and Applied Mathematics 121 (2000) 247-296

.. [BZ1992] \J. Borwein and I. Zucker, "Fast evaluation of the gamma function for small rational fractions using complete elliptic integrals of the first kind", IMA Journal of Numerical Analysis 12 (1992) 519-526

.. [Bog2012] \I. Bogaert, B. Michiels and J. Fostier, "O(1) computation of Legendre polynomials and Gauss-Legendre nodes and weights for parallel computing", SIAM Journal on Scientific Computing 34:3 (2012), C83-C101

.. [Bor1987] \P. Borwein, "Reduced complexity evaluation of hypergeometric functions", Journal of Approximation Theory 50:3 (1987)

.. [Bor2000] \P. Borwein, "An Efficient Algorithm for the Riemann Zeta Function", Constructive experimental and nonlinear analysis, CMS Conference Proc. 27 (2000) 29-34, http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf

.. [BM1980] \R. P. Brent and E. M. McMillan, "Some new algorithms for high-precision computation of Euler's constant", Mathematics of Computation 34 (1980) 305-312.

.. [Bre1978] \R. P. Brent, "A Fortran multiple-precision arithmetic package", ACM Transactions on Mathematical Software, 4(1):57–70, March 1978.

.. [Bre1979] \R. P. Brent, "On the Zeros of the Riemann Zeta Function in the Critical Strip", Mathematics of Computation 33 (1979), 1361-1372, https://doi.org/10.1090/S0025-5718-1979-0537983-2

.. [Bre2010] \R. P. Brent, "Ramanujan and Euler's Constant", http://wwwmaths.anu.edu.au/~brent/pd/Euler_CARMA_10.pdf

.. [BJ2013] \R. P. Brent and F. Johansson, "A bound for the error term in the Brent-McMillan algorithm", preprint (2013), http://arxiv.org/abs/1312.0039

.. [BZ2011] \R. P. Brent and P. Zimmermann, *Modern Computer Arithmetic*, Cambridge University Press (2011), http://www.loria.fr/~zimmerma/mca/pub226.html

.. [Car1995] \B. C. Carlson, "Numerical computation of real or complex elliptic integrals". Numerical Algorithms, 10(1):13-26 (1995).

.. [CP2005] \R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective*, second edition, Springer (2005).

.. [CGHJK1996] \R. M. Corless, G. H. Gonnet, D. E. Hare, D. J. Jeffrey and D. E. Knuth, "On the Lambert W function", Advances in Computational Mathematics, 5(1) (1996), 329-359

.. [Dup2006] \R. Dupont. "Moyenne arithmético-géométrique, suites de Borchardt et applications." These de doctorat, École polytechnique, Palaiseau (2006). http://http://www.lix.polytechnique.fr/Labo/Regis.Dupont/these_soutenance.pdf

.. [DYF1999] \A. Dzieciol, S. Yngve and P. O. Fröman, "Coulomb wave functions with complex values of the variable and the parameters", J. Math. Phys. 40, 6145 (1999), https://doi.org/10.1063/1.533083

.. [EHJ2016] \A. Enge, W. Hart and F. Johansson, "Short addition sequences for theta functions", preprint (2016), https://arxiv.org/abs/1608.06810

.. [EM2004] \O. Espinosa and V. Moll, "A generalized polygamma function", Integral Transforms and Special Functions (2004), 101-115.

.. [Fil1992] \S. Fillebrown, "Faster Computation of Bernoulli Numbers", Journal of Algorithms 13 (1992) 431-445

.. [Gas2018]  \D. Gaspard, "Connection formulas between Coulomb wave functions" (2018), https://arxiv.org/abs/1804.10976

.. [GG2003] \J. von zur Gathen and J. Gerhard, *Modern Computer Algebra*, second edition, Cambridge University Press (2003)

.. [GVL1996] \G. H. Golub and C. F. Van Loan, *Matrix Computations*, third edition, Johns Hopkins University Press (1996).

.. [GS2003] \X. Gourdon and P. Sebah, "Numerical evaluation of the Riemann Zeta-function" (2003), http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.pdf

.. [HS1967] \E. Hansen and R. Smith, "Interval Arithmetic in Matrix Computations, Part II", SIAM Journal of Numerical Analysis, 4(1):1-9 (1967). https://doi.org/10.1137/0704001

.. [HZ2004] \G. Hanrot and P. Zimmermann, "Newton Iteration Revisited" (2004), http://www.loria.fr/~zimmerma/papers/fastnewton.ps.gz

.. [Har2010] \D. Harvey,  "A multimodular algorithm for computing Bernoulli numbers" (2010), Mathematics of Computation 79.272: 2361-2370

.. [Hoe2009] \J. van der Hoeven, "Ball arithmetic", Technical Report, HAL 00432152 (2009), http://www.texmacs.org/joris/ball/ball-abs.html

.. [Hoe2001] \J. van der Hoeven. "Fast evaluation of holonomic functions near and in regular singularities", Journal of Symbolic Computation, 31(6):717-743 (2001).

.. [HM2017] \J. van der Hoeven and B. Mourrain. "Efficient certification of numeric solutions to eigenproblems", MACIS 2017, 81-94, (2017), https://hal.archives-ouvertes.fr/hal-01579079

.. [JB2018] \F. Johansson and I. Blagouchine. "Computing Stieltjes constants using complex integration", preprint (2018), https://arxiv.org/abs/1804.01679

.. [Joh2012] \F. Johansson, "Efficient implementation of the Hardy-Ramanujan-Rademacher formula", LMS Journal of Computation and Mathematics, Volume 15 (2012), 341-359, http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=8710297

.. [Joh2013] \F. Johansson, "Rigorous high-precision computation of the Hurwitz zeta function and its derivatives", Numerical Algorithms, http://arxiv.org/abs/1309.2877 http://dx.doi.org/10.1007/s11075-014-9893-1

.. [Joh2014a] \F. Johansson, *Fast and rigorous computation of special functions to high precision*, PhD thesis, RISC, Johannes Kepler University, Linz, 2014. http://fredrikj.net/thesis/

.. [Joh2014b] \F. Johansson, "Evaluating parametric holonomic sequences using rectangular splitting", ISSAC 2014, 256-263. http://dx.doi.org/10.1145/2608628.2608629

.. [Joh2014c] \F. Johansson, "Efficient implementation of elementary functions in the medium-precision range", http://arxiv.org/abs/1410.7176

.. [Joh2015] \F. Johansson, "Computing Bell numbers", http://fredrikj.net/blog/2015/08/computing-bell-numbers/

.. [Joh2016] \F. Johansson, "Computing hypergeometric functions rigorously", preprint (2016), https://arxiv.org/abs/1606.06977

.. [Joh2017a] \F. Johansson. "Arb: efficient arbitrary-precision midpoint-radius interval arithmetic", IEEE Transactions on Computers, 66(8):1281-1292 (2017). https://doi.org/10.1109/TC.2017.2690633

.. [Joh2017b] \F. Johansson, "Computing the Lambert W function in arbitrary-precision complex interval arithmetic", preprint (2017), https://arxiv.org/abs/1705.03266

.. [Joh2018a] \F. Johansson, "Numerical integration in arbitrary-precision ball arithmetic", preprint (2018), https://arxiv.org/abs/1802.07942

.. [Joh2018b] \F. Johansson and others, "mpmath: a Python library for arbitrary-precision floating-point arithmetic (version 1.1.0)", December 2018. http://mpmath.org/

.. [JM2018] \F. Johansson and M. Mezzarobba, "Fast and rigorous arbitrary-precision computation of Gauss-Legendre quadrature nodes and weights", preprint (2018), https://arxiv.org/abs/1802.03948

.. [Kar1998] \E. A. Karatsuba, "Fast evaluation of the Hurwitz zeta function and Dirichlet L-series", Problems of Information Transmission 34:4 (1998), 342-353, http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=ppi&paperid=425&option_lang=eng

.. [Kob2010] \A. Kobel, "Certified Complex Numerical Root Finding", Seminar on Computational Geometry and Geometric Computing (2010), http://www.mpi-inf.mpg.de/departments/d1/teaching/ss10/Seminar_CGGC/Slides/02_Kobel_NRS.pdf

.. [Kri2013] \A. Krishnamoorthy and D. Menon, "Matrix Inversion Using Cholesky Decomposition" Proc. of the International Conference on Signal Processing Algorithms, Architectures, Arrangements, and Applications (SPA-2013), pp. 70-72, 2013.

.. [Leh1970] \R. S. Lehman, "On the Distribution of Zeros of the Riemann Zeta-Function", Proc. of the London Mathematical Society 20(3) (1970), 303-320, https://doi.org/10.1112/plms/s3-20.2.303

.. [Mic2007] \N. Michel, "Precise Coulomb wave functions for a wide range of complex l, eta and z", Computer Physics Communications, Volume 176, Issue 3, (2007), 232-249, https://doi.org/10.1016/j.cpc.2006.10.004

.. [Miy2010] \S. Miyajima, "Fast enclosure for all eigenvalues in generalized eigenvalue problems", Journal of Computational and Applied Mathematics, 233 (2010), 2994-3004, https://dx.doi.org/10.1016/j.cam.2009.11.048

.. [MPFR2012] The MPFR team, "MPFR Algorithms" (2012), http://www.mpfr.org/algo.html

.. [NIST2012] National Institute of Standards and Technology, *Digital Library of Mathematical Functions* (2012), http://dlmf.nist.gov/

.. [Olv1997] \F. Olver, *Asymptotics and special functions*, AKP Classics, AK Peters Ltd., Wellesley, MA, 1997. Reprint of the 1974 original.

.. [Rad1973] \H. Rademacher, *Topics in analytic number theory*, Springer, 1973.

.. [Pet1999] \K. Petras, "On the computation of the Gauss-Legendre quadrature formula with a given precision", Journal of Computational and Applied Mathematics 112 (1999), 253-267

.. [Pla2011] \D. J. Platt, "Computing degree 1 L-functions rigorously", Ph.D. Thesis, University of Bristol (2011), https://people.maths.bris.ac.uk/~madjp/thesis5.pdf

.. [Pla2017] \D. J. Platt, "Isolating some non-trivial zeros of zeta", Mathematics of Computation 86 (2017), 2449-2467, https://doi.org/10.1090/mcom/3198 

.. [PP2010] \K. H. Pilehrood and T. H. Pilehrood. "Series acceleration formulas for beta values", Discrete Mathematics and Theoretical Computer Science, DMTCS, 12 (2) (2010), 223-236, https://hal.inria.fr/hal-00990465/

.. [PS1973] \M. S. Paterson and L. J. Stockmeyer, "On the number of nonscalar multiplications necessary to evaluate polynomials", SIAM J. Comput (1973)

.. [PS1991] \G. Pittaluga and L. Sacripante, "Inequalities for the zeros of the Airy functions", SIAM J. Math. Anal. 22:1 (1991), 260-267.

.. [Rum2010] \S. M. Rump, "Verification methods: Rigorous results using floating-point arithmetic", Acta Numerica 19 (2010), 287-449.

.. [Smi2001] \D. M. Smith, "Algorithm: Fortran 90 Software for Floating-Point Multiple Precision Arithmetic, Gamma and Related Functions", Transactions on Mathematical Software 27 (2001) 377-387, http://myweb.lmu.edu/dmsmith/toms2001.pdf

.. [Tak2000] \D. Takahashi, "A fast algorithm for computing large Fibonacci numbers", Information Processing Letters 75 (2000) 243-246, http://www.ii.uni.wroc.pl/~lorys/IPL/article75-6-1.pdf

.. [Tre2008] \L. N. Trefethen, "Is Gauss Quadrature Better than Clenshaw-Curtis?", SIAM Review, 50:1 (2008), 67-87, https://doi.org/10.1137/060659831

.. [Tru2011] \T. S. Trudgian, "Improvements to Turing's method", Mathematics of Computation 80 (2011), 2259-2279, https://doi.org/10.1090/S0025-5718-2011-02470-1 

.. [Tru2014] \T. S. Trudgian, "An improved upper bound for the argument of the Riemann zeta-function on the critical line II", Journal of Number Theory 134 (2014), 280-292, https://doi.org/10.1016/j.jnt.2013.07.017

.. [Tur1953] \A. M. Turing, "Some Calculations of the Riemann Zeta-Function", Proc. of the London Mathematical Society 3(3) (1953), 99-117, https://doi.org/10.1112/plms/s3-3.1.99

.. [Wei2000] \A. Weilert, "(1+i)-ary GCD computation in Z[i] as an analogue to the binary GCD algorithm", Journal of Symbolic Computation 30.5 (2000): 605-617, https://doi.org/10.1006/jsco.2000.0422
