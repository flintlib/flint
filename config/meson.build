# -----------------------------------------------------------------------------
#   FLINT configuration
# -----------------------------------------------------------------------------

# Data set on cfg_data is used to generate flint-config.h
# Data set on cfg_data_internal (below) is also used to generate config.h
cfg_data = configuration_data()

# -----------------------------------------------------------------------------
#   Project options
# -----------------------------------------------------------------------------

libtool_lock_opt = get_option('libtool_lock')
pthread_opt = get_option('pthread')
reentrant_opt = get_option('reentrant')
thread_safe_opt = get_option('thread_safe')
abi_opt = get_option('abi')
assert_opt = get_option('assert')
pretty_tests_opt = get_option('pretty_tests')
gmp_internals_opt = get_option('gmp_internals')
assembly_opt = get_option('assembly')
avx2_opt = get_option('avx2')
avx512_opt = get_option('avx512')
fft_small_opt = get_option('fft_small')
blas_opt = get_option('blas')
gc_opt = get_option('gc')
ntl_opt = get_option('ntl')
coverage_opt = get_option('b_coverage')
python = import('python').find_installation()

flint_abi_c_args = []
if abi_opt != 'auto'
  abi_flag = '-m' + abi_opt
  if not cc.has_argument(abi_flag)
    error(meson.get_compiler('c').get_id() + ' does not support ' + abi_flag)
  endif
  flint_abi_c_args += abi_flag
endif

# -----------------------------------------------------------------------------
#   libm
# -----------------------------------------------------------------------------

if not cc.has_function('atan2', dependencies: m_dep)
  error('libm is required')
endif

# -----------------------------------------------------------------------------
#   POSIX threads
# -----------------------------------------------------------------------------

pthread_dep = declare_dependency()
have_pthread = false

if pthread_opt.enabled()
  pthread_candidates = []

  if host_machine.system() == 'windows'
    pthread_candidates += [
      dependency('unofficial-pthreads', required: false, method: 'cmake', modules: ['unofficial::pthreads::pthread']),
      dependency('PThreads4W', required: false, method: 'cmake', modules: ['PThreads4W::PThreads4W']),
      dependency('PThreads4W', required: false, method: 'cmake', modules: ['PThreads4W::pthread']),
      dependency('PThreads4W', required: false, method: 'cmake', modules: ['PThreads4W::pthreadVC3']),
      dependency('pthreads4w', required: false),
    ]
  endif

  pthread_candidates += [
    dependency('pthread', required: false),
    dependency('pthreads', required: false),
    dependency('threads', required: false),
  ]

  foreach dep : pthread_candidates
    if dep.found() and cc.has_header('pthread.h', dependencies: dep) and cc.has_function('pthread_create', dependencies: dep)
      pthread_dep = dep
      have_pthread = true
      break
    endif
  endforeach

  if not have_pthread
    foreach libname : ['pthread', 'pthreadVC3', 'pthreadVCE3', 'pthreadVSE3']
      libpthread_dep = cc.find_library(libname, required: false)
      if libpthread_dep.found() and cc.has_header('pthread.h', dependencies: libpthread_dep) and cc.has_function('pthread_create', dependencies: libpthread_dep)
        pthread_dep = libpthread_dep
        have_pthread = true
        break
      endif
    endforeach
  endif

  if not have_pthread
    error('POSIX threads was not found.')
  endif
endif

# -----------------------------------------------------------------------------
#   BLAS
# -----------------------------------------------------------------------------

blas_dep = declare_dependency()
have_blas = false

if not blas_opt.disabled()
  blas_candidates = [
    dependency('cblas', required: false),
    dependency('openblas', required: false),
  ]

  foreach dep : blas_candidates
    if dep.found()
      dep_has_cblas_h = cc.has_header('cblas.h', dependencies: dep)
      dep_has_cblas_dgemm = dep_has_cblas_h and cc.has_function(
        'cblas_dgemm',
        prefix: '#include <cblas.h>',
        dependencies: dep)
      if dep_has_cblas_dgemm
        blas_dep = dep
        have_blas = true
        break
      endif
    endif
  endforeach

  if not have_blas
    foreach libname : ['cblas', 'openblas', 'blas']
      libdep = cc.find_library(libname, required: false)
      if libdep.found()
        lib_has_cblas_h = cc.has_header('cblas.h')
        lib_has_cblas_dgemm = lib_has_cblas_h and cc.has_function(
          'cblas_dgemm',
          prefix: '#include <cblas.h>',
          dependencies: libdep)
        if lib_has_cblas_dgemm
          blas_dep = libdep
          have_blas = true
          break
        endif
      endif
    endforeach
  endif

  if blas_opt.enabled() and not have_blas
    if not cc.has_header('cblas.h')
      error('Could not find cblas.h')
    else
      error('BLAS library was not found.')
    endif
  endif
endif

# -----------------------------------------------------------------------------
#   NTL
# -----------------------------------------------------------------------------

ntl_dep = declare_dependency()
have_ntl = false

if not ntl_opt.disabled()
  have_cpp = add_languages('cpp', native: false, required: false)

  if have_cpp
    cpp = meson.get_compiler('cpp')
    add_project_arguments(flint_abi_c_args, language: 'cpp')
    add_project_link_arguments(flint_abi_c_args, language: 'cpp')

    ntl_candidates = [
      dependency('ntl', required: false, static: true),
      dependency('ntl', required: false),
    ]

    foreach dep : ntl_candidates
      if dep.found()
        dep_has_ntl_h = cpp.has_header('NTL/ZZ.h', dependencies: dep)
        dep_links_ntl = dep_has_ntl_h and cpp.links(
          '''
          #include <iostream>
          #include <NTL/ZZ.h>

          int main() {
            NTL::ZZ a, b, c;
            std::cin >> a;
            std::cin >> b;
            c = (a + 1) * (b + 1);
            std::cout << c << "\\n";
            return 0;
          }
          ''',
          dependencies: dep,
          name: 'if linking with NTL works')
        if dep_links_ntl
          ntl_dep = dep
          have_ntl = true
          break
        endif
      endif
    endforeach

    if not have_ntl
      libntl_dep = cpp.find_library('ntl', required: false)
      if libntl_dep.found() and cpp.has_header('NTL/ZZ.h')
        libntl_full_dep = declare_dependency(
          dependencies: [libntl_dep, gmp_link_dep],
        )
        lib_links_ntl = cpp.links(
          '''
          #include <iostream>
          #include <NTL/ZZ.h>

          int main() {
            NTL::ZZ a, b, c;
            std::cin >> a;
            std::cin >> b;
            c = (a + 1) * (b + 1);
            std::cout << c << "\\n";
            return 0;
          }
          ''',
          dependencies: libntl_full_dep,
          name: 'if linking with NTL works')
        if lib_links_ntl
          ntl_dep = libntl_full_dep
          have_ntl = true
        endif
      endif
    endif

    if ntl_opt.enabled() and not have_ntl
      if not cpp.has_header('NTL/ZZ.h')
        error('Could not find NTL/ZZ.h')
      else
        error('NTL library was not found.')
      endif
    endif
  elif ntl_opt.enabled()
    error('A C++ compiler is required for NTL support.')
  endif
endif

# -----------------------------------------------------------------------------
#   Boehm GC
# -----------------------------------------------------------------------------

gc_dep = declare_dependency()
have_gc = false

if not gc_opt.disabled()
  gc_candidates = [
    dependency('bdw-gc', required: false),
    dependency('gc', required: false),
  ]

  foreach dep : gc_candidates
    if dep.found()
      dep_has_gc_h = cc.has_header('gc.h', dependencies: dep)
      if dep_has_gc_h
        gc_dep = dep
        have_gc = true
        break
      endif
    endif
  endforeach

  if not have_gc
    libgc_dep = cc.find_library('gc', required: false)
    if libgc_dep.found() and cc.has_header('gc.h')
      gc_dep = libgc_dep
      have_gc = true
    endif
  endif

  if gc_opt.enabled() and not have_gc
    if not cc.has_header('gc.h')
      error('Could not find gc.h')
    else
      error('GC library was not found.')
    endif
  endif
endif

# -----------------------------------------------------------------------------
#   fmpz memory management
# -----------------------------------------------------------------------------

if gc_opt.enabled()
  if thread_safe_opt.enabled()
    error('GC does not support TLS.')
  else
    fmpz_c_in = 'fmpz_gc.c'
  endif
else
  if reentrant_opt.enabled()
    fmpz_c_in = 'fmpz_reentrant.c'
  else
    fmpz_c_in = 'fmpz_single.c'
  endif
endif

# -----------------------------------------------------------------------------
#   GMP
# -----------------------------------------------------------------------------

if not cc.has_header('gmp.h', dependencies: gmp_dep)
  error('Could not find gmp.h')
endif

if not cc.has_header('mpfr.h', dependencies: mpfr_dep)
  error('Could not find mpfr.h')
endif

gmp_long_long_limb = cc.compiles('''
#include <gmp.h>
#if !defined(_LONG_LONG_LIMB)
# error mp_limb_t != unsigned long long int
#endif
''', dependencies: gmp_dep)

if gmp_long_long_limb
  cfg_data.set('FLINT_LONG_LONG', 1,
    description: 'Define to use long long limbs')
endif

if abi_opt != 'auto'
  flint_bits = abi_opt == '64' ? 64 : 32
else
gmp_64_bit_limb = cc.compiles('''
#include <gmp.h>
#if GMP_LIMB_BITS != 64
# error GMP_LIMB_BITS != 64
#endif
''', args: flint_abi_c_args, dependencies: gmp_dep, name: 'GMP_LIMB_BITS is 64')

flint_bits = gmp_64_bit_limb ? 64 : 32
endif

cfg_data.set('FLINT_BITS', flint_bits,
  description: 'Define according to the ABI FLINT was compiled with')

# We need to check with cc.links rather than cc.has_function etc because we
# are checking if we can link symbols not in the headers.

gmp_test_code = '''
#include <gmp.h>

@0@

int main() {
  @1@
  return 0;
}
'''

gmp_abi_error_suffix = abi_opt == 'auto' ? '' : ' for ABI ' + abi_opt + ' (check that GMP/MPFR are built for this ABI)'

gmp_base_tests = {
  'mpz_init': ['', 'mpz_t x; __gmpz_init(x);'],
  'mpn_mul_basecase': [
    'void __gmpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);',
    'mp_limb_t a[2], b[2], p[4]; __gmpn_mul_basecase(p, a, (mp_size_t) 1, b, (mp_size_t) 1);',
  ],
}

foreach name, snippets : gmp_base_tests
  _code = gmp_test_code.format(snippets[0], snippets[1])
  if not cc.links(_code, args: flint_abi_c_args, dependencies: gmp_dep, name: name + ' links')
    message('This code did not link with -lgmp: ', _code)
    error('GMP does not have ' + name + gmp_abi_error_suffix)
  endif
endforeach

if gmp_internals_opt.allowed()

  gmp_internals_required = {
    'mpn_gcd_11': ['', 'mp_limb_t res, a=0, b=0; res = __gmpn_gcd_11(a, b);'],
    'mpn_div_q': [
      'void __gmpn_div_q(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, mp_ptr);',
      'mp_limb_t q[2], n[2], d[2], t[2]; __gmpn_div_q(q, t, (mp_size_t) 1, n, (mp_size_t) 1, d);'
    ],
  }

  foreach name, snippets : gmp_internals_required
    _code = gmp_test_code.format(snippets[0], snippets[1])
    if not cc.links(_code, args: flint_abi_c_args, dependencies: gmp_dep, name: name + ' links')
      message('This code did not link with -lgmp: ', _code)
      # Match autotools: GMP internals are enabled by default and these symbols
      # are required unless -Dgmp_internals=disabled is used.
      error('GMP does not have ' + name + gmp_abi_error_suffix)
    endif
  endforeach

  gmp_internals_optional = {
    'mpn_add_n_sub_n': [
      'void __gmpn_add_n_sub_n(mp_ptr, mp_ptr, mp_ptr, mp_size_t);',
      'mp_limb_t a[2], b[2], c[2], d[2]; __gmpn_add_n_sub_n(a, b, c, (mp_size_t) 1);',
    ],
    'mpn_add_nc': [
      'mp_limb_t __gmpn_add_nc(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_add_nc(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_addlsh1_n': [
      'mp_limb_t __gmpn_addlsh1_n(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_addlsh1_n(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_addlsh1_n_ip1': [
      'mp_limb_t __gmpn_addlsh1_n_ip1(mp_ptr, mp_srcptr, mp_size_t);',
      'mp_limb_t a[2], b[2]; __gmpn_addlsh1_n_ip1(a, b, (mp_size_t) 1);',
    ],
    'mpn_addmul_2': [
      'void __gmpn_addmul_2(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_addmul_2(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_modexact_1_odd': [
      'mp_limb_t __gmpn_modexact_1_odd(mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2]; __gmpn_modexact_1_odd(a, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_rsh1add_n': [
      'mp_limb_t __gmpn_rsh1add_n(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_rsh1add_n(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_rsh1sub_n': [
      'mp_limb_t __gmpn_rsh1sub_n(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_rsh1sub_n(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
    'mpn_sub_nc': [
      'mp_limb_t __gmpn_sub_nc(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);',
      'mp_limb_t a[2], b[2]; __gmpn_sub_nc(a, b, (mp_size_t) 1, (mp_limb_t) 0);',
    ],
  }

  foreach name, snippets : gmp_internals_optional
    _code = gmp_test_code.format(snippets[0], snippets[1])
    if cc.links(_code, args: flint_abi_c_args, dependencies: gmp_dep, name: name + ' links')
      cfg_data.set('FLINT_HAVE_NATIVE_' + name, 1,
        description: 'Define if GMP has ' + name)
    endif
  endforeach

endif

# -----------------------------------------------------------------------------
#    CPU detection
# -----------------------------------------------------------------------------

# We check here for gcc or clang because the config.guess script probably does
# not work for other compilers. This check will also allow MinGW and cygwin on
# Windows but not MSVC.

# The config.guess script inspects the current CPU so should not be used for a
# cross build. Maybe there should be a way to configure the exact CPU as a
# project build option.

if not meson.is_cross_build() and cc.get_id() in ['gcc', 'clang']
  config_guess_result = run_command('./config.guess', check: true)
  target_triple = config_guess_result.stdout().strip()
  exact_cpu = target_triple.split('-')[0]
else
  # This could be 'aarch64' or 'x86_64'
  # https://mesonbuild.com/Reference-tables.html#cpu-families
  exact_cpu = host_machine.cpu_family()
endif

message('EXACT CPU:', exact_cpu)

################################################################################
# architecture specifics
################################################################################

gcc_cflags = '-O3 -pedantic -std=c11'
gcc_warnings = '-Werror=implicit-function-declaration -Werror=newline-eof -Wno-stringop-overread -Wno-stringop-overflow -Wno-deprecated-declarations'
gcc_cflags_optlist = ''
gcc_cflags_arch = ''
gcc_cflags_tune = ''
flint_warning_c_args = cc.get_supported_arguments(gcc_warnings.split())
flint_msvc_c_args = []
if cc.get_argument_syntax() == 'msvc'
  flint_msvc_c_args += cc.get_supported_arguments('/wd4146')
  if pthread_opt.enabled()
    flint_msvc_c_args += cc.get_supported_arguments('/experimental:c11atomics')
  endif
endif
flint_cpu_c_args = []

# We only try to provide specifics for those systems that currently supports
# our assembly routines. If more combinations are wished for than what is
# specified, please open up an issue at
# <https://github.com/flintlib/flint/issues/> and we will consider it.
#
# For these systems we aim to provide:
#
#   gcc_cflags     flags for GCC-compatible compilers
#
#   asm_path       directory for assembly, relative to src/mpn_extras
#   param_path     directory for flint-mparam.h, relative to src/mpn_extras
#
# For x86_64 systems, we can also set:
#
#   have_avx512    system has AVX512F and AVX512DQ (we assume this implies have_avx2)
#   have_avx2      system has AVX2 and FMA

# x86_64 CPUs that have AVX2:
X86_64_ADX_CPU_NAMES = [
  'zen1',
  'zen2',
  'zen3',
  'zen4',
  'zen5',
  'coreibwl',
  'broadwell',
  'skylake',
  'skylake_server',
  'cannonlake',
  'kabylake',
  'icelake',
  'icelake_server',
  'rocketlake',
  'tigerlake',
  'alderlake',
  'raptorlake',
  'meteorlake',
  'knightslanding',
  'sapphirerapids',
  'cometlake',
]

if host_machine.cpu_family() == 'aarch64'

  gcc_cflags  =  gcc_cflags
  gcc_cflags_optlist = 'arch tune'
  asm_path = 'arm64'
  param_path = 'arm64'
  flint_cv_have_fft_small_arm_i = 'yes'
  flint_know_strong_order = 'no'

  # NOTE: Cortex values where taken from
  # https://developer.arm.com/Processors/Cortex-XXX

  if exact_cpu in ['armcortexa35', 'armcortexa35neon']
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = '-mtune=cortex-a35'

  elif exact_cpu in ['armcortexa53', 'armcortexa53neon']
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = '-mtune=cortex-a53'

  elif exact_cpu in ['armcortexa55', 'armcortexa55neon']
    gcc_cflags_arch = '-march=armv8.2-a'
    gcc_cflags_tune = '-mtune=cortex-a55'

  elif exact_cpu in ['armcortexa57', 'armcortexa57neon']
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = '-mtune=cortex-a57'

  elif exact_cpu in ['armcortexa72', 'armcortexa72neon']
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = '-mtune=cortex-a72'

  elif exact_cpu in ['armcortexa73', 'armcortexa73neon']
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = '-mtune=cortex-a73'

  elif exact_cpu in ['armcortexa75', 'armcortexa75neon']
    gcc_cflags_arch = '-march=armv8.2-a'
    gcc_cflags_tune = '-mtune=cortex-a75'

  elif exact_cpu in ['armcortexa76', 'armcortexa76neon']
    gcc_cflags_arch = '-march=armv8.2-a'
    gcc_cflags_tune = '-mtune=cortex-a76'

  elif exact_cpu in ['armcortexa77', 'armcortexa77neon']
    gcc_cflags_arch = '-march=armv8.2-a'
    gcc_cflags_tune = '-mtune=cortex-a77'

  elif exact_cpu in ['armcortexa65', 'armcortexa65neon']
    gcc_cflags_arch = '-march=armv8.2-a'
    gcc_cflags_tune = '-mtune=cortex-a65'

  elif exact_cpu in ['armcortexa34', 'armcortexa34neon']
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = '-mtune=cortex-a34'

  elif exact_cpu in ['armcortexa78', 'armcortexa78neon']
    gcc_cflags_arch = '-march=armv8.2-a'
    gcc_cflags_tune = '-mtune=cortex-a78'

  elif exact_cpu in ['armexynosm1']
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = '-mtune=exynos-m1'

  elif exact_cpu in ['armthunderx']
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = '-mtune=thunderx'

  elif exact_cpu in ['armxgene1']
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = '-mtune=xgene1'

  elif exact_cpu.startswith('applem1')
    gcc_cflags_arch = '-march=armv8.5-a'
    gcc_cflags_tune = ''

  elif exact_cpu.startswith('applem2')
    gcc_cflags_arch = '-march=armv8.6-a'
    gcc_cflags_tune = ''

  elif exact_cpu.startswith('applem3')
    gcc_cflags_arch = '-march=armv8.6-a'
    gcc_cflags_tune = ''

  elif exact_cpu.startswith('applem4')
    gcc_cflags_arch = '-march=armv8.7-a'
    gcc_cflags_tune = ''

  elif exact_cpu.startswith('aarch64') or exact_cpu.startswith('armv8')
    gcc_cflags_arch = '-march=armv8-a'
    gcc_cflags_tune = ''

  else
    error_msg = '''
    Error in deciding flags for @0@.
    Please report at <https://github.com/flintlib/flint/issues/>
    '''
    error(error_msg.format(exact_cpu))

  endif

elif host_machine.cpu_family() == 'x86_64' and exact_cpu in X86_64_ADX_CPU_NAMES

  gcc_cflags = gcc_cflags
  gcc_cflags_optlist = 'arch'
  asm_path = 'x86_64/broadwell'
  param_path = 'x86_64/broadwell'

  have_avx512 = 'no'
  have_avx2 = 'yes'
  flint_cv_have_fft_small_x86_i = 'yes'
  flint_know_strong_order = 'yes'

  if exact_cpu == 'zen'
    gcc_cflags_arch = '-march=znver1'
    param_path = 'x86_64/zen3'

  elif exact_cpu == 'zen2'
    gcc_cflags_arch = '-march=znver2'
    param_path = 'x86_64/zen3'

  elif exact_cpu == 'zen3'
    gcc_cflags_arch = '-march=znver3'
    param_path = 'x86_64/zen3'

  elif exact_cpu == 'zen4'
    gcc_cflags_arch = '-march=znver4'
    param_path = 'x86_64/zen4'
    have_avx512 = 'yes'

  elif exact_cpu == 'zen5'
    gcc_cflags_arch = '-march=znver5 -march=znver4'
    param_path = 'x86_64/zen4'
    have_avx512 = 'yes'

  elif exact_cpu in ['coreibwl', 'broadwell']
    gcc_cflags_arch = '-march=broadwell'
    param_path = 'x86_64/broadwell'

  elif exact_cpu == 'skylake'
    gcc_cflags_arch = '-march=skylake'
    param_path = 'x86_64/skylake'

  elif exact_cpu == 'skylake_server'
    gcc_cflags_arch = '-march=skylake-avx512'
    param_path = 'x86_64/skylake'
    have_avx512 = 'yes'

  elif exact_cpu == 'cannonlake'
    gcc_cflags_arch = '-march=cannonlake'
    param_path = 'x86_64/icelake'
    have_avx512 = 'yes'

  elif exact_cpu == 'kabylake'
    gcc_cflags_arch = '-march=skylake'
    param_path = 'x86_64/skylake'

  elif exact_cpu == 'icelake'
    gcc_cflags_arch = '-march=icelake-client'
    param_path = 'x86_64/icelake'
    have_avx512 = 'yes'

  elif exact_cpu == 'icelake_server'
    gcc_cflags_arch = '-march=icelake-server'
    param_path = 'x86_64/icelake'
    have_avx512 = 'yes'

  elif exact_cpu == 'rocketlake'
    gcc_cflags_arch = '-march=rocketlake'
    param_path = 'x86_64/icelake'
    have_avx512 = 'yes'

  elif exact_cpu == 'tigerlake'
    gcc_cflags_arch = '-march=tigerlake'
    param_path = 'x86_64/tigerlake'
    have_avx512 = 'yes'

  elif exact_cpu == 'alderlake'
    gcc_cflags_arch = '-march=alderlake'
    param_path = 'x86_64/icelake'

  elif exact_cpu == 'raptorlake'
    gcc_cflags_arch = '-march=raptorlake -march=alderlake'
    param_path = 'x86_64/raptorlake'

  elif exact_cpu == 'meteorlake'
    gcc_cflags_arch = '-march=meteorlake -march=alderlake'
    param_path = 'x86_64/meteorlake'

  elif exact_cpu == 'knightslanding'
    gcc_cflags_arch = '-march=knl'
    param_path = 'x86_64/skylake'

  elif exact_cpu == 'sapphirerapids'
    gcc_cflags_arch = '-march=sapphirerapids'
    param_path = 'x86_64/icelake'
    have_avx512 = 'yes'

  elif exact_cpu == 'cometlake'
    gcc_cflags_arch = '-march=kabylake'
    param_path = 'x86_64/skylake'

  else
    error_msg = '''
    Error in deciding flags for @0@.
    Please report at <https://github.com/flintlib/flint/issues/>
    '''
    error(error_msg.format(exact_cpu))

  endif

  # Disable assembly on Windows
  if host_machine.system() == 'windows'
    flint_nonstd_abi = 'yes'
    asm_path = ''
    param_path = 'x86_64'
  endif

elif host_machine.cpu_family() == 'x86_64'

  # NOTE: We do not care if the user specifies noavx!
  gcc_cflags = gcc_cflags
  gcc_cflags_optlist = 'arch'
  asm_path = ''
  param_path = 'x86_64'

  have_avx512 = 'no'
  have_avx2 = 'no'
  flint_know_strong_order = 'yes'

  if exact_cpu in ['coreiibr', 'ivybridge']
    gcc_cflags_arch = '-march=ivybridge'

  elif exact_cpu in ['coreihwl', 'haswell']
    gcc_cflags_arch = '-march=haswell'

  elif exact_cpu in ['piledriver', 'bd2']
    gcc_cflags_arch = '-march=bdver2'

  elif exact_cpu in ['steamroller', 'bd3']
    gcc_cflags_arch = '-march=bdver3'

  elif exact_cpu in ['excavator', 'bd4']
    gcc_cflags_arch = '-march=bdver4'
    have_avx2 = 'yes'

  elif exact_cpu in ['x86_64v3']
    gcc_cflags_arch = '-march=x86-64-v3'
    have_avx2 = 'yes'

  elif exact_cpu in ['x86_64v4']
    gcc_cflags_arch = '-march=x86-64-v4'
    have_avx512 = 'yes'

  else
    flint_know_strong_order = 'no'

  endif

  if have_avx512 == 'yes'
    have_avx2 = 'yes'
  endif

  if have_avx2 == 'yes'
    flint_cv_have_fft_small_x86_i = 'yes'
  else
    flint_cv_have_fft_small_x86_i = 'no'
  endif

else
  # Generic case for unrecognised CPU
  asm_path = ''
  param_path = 'generic'
  gcc_cflags = gcc_cflags

  have_avx512 = 'no'
  have_avx2 = 'no'
  flint_cv_have_fft_small_x86_i = 'no'
  flint_know_strong_order = 'no'

endif

# Probe and collect CPU-specific compiler flags. For architecture flags, the
# autotools build treats a whitespace-separated list as ordered fallbacks.
foreach opt : gcc_cflags_optlist.split()
  candidate_flags = get_variable('gcc_cflags_' + opt, '').split()
  found_single = false

  foreach flag : candidate_flags
    if found_single
      break
    endif

    if cc.has_argument(flag)
      flint_cpu_c_args += flag
      if opt == 'arch'
        found_single = true
      endif
    elif opt == 'arch'
      flint_cv_have_fft_small_x86_i = 'no'
      flint_cv_have_fft_small_arm_i = 'no'
    endif
  endforeach
endforeach

if avx2_opt.enabled()
  if host_machine.cpu_family() != 'x86_64'
    error('AVX2 is only supported for x86_64 builds')
  endif
  if not cc.has_argument('-mavx2')
    error(meson.get_compiler('c').get_id() + ' does not support -mavx2')
  endif
  flint_cpu_c_args += '-mavx2'
  if cc.has_argument('-mfma')
    flint_cpu_c_args += '-mfma'
  endif
  flint_cv_have_fft_small_x86_i = 'yes'
endif

if avx512_opt.enabled()
  if host_machine.cpu_family() != 'x86_64'
    error('AVX512 is only supported for x86_64 builds')
  endif
  if not cc.has_argument('-mavx512f')
    error(meson.get_compiler('c').get_id() + ' does not support -mavx512f')
  endif
  flint_cpu_c_args += '-mavx512f'
endif

message('CPU-specific C arguments:', ' '.join(flint_cpu_c_args))

# -----------------------------------------------------------------------------
#   mparam
# -----------------------------------------------------------------------------

flint_mparam = param_path

# -----------------------------------------------------------------------------
#    fft_small
# -----------------------------------------------------------------------------

fft_small_opt = get_option('fft_small')

fft_small_arm_code = '''
#include <arm_neon.h>
#if !(defined(__GNUC__) && defined(__ARM_NEON))
# if !(defined(_MSC_VER) && defined(_M_ARM64))
#  error
error
# endif
#endif
int main(){return 0;}
'''

fft_small_x86_code = '''
#if defined(__GNUC__)
# include <x86intrin.h>
#elif defined(_MSC_VER)
# include <intrin.h>
#else
# error
error
#endif

#if !defined(__AVX2__)
# error
error
#endif
int main(){return 0;}
'''

if host_machine.cpu_family() == 'aarch64'
  fft_small_supported = (
    flint_cv_have_fft_small_arm_i == 'yes'
    and cc.compiles(fft_small_arm_code,
      args: flint_cpu_c_args,
      dependencies: gmp_dep)
  )
else
  fft_small_supported = (
    flint_cv_have_fft_small_x86_i == 'yes'
    and cc.compiles(fft_small_x86_code,
      args: flint_cpu_c_args,
      dependencies: gmp_dep)
  )
endif

# Error if fft_small enabled but not possible, otherwise disable if auto
fft_small_opt = fft_small_opt.require(fft_small_supported,
    error_message: 'host CPU does not support AVX/NEON. Set fft_small to disabled or auto')

# If fft_small was auto we enable it if possible
fft_small_opt = fft_small_opt.enable_auto_if(fft_small_supported)
have_fft_small = fft_small_opt.enabled()

message('FFT_SMALL: ', have_fft_small ? 'enabled' : 'disabled')

# -----------------------------------------------------------------------------
#   Assembly
# -----------------------------------------------------------------------------

asm_text = '.text'
asm_data = '.data'
asm_label_suffix = ':'
asm_globl = '.globl'
asm_globl_attr = ''
asm_gsym_prefix = ''
asm_rodata = '.section .rodata'
asm_type = ''
asm_size = ''
asm_lsym_prefix = 'L'
asm_align_logarithmic = 'yes'
asm_align_fill_0x90 = 'no'
asm_have_coff_type = 'no'

if asm_path != '' and flint_bits == 64
  assembly_cpu_supported = true
  assembly_dir = 'mpn_extras/' + asm_path
  if host_machine.cpu_family() == 'aarch64'
    FLINT_HAVE_ASSEMBLY = 'armv8'
  elif host_machine.cpu_family() == 'x86_64'
    FLINT_HAVE_ASSEMBLY = 'x86_64_adx'
  else
    error('Assembly only supported for x86_64 or aarch64')
  endif
else
  assembly_cpu_supported = false
  FLINT_HAVE_ASSEMBLY = ''
endif

assembly_opt = assembly_opt.require(host_machine.system() != 'windows',
    error_message: 'No assembly support for Windows')

assembly_opt = assembly_opt.require(assembly_cpu_supported,
    error_message: 'No assembly support for host CPU')

if host_machine.cpu_family() == 'aarch64' or host_machine.cpu_family() == 'x86_64'
  asm_detect = run_command(
    python,
    files('check_assembly.py'),
    [
      '--host-cpu', host_machine.cpu_family(),
      '--host-system', host_machine.system(),
      '--',
    ] + cc.cmd_array() + flint_abi_c_args + flint_cpu_c_args,
    check: true,
  )

  asm_detect_lines = asm_detect.stdout().strip().split('\n')
  foreach line : asm_detect_lines
    if line.contains('=')
      asm_pair = line.split('=')
      asm_key = asm_pair[0]
      asm_value = ''
      if asm_pair.length() > 1
        asm_value = asm_pair[1]
      endif
      if asm_key == 'TEXT'
        asm_text = asm_value
      elif asm_key == 'DATA'
        asm_data = asm_value
      elif asm_key == 'LABEL_SUFFIX'
        asm_label_suffix = asm_value
      elif asm_key == 'GLOBL'
        asm_globl = asm_value
      elif asm_key == 'GLOBL_ATTR'
        asm_globl_attr = asm_value
      elif asm_key == 'GSYM_PREFIX'
        asm_gsym_prefix = asm_value
      elif asm_key == 'RODATA'
        asm_rodata = asm_value
      elif asm_key == 'TYPE'
        asm_type = asm_value
      elif asm_key == 'SIZE'
        asm_size = asm_value
      elif asm_key == 'LSYM_PREFIX'
        asm_lsym_prefix = asm_value
      elif asm_key == 'ALIGN_LOGARITHMIC'
        asm_align_logarithmic = asm_value
      elif asm_key == 'ALIGN_FILL_0x90'
        asm_align_fill_0x90 = asm_value
      elif asm_key == 'HAVE_COFF_TYPE'
        asm_have_coff_type = asm_value
      endif
    endif
  endforeach
endif

asm_m4_includes = ''
if host_machine.cpu_family() == 'aarch64'
  if host_machine.system() == 'darwin'
    asm_m4_includes = 'include(' + meson.project_source_root() + '/src/mpn_extras/arm64/darwin.m4)\n'
  else
    asm_m4_includes = 'include(' + meson.project_source_root() + '/src/mpn_extras/arm64/arm64-defs.m4)\n'
  endif
elif host_machine.cpu_family() == 'x86_64'
  asm_m4_includes = 'include(' + meson.project_source_root() + '/src/mpn_extras/x86_64/x86_64-defs.m4)\n'
  if host_machine.system() == 'darwin'
    asm_m4_includes += 'include(' + meson.project_source_root() + '/src/mpn_extras/x86_64/darwin.m4)\n'
  endif
endif

asm_cfg = configuration_data()
asm_cfg.set('SOURCE_DIR', meson.project_source_root())
asm_cfg.set('TEXT', asm_text)
asm_cfg.set('DATA', asm_data)
asm_cfg.set('LABEL_SUFFIX', asm_label_suffix)
asm_cfg.set('GLOBL', asm_globl)
asm_cfg.set('GLOBL_ATTR', asm_globl_attr)
asm_cfg.set('GSYM_PREFIX', asm_gsym_prefix)
asm_cfg.set('RODATA', asm_rodata)
asm_cfg.set('TYPE', asm_type)
asm_cfg.set('SIZE', asm_size)
asm_cfg.set('LSYM_PREFIX', asm_lsym_prefix)
asm_cfg.set('ALIGN_LOGARITHMIC', asm_align_logarithmic)
asm_cfg.set('ALIGN_FILL_0x90', asm_align_fill_0x90)
asm_cfg.set('HAVE_COFF_TYPE', asm_have_coff_type)
asm_cfg.set('ASM_M4_INCLUDES', asm_m4_includes)

assembly_opt = assembly_opt.enable_auto_if(true)
have_assembly = assembly_opt.enabled()

message('ASSEMBLY: ', have_assembly ? 'enabled' : 'disabled')
if have_assembly
  message('ASM TEXT: "@0@"'.format(asm_text))
  message('ASM DATA: "@0@"'.format(asm_data))
  message('ASM LABEL_SUFFIX: "@0@"'.format(asm_label_suffix))
  message('ASM GLOBL: "@0@"'.format(asm_globl))
  message('ASM GLOBL_ATTR: "@0@"'.format(asm_globl_attr))
  message('ASM GSYM_PREFIX: "@0@"'.format(asm_gsym_prefix))
  message('ASM RODATA: "@0@"'.format(asm_rodata))
  message('ASM TYPE: "@0@"'.format(asm_type))
  message('ASM SIZE: "@0@"'.format(asm_size))
  message('ASM LSYM_PREFIX: "@0@"'.format(asm_lsym_prefix))
  message('ASM ALIGN_LOGARITHMIC: "@0@"'.format(asm_align_logarithmic))
  message('ASM ALIGN_FILL_0x90: "@0@"'.format(asm_align_fill_0x90))
  message('ASM HAVE_COFF_TYPE: "@0@"'.format(asm_have_coff_type))
endif

# -----------------------------------------------------------------------------
#   cpu_set_t support
# -----------------------------------------------------------------------------

have_cpu_set_t = cc.compiles('''
#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>

int main() {
  cpu_set_t s;
  CPU_ZERO(&s);
  pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &s);
  return 0;
}
''')

# -----------------------------------------------------------------------------
#  external configuration for flint_config.h
# -----------------------------------------------------------------------------

if host_machine.endian() == 'big'
  cfg_data.set('FLINT_BIG_ENDIAN', 1,
    description: 'Define if system is big endian.')
endif

if have_assembly and FLINT_HAVE_ASSEMBLY == 'armv8'
  cfg_data.set('FLINT_HAVE_ASSEMBLY_armv8', 1,
    description: 'Define if Arm v8 assembly is available')
elif have_assembly and FLINT_HAVE_ASSEMBLY == 'x86_64_adx'
  cfg_data.set('FLINT_HAVE_ASSEMBLY_x86_64_adx', 1,
    description: 'Define if x86_64 ADX assembly is available')
endif

if have_assembly
  cfg_data.set('FLINT_WANT_ASSEMBLY', 1,
    description: 'Define to enable inline assembly.')
endif

if have_fft_small
  cfg_data.set('FLINT_HAVE_FFT_SMALL', 1,
    description: 'Define to use the fft_small module')
endif

if flint_know_strong_order == 'yes'
  cfg_data.set('FLINT_KNOW_STRONG_ORDER', 1,
    description: 'Define if system is strongly ordered')
endif

if reentrant_opt.enabled()
  cfg_data.set('FLINT_REENTRANT', 1,
    description: 'Define to enable reentrant.')
endif

cfg_data.set('FLINT_UNROLL_LOOPS', 1,
  description: 'Define to locally unroll some loops')

if blas_opt.enabled()
  cfg_data.set('FLINT_USES_BLAS', 1,
    description: 'Define to enable BLAS.')
endif

if have_cpu_set_t
  cfg_data.set('FLINT_USES_CPUSET', 1,
    description: 'Define if system has cpu_set_t')
endif

if gc_opt.enabled()
  cfg_data.set('FLINT_USES_GC', 1,
    description: 'Define to enable the Boehm-Demers-Weise garbage collector.')
endif

if pthread_opt.enabled()
  cfg_data.set('FLINT_USES_PTHREAD', 1,
    description: 'Define to enable the use of pthread.')
endif

if thread_safe_opt.enabled()
  cfg_data.set('FLINT_USES_TLS', 1,
    description: 'Define to enable thread-local storage.')
endif

if assert_opt.enabled()
  cfg_data.set('FLINT_WANT_ASSERT', 1,
    description: 'Define to enable use of asserts.')
endif

if gmp_internals_opt.allowed()
  cfg_data.set('FLINT_WANT_GMP_INTERNALS', 1,
    description: 'Define to enable use of GMP internals.')
endif

if pretty_tests_opt.enabled()
  cfg_data.set('FLINT_WANT_PRETTY_TESTS', 1,
    description: 'Define to enable pretty printing for tests.')
endif

if cc.get_argument_syntax() == 'msvc' and get_option('default_library') == 'static'
  cfg_data.set('FLINT_STATIC_BUILD', 1,
    description: 'Define when building FLINT statically with MSVC.')
endif

# -----------------------------------------------------------------------------
#  compiler pragmas
# -----------------------------------------------------------------------------

pragma_compiler = ''
if cc.get_id() == 'clang'
  pragma_compiler = 'clang'
elif cc.get_id() == 'gcc'
  pragma_compiler = 'GCC'
endif

pragma_checks = {
  'DIAGNOSTIC_PUSH': [
    'diagnostic push',
    '',
  ],
  'DIAGNOSTIC_POP': [
    'diagnostic pop',
    'diagnostic push',
  ],
  'DIAGNOSTIC_IGNORE_INCOMPATIBLE_FUNCTION_POINTER_TYPES': [
    'diagnostic ignored \\"-Wincompatible-function-pointer-types\\"',
    '',
  ],
  'DIAGNOSTIC_IGNORE_DISCARDED_QUALIFIERS': [
    'diagnostic ignored \\"-Wdiscarded-qualifiers\\"',
    '',
  ],
  'DIAGNOSTIC_IGNORE_FORMAT': [
    'diagnostic ignored \\"-Wformat\\"',
    '',
  ],
  'DIAGNOSTIC_IGNORE_DANGLING_POINTER': [
    'diagnostic ignored \\"-Wdangling-pointer\\"',
    '',
  ],
  'DIAGNOSTIC_IGNORE_CAST_FUNCTION_TYPE': [
    'diagnostic ignored \\"-Wcast-function-type\\"',
    '',
  ],
  'DIAGNOSTIC_IGNORE_OVERLENGTH_STRINGS': [
    'diagnostic ignored \\"-Woverlength-strings\\"',
    '',
  ],
  'DIAGNOSTIC_IGNORE_UNUSED_VARIABLE': [
    'diagnostic ignored \\"-Wunused-variable\\"',
    '',
  ],
  'DIAGNOSTIC_IGNORE_MAYBE_UNINITIALIZED': [
    'diagnostic ignored \\"-Wmaybe-uninitialized\\"',
    '',
  ],
  'PUSH_OPTIONS': [
    'push_options',
    '',
  ],
  'POP_OPTIONS': [
    'pop_options',
    'push_options',
  ],
  'OPTIMIZE_O2': [
    'optimize (\\"O2\\")',
    '',
  ],
  'OPTIMIZE_OSIZE': [
    'optimize (\\"Os\\")',
    '',
  ],
  'OPTIMIZE_UNROLL_LOOPS': [
    'optimize (\\"unroll-loops\\")',
    '',
  ],
}

foreach macro_name, pragma_info : pragma_checks
  cfg_data.set(macro_name, '',
    description: 'Define to the corresponding compiler pragma if supported.')

  if pragma_compiler != ''
    pragma = pragma_info[0]
    pre_pragma = pragma_info[1]
    pragma_value = '_Pragma("@0@ @1@")'.format(pragma_compiler, pragma)
    pragma_test = ''

    if pre_pragma != ''
      pragma_test += '_Pragma("@0@ @1@")\n'.format(pragma_compiler, pre_pragma)
    endif

    pragma_test += pragma_value + '\nint main(void) { return 0; }\n'

    if cc.compiles(pragma_test,
        args: ['-Werror', '-Wunknown-pragmas'],
        name: macro_name + ' pragma')
      cfg_data.set(macro_name, pragma_value,
        description: 'Define to the corresponding compiler pragma if supported.')
    endif
  endif
endforeach

# -----------------------------------------------------------------------------
#  internal configuration for config.h
# -----------------------------------------------------------------------------

cfg_data_internal = cfg_data

#cfg_data_internal.set('AC_APPLE_UNIVERSAL_BUILD', 1)

if coverage_opt
  cfg_data_internal.set('FLINT_COVERAGE', 1)
endif

if cc.has_header_symbol('stdlib.h', 'aligned_alloc')
  cfg_data_internal.set('HAVE_ALIGNED_ALLOC', 1,
    description: 'Define to 1 if you have the \'aligned_alloc\' function.')
endif

if cc.has_function('_aligned_malloc', prefix: '#include <malloc.h>')
  cfg_data_internal.set('HAVE__ALIGNED_MALLOC', 1,
    description: 'Define to 1 if you have the \'_aligned_malloc\' function.')
endif

mandatory_headers = [
  'stdarg.h',
  'math.h',
  'float.h',
  'errno.h',
]

foreach h : mandatory_headers
  if not cc.has_header(h)
    error('Could not find a mandatory header: ' + h)
  endif
endforeach

if pthread_opt.enabled() and not have_pthread and not cc.has_header('pthread.h')
  error('Could not find pthread.h')
endif

if assert_opt.enabled() and not cc.has_header('assert.h')
  error('Could not find assert.h')
endif

ac_headers_check = [
  'alloca.h',
  'arm_neon.h',
  'dlfcn.h',
  'errno.h',
  'fenv.h',
  'float.h',
  'inttypes.h',
  'malloc.h',
  'math.h',
  'pthread_np.h',
  'stdarg.h',
  'stdint.h',
  'stdio.h',
  'stdlib.h',
  'string.h',
  'strings.h',
  'sys/param.h',
  'sys/stat.h',
  'sys/types.h',
  'unistd.h',
  'windows.h',
]

foreach h : ac_headers_check
  if cc.has_header(h)
    cfg_data_internal.set('HAVE_' + h.underscorify().to_upper(), 1,
      description: f'Define to 1 if you have the <@h@> header file.')
  endif
endforeach

# XXX: Need better checking for LSYM_PREFIX and related assembly options
cfg_data_internal.set_quoted('LSYM_PREFIX', asm_lsym_prefix,
  description: 'Assembler local label prefix')

cfg_data_internal.set_quoted('LT_OBJDIR', '.libs/',
  description: 'Define to the sub-directory where libtool stores uninstalled libraries.')
cfg_data_internal.set_quoted('PACKAGE_BUGREPORT',
  'https://github.com/flintlib/flint/issues/',
  description: 'Define to the address where bug reports for this package should be sent.')
cfg_data_internal.set_quoted('PACKAGE_NAME', 'FLINT',
  description: 'Define to the full name of this package.')
cfg_data_internal.set_quoted('PACKAGE_STRING', 'FLINT ' + FLINT_VERSION_FULL,
  description: 'Define to the full name and version of this package.')
cfg_data_internal.set_quoted('PACKAGE_TARNAME', 'flint',
  description: 'Define to the one symbol short name of this package.')
cfg_data_internal.set_quoted('PACKAGE_URL', 'https://flintlib.org/',
  description: 'Define to the home page for this package.')
cfg_data_internal.set_quoted('PACKAGE_VERSION', FLINT_VERSION_FULL,
  description: 'Define to the version of this package.')

cfg_data_internal.set('STDC_HEADERS', 1,
  description: '''Define to 1 if all of the C89 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it.''')
