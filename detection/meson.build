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

# -----------------------------------------------------------------------------
#   Unsupported options
# -----------------------------------------------------------------------------

if avx2_opt.enabled()
  error('AVX2 support is not yet implemented')
elif avx512_opt.enabled()
  error('AVX512 support is not yet implemented')
elif blas_opt.enabled()
  error('BLAS support is not yet implemented')
elif gc_opt.enabled()
  error('Garbage collection is not yet implemented')
elif ntl_opt.enabled()
  error('NTL support is not yet implemented')
endif

# -----------------------------------------------------------------------------
#   fmpz memory management
# -----------------------------------------------------------------------------

if gc_opt.allowed()
  if thread_safe_opt.allowed()
    error('Garbage collection is not thread safe')
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
#   libm
# -----------------------------------------------------------------------------

if not cc.has_function('atan2', dependencies: m_dep)
  error('libm is required')
endif

# -----------------------------------------------------------------------------
#   GMP
# -----------------------------------------------------------------------------

gmp_long_long_limb = cc.compiles('''
#include <gmp.h>
#if !defined(_LONG_LONG_LIMB)
# error mp_limb_t != unsigned long long int
#endif
''')

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

gmp_base_tests = {
  'mpz_init': ['', 'mpz_t x; __gmpz_init(x);'],
  'mpn_mul_basecase': [
    'void __gmpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);',
    'mp_limb_t a[2], b[2], p[4]; __gmpn_mul_basecase(p, a, (mp_size_t) 1, b, (mp_size_t) 1);',
  ],
}

foreach name, snippets : gmp_base_tests
  _code = gmp_test_code.format(snippets[0], snippets[1])
  if not cc.links(_code, dependencies: gmp_dep, name: name + ' links')
    message('This code did not link with -lgmp: ', _code)
    error('GMP does not have ' + name)
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
    if not cc.links(_code, dependencies: gmp_dep, name: name + ' links')
      message('This code did not link with -lgmp: ', _code)
      # Maybe disable gmp_internals rather than error?
      error('GMP does not have ' + name)
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
    if cc.links(_code, dependencies: gmp_dep, name: name + ' links')
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
  config_guess_result = run_command('../config/config.guess', check: true)
  target_triple = config_guess_result.stdout().strip()
  exact_cpu = target_triple.split('-')[0]
  message('EXACT CPU:', exact_cpu)
endif

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

# XXX: Do we need to check ARM64_PATTERN/X86_64_PATTERN?
fft_small_supported = (
  cc.compiles(fft_small_arm_code, dependencies: gmp_dep)
  or cc.compiles(fft_small_x86_code, dependencies: gmp_dep)
)

# Error if fft_small enabled but not possible, otherwise disable if auto
fft_small_opt = fft_small_opt.require(fft_small_supported,
    error_message: 'host CPU does not support AVX/NEON. Set fft_small to disabled or auto')

# If fft_small was auto we enable it if possible
fft_small_opt = fft_small_opt.enable_auto_if(fft_small_supported)
have_fft_small = fft_small_opt.enabled()

message('FFT_SMALL: ', have_fft_small ? 'enabled' : 'disabled')

# -----------------------------------------------------------------------------
#   assembly
# -----------------------------------------------------------------------------

assembly_opt = get_option('assembly')

# Might be better to check for known good operating systems rather than
# exclude known bad ones.
assembly_opt = assembly_opt.require(host_machine.system() != 'windows',
    error_message: 'No assembly support for Windows')

# Checking for assembly support needs to be more complicated than this. ARM
# needs to be v8 and not sure exactly what to check for x86_64. For now
# the checks below would attempt to build assembly on some systems that do not
# support it. The feature can be disabled with -Dassembly=disabled to build on
# those systems.
if host_machine.cpu_family() == 'aarch64'
  assembly_cpu_supported = true
  assembly_dir = 'mpn_extras/arm64'
  FLINT_HAVE_ASSEMBLY = 'armv8'
elif host_machine.cpu_family() == 'x86_64'
  assembly_cpu_supported = true
  assembly_dir = 'mpn_extras/x86_64/broadwell'
  FLINT_HAVE_ASSEMBLY = 'x86_64_adx'
else
  assembly_cpu_supported = false
  FLINT_HAVE_ASSEMBLY = ''
endif

assembly_opt = assembly_opt.require(assembly_cpu_supported,
    error_message: 'No assembly support for host CPU')

assembly_opt = assembly_opt.enable_auto_if(true)
have_assembly = assembly_opt.enabled()

message('ASSEMBLY: ', have_assembly ? 'enabled' : 'disabled')

# XXX: Need better detection of ASM language support
if have_assembly
  LSYM_PREFIX = 'L'
endif

# -----------------------------------------------------------------------------
#   mparam
# -----------------------------------------------------------------------------

flint_mparam = 'generic'

if host_machine.system() != 'windows'
  if host_machine.cpu_family() == 'aarch64'
    if host_machine.system() == 'darwin'
      flint_mparam = 'arm64/applem1'
    else
      flint_mparam = 'arm64'
    endif
  elif host_machine.cpu_family() == 'x86_64'
    # XXX: Need better CPU detection
    if have_assembly
      flint_mparam = 'x86_64/broadwell'
    else
      flint_mparam = 'generic'
    endif
  endif
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

if have_fft_small
  cfg_data.set('FLINT_HAVE_FFT_SMALL', 1,
    description: 'Define to use the fft_small module')
endif

# XXX: Need better CPU detection
# cfg_data.set('FLINT_KNOW_STRONG_ORDER', 1,
#   description: 'Define if system is strongly ordered')

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

# -----------------------------------------------------------------------------
#  internal configuration for config.h
# -----------------------------------------------------------------------------

cfg_data_internal = cfg_data

#cfg_data_internal.set('AC_APPLE_UNIVERSAL_BUILD', 1)

# XXX: Change this when coverage options are added
#cfg_data_internal.set('FLINT_COVERAGE', 1)

if cc.has_header_symbol('stdlib.h', 'aligned_alloc')
  cfg_data_internal.set('HAVE_ALIGNED_ALLOC', 1,
    description: 'Define to 1 if you have the \'aligned_alloc\' function.')
endif

ac_headers_check = [
  'alloca.h',
  'arm_neon.h',
  'dlfcn.h',
  'errno.h',
  'fenv.h',
  'float.h',
  'immintrin.h', # XXX: has_header is not equivalent to autotools check...
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

# XXX: Do we need to check for _alligned_malloc?

# XXX: Need better checking for LSYM_PREFIX and related assembly options
cfg_data_internal.set_quoted('LSYM_PREFIX', 'L',
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

# XXX: set WORDS_BIGENDIAN or undef inline here?
