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
coverage_opt = get_option('b_coverage')

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
gcc_warnings = '-Werror=implicit-function-declaration -Wall -Wno-stringop-overread -Wno-stringop-overflow'

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
    param_path = 'x86_64/zen3'
    have_avx512 = 'yes'

  elif exact_cpu in ['coreibwl', 'broadwell']
    gcc_cflags_arch = '-march=broadwell'
    param_path = 'x86_64/skylake'

  elif exact_cpu == 'skylake'
    gcc_cflags_arch = '-march=skylake'
    param_path = 'x86_64/skylake'

  elif exact_cpu == 'skylake_server'
    gcc_cflags_arch = '-march=skylake-avx512'
    param_path = 'x86_64/skylake'
    have_avx512 = 'yes'

  elif exact_cpu == 'cannonlake'
    gcc_cflags_arch = '-march=cannonlake'
    param_path = 'x86_64/skylake'
    have_avx512 = 'yes'

  elif exact_cpu == 'kabylake'
    gcc_cflags_arch = '-march=skylake'
    param_path = 'x86_64/skylake'

  elif exact_cpu == 'icelake'
    gcc_cflags_arch = '-march=icelake-client'
    param_path = 'x86_64/skylake'
    have_avx512 = 'yes'

  elif exact_cpu == 'icelake_server'
    gcc_cflags_arch = '-march=icelake-server'
    param_path = 'x86_64/skylake'
    have_avx512 = 'yes'

  elif exact_cpu == 'rocketlake'
    gcc_cflags_arch = '-march=rocketlake'
    param_path = 'x86_64/skylake'
    have_avx512 = 'yes'

  elif exact_cpu == 'tigerlake'
    gcc_cflags_arch = '-march=tigerlake'
    param_path = 'x86_64/skylake'
    have_avx512 = 'yes'

  elif exact_cpu == 'alderlake'
    gcc_cflags_arch = '-march=alderlake'
    param_path = 'x86_64/skylake'

  elif exact_cpu == 'raptorlake'
    gcc_cflags_arch = '-march=alderlake'
    param_path = 'x86_64/skylake'

  elif exact_cpu == 'knightslanding'
    gcc_cflags_arch = '-march=knl'
    param_path = 'x86_64/skylake'

  elif exact_cpu == 'sapphirerapids'
    gcc_cflags_arch = '-march=sapphirerapids'
    param_path = 'x86_64/skylake'
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
  gcc_cflags = '$gcc_cflags'
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
    and cc.compiles(fft_small_arm_code, dependencies: gmp_dep)
  )
else
  fft_small_supported = (
    flint_cv_have_fft_small_x86_i == 'yes'
    and cc.compiles(fft_small_x86_code, dependencies: gmp_dep)
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

if asm_path != ''
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

assembly_opt = assembly_opt.enable_auto_if(true)
have_assembly = assembly_opt.enabled()

message('ASSEMBLY: ', have_assembly ? 'enabled' : 'disabled')

# XXX: Need better detection of ASM language support
if have_assembly
  LSYM_PREFIX = 'L'
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
