dnl
dnl Copyright (C) 2024 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

dnl TODO:
dnl * Do stuff in between to avoid latency penalties.
dnl * Convert %rX registers using 32-bit operations to %rax, ..., %rbp
dnl   registers to save a few bytes.

define(`rp',       `%rdi')
define(`ap',       `%rsi')
define(`bp_param', `%rdx')
define(`bp',	   `%r8')

define(`s0', `%rcx')
define(`s1', `%r9')
define(`s2', `%r10')
define(`s3', `%r11')
define(`s4', `%rbx')
define(`s5', `%rbp')
define(`s6', `%r12')
define(`s7', `%r13')
define(`s8', `%r14')
define(`s9', `%r15')

define(`sx',  `%rax')

	TEXT

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_1)
	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), sx, s0
	mov	s0, 0*8(rp)
	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_2)
	mov	1*8(bp_param), bp
	mov	0*8(bp_param), %rdx
	mulx	0*8(ap), sx, sx
	mulx	1*8(ap), s0, s1
	mov	bp, %rdx
	mulx	0*8(ap), s2, s3
	mulx	1*8(ap), %rdx, bp
	C (x, 0, 2), (1, 3, rdx), bp

	add	s0, sx
	adc	s3, s1
	adc	$0, bp
	C (x, 2), (1, rdx), bp

	add	s2, sx
	adc	%rdx, s1
	adc	$0, bp
	C x, 1, bp

	mov	s1, 0*8(rp)
	mov	bp, 1*8(rp)
	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_3)
	push	%rbx

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	1*8(%rsi), %r8, %rax
	mulx	2*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbx
	mulx	1*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	2*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbx
	adcx	%r8, %rax
	adcx	%rbx, %r10
	mulx	1*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adox	%rbx, %r11
	mulx	2*8(%rsi), %r8, %rbx
	adcx	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_4)
	push	%rbx
	push	%rbp

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	2*8(%rsi), %r8, %rax
	mulx	3*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %rbx
	mulx	2*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	3*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbp
	mulx	1*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	2*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	3*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbp
	adcx	%r8, %rax
	adcx	%rbp, %r10
	mulx	1*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adox	%rbp, %r11
	mulx	2*8(%rsi), %r8, %rbp
	adcx	%r8, %r11
	adcx	%rbp, %rbx
	mulx	3*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_5)
	push	%rbx
	push	%rbp
	push	%r12

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	3*8(%rsi), %r8, %rax
	mulx	4*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %rbx
	mulx	3*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	4*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %rbp
	mulx	2*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	3*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	4*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r12
	mulx	1*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	2*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	3*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	4*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r12
	adcx	%r8, %rax
	adcx	%r12, %r10
	mulx	1*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adox	%r12, %r11
	mulx	2*8(%rsi), %r8, %r12
	adcx	%r8, %r11
	adcx	%r12, %rbx
	mulx	3*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adox	%r12, %rbp
	mulx	4*8(%rsi), %r8, %r12
	adcx	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_6)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	4*8(%rsi), %r8, %rax
	mulx	5*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %rbx
	mulx	4*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	5*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %rbp
	mulx	3*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	4*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	5*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r12
	mulx	2*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	3*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	4*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r13
	mulx	1*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	2*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	3*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r13
	adcx	%r8, %rax
	adcx	%r13, %r10
	mulx	1*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adox	%r13, %r11
	mulx	2*8(%rsi), %r8, %r13
	adcx	%r8, %r11
	adcx	%r13, %rbx
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adox	%r13, %rbp
	mulx	4*8(%rsi), %r8, %r13
	adcx	%r8, %rbp
	adcx	%r13, %r12
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_7)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	5*8(%rsi), %r8, %rax
	mulx	6*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %rbx
	mulx	5*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	6*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %rbp
	mulx	4*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	5*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	6*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r12
	mulx	3*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	4*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r13
	mulx	2*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	3*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r14
	mulx	1*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	2*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r14
	adcx	%r8, %rax
	adcx	%r14, %r10
	mulx	1*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adox	%r14, %r11
	mulx	2*8(%rsi), %r8, %r14
	adcx	%r8, %r11
	adcx	%r14, %rbx
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adox	%r14, %rbp
	mulx	4*8(%rsi), %r8, %r14
	adcx	%r8, %rbp
	adcx	%r14, %r12
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adox	%r14, %r13
	mulx	6*8(%rsi), %r8, %r14
	adcx	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)
	mov	%r14, 6*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_8)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	6*8(%rsi), %r8, %rax
	mulx	7*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	5*8(%rsi), %r8, %rbx
	mulx	6*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	7*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %rbp
	mulx	5*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	6*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	7*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %r12
	mulx	4*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	5*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	7*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r13
	mulx	3*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	4*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	7*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r14
	mulx	2*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	3*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	7*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r15
	mulx	1*8(%rsi), %r8, %r14
	adcx	%r15, %rax
	adox	%r8, %rax
	adcx	%r14, %r10
	mulx	2*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adcx	%r14, %r11
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %r11
	adcx	%r14, %rbx
	mulx	4*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adcx	%r14, %rbp
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %rbp
	adcx	%r14, %r12
	mulx	6*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adcx	%r14, %r13
	mulx	7*8(%rsi), %r8, %r14
	adox	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	7*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %r15
	adcx	%r8, %rax
	adcx	%r15, %r10
	mulx	1*8(%rsi), %r8, %r15
	adox	%r8, %r10
	adox	%r15, %r11
	mulx	2*8(%rsi), %r8, %r15
	adcx	%r8, %r11
	adcx	%r15, %rbx
	mulx	3*8(%rsi), %r8, %r15
	adox	%r8, %rbx
	adox	%r15, %rbp
	mulx	4*8(%rsi), %r8, %r15
	adcx	%r8, %rbp
	adcx	%r15, %r12
	mulx	5*8(%rsi), %r8, %r15
	adox	%r8, %r12
	adox	%r15, %r13
	mulx	6*8(%rsi), %r8, %r15
	adcx	%r8, %r13
	adcx	%r15, %r14
	mulx	7*8(%rsi), %r8, %r15
	adox	%r8, %r14
	adcx	%r9, %r15
	adox	%r9, %r15

	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)
	mov	%rbp, 3*8(%rdi)
	mov	%r12, 4*8(%rdi)
	mov	%r13, 5*8(%rdi)
	mov	%r14, 6*8(%rdi)
	mov	%r15, 7*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mulhigh_9)
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	push	%rdi

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	7*8(%rsi), %r8, %rax
	mulx	8*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	6*8(%rsi), %r8, %rbx
	mulx	7*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	8*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	5*8(%rsi), %r8, %rbp
	mulx	6*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	7*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	8*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %r12
	mulx	5*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	7*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	8*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %r13
	mulx	4*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	7*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	8*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r14
	mulx	3*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	7*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	8*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r15
	mulx	2*8(%rsi), %r8, %r14
	adcx	%r15, %rax
	adox	%r8, %rax
	adcx	%r14, %r10
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adcx	%r14, %r11
	mulx	4*8(%rsi), %r8, %r14
	adox	%r8, %r11
	adcx	%r14, %rbx
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adcx	%r14, %rbp
	mulx	6*8(%rsi), %r8, %r14
	adox	%r8, %rbp
	adcx	%r14, %r12
	mulx	7*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adcx	%r14, %r13
	mulx	8*8(%rsi), %r8, %r14
	adox	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	7*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	mulx	1*8(%rsi), %r8, %r15
	adcx	%rdi, %rax
	adox	%r8, %rax
	adcx	%r15, %r10
	mulx	2*8(%rsi), %r8, %r15
	adox	%r8, %r10
	adcx	%r15, %r11
	mulx	3*8(%rsi), %r8, %r15
	adox	%r8, %r11
	adcx	%r15, %rbx
	mulx	4*8(%rsi), %r8, %r15
	adox	%r8, %rbx
	adcx	%r15, %rbp
	mulx	5*8(%rsi), %r8, %r15
	adox	%r8, %rbp
	adcx	%r15, %r12
	mulx	6*8(%rsi), %r8, %r15
	adox	%r8, %r12
	adcx	%r15, %r13
	mulx	7*8(%rsi), %r8, %r15
	adox	%r8, %r13
	adcx	%r15, %r14
	mulx	8*8(%rsi), %r8, %r15
	adox	%r8, %r14
	adcx	%r9, %r15
	adox	%r9, %r15

	mov	8*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	adcx	%r8, %rax
	adcx	%rdi, %r10
	mulx	1*8(%rsi), %r8, %rdi
	adox	%r8, %r10
	adox	%rdi, %r11
	mulx	2*8(%rsi), %r8, %rdi
	adcx	%r8, %r11
	adcx	%rdi, %rbx
	mulx	3*8(%rsi), %r8, %rdi
	adox	%r8, %rbx
	adox	%rdi, %rbp
	mulx	4*8(%rsi), %r8, %rdi
	adcx	%r8, %rbp
	adcx	%rdi, %r12
	mulx	5*8(%rsi), %r8, %rdi
	adox	%r8, %r12
	adox	%rdi, %r13
	mulx	6*8(%rsi), %r8, %rdi
	adcx	%r8, %r13
	adcx	%rdi, %r14
	mulx	7*8(%rsi), %r8, %rdi
	adox	%r8, %r14
	adox	%rdi, %r15
	mulx	8*8(%rsi), %r8, %rdi
	adcx	%r8, %r15
	adcx	%r9, %rdi
	pop	%r8
	adox	%r9, %rdi

	mov	%r10, 0*8(%r8)
	mov	%r11, 1*8(%r8)
	mov	%rbx, 2*8(%r8)
	mov	%rbp, 3*8(%r8)
	mov	%r12, 4*8(%r8)
	mov	%r13, 5*8(%r8)
	mov	%r14, 6*8(%r8)
	mov	%r15, 7*8(%r8)
	mov	%rdi, 8*8(%r8)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()
