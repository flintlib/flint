dnl
dnl Copyright (C) 2023 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

define(`rp', `%rdi')
define(`ap', `%rsi')

define(`s0', `%rcx')
define(`s1', `%r8')
define(`s2', `%r9')
define(`s3', `%r10')
define(`s4', `%r11')
define(`s5', `%rbx')
define(`s6', `%rbp')
define(`s7', `%r12')
define(`s8', `%r13')
define(`s9', `%r14')
define(`s10', `%r15')

define(`sx', `%rax')

	TEXT

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_1)
	mov	0*8(ap), %rdx
	mulx	%rdx, s0, sx
	mov	s0, 0*8(rp)
	mov	sx, 1*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_2)
	mov	0*8(ap), %rdx
	mulx	1*8(ap), s1, s2		C a0 a1
	mulx	%rdx, s3, s4		C a0^2
	add	s1, s1
	adc	s2, s2
	mov	1*8(ap), %rdx
	mulx	%rdx, s0, sx		C a1^2
	mov	s3, 0*8(rp)
	adc	$0, sx
	add	s4, s1
	adc	s0, s2
	adc	$0, sx
	mov	s1, 1*8(rp)
	mov	s2, 2*8(rp)
	mov	sx, 3*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_3)
	mov	0*8(ap), %rdx
	mulx	1*8(ap), sx, s0		C a0 a1
	mulx	2*8(ap), s1, s2		C a0 a2

	add	s0, s1
	mov	1*8(ap), %rdx
	mulx	2*8(ap), s3, s4		C a1 a2

	mov	$0, R32(s0)
	adc	s2, s3
	adc	s0, s4

	add	sx, sx
	adc	s1, s1
	adc	s3, s3
	adc	s4, s4

	mov	0*8(ap), %rdx
	mulx	%rdx, %rdx, s2		C a0^2
	adc	R32(s0), R32(s0)
	mov	%rdx, 0*8(rp)
	add	s2, sx
	mov	sx, 1*8(rp)

	mov	1*8(ap), %rdx
	mulx	%rdx, s2, sx		C a1^2
	adc	s2, s1
	adc	sx, s3

	mov	2*8(ap), %rdx
	mulx	%rdx, s2, sx		C a2^2
	mov	s1, 2*8(rp)
	mov	s3, 3*8(rp)
	adc	s4, s2
	adc	s0, sx
	mov	s2, 4*8(rp)
	mov	sx, 5*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_4)
	mov	0*8(ap), %rdx
	push	s5
	push	s6
	xor	R32(s0), R32(s0)

	mulx	1*8(ap), s1, s4		C a0 a1
	mulx	2*8(ap), s3, s2		C a0 a2
	mulx	3*8(ap), s5, s6		C a0 a3
	mov	1*8(ap), %rdx
	adox	s4, s3
	adox	s2, s5
	mulx	2*8(ap), s4, s2		C a1 a2
	adcx	s4, s5
	adcx	s2, s6
	mulx	3*8(ap), s4, s2		C a1 a3
	mov	2*8(ap), %rdx
	adox	s4, s6
	adox	s0, s2
	mulx	3*8(ap), %rdx, s4	C a1 a3
	adc	%rdx, s2
	adc	s0, s4
	mov	0*8(ap), %rdx
	add	s1, s1
	adc	s3, s3
	adc	s5, s5
	adc	s6, s6
	adc	s2, s2
	adc	s4, s4
	setc	R8(s0)

	mulx	%rdx, %rdx, sx		C a0^2
	mov	%rdx, 0*8(rp)
	add	sx, s1
	mov	1*8(ap), %rdx
	mov	s1, 1*8(rp)
	mulx	%rdx, %rdx, sx		C a1^2
	adc	%rdx, s3
	adc	sx, s5
	mov	2*8(ap), %rdx
	mov	s3, 2*8(rp)
	mov	s5, 3*8(rp)
	mulx	%rdx, s1, sx		C a2^2
	adc	s1, s6
	adc	sx, s2
	mov	3*8(ap), %rdx
	mov	s6, 4*8(rp)
	mov	s2, 5*8(rp)
	mulx	%rdx, s1, sx		C a3^2
	adc	s1, s4
	adc	s0, sx
	mov	s4, 6*8(rp)
	mov	sx, 7*8(rp)

	pop	s6
	pop	s5

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_5)
	mov	0*8(ap), %rdx
	push	s5
	push	s6
	push	s7
	xor	R32(s0), R32(s0)

	mulx	1*8(ap), sx, s2		C a0 a1
	mulx	2*8(ap), s3, s4		C a0 a2
	mulx	3*8(ap), s5, s1		C a0 a3
	mulx	4*8(ap), s6, s7		C a0 a4
	adcx	s2, s3
	adcx	s4, s5
	adcx	s1, s6
	C x, 3, 5, 6, 7

	mov	1*8(ap), %rdx
	mulx	2*8(ap), s2, s4		C a1 a2
	adcx	s0, s7
	adcx	s2, s5
	mulx	3*8(ap), s2, s1		C a1 a3
	adcx	s4, s6
	adox	s2, s6
	mulx	4*8(ap), s2, s4		C a1 a4
	adox	s1, s7
	adcx	s2, s7
	adcx	s0, s4
	C x, 3, 5, 6, 7, 4

	mov	2*8(ap), %rdx
	mulx	3*8(ap), s1, s2		C a2 a3
	adcx	s1, s7
	adcx	s2, s4
	mulx	4*8(ap), s1, s2		C a2 a4
	adox	s1, s4
	C x, 3, 5, 6, 7, 4, 2

	mov	3*8(ap), %rdx
	mulx	4*8(ap), %rdx, s1	C a3 a4
	adox	s0, s2
	adc	%rdx, s2
	adc	s0, s1
	C x, 3, 5, 6, 7, 4, 2, 1

	mov	0*8(ap), %rdx
	add	sx, sx
	adc	s3, s3
	mov	sx, -1*8(%rsp)
	adc	s5, s5
	adc	s6, s6
	mulx	%rdx, %rdx, sx		C a0^2
	adc	s7, s7
	adc	s4, s4
	adc	s2, s2
	adc	s1, s1
	setc	R8(s0)
	C x, 3, 5, 6, 7, 4, 2, 1, 0

	mov	%rdx, 0*8(rp)
	add	-1*8(%rsp), sx
	mov	sx, 1*8(rp)

	mov	1*8(ap), %rdx
	mulx	%rdx, %rdx, sx		C a1^2
	adc	%rdx, s3
	adc	sx, s5
	mov	s3, 2*8(rp)

	mov	2*8(ap), %rdx
	mulx	%rdx, s3, sx		C a2^2
	mov	s5, 3*8(rp)
	adc	s3, s6
	adc	sx, s7
	mov	s6, 4*8(rp)

	mov	3*8(ap), %rdx
	mulx	%rdx, s3, sx		C a3^2
	mov	s7, 5*8(rp)
	adc	s3, s4
	adc	sx, s2
	mov	s4, 6*8(rp)

	mov	4*8(ap), %rdx
	mulx	%rdx, s3, sx		C a4^2
	mov	s2, 7*8(rp)
	adc	s1, s3
	adc	s0, sx
	mov	s3, 8*8(rp)
	mov	sx, 9*8(rp)

	pop	s7
	pop	s6
	pop	s5

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_6)
	mov	0*8(ap), %rdx
	push	s5
	push	s6
	push	s7
	push	s8
	push	s9
	xor	R32(s0), R32(s0)

	mulx	1*8(ap), sx, s2		C a0 a1
	mulx	2*8(ap), s3, s4		C a0 a2
	mulx	3*8(ap), s5, s1		C a0 a3
	mulx	4*8(ap), s6, s9		C a0 a4
	mulx	5*8(ap), s7, s8		C a0 a5
	adox	s2, s3
	adox	s4, s5
	adox	s1, s6
	adox	s9, s7
	adox	s0, s8
	C x, 3, 5, 6, 7, 8

	mov	1*8(ap), %rdx
	mulx	2*8(ap), s2, s4		C a1 a2
	mulx	3*8(ap), s1, s9		C a1 a3
	adcx	s2, s5
	adcx	s4, s6
	adox	s1, s6
	adox	s9, s7
	mulx	4*8(ap), s2, s4		C a1 a4
	mulx	5*8(ap), s1, s9		C a1 a5
	adcx	s2, s7
	adcx	s4, s8
	adox	s1, s8
	C x, 3, 5, 6, 7, 8, 9

	mov	2*8(ap), %rdx
	mulx	3*8(ap), s1, s2		C a2 a3
	adcx	s0, s9
	adox	s0, s9
	adcx	s1, s7
	adcx	s2, s8
	mulx	4*8(ap), s4, s1		C a2 a4
	mulx	5*8(ap), %rdx, s2	C a2 a5
	adox	s4, s8
	adox	s1, s9
	adcx	%rdx, s9
	C x, 3, 5, 6, 7, 8, 9, 2

	mov	3*8(ap), %rdx
	mulx	4*8(ap), s1, s4		C a3 a4
	adcx	s0, s2
	adcx	s1, s9
	mulx	5*8(ap), %rdx, s1	C a3 a5
	adcx	s4, s2
	adox	%rdx, s2
	adox	s0, s1
	C x, 3, 5, 6, 7, 8, 9, 2, 1

	mov	4*8(ap), %rdx
	mulx	5*8(ap), %rdx, s4	C a4 a5
	adcx	%rdx, s1
	adcx	s0, s4
	C x, 3, 5, 6, 7, 8, 9, 2, 1, 4

	add	sx, sx
	adc	s3, s3
	mov	0*8(ap), %rdx
	mov	sx, -1*8(%rsp)
	adc	s5, s5
	adc	s6, s6
	mulx	%rdx, %rdx, sx
	adc	s7, s7
	adc	s8, s8
	adc	s9, s9
	adc	s2, s2
	adc	s1, s1
	adc	s4, s4
	setc	R8(s0)
	C x, 3, 5, 6, 7, 8, 9, 2, 1, 4, 0

	mov	%rdx, 0*8(rp)
	add	-1*8(%rsp), sx
	mov	sx, 1*8(rp)

	mov	1*8(ap), %rdx
	mulx	%rdx, %rdx, sx		C a1^2
	adc	%rdx, s3
	adc	sx, s5
	mov	s3, 2*8(rp)

	mov	2*8(ap), %rdx
	mulx	%rdx, s3, sx		C a2^2
	mov	s5, 3*8(rp)
	adc	s3, s6
	adc	sx, s7
	mov	s6, 4*8(rp)

	mov	3*8(ap), %rdx
	mulx	%rdx, s3, sx		C a3^2
	mov	s7, 5*8(rp)
	adc	s3, s8
	adc	sx, s9
	mov	s8, 6*8(rp)

	mov	4*8(ap), %rdx
	mulx	%rdx, s3, sx		C a4^2
	mov	s9, 7*8(rp)
	adc	s3, s2
	adc	sx, s1
	mov	s2, 8*8(rp)

	mov	5*8(ap), %rdx
	mulx	%rdx, s3, sx		C a5^2
	mov	s1, 9*8(rp)
	adc	s4, s3
	adc	s0, sx
	mov	s3, 10*8(rp)
	mov	sx, 11*8(rp)

	pop	s9
	pop	s8
	pop	s7
	pop	s6
	pop	s5

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqr_7)
	mov	0*8(ap), %rdx
	push	s5
	push	s6
	push	s7
	push	s8
	push	s9
	push	s10
	xor	R32(s0), R32(s0)

	mulx	1*8(ap), sx, s1		C a0 a1
	mulx	2*8(ap), s3, s2		C a0 a2
	mulx	3*8(ap), s5, s4		C a0 a3
	mulx	4*8(ap), s6, s10	C a0 a4
	mov	sx, -2*8(%rsp)
	adox	s1, s3
	adox	s2, s5
	mov	s3, -1*8(%rsp)
	mulx	5*8(ap), s7, s1		C a0 a5
	mulx	6*8(ap), s8, s9		C a0 a6
	adox	s4, s6
	adox	s10, s7
	adox	s1, s8
	adox	s0, s9
	C (-, -,) 5, 6, 7, 8, 9
	C x, 1, 2, 3, 4, 10

	mov	1*8(ap), %rdx
	mulx	2*8(ap), sx, s1		C a1 a2
	mulx	3*8(ap), s2, s3		C a1 a3
	adox	sx, s5
	adox	s1, s6
	adcx	s2, s6
	adcx	s3, s7
	mulx	4*8(ap), sx, s1		C a1 a4
	mulx	5*8(ap), s2, s3		C a1 a5
	mulx	6*8(ap), s4, s10	C a1 a6
	adox	sx, s7
	adox	s1, s8
	adcx	s2, s8
	adcx	s3, s9
	adox	s4, s9
	adox	s0, s10
	adcx	s0, s10
	C (-, -,) 5, 6, 7, 8, 9, 10
	C x, 1, 2, 3, 4

	mov	2*8(ap), %rdx
	mulx	3*8(ap), sx, s1		C a2 a3
	mulx	4*8(ap), s2, s3		C a2 a4
	adox	sx, s7
	adox	s1, s8
	adcx	s2, s8
	adcx	s3, s9
	mulx	5*8(ap), sx, s1		C a2 a5
	mulx	6*8(ap), s2, s3		C a2 a6
	adox	sx, s9
	adox	s1, s10
	adcx	s2, s10
	adcx	s0, s3
	adox	s0, s3
	C (-, -,) 5, 6, 7, 8, 9, 10, 3
	C x, 1, 2, 4

	mov	3*8(ap), %rdx
	mulx	4*8(ap), sx, s1		C a3 a4
	mulx	5*8(ap), s2, s4		C a3 a5
	adox	sx, s9
	adox	s1, s10
	mulx	6*8(ap), sx, s1		C a3 a6
	adcx	s2, s10
	adcx	s4, s3
	adox	sx, s3
	adox	s0, s1
	C (-, -,) 5, 6, 7, 8, 9, 10, 3, 1
	C x, 2, 4

	mov	4*8(ap), %rdx
	mulx	5*8(ap), sx, s2		C a4 a5
	mulx	6*8(ap), %rdx, s4	C a4 a6
	adox	sx, s3
	adox	s2, s1
	adcx	%rdx, s1
	adox	s0, s4
	C (-, -,) 5, 6, 7, 8, 9, 10, 3, 1, 4
	C x, 2

	mov	5*8(ap), %rdx
	mulx	6*8(ap), sx, s2		C a5 a6
	adc	sx, s4
	adc	s0, s2
	C (-, -,) 5, 6, 7, 8, 9, 10, 3, 1, 4, 2
	C x

	shlq	-2*8(%rsp)
	mov	-1*8(%rsp), sx
	mov	0*8(ap), %rdx
	adc	sx, sx
	adc	s5, s5
	mov	sx, -1*8(%rsp)
	adc	s6, s6
	adc	s7, s7
	mulx	%rdx, %rdx, sx
	adc	s8, s8
	adc	s9, s9
	adc	s10, s10
	adc	s3, s3
	adc	s1, s1
	adc	s4, s4
	adc	s2, s2
	setc	R8(s0)
	C (-, -,) 5, 6, 7, 8, 9, 10, 3, 1, 4, 2, 0
	C x

	add	-2*8(%rsp), sx
	mov	%rdx, 0*8(rp)
	mov	sx, 1*8(rp)

	mov	1*8(ap), %rdx
	mulx	%rdx, %rdx, sx		C a1^2
	adc	-1*8(%rsp), %rdx
	mov	%rdx, 2*8(rp)
	adc	sx, s5

	mov	2*8(ap), %rdx
	mulx	%rdx, %rdx, sx		C a2^2
	mov	s5, 3*8(rp)
	adc	%rdx, s6
	adc	sx, s7
	mov	s6, 4*8(rp)

	mov	3*8(ap), %rdx
	mulx	%rdx, s5, sx		C a3^2
	mov	s7, 5*8(rp)
	adc	s5, s8
	adc	sx, s9
	mov	s8, 6*8(rp)

	mov	4*8(ap), %rdx
	mulx	%rdx, s5, sx		C a4^2
	mov	s9, 7*8(rp)
	adc	s5, s10
	adc	sx, s3
	mov	s10, 8*8(rp)

	mov	5*8(ap), %rdx
	mulx	%rdx, s5, sx		C a5^2
	mov	s3, 9*8(rp)
	adc	s5, s1
	adc	sx, s4
	mov	s1, 10*8(rp)

	mov	6*8(ap), %rdx
	mulx	%rdx, s3, sx		C a6^2
	mov	s4, 11*8(rp)
	adc	s2, s3
	adc	s0, sx
	mov	s3, 12*8(rp)
	mov	sx, 13*8(rp)

	pop	s10
	pop	s9
	pop	s8
	pop	s7
	pop	s6
	pop	s5

	ret
EPILOGUE()
