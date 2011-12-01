/*
 * Matrix op worker threads
 *
 * Copyright (C) 2011 Steve Kerrison <github@stevekerrison.com>
 *
 * This software is freely distributable under a derivative of the
 * University of Illinois/NCSA Open Source License posted in
 * LICENSE.txt and at <http://github.xcore.com/>
 */

#include "matrix_worker.h"
#include <stdlib.h>

/* If rand() isn't 32-bit let's fix it. */
#if RAND_MAX < 0xffff
#define myrand() (((rand() << 24) & 0xff000000) | ((rand() << 16) & 0x00ff000)) \
	((rand() << 8) & 0x0000ff00)) | (rand() & 0x000000ff))
#elif RAND_MAX < 0xffffffff
#define myrand() (((rand() << 16) & 0xffff0000) | (rand() & 0xffff))
#else
#define myrand() rand()
#endif

void matrix_mul_worker(int ptA, int ptDimA, int ptB, int ptDimB, int ptC,
	short startC, short lenC, int ptOps)
{
	int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps;
	short *dimA = (short *)ptDimA, *dimB = (short *)ptDimB;
	int col = startC % dimB[1], bstep = dimB[1], row = startC / dimB[1], elim = dimA[1];
	int result_hi = 0; //We throw this guy away.
	for (int dst = startC; dst < startC + lenC; dst += 1)
	{
		int result = 0;
		for (int e = 0; e < elim; e += 1)
		{
			asm("maccu %1, %0, %4, %5" :
				"=r"(result), "=r"(result_hi) :
				"0"(result),"1"(result_hi),
					"r"(A[row * elim + e]),
					"r"(B[bstep * e + col]));
		}
		C[dst] = result;
		col += 1;
		if (col == bstep)
		{
			col = 0;
			row += 1;
		}
	}
	*ops = lenC * (dimA[1] + 1);
	return;
}

void matrix_arr_worker_add(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len)
{
		int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = A[base] + B[base];
	}
	*ops = len;
	return;
}

void matrix_arr_worker_sub(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len)
{
		int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = A[base] - B[base];
	}
	*ops = len;
	return;
}

void matrix_arr_worker_mul(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len)
{
		int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = A[base] * B[base];
	}
	*ops = len;
	return;
}

void matrix_arr_worker_div(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len)
{
		int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = A[base] / B[base];
	}
	*ops = len;
	return;
}

void matrix_arr_worker_udiv(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len)
{
		int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = (unsigned)A[base] / (unsigned)B[base];
	}
	*ops = len;
	return;
}

void matrix_sca_worker_add(int ptA, int S, int ptC,
	int ptOps, short offset, short len)
{
	int *A = (int *)ptA, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = A[base] + S;
	}
	*ops = len;
	return;
}

void matrix_sca_worker_sub(int ptA, int S, int ptC,
	int ptOps, short offset, short len)
{
	int *A = (int *)ptA, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = A[base] - S;
	}
	*ops = len;
	return;
}

void matrix_sca_worker_mul(int ptA, int S, int ptC,
	int ptOps, short offset, short len)
{
	int *A = (int *)ptA, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = A[base] * S;
	}
	*ops = len;
	return;
}

void matrix_sca_worker_div(int ptA, int S, int ptC,
	int ptOps, short offset, short len)
{
	int *A = (int *)ptA, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = A[base] / S;
	}
	*ops = len;
	return;
}

void matrix_sca_worker_udiv(int ptA, int S, int ptC,
	int ptOps, short offset, short len)
{
	int *A = (int *)ptA, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = (unsigned)A[base] / (unsigned)S;
	}
	*ops = len;
	return;
}

void matrix_sca_worker_rand(int ptC, int ptOps, short offset, short len)
{
	int *C = (int *)ptC, *ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = myrand();
	}
	*ops = len;
	return;
}

void matrix_sca_worker_set(int S, int ptC, int ptOps, short offset, short len)
{
	int *C = (int *)ptC, *ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = S;
	}
	*ops = len;
	return;
}

void matrix_sca_worker_shr(int ptA, int S, int ptC, int ptOps, short offset, short len)
{
	int *A = (int *)ptA, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = (unsigned)A[base] >> S;
	}
	*ops = len;
	return;
}

void matrix_sca_worker_ashr(int ptA, int S, int ptC, int ptOps, short offset, short len)
{
	int *A = (int *)ptA, *C = (int *)ptC,
		*ops = (int *)ptOps, base;
	for (base = offset; base < offset + len; base += 1)
	{
		C[base] = A[base] >> S;
	}
	*ops = len;
	return;
}
