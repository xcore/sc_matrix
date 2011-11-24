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

void matrix_mul_worker_new(int ptA, int ptDimA, int ptB, int ptDimB, int ptC,
	short startC, short lenC, int ptOps)
{
	int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps;
	short *dimA = (short *)ptDimA, *dimB = (short *)ptDimB;
	int col = startC % dimB[1], row = startC / dimB[1];
	for (int dst = startC; dst < startC + lenC; dst += 1)
	{
		C[dst] = 0;
		for (int e = 0; e < dimA[1]; e += 1)
		{
			C[dst] += A[row * dimA[1] + e] * B[dimB[1] * e + col];
		}
		col += 1;
		if (col == dimB[1])
		{
			col = 0;
			row += 1;
		}
	}
	/*for (int c = colB; c < colB + lenB; c+= 1)
	{
		for (int r = rowA; r < rowA + lenA; r += 1)
		{
			int dst = r * dimB[1] + c;
			C[dst] = 0;
			for (int e = 0; e < dimA[1]; e += 1)
			{
				C[dst] += A[r * dimA[1] + e] * B[c + dimB[1] * e];
			}
		}
	}*/
	*ops = lenC * (dimA[1] + 1);
	return;
}

void matrix_mul_worker(int ptA, int ptDimA, int ptB, int ptDimB, int ptC,
	int ptOps, char nThreads, char offset)
{
	int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps;
	short *dimA = (short *)ptDimA, *dimB = (short *)ptDimB;
	int p,r,c,lops = 0;
	ops += offset;
	int rlim = dimA[0], clim = dimB[1], plim = dimA[1];
	for (r = 0; r < rlim; r++)
	{
		for (c = offset; c < clim; c += nThreads)
		{
			C[r * rlim + c] = 0;
			for (p = 0; p < plim; p += 1)
			{
				C[r * rlim + c] += A[r * plim + p] * B[p * clim + c];
				lops += 1;
			}
		}
	}
	*ops = lops;
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

