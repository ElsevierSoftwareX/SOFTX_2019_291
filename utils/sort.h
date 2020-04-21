/**
* @file sort.h
* @brief Binary sorting algorithms.
*
* @author Manfred Liebmann
* @version 
* @date 2016-12-13
*/


/**
* Private procedure, see binary_sort.
*/
template<class T> inline
void _binary_sort(T *inP, T *inQ, T s)
{
	T *P = inP, *Q = inQ;
	T p=0, q=0;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		while(((p = P[0]) & s) == 0) P++;
		while(((q = Q[-1]) & s) != 0) Q--;
	}
	s >>= 1;
	if(s)
	{
		if(inQ - P > 1) _binary_sort(P, inQ, s);
		if(Q - inP > 1) _binary_sort(inP, Q, s);
	}
}

/**
* Private procedure, see binary_sort_copy.
*/
template<class T, class S> inline
void _binary_sort_copy(T *inP, T *inQ, S *inU, S* inV, T s)
{
	T *P = inP, *Q = inQ;
	S *U = inU, *V = inV;
	T p=0, q=0;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		U++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		V--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		S u, v;
		u = *U;
		v = *--V;
		*U++ = v;
		*V = u;
		while(((p = P[0]) & s) == 0) P++, U++;
		while(((q = Q[-1]) & s) != 0) Q--, V--;
	}
	s >>= 1;
	if(s)
	{
		if(inQ - P > 1) _binary_sort_copy(P, inQ, U, inV, s);
		if(Q - inP > 1) _binary_sort_copy(inP, Q, inU, V, s);
	}
}

/**
* Private procedure, see binary_sort_copy_copy.
*/
template<class T, class S, class R> inline
void _binary_sort_copy_copy(T *inP, T *inQ, S *inA, S *inB, R *inU, R *inV, T s)
{
	T *P = inP, *Q = inQ;
	S *A = inA, *B = inB;
	R *U = inU, *V = inV;
	T p=0, q=0;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		A++;
		U++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		B--;
		V--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		S u, v;
		u = *U;
		v = *--V;
		*U++ = v;
		*V = u;
		S a, b;
		a = *A;
		b = *--B;
		*A++ = b;
		*B = a;
		while(((p = P[0]) & s) == 0) P++, A++, U++;
		while(((q = Q[-1]) & s) != 0) Q--, B--, V--;
	}
	s >>= 1;
	if(s)
	{
		if(inQ - P > 1) _binary_sort_copy_copy(P, inQ, A, inB, U, inV, s);
		if(Q - inP > 1) _binary_sort_copy_copy(inP, Q, inA, B, inU, V, s);
	}
}

/**
* Private procedure, see binary_sort_sort.
*/
template<class T> inline
void _binary_sort_sort(T *inP, T *inQ, T *inA, T *inB, T s, T t)
{
	T *P = inP, *Q = inQ, *A = inA, *B = inB;
	T p=0, q=0;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		A++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		B--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		T a, b;
		a = *A;
		b = *--B;
		*A++ = b;
		*B = a;
		while(((p = P[0]) & s) == 0) P++, A++;
		while(((q = Q[-1]) & s) != 0) Q--, B--;
	}
	s >>= 1;
	if(s)
	{
		if(inQ - P > 1) _binary_sort_sort(P, inQ, A, inB, s, t);
		if(Q - inP > 1) _binary_sort_sort(inP, Q, inA, B, s, t);
	}
	else if(t)
	{
		if(inQ - P > 1) _binary_sort_sort(A, inB, P, inQ, t, s);
		if(Q - inP > 1) _binary_sort_sort(inA, B, inP, Q, t, s);
	}
}

/**
* Private procedure, see binary_sort_sort_copy.
*/
template<class T, class S> inline
void _binary_sort_sort_copy(T *inP, T *inQ, T *inA, T *inB, S *inU, S* inV, T s, T t)
{
	T *P = inP, *Q = inQ, *A = inA, *B = inB;
	S *U = inU, *V = inV;
	T p=0, q=0;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		A++;
		U++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		B--;
		V--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		T a, b;
		a = *A;
		b = *--B;
		*A++ = b;
		*B = a;
		S u, v;
		u = *U;
		v = *--V;
		*U++ = v;
		*V = u;
		while(((p = P[0]) & s) == 0) P++, A++, U++;
		while(((q = Q[-1]) & s) != 0) Q--, B--, V--;
	}
	s >>= 1;
	if(s)
	{
		if(inQ - P > 1) _binary_sort_sort_copy(P, inQ, A, inB, U, inV, s, t);
		if(Q - inP > 1) _binary_sort_sort_copy(inP, Q, inA, B, inU, V, s, t);
	}
	else if(t)
	{
		if(inQ - P > 1) _binary_sort_sort_copy(A, inB, P, inQ, U, inV, t, s);
		if(Q - inP > 1) _binary_sort_sort_copy(inA, B, inP, Q, inU, V, t, s);
	}
}

/**
* Private procedure, see fractal_sort_sort_copy.
*/
template<class T, class S> inline
void _fractal_sort_sort_copy(T *inP, T *inQ, T *inA, T *inB, S *inU, S* inV, T s, T t)
{
	T *P = inP, *Q = inQ, *A = inA, *B = inB;
	S *U = inU, *V = inV;
	T p=0, q=0;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		A++;
		U++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		B--;
		V--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		T a, b;
		a = *A;
		b = *--B;
		*A++ = b;
		*B = a;
		S u, v;
		u = *U;
		v = *--V;
		*U++ = v;
		*V = u;
		while(((p = P[0]) & s) == 0) P++, A++, U++;
		while(((q = Q[-1]) & s) != 0) Q--, B--, V--;
	}
	s >>= 1;
	if(t != 0)
	{
		if(inQ - P > 1) _fractal_sort_sort_copy(A, inB, P, inQ, U, inV, t, s);
		if(Q - inP > 1) _fractal_sort_sort_copy(inA, B, inP, Q, inU, V, t, s);
	}
}

/**
* Private procedure, see binary_sort.
*/
template<class T> inline
T _binary_log(const T *P, const T *Q)
{
	T s = 0;
	while(P != Q)
	{
		s |= *P++;
	}
	T t = ~0;
	while(s & t)
	{
		s &= t;
		t <<= 1;
	}
//	if(s < 0)
//	{
//		cout << "SORT ERROR!" << endl;
//		s = 0;
//	}
	return s;
}


/**
* The binary_sort procedure sorts a vector of nonnegative integers in place in ascending order.
* \param inV Input: Vector of nonnegative integers to sort. Output: Sorted vector in ascending order.
*/
template<class T> inline
void binary_sort(mt_vector<T> &inV)
{
	if(inV.size() < 2) return;
	_binary_sort(&inV[0], &inV[0]+inV.size(), _binary_log(&inV[0], &inV[0]+inV.size()));
}

/**
* The binary_sort_copy procedure partially sorts pairs in place in ascending order only looking at the first argument.
* \param inV Input: Vector of nonnegative integers to sort, first argument. Output: Partially sorted vector in ascending order, first argument.
* \param inW Input: Vector of arbitrary type to sort, second argument. Output: Partially sorted vector in ascending order, second argument.
*/
template<class T, class S> inline
void binary_sort_copy(mt_vector<T> &inV, mt_vector<S> &inW)
{
	if(inV.size() < 2) return;
	_binary_sort_copy(&inV[0], &inV[0]+inV.size(), &inW[0], &inW[0]+inW.size(), _binary_log(&inV[0], &inV[0]+inV.size()));
}

/**
* The binary_sort_copy_copy procedure partially sorts triples in place in ascending order only looking at the first argument.
* \param inV Input: Vector of nonnegative integers to sort, first argument. Output: Partially sorted vector in ascending order, first argument.
* \param inW Input: Vector of arbitrary type to sort, second argument. Output: Partially sorted vector in ascending order, second argument.
* \param inA Input: Vector of arbitrary type to sort, third argument. Output: Partially sorted vector in ascending order, second argument.
*/
template<class T, class S, class R> inline
void binary_sort_copy_copy(mt_vector<T> &inV, mt_vector<S> &inW, mt_vector<R> &inA)
{
	if(inV.size() < 2) return;
	_binary_sort_copy_copy(&inV[0], &inV[0]+inV.size(), &inW[0], &inW[0]+inW.size(), &inA[0], &inA[0]+inA.size(), _binary_log(&inV[0], &inV[0]+inV.size()));
}

/**
* The binary_sort_sort procedure sorts pairs of nonnegative integers in place in lexicographic order looking at both arguments.
* \param inV Input: Vector of nonnegative integers to sort, first argument. Output: Sorted vector in ascending order, first argument.
* \param inW Input: Vector of nonnegative integers to sort, second argument. Output: Sorted vector in ascending order, second argument.
*/
template<class T> inline
void binary_sort_sort(mt_vector<T> &inV, mt_vector<T> &inW)
{
	if(inV.size() < 2) return;
	_binary_sort_sort(&inV[0], &inV[0]+inV.size(), &inW[0], &inW[0]+inW.size(), 
		_binary_log(&inV[0], &inV[0]+inV.size()), _binary_log(&inW[0], &inW[0]+inW.size()));
}

/**
* The binary_sort_sort_copy procedure partially sorts triples in place in lexicographic order only looking at the first and second argument.
* \param inV Input: Vector of nonnegative integers to sort, first argument. Output: Partially sorted vector in ascending order, first argument.
* \param inW Input: Vector of nonnegative integers to sort, second argument. Output: Partially sorted vector in ascending order, second argument.
* \param inA Input: Vector of arbitrary type to sort, third argument. Output: Partially sorted vector in ascending order, third argument.
*/
template<class T, class S> inline
void binary_sort_sort_copy(mt_vector<T> &inV, mt_vector<T> &inW, mt_vector<S> &inA)
{
	if(inV.size() < 2) return;
	_binary_sort_sort_copy(&inV[0], &inV[0]+inV.size(), &inW[0], &inW[0]+inW.size(), &inA[0], &inA[0]+inA.size(), 
		_binary_log(&inV[0], &inV[0]+inV.size()), _binary_log(&inW[0], &inW[0]+inW.size()));
}

/**
* The fractal_sort_sort_copy procedure partially sorts triples in place in fractal order only looking at the first and second argument.
* \param inV Input: Vector of nonnegative integers to sort, first argument. Output: Partially sorted vector in ascending order, first argument.
* \param inW Input: Vector of nonnegative integers to sort, second argument. Output: Partially sorted vector in ascending order, second argument.
* \param inA Input: Vector of arbitrary type to sort, third argument. Output: Partially sorted vector in ascending order, third argument.
*/
template<class T, class S> inline
void fractal_sort_sort_copy(mt_vector<T> &inV, mt_vector<T> &inW, mt_vector<S> &inA)
{
	if(inV.size() < 2) return;
	T log = max(_binary_log(&inV[0], &inV[0]+inV.size()), _binary_log(&inW[0], &inW[0]+inW.size()));
	_fractal_sort_sort_copy(&inV[0], &inV[0]+inV.size(), &inW[0], &inW[0]+inW.size(), &inA[0], &inA[0]+inA.size(), log, log);
}

/**
* The bucket_sort_count procedure counts how often an index appears in the first vector and stores the counts in the second.
* \param inU Input: Vector of indices.
* \param inW Input: Vector initialized to zero. Output: Sorted Vector of index counts.
*/
template<class T, class S> inline
void bucket_sort_count(const mt_vector<T> &inU, mt_vector<S> &inW)
{
	const T *U = &inU[0], *V = &inU[0] + inU.size();
	S *W = &inW[0];

	while(U != V)
	{
		W[*U++]++;
	}
}

/**
* The bucket_sort_size function calculates the sum over all elements in the vector.
* \param inU Input: Vector of counts.
* \return Output: Sum over all elements in the vector.
*/
template<class T> inline
T bucket_sort_size(const mt_vector<T> &inU)
{
  if(inU.size() == 0) return 0;
  const T *U = &inU[0], *V = &inU[0] + inU.size();

  T t = 0;
  while(U != V)
  {
    t += *U++;
  }
  return t;
}

/**
* The bucket_sort_offset procedure calculates the displacements from a vector of counts.
* \param inU Input: Vector of counts.
* \param inW Output: Vector of displacements.
*/
template<class T> inline
void bucket_sort_offset(const mt_vector<T> &inU, mt_vector<T> &inW)
{
  if(inU.size() == 0) return;
  inW.resize(inU.size());

  const T *U = &inU[0], *V = &inU[0] + inU.size();
  T *W = &inW[0];

  T t = 0;
  while(U != V)
  {
    T s = *U++;
    *W++ = t;
    t += s;
  }
}

/**
* The bucket_sort_copy procedure partially sorts pairs in ascending order only looking at the first argument and with partial output.
* \param inU Input: Vector to sort, first argument.
* \param inA Input: Vector to sort, second argument.
* \param inB Output: Sorted vector in ascending order, second argument.
* \param inW Input: Vector of displacements. Output: Vector of final displacements.
* \param _n Input: Size of second argument.
*/
template<class T, class S, class C> inline
void bucket_sort_copy(const mt_vector<T> &inU, const mt_vector<C> &inA,  mt_vector<C> &inB, mt_vector<S> inW, int _n)
{
	const T *U = &inU[0], *V = &inU[0] + inU.size();
	const C *A = &inA[0];
	C *B = &inB[0];
	S *W = &inW[0];

	while(U != V)
	{
		T t = *U++;
		S s = W[t];
		W[t]++;
		B = &inB[_n * s];
		for(int i = 0; i < _n; i++)
		{
			*B++ = *A++;
		}
	}
}

/**
* The unique_resize procedure calculates the set union of an in ascending order sorted vector.
* \param inP Input: Vector sorted in ascending order. Output: Vector sorted in ascending order with only unique elements.
*/
template<class T> inline
void unique_resize(mt_vector<T> &inP)
{
	if(inP.size() < 2) return;
	T *P = &inP[0], *Q = &inP[0] + inP.size();

	if(P != Q)
	{
		T* R = P;
		++P;
		while(P != Q)
		{
			if ((*R != *P))
			{
				*++R = *P;
			}
			++P;
		}
		++R;
		inP.resize(int(R - &inP[0]));
	}
}

/**
* The unique_resize procedure calculates the set union of an in ascending order sorted vector.
* \param inP Input: Vector sorted in ascending order. Output: Vector sorted in ascending order with only unique elements.
* \param inU Input: Vector sorted in ascending order. Output: Vector sorted in ascending order with only unique elements.
*/
template<class T> inline
void unique_resize(mt_vector<T> &inP, mt_vector<T> &inU)
{
	if(inP.size() < 2) return;
	T *P = &inP[0], *Q = &inP[0] + inP.size();
	T *U = &inU[0];

	if(P != Q)
	{
		T* R = P;
		T* W = U;
		++P; ++U;
		while(P != Q)
		{
			if ((*R != *P) || (*W != *U))
			{
				*++R = *P;
				*++W = *U;
			}
			++P; ++U;
		}
		++R;
		++W;
		inP.resize(int(R - &inP[0]));
		inU.resize(int(W - &inU[0]));
	}
}

/**
* The unique_accumulate procedure calculates a partial set union of pairs only looking at the first argument and accumulates the values of the second argument of partially matching pairs.
* \param inP Input: Vector sorted in ascending order, first argument. Output: Vector sorted in ascending order with only unique elements, first argument.
* \param inA Input: Vector to accumulate, second argument. Output: Vector of accumulated values, second argument.
*/
template<class T, class S> inline
void unique_accumulate(mt_vector<T> &inP, mt_vector<S> &inA)
{
	if(inP.size() < 2) return;
	T *P = &inP[0], *Q = &inP[0] + inP.size();
	S *A = &inA[0];

	if(P != Q)
	{
		T* R = P;
		S* C = A;
		++P; ++A;
		while(P != Q)
		{
			if (*R == *P)
			{
				*C += *A;
			}
			else
			{
				*++R = *P;
				*++C = *A;
			}
			++P; ++A;
		}
		++R;
		++C;
		inP.resize(int(R - &inP[0]));
		inA.resize(int(C - &inA[0]));
	}
}

/**
* The unique_accumulate procedure calculates a partial set union of triples only looking at the first and second argument and accumulates the values of the third argument of partially matching triples.
* \param inP Input: Vector sorted in lexicographic order, first argument. Output: Vector sorted in lexicographic order with only unique elements, first argument.
* \param inU Input: Vector sorted in lexicographic order, second argument. Output: Vector sorted in lexicographic order with only unique elements, second argument.
* \param inA Input: Vector to accumulate, third argument. Output: Vector of accumulated values, third argument.
*/
template<class T, class S> inline
void unique_accumulate(mt_vector<T> &inP, mt_vector<T> &inU, mt_vector<S> &inA)
{
	if(inP.size() < 2) return;
	T *P = &inP[0], *Q = &inP[0] + inP.size();
	T *U = &inU[0];
	S *A = &inA[0];

	if(P != Q)
	{
		T* R = P;
		T* W = U;
		S* C = A;
		++P; ++U; ++A;
		while(P != Q)
		{
			if ((*R == *P) && (*W == *U))
			{
				*C += *A;
			}
			else
			{
				*++R = *P;
				*++W = *U;
				*++C = *A;
			}
			++P; ++U; ++A;
		}
		++R;
		++W;
		++C;
		inP.resize(int(R - &inP[0]));
		inU.resize(int(W - &inU[0]));
		inA.resize(int(C - &inA[0]));
	}
}
/**
* The global_intersection procedure calculates the intersection of multiple index sets.
* \param _size Input: Total number of processes.
* \param _rank Input: Rank of the current process.
* \param inP Input: Vector of global node numbers of the current process to intersect.
* \param inU Input: Vector of local node numbers of the current process.
* \param inA Input: Vector of global node numbers for all processes to intersect. Output: Vector of intersecting sets with local node numbers for all processes.
* \param inC Input: Vector of initial counts. Output: Vector of final counts for intersecting sets.
* \param inD Input: Vector of initial displacements. Output: Vector of final displacements for intersecting sets.
*/
template<class T> inline
void global_intersection(const int _size, const int _rank, const mt_vector<T> &inP, const mt_vector<T> &inU, mt_vector<T> &inA, mt_vector<T> &inC, mt_vector<T> &inD)
{
	const T *P = &inP[0], *U = &inU[0];
	T *A = &inA[0], *C = &inC[0], *D = &inD[0];

	int j = 0, k, l;
	int m = (int)inP.size(), n;
	for(int i = 0; i < _size; i++)
	{
		k = 0, l = D[i], n = l + C[i];

		D[i] = j;
		if(_rank != i)
		{
			while((k < m) && (l < n))
			{
				if (P[k] < A[l])
				{
					k++;
				}
				else if (A[l] < P[k])
				{
					l++;
				}
				else
				{
					A[j] = U[k];
					j++;
					k++;
					l++;
				}
			}
		}
		C[i] = j - D[i];
	}
	inA.resize(j);
}
