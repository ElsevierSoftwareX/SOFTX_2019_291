/**
* @file dense_mat.hpp
* @brief Dense matrix class and associated funcs.
* @author Aurel Neic
* @version 
* @date 2017-10-31
*/
#ifndef DENSE_MAT_H
#define DENSE_MAT_H

#include <utility>
#include <assert.h>

/**
* @brief Dense matrix class.
*
* @tparam S  Floating point type.
*/
template<class S>
class dmat
{
  private:
  S*     _data;   ///< the data array
  short  _rows;   ///< number of rows
  short  _cols;   ///< number of cols
  short* _ps;     ///< row permutation for LU decomp

  public:
  /// empty constructor
  dmat() : _data(NULL), _rows(0), _cols(0), _ps(NULL)
  {}
  /// constructor that initializes the dimensions
  dmat(const short irows, const short icols) : _data(NULL), _rows(0), _cols(0), _ps(NULL)
  {
    set_size(irows, icols);
  }
  /// constructor that deep-copies a given dmat
  dmat(const dmat<S> & m) : _data(NULL), _rows(0), _cols(0), _ps(NULL)
  {
    this->assign(m);
  }
  /// destructor
  ~dmat()
  {
    if(_data)
      delete [] _data;
    if(_ps)
      delete [] _ps;
  }
  /// set the matrix dimensions
  inline void set_size(const short irows, const short icols)
  {
    if(irows*icols > _rows*_cols) {
      if(_data) delete [] _data;
      _data = new S[irows*icols];
    }
    _rows = irows;
    _cols = icols;
  }

  inline dmat<S>& operator= (dmat<S>&& m)
  {
    std::swap(_data, m._data);
    std::swap(_rows, m._rows);
    std::swap(_cols, m._cols);
    return *this;
  }

  inline dmat<S>& operator= (const dmat<S>& m)
  {
    this->set_size(m._rows, m._cols);
    for(int idx=0; idx < _rows*_cols; idx++)
      _data[idx] = m._data[idx];

    return *this;
  }
  /// [] operator returns the pointer to the beginning of a given row
  inline const S* operator[] (short ridx) const
  {
    return (_data +(ridx*_cols));
  }
  /// [] operator returns the pointer to the beginning of a given row
  inline S* operator[] (short ridx)
  {
    return (_data +(ridx*_cols));
  }
 
  inline void operator+= (const dmat<S>& m)
  {
    if ((_rows == 0)||(_cols == 0)) assign(m._rows, m._cols, S());
    assert (_rows == m._rows && _cols == m._cols);
    for(int i=0; i<_rows*_cols; i++) _data[i] += m._data[i];
  }

  inline void operator-= (const dmat<S>& m)
  {
    if ((_rows == 0)||(_cols == 0)) assign(m._rows, m._cols, S());
    assert (_rows == m._rows && _cols == m._cols);
    for(int i=0; i<_rows*_cols; i++) _data[i] -= m._data[i];
  }

  inline void operator*= (const S s)
  {
    for(int i=0; i<_rows*_cols; i++) _data[i] *= s;
  }

  inline void operator/= (const S s)
  {
    for(int i=0; i<_rows*_cols; i++) _data[i] /= s;
  }

  /// copy a mtrix.
  inline void assign(const dmat<S> & m)
  {
    set_size(m.rows(), m.cols());
    const S* ele = m[0];

    for(short i=0; i<_rows*_cols; i++) _data[i] = ele[i];
  }
  /// set all entries to a value
  inline void assign(const S v)
  {
    for(short i=0; i<_rows*_cols; i++) _data[i] = v;
  }
  /// set all entries to a value
  inline void assign(const S *v)
  {
    for(short i=0; i<_rows*_cols; i++) _data[i] = v[i];
  }
  /// resize and set all entries to a value
  inline void assign(const short irows, const short icols, const S v)
  {
    this->set_size(irows, icols);
    this->assign(v);
  }
  /// resize and set all entries to a value
  inline void assign(const short irows, const short icols, const S *v)
  {
    this->set_size(irows, icols);
    this->assign(v);
  }

  inline short rows() const
  {
    return _rows;
  }
  inline short cols() const
  {
    return _cols;
  }
  /// mat-mat multiplication
  inline void mult(const dmat<S> & in, dmat<S> & out) const
  {
    out.set_size(_rows, in.cols());
    out.assign(S(0));

    for(short i=0; i<_rows; i++) {
      for(short j=0; j<in.cols(); j++) {
        for(short k=0; k<_cols; k++)
          out[i][j] += _data[i*_cols+k] * in[k][j];
      }
    }
  }
  /// mat-vec multiplication
  inline void mult(const S* in, S* out)
  {
    for(short i=0; i<_rows; i++) {
      out[i] = S(0);
      for(short j=0; j<_cols; j++)
        out[i] += _data[i*_cols+j] * in[j];
    }
  }

  inline void transpose()
  {
    if(_cols == _rows) {
      // for quadratic matrices we dont need to resize
      // therefore we swap the entries in-place
      for(short i=0; i<_rows; i++) {
        for(short j=i+1; j<_cols; j++)
        {
          S up = _data[i*_cols+j];
          S lo = _data[j*_cols+i];
          _data[j*_cols+i] = up;
          _data[i*_cols+j] = lo;
        }
      }
    }
    else {
      // since we need to resize, we create a local copy
      // and then resize and copy
      dmat<S> t(*this);
      set_size(t.cols(), t.rows());
      for(short i=0; i<t.rows(); i++) {
        for(short j=0; j<t.cols(); j++)
          _data[j*_cols+i] = t[i][j];
      }
    }
  }

  inline void disp(const char* name)
  {
    printf("\n%s\n", name);
    for(short i=0; i<_rows; i++) {
      for(short j=0; j<_cols; j++)
        printf("  %.2f ", _data[i*_cols+j]);
      printf("\n");
    }
  }

  inline bool lu_decomp()
  {
    assert(_rows == _cols);
    const short n = _rows;
    dmat<S> & lu = *this;
    _ps = new short[n];

    S* scales = new S[n];
    S pivot, biggest, mlt, tempf;
    short pivotindex = 0;
    short i, j, k;
    // double d = 1.0; // No row interchanges yet.

    // For each row.
    for (i = 0; i < n; i++) {
      // Find the largest element in each row for row equilibration
      biggest = 0.0;
      for (j = 0; j < n; j++)
        if (biggest < (tempf = fabs(lu[i][j])))
          biggest  = tempf;
      if (biggest != 0.0)
        scales[i] = 1.0 / biggest;
      else {
        scales[i] = 0.0;
        return false; // Zero row: singular matrix.
      }
      _ps[i] = i; // Initialize pivot sequence.
    }

    // For each column.
    for (k = 0; k < n - 1; k++) {
      // Find the largest element in each column to pivot around.
      biggest = 0.0;
      for (i = k; i < n; i++) {
        if (biggest < (tempf = fabs(lu[_ps[i]][k]) * scales[_ps[i]])) {
          biggest = tempf;
          pivotindex = i;
        }
      }

      if (biggest == 0.0)
        return false; // Zero column: singular matrix.

      // Update pivot sequence.
      if (pivotindex != k) {
        j = _ps[k];
        _ps[k] = _ps[pivotindex];
        _ps[pivotindex] = j;
        // d = -(d); // ...and change the parity of d.
      }

      // Pivot, eliminating an extra variable each time
      pivot = lu[_ps[k]][k];
      for (i = k + 1; i < n; i++) {
        lu[_ps[i]][k] = mlt = lu[_ps[i]][k] / pivot;

        if (mlt != 0.0) {
          for (j = k + 1; j < n; j++)
            lu[_ps[i]][j] -= mlt * lu[_ps[k]][j];
        }
      }
    }

    delete [] scales;

    // (lu[ps[n + N - 1]][n + N - 1] == 0.0) ==> A is singular.
    return lu[_ps[n - 1]][n - 1] != 0.0;
  }

  void lu_solve(S* rhs)
  {
    short i, j;
    const short n = _rows;
    dmat<S> & lu = *this;

    S *X = new S[n], dot_prod;

    for (i = 0; i < n; i++) X[i] = 0.0;

    // Vector reduction using U triangular matrix.
    for (i = 0; i < n; i++) {
      dot_prod = 0.0;
      for (j = 0; j < i; j++) dot_prod += lu[_ps[i]][j] * X[j];
      X[i] = rhs[_ps[i]] - dot_prod;
    }

    // Back substitution, in L triangular matrix.
    for (i = n - 1; i >= 0; i--) {
      dot_prod = 0.0;
      for (j = i + 1; j < n; j++) dot_prod += lu[_ps[i]][j] * X[j];
      X[i] = (X[i] - dot_prod) / lu[_ps[i]][i];
    }

    for (i = 0; i < n; i++) rhs[i] = X[i];

    delete [] X;
  }
};

template<class S>
dmat<S> operator* (const dmat<S> & a, const dmat<S> & b)
{
  dmat<S> r;
  a.mult(b, r);
  return r;
}

template<class S>
S* operator* (const dmat<S> & a, const S* v)
{
  S* r = new S[a.rows()];
  a.mult(v, r);
  return r;
}

template<class S>
dmat<S> operator* (const dmat<S> & a, const S v)
{
  dmat<S> r(a);
  r *= v;
  return r;
}

template<class S>
dmat<S> operator/ (const dmat<S> & a, const S v)
{
  dmat<S> r(a);
  r /= v;
  return r;
}

template<class S>
dmat<S> operator+ (const dmat<S> & a, const dmat<S> & b)
{
  dmat<S> r(a);
  r += b;
  return r;
}

template<class S>
dmat<S> operator- (const dmat<S> & a, const dmat<S> & b)
{
  dmat<S> r(a);
  r -= b;
  return r;
}

template<class S>
void invert_3x3(dmat<S> & m, S & det)
{
  assert(m.rows() == 3 && m.cols() == 3);

  // save block entries in temp variables
  S* ele = m[0];
  S e0 = ele[0], e1 = ele[1], e2 = ele[2], e3 = ele[3];
  S e4 = ele[4], e5 = ele[5], e6 = ele[6], e7 = ele[7], e8 = ele[8];
  // compute determinant
  det = e0*e4*e8 + e1*e5*e6 + e2*e3*e7 - e2*e4*e6 - e1*e3*e8 - e0*e5*e7;
  S idet = 1.0 / det;

  // invert block
  ele[0] = (e4*e8 - e5*e7) * idet;
  ele[1] = (e2*e7 - e1*e8) * idet;
  ele[2] = (e1*e5 - e2*e4) * idet;
  ele[3] = (e5*e6 - e3*e8) * idet;
  ele[4] = (e0*e8 - e2*e6) * idet;
  ele[5] = (e2*e3 - e0*e5) * idet;
  ele[6] = (e3*e7 - e4*e6) * idet;
  ele[7] = (e1*e6 - e0*e7) * idet;
  ele[8] = (e0*e4 - e1*e3) * idet;
}

template<class S>
void invert_2x2(dmat<S> & m, S & det)
{
  assert(m.rows() == 2 && m.cols() == 2);

  // save block entries in temp variables
  S* ele = m[0];
  S e0 = ele[0], e1 = ele[1], e2 = ele[2], e3 = ele[3];
  // compute determinant
  det = e0*e3 - e1*e2;
  S idet = 1.0 / det;

  // invert block
  ele[0] =  e3 * idet;
  ele[1] = -e1 * idet;
  ele[2] = -e2 * idet;
  ele[3] =  e0 * idet;
}

template<class S, class V>
void array_to_tensors(const mt_vector<S> & arr, mt_vector<dmat<V> > & m)
{
  m.resize(arr.size()/9);
  for(size_t i=0; i<m.size(); i++)
  {
    m[i].assign(3, 3, &arr[i*9]);
  }
}

template<class S, class V>
void tensors_to_array(const mt_vector<dmat<S> > & m, mt_vector<V> & arr)
{
  arr.resize(m.size()*9);
  for(size_t i=0; i<m.size(); i++)
  {

    arr[i*9+0] = m[i][0][0]; arr[i*9+1] = m[i][0][1]; arr[i*9+2] = m[i][0][2];
    arr[i*9+3] = m[i][1][0]; arr[i*9+4] = m[i][1][1]; arr[i*9+5] = m[i][1][2];
    arr[i*9+6] = m[i][2][0]; arr[i*9+7] = m[i][2][1]; arr[i*9+8] = m[i][2][2];
  }
}

#endif
