/**
* @file mt_vector.h
* @brief Vector class similar to std::vector.
*
* This implementation was inspired by the toolbox_vector implementation by Manfred Liebmann.
*
* @author Aurel Neic
* @version
* @date 2016-12-13
*/

#ifndef _TOOLBOX_VECTOR
#define _TOOLBOX_VECTOR

#include <assert.h>
#include <iterator>


template <class T>
class mt_vector
{
public:
  mt_vector(): _data(NULL), _capacity(0), _size(0)
  {}

  mt_vector(size_t n): _data(NULL), _capacity(0), _size(0)
  {
    this->resize(n);
  }

  mt_vector(size_t n, const T val): _data(NULL), _capacity(0), _size(0)
  {
    this->assign(n, val);
  }

  mt_vector(const mt_vector<T> &vec): _data(NULL), _capacity(0), _size(0)
  {
    this->assign(vec.begin(), vec.end());
  }

  virtual ~mt_vector()
  {
    delete [] _data;
  }

  const T& operator[](size_t i) const
  {
    return _data[i];
  }

  T& operator[](size_t i)
  {
    return _data[i];
  }

  void operator= (const mt_vector<T> & vec)
  {
    this->assign(vec.begin(), vec.end());
  }

  T* data()
  {
    return _data;
  }

  const T* data() const
  {
    return _data;
  }

  size_t size() const
  {
    return _size;
  }
  size_t capacity() const
  {
    return _capacity;
  }

  const T* begin() const
  {
    return _data;
  }

  T* begin()
  {
    return _data;
  }

  const T* end() const
  {
    return _data + _size;
  }

  T* end()
  {
    return _data + _size;
  }

  #if 0
  void assign(const T* s, const T* e) {
    size_t n = (e - s);
    if (_capacity < n) {
      delete [] _data;
      _data = new T[n];
      _capacity = n;
    }
    for (size_t i = 0; i < n; i++) _data[i] = s[i];
    _size = n;
  }
  #endif

  template<class InputIterator>
  void assign(InputIterator s, InputIterator e)
  {
    long int n = std::distance(s, e);
    assert(n >= 0);

    // realloc if necessary
    if ( ((long int)_capacity) < n) {
      delete [] _data;
      _data = new T[n];
      _capacity = n;
    }

    // copy data
    size_t idx = 0;
    while(s != e) {
      _data[idx++] = *s;
      ++s;
    }

    _size = n;
  }

  void assign(size_t n, T val = T()) {
    if (_capacity < n) {
      if(_size > 0) delete [] _data;
      _data = new T[n];
      _capacity = n;
    }
    for (size_t i = 0; i < n; i++) _data[i] = val;
    _size = n;
  }

  void assign(size_t n, T* array, bool del = true) {
    if (del) delete [] _data;
    _data = array;
    _size = n;
    _capacity = n;
  }

  void resize(size_t n)
  {
    if(_capacity < n)
    {
      if(_size > 0) {
        T* buf = _data;
        _data = new T[n];
        for (size_t i = 0; i < _size; i++) _data[i] = buf[i];
        delete [] buf;
      }
      else {
        _data = new T[n];
      }
      _capacity = n;
    }
    _size = n;
  }

  void resize(size_t n, const T val)
  {
    size_t oldsize = this->size();
    this->resize(n);

    if(n > oldsize) {
      T*     out = _data + oldsize;
      size_t num = n - oldsize;

      for (size_t i = 0; i < num; i++) out[i] = val;
    }
  }

  void reserve(size_t n)
  {
    size_t osize = _size;
    this->resize(n);
    this->resize(osize);
  }

  void zero()
  {
    if(_size != 0) memset(_data, 0, _size * sizeof(T));
  }

  void reallocate()
  {
    _capacity = _size;
    T* buf = _data;
    _data = new T[_size];
    for (size_t i = 0; i < _size; i++) _data[i] = buf[i];
    delete [] buf;
  }

  void append(const T* s, const T* e)
  {
    size_t addsize = (e - s);
    size_t oldsize = _size;
    this->resize(oldsize + addsize);
    for (size_t i = 0; i < addsize; i++) _data[oldsize + i] = s[i];
  }

  T& push_back(const T & val)
  {
    if(_capacity > _size + 1)
    {
      _data[_size] = val;
      _size++;
    }
    else {
      // if current capacity is not larger than size + 1 we reallocate, but keep size to
      // old value plus one.
      size_t newcap = _capacity > 0 ? _capacity * 2 : 1;
      size_t cursize = _size;
      this->resize(newcap);
      _data[cursize] = val;
      _size = cursize + 1;
    }

    return _data[_size-1];
  }

  size_t write(FILE* fd) const
  {
    size_t items_written = fwrite(&_size, sizeof(size_t), 1, fd);
    if(items_written != 1) {
      fprintf(stderr, "%s error: Cannot write vector size to provided FD!\n", __func__);
      return items_written;
    }

    items_written = fwrite(_data, sizeof(T), _size, fd);

    if(items_written != _size)
      fprintf(stderr, "%s error: Cannot write vector data to provided FD!\n", __func__);

    return items_written;
  }

  size_t read(FILE* fd)
  {
    size_t newsize, items_read = 0;
    items_read = fread(&newsize, sizeof(size_t), 1, fd);
    if(items_read != 1) {
      fprintf(stderr, "%s error: Cannot read vector size from provided FD!\n", __func__);
      return items_read;
    }

    this->resize(newsize);

    items_read = fread(_data, sizeof(T), newsize, fd);
    if(items_read != newsize)
      fprintf(stderr, "%s error: Cannot read vector data from provided FD!\n", __func__);

    return items_read;
  }

private:
  T* _data;
  size_t _capacity;
  size_t _size;
};

/**
* @brief Class for storing a set of positive indices in a similar way to std::set, but using a
*        bool vector to be much faster for insert, erase and count.
*
* insert, erase and count are very very fast, but the class offers no iterator access.
* A vector of currently inserted entries can be retrieved via entries().
*
*/
class mt_mask {

  private:
  mt_vector<bool> _mask;      ///< the boolean mask holding which indices have been inserted
  mt_vector<size_t> _entries; ///< a buffer that gets populated when entries() is called

  /// check whether _mask is big enough to hold the new index
  inline bool can_access(size_t idx) const {
    return (idx < _mask.size());
  }

  public:
  mt_mask() {}

  /// resize _mask
  inline void resize(size_t s) {
    _mask.resize(s, false);
  }

  /// constructor with pre-allocation of _mask
  mt_mask(size_t s) {
    resize(s);
  }
  /// insert an index
  inline void insert(size_t idx) {
    if(!can_access(idx))
      resize(idx+1);

    _mask[idx] = true;
  }
  /// insert an index range
  template<typename InputIterator>
  inline void insert(InputIterator s, InputIterator e) {
    while(s != e) {
      insert(size_t(*s));
      ++s;
    }
  }
  /// check whether an index is present
  inline int count(size_t idx) const
  {
    if(!can_access(idx)) return 0;
    else                 return int(_mask[idx]);
  }
  /// erase an index
  inline void erase(size_t idx) {
    if(can_access(idx))
      _mask[idx] = false;
  }
  /// erase an index range
  template<typename InputIterator>
  inline void erase(InputIterator s, InputIterator e) {
    while(s != e) {
      erase(size_t(*s));
      ++s;
    }
  }
  /// clear the mask
  inline void clear() {
    _mask.assign(_mask.size(), false);
  }
  /// get a const reference to the inserted indices
  inline const mt_vector<size_t> & entries()
  {
    _entries.resize(_mask.size());
    size_t widx=0;
    for(size_t ridx=0; ridx<_mask.size(); ridx++)
      if(_mask[ridx]) _entries[widx++] = ridx;

    _entries.resize(widx);
    return _entries;
  }
};

template<class V>
V sum(const mt_vector<V> & vec)
{
  V ret = V();
  for(const V & e : vec) ret += e;

  return ret;
}

#endif

