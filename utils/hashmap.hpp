/**
* @file mt_hashmap.h
* @brief Classes similar to unordered_set and unordered_map, but better performance.
*
* This code is based on hashlib.h by Clifford Wolf <clifford@clifford.at>.
*
* @author Aurel Neic
* @version
* @date 2017-04-07
*/

#ifndef _MT_HASHMAP_H
#define _MT_HASHMAP_H

#include <cstdint>
#include <climits>

#include <stdexcept>
#include <algorithm>
#include <string>
#include <vector>

// #define NDEBUG

typedef unsigned long int hm_uint;
typedef long int hm_int;

namespace hashmap {

const hm_int hashtable_size_trigger = 2;
const hm_int hashtable_size_factor = 3;

// The XOR version of DJB2
inline hm_uint mkhash(hm_uint a, hm_uint b) {
  return ((a << 5) + a) ^ b;
}

// traditionally 5381 is used as starting value for the djb2 hash
const hm_uint mkhash_init = 5381;

// The ADD version of DJB2
// (use this version for cache locality in b)
inline hm_uint mkhash_add(hm_uint a, hm_uint b) {
  return ((a << 5) + a) + b;
}

inline hm_uint mkhash_xorshift(hm_uint a) {
  if (sizeof(a) == 4) {
    a ^= a << 13;
    a ^= a >> 17;
    a ^= a << 5;
  } else if (sizeof(a) == 8) {
    a ^= a << 13;
    a ^= a >> 7;
    a ^= a << 17;
  } else
    throw std::runtime_error("mkhash_xorshift() only implemented for 32 bit and 64 bit ints");
  return a;
}

// ================== Hashing structs ==============================

/**
* @brief Base hashing class
*
* @tparam T Datatype with a .hash() member.
*/
template<typename T>
struct hash_ops
{
  static inline bool cmp(const T &a, const T &b) {
    return a == b;
  }
  static inline hm_uint hash(const T &a) {
    return a.hash();
  }
};

struct hash_int_ops {
  template<typename T>
  static inline bool cmp(T a, T b) {
    return a == b;
  }
};

template<> struct hash_ops<int32_t> : hash_int_ops
{
  static inline hm_uint hash(int32_t a) {
    return a;
  }
};
template<> struct hash_ops<int64_t> : hash_int_ops
{
  static inline hm_uint hash(int64_t a) {
    return mkhash((hm_uint)(a), (hm_uint)(a >> 32));
  }
};

// in the case that long is not an int64 we define an additional long hasher
#if defined(__APPLE__) || defined(WINBUILD)
template<> struct hash_ops<long> : hash_int_ops
{
  static inline hm_uint hash(long a) {
    return (hm_uint)a;
  }
};
#endif

template<> struct hash_ops<std::string> {
  static inline bool cmp(const std::string &a, const std::string &b) {
    return a == b;
  }
  static inline hm_uint hash(const std::string &a) {
    hm_uint v = 0;
    for (auto c : a)
      v = mkhash(v, c);
    return v;
  }
};

template<typename P, typename Q> struct hash_ops<std::pair<P, Q>> {
  static inline bool cmp(std::pair<P, Q> a, std::pair<P, Q> b) {
    return a == b;
  }
  static inline hm_uint hash(std::pair<P, Q> a) {
    return mkhash(hash_ops<P>::hash(a.first), hash_ops<Q>::hash(a.second));
  }
};

template<typename... T> struct hash_ops<std::tuple<T...>> {
  static inline bool cmp(std::tuple<T...> a, std::tuple<T...> b) {
    return a == b;
  }
  template<size_t I = 0>
  static inline typename std::enable_if<I == sizeof...(T), hm_uint>::type hash(std::tuple<T...>) {
    return mkhash_init;
  }
  template<size_t I = 0>
  static inline typename std::enable_if<I != sizeof...(T), hm_uint>::type hash(std::tuple<T...> a) {
    typedef hash_ops<typename std::tuple_element<I, std::tuple<T...>>::type> element_ops_t;
    return mkhash(hash<I+1>(a), element_ops_t::hash(std::get<I>(a)));
  }
};

template<typename T> struct hash_ops<std::vector<T>> {
  static inline bool cmp(std::vector<T> a, std::vector<T> b) {
    return a == b;
  }
  static inline hm_uint hash(std::vector<T> a) {
    hm_uint h = mkhash_init;
    for (auto k : a)
      h = mkhash(h, hash_ops<T>::hash(k));
    return h;
  }
};

struct hash_cstr_ops {
  static inline bool cmp(const char *a, const char *b) {
    for (int i = 0; a[i] || b[i]; i++)
      if (a[i] != b[i])
        return false;
    return true;
  }
  static inline hm_uint hash(const char *a) {
    hm_uint hash = mkhash_init;
    while (*a)
      hash = mkhash(hash, *(a++));
    return hash;
  }
};

struct hash_ptr_ops {
  static inline bool cmp(const void *a, const void *b) {
    return a == b;
  }
  static inline hm_uint hash(const void *a) {
    return *((hm_uint*)(a));
  }
};

struct hash_obj_ops {
  static inline bool cmp(const void *a, const void *b) {
    return a == b;
  }
  template<typename T>
  static inline hm_uint hash(const T *a) {
    return a ? a->hash() : 0;
  }
};

template<typename T>
inline hm_uint mkhash(const T &v) {
  return hash_ops<T>().hash(v);
}

inline hm_int hashtable_size(hm_int min_size)
{
  static std::vector<hm_int> zero_and_some_primes = {
    0, 23, 29, 37, 47, 59, 79, 101, 127, 163, 211, 269, 337, 431, 541, 677,
    853, 1069, 1361, 1709, 2137, 2677, 3347, 4201, 5261, 6577, 8231, 10289,
    12889, 16127, 20161, 25219, 31531, 39419, 49277, 61603, 77017, 96281,
    120371, 150473, 188107, 235159, 293957, 367453, 459317, 574157, 717697,
    897133, 1121423, 1401791, 1752239, 2190299, 2737937, 3422429, 4278037,
    5347553, 6684443, 8355563, 10444457, 13055587, 16319519, 20399411,
    25499291, 31874149, 39842687, 49803361, 62254207, 77817767, 97272239,
    121590311, 151987889, 189984863, 237481091, 296851369, 371064217
  };

  for (auto p : zero_and_some_primes)
    if (p >= min_size) return p;

  if (sizeof(hm_int) == 4)
    throw std::length_error("hash table exceeded maximum size. use a ILP64 abi for larger tables.");

  for (auto p : zero_and_some_primes)
    if (100129 * p > min_size) return 100129 * p;

  throw std::length_error("hash table exceeded maximum size.");
}

template<typename K, typename T, typename OPS = hash_ops<K>> class unordered_map;
template<typename K, hm_int offset = 0, typename OPS = hash_ops<K>> class idict;
template<typename K, typename OPS = hash_ops<K>> class unordered_set;
template<typename K, typename OPS = hash_ops<K>> class mfp;

/**
* @brief
*
* @tparam K
* @tparam T
* @tparam OPS
*/
template<typename K, typename T, typename OPS>
class unordered_map
{
  /// internel entry type
  struct entry_t
  {
    std::pair<K, T> udata;
    hm_int next;

    entry_t() { }
    entry_t(const std::pair<K, T> & idata, hm_int inext) : udata(idata), next(inext) { }
    entry_t(std::pair<K, T> && idata, hm_int inext) : udata(std::move(idata)), next(inext) { }
  };

  /// the hashtable
  std::vector<hm_int> hashtable;
  /// the stored entries
  std::vector<entry_t> entries;
  /// the hash generator
  OPS ops;

  #ifdef NDEBUG
  static inline void do_assert(bool) { }
  #else
  static inline void do_assert(bool cond) {
    if (!cond) throw std::runtime_error("unordered_map<> assert failed.");
  }
  #endif


  /**
  * @brief Generate a hash from a key.
  *
  * Hash in always in the range of the current hash table size.
  *
  * @param key The key.
  *
  * @return The hash.
  */
  hm_int do_hash(const K &key) const
  {
    hm_uint hash = 0;
    if (!hashtable.empty())
      hash = ops.hash(key) % (hm_uint)(hashtable.size());
    return hash;
  }

  /**
  * @brief Resize the hashtable and compute new hashes.
  */
  void do_rehash()
  {
    hashtable.clear();
    hashtable.resize(hashtable_size(entries.capacity() * hashtable_size_factor), -1);

    for (hm_int i = 0; i < hm_int(entries.size()); i++) {
      do_assert(-1 <= entries[i].next && entries[i].next < hm_int(entries.size()));
      hm_int hash = do_hash(entries[i].udata.first);
      entries[i].next = hashtable[hash];
      hashtable[hash] = i;
    }
  }

  /**
  * @brief Remove an entry.
  *
  * @param index The index in the entries vector.
  * @param hash  The hash -> index w.r.t. hashtable.
  *
  * @return
  */
  hm_int do_erase(hm_int index, hm_int hash)
  {
    do_assert(index < hm_int(entries.size()));
    if (hashtable.empty() || index < 0)
      return 0;

    hm_int k = hashtable[hash];
    do_assert(0 <= k && k < hm_int(entries.size()));

    if (k == index) {
      hashtable[hash] = entries[index].next;
    } else {
      while (entries[k].next != index) {
        k = entries[k].next;
        do_assert(0 <= k && k < hm_int(entries.size()));
      }
      entries[k].next = entries[index].next;
    }

    hm_int back_idx = entries.size()-1;

    if (index != back_idx)
    {
      hm_int back_hash = do_hash(entries[back_idx].udata.first);

      k = hashtable[back_hash];
      do_assert(0 <= k && k < hm_int(entries.size()));

      if (k == back_idx) {
        hashtable[back_hash] = index;
      } else {
        while (entries[k].next != back_idx) {
          k = entries[k].next;
          do_assert(0 <= k && k < hm_int(entries.size()));
        }
        entries[k].next = index;
      }

      entries[index] = std::move(entries[back_idx]);
    }

    entries.pop_back();

    if (entries.empty())
      hashtable.clear();

    return 1;
  }

  /**
  * @brief Return hash and index for a key.
  *
  * @param [in]  key   The key value.
  * @param [out] hash  The hash, i.e. the hashtable index.
  *
  * @return The enties index.
  */
  hm_int do_lookup(const K &key, hm_int &hash) const
  {
    if (hashtable.empty())
      return -1;

    if (entries.size() * hashtable_size_trigger > hashtable.size()) {
      ((unordered_map*)this)->do_rehash();
      hash = do_hash(key);
    }

    hm_int index = hashtable[hash];

    while (index >= 0 && !ops.cmp(entries[index].udata.first, key)) {
      index = entries[index].next;
      do_assert(-1 <= index && index < hm_int(entries.size()));
    }

    return index;
  }

  /**
  * @brief Insert a pair consisting of a key and a default (empty) value.
  *
  * New values are always added to the end of the entries array.
  * Their index is inserted into the hashtable at the hash index.
  *
  * @param [in]  key   The key.
  * @param [out] hash  The hash index.
  *
  * @return Insert index.
  */
  hm_int do_insert(const K &key, hm_int &hash)
  {
    if (hashtable.empty()) {
      entries.push_back(entry_t(std::pair<K, T>(key, T()), -1));
      do_rehash();
      hash = do_hash(key);
    } else {
      entries.push_back(entry_t(std::pair<K, T>(key, T()), hashtable[hash]));
      hashtable[hash] = entries.size() - 1;
    }
    return entries.size() - 1;
  }

  /**
  * @brief Insert a key-value pair.
  *
  * New values are always added to the end of the entries array.
  * Their index is inserted into the hashtable at the hash index.
  *
  * @param [in]  key   The key-value pair.
  * @param [out] hash  The hash index.
  *
  * @return Insert index.
  */
  hm_int do_insert(const std::pair<K, T> &value, hm_int &hash)
  {
    if (hashtable.empty()) {
      entries.push_back(entry_t(value, -1));
      do_rehash();
      hash = do_hash(value.first);
    } else {
      entries.push_back(entry_t(value, hashtable[hash]));
      hashtable[hash] = entries.size() - 1;
    }
    return entries.size() - 1;
  }



  public:

  /// Const iterator.
  class const_iterator : public std::iterator<std::forward_iterator_tag, std::pair<K, T>>
  {
    friend class unordered_map;
    protected:
    const unordered_map *ptr;
    hm_int index;
    const_iterator(const unordered_map *iptr, hm_int iindex) : ptr(iptr), index(iindex) { }
    public:
    const_iterator() { }
    const_iterator & operator++()  { index--; return *this; }  // pre-increment
    const_iterator operator++(int) { index--; return *this; }  // post-increment
    const_iterator & operator--()  { index++; return *this; }  // pre-decrement
    const_iterator operator--(int) { index++; return *this; }  // post-decrement
    bool operator<(const const_iterator &other) const { return index > other.index; }
    bool operator==(const const_iterator &other) const { return index == other.index; }
    bool operator!=(const const_iterator &other) const { return index != other.index; }
    const std::pair<K, T> &operator*() const { return ptr->entries[index].udata; }
    const std::pair<K, T> *operator->() const { return &ptr->entries[index].udata; }
  };

  /// Default iterator.
  class iterator : public std::iterator<std::forward_iterator_tag, std::pair<K, T>>
  {
    friend class unordered_map;
    protected:
    unordered_map *ptr;
    hm_int index;
    iterator(unordered_map *iptr, hm_int iindex) : ptr(iptr), index(iindex) { }
    public:
    iterator() { }
    iterator & operator++()  { index--; return *this; }  // pre-increment
    iterator operator++(int) { index--; return *this; }  // post-increment
    iterator & operator--()  { index++; return *this; }  // pre-decrement
    iterator operator--(int) { index++; return *this; }  // post-decrement
    bool operator<(const iterator &other) const { return index > other.index; }
    bool operator==(const iterator &other) const { return index == other.index; }
    bool operator!=(const iterator &other) const { return index != other.index; }
    std::pair<K, T> &operator*() { return ptr->entries[index].udata; }
    std::pair<K, T> *operator->() { return &ptr->entries[index].udata; }
    const std::pair<K, T> &operator*() const { return ptr->entries[index].udata; }
    const std::pair<K, T> *operator->() const { return &ptr->entries[index].udata; }
    operator const_iterator() const { return const_iterator(ptr, index); }
  };

  /// Empty constructor.
  unordered_map()
  {}

  /// Construct from another map.
  unordered_map(const unordered_map &other)
  {
    entries = other.entries;
    do_rehash();
  }

  /// Construct map form another map.
  unordered_map(unordered_map &&other)
  {
    swap(other);
  }

  unordered_map &operator=(const unordered_map &other) {
    entries = other.entries;
    do_rehash();
    return *this;
  }

  unordered_map &operator=(unordered_map &&other) {
    clear();
    swap(other);
    return *this;
  }

  unordered_map(const std::initializer_list<std::pair<K, T>> &list)
  {
    for (auto &it : list)
      insert(it);
  }

  /// Construct from Iterator range.
  template<class InputIterator>
  unordered_map(InputIterator first, InputIterator last)
  {
    insert(first, last);
  }
  /// Insert Iterator range.
  template<class InputIterator>
  void insert(InputIterator first, InputIterator last)
  {
    for (; first != last; ++first)
      insert(*first);
  }
  /// User insert as key lookup.
  std::pair<iterator, bool> insert(const K &key)
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    if (i >= 0)
      return std::pair<iterator, bool>(iterator(this, i), false);
    i = do_insert(key, hash);
    return std::pair<iterator, bool>(iterator(this, i), true);
  }
  /// Insert key-value pair.
  std::pair<iterator, bool> insert(const std::pair<K, T> &value)
  {
    hm_int hash = do_hash(value.first);
    hm_int i = do_lookup(value.first, hash);
    if (i >= 0)
      return std::pair<iterator, bool>(iterator(this, i), false);
    i = do_insert(value, hash);
    return std::pair<iterator, bool>(iterator(this, i), true);
  }
  /// Erase by key.
  hm_int erase(const K &key)
  {
    hm_int hash = do_hash(key);
    hm_int index = do_lookup(key, hash);
    return do_erase(index, hash);
  }
  /// Erase by iterator.
  iterator erase(iterator it)
  {
    hm_int hash = do_hash(it->first);
    do_erase(it.index, hash);
    return ++it;
  }
  /// Check if key exists.
  hm_int count(const K &key) const
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    return i < 0 ? 0 : 1;
  }
  /// Check if key exists and matches iterator.
  hm_int count(const K &key, const_iterator it) const
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    return i < 0 || i > it.index ? 0 : 1;
  }
  /// Search for key. Return iterator.
  iterator find(const K &key)
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    if (i < 0)
      return end();
    return iterator(this, i);
  }

  /// Search for key. Return const_iterator.
  const_iterator find(const K &key) const
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    if (i < 0)
      return end();
    return const_iterator(this, i);
  }

  /// Data access by key.
  T& at(const K &key)
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    if (i < 0)
      throw std::out_of_range("unordered_map::at()");
    return entries[i].udata.second;
  }

  /// Const data access by key.
  const T& at(const K &key) const
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    if (i < 0)
      throw std::out_of_range("unordered_map::at()");
    return entries[i].udata.second;
  }

  /// Return data if existent or default value.
  T at(const K &key, const T &defval) const
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    if (i < 0)
      return defval;
    return entries[i].udata.second;
  }

  /// Data access or empty insert.
  T& operator[](const K &key)
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    if (i < 0)
      i = do_insert(std::pair<K, T>(key, T()), hash);
    return entries[i].udata.second;
  }

  /**
  * @brief Sort data entries.
  *
  * @tparam Compare
  * @param comp
  */
  template<typename Compare = std::less<K>>
  void sort(Compare comp = Compare())
  {
    std::sort(entries.begin(), entries.end(), [comp](const entry_t &a, const entry_t &b){ return comp(b.udata.first, a.udata.first); });
    do_rehash();
  }

  void swap(unordered_map &other)
  {
    hashtable.swap(other.hashtable);
    entries.swap(other.entries);
  }

  bool operator==(const unordered_map &other) const {
    if (size() != other.size())
      return false;
    for (auto &it : entries) {
      auto oit = other.find(it.udata.first);
      if (oit == other.end() || !(oit->second == it.udata.second))
        return false;
    }
    return true;
  }

  bool operator!=(const unordered_map &other) const {
    return !operator==(other);
  }

  void reserve(size_t n) { entries.reserve(n); }
  size_t size() const { return entries.size(); }
  bool empty() const { return entries.empty(); }
  void clear() { hashtable.clear(); entries.clear(); }

  iterator begin() { return iterator(this, hm_int(entries.size())-1); }
  iterator end() { return iterator(nullptr, -1); }

  const_iterator begin() const { return const_iterator(this, hm_int(entries.size())-1); }
  const_iterator end() const { return const_iterator(nullptr, -1); }
};

/**
* @brief Custom unordered_set implementation.
*
* @tparam K    Key type.
* @tparam OPS  Hashing and comparison class.
*/
template<typename K, typename OPS>
class unordered_set
{
  template<typename, hm_int, typename> friend class idict;

  protected:
  /// internel entry type
  struct entry_t
  {
    K udata;
    hm_int next;

    entry_t() { }
    entry_t(const K &idata, hm_int inext) : udata(idata), next(inext) { }
  };

  /// the hashtable
  std::vector<hm_int> hashtable;
  /// the stored entries
  std::vector<entry_t> entries;
  /// the hash generator
  OPS ops;

  #ifdef NDEBUG
  static inline void do_assert(bool) { }
  #else
  static inline void do_assert(bool cond) {
    if (!cond) throw std::runtime_error("unordered_set<> assert failed.");
  }
  #endif

  /**
  * @brief Generate a hash from a key.
  *
  * Hash in always in the range of the current hash table size.
  *
  * @param key The key.
  *
  * @return The hash.
  */
  hm_int do_hash(const K &key) const
  {
    hm_uint hash = 0;
    if (!hashtable.empty())
      hash = ops.hash(key) % (hm_uint)(hashtable.size());
    return hash;
  }

  /**
  * @brief Resize the hashtable and compute new hashes.
  */
  void do_rehash()
  {
    hashtable.clear();
    hashtable.resize(hashtable_size(entries.capacity() * hashtable_size_factor), -1);

    for (hm_int i = 0; i < hm_int(entries.size()); i++) {
      do_assert(-1 <= entries[i].next && entries[i].next < hm_int(entries.size()));
      hm_int hash = do_hash(entries[i].udata);
      entries[i].next = hashtable[hash];
      hashtable[hash] = i;
    }
  }

  /**
  * @brief Remove an entry.
  *
  * @param index The index in the entries vector.
  * @param hash  The hash -> index w.r.t. hashtable.
  *
  * @return
  */
  hm_int do_erase(hm_int index, hm_int hash)
  {
    do_assert(index < hm_int(entries.size()));
    if (hashtable.empty() || index < 0)
      return 0;

    hm_int k = hashtable[hash];
    if (k == index) {
      hashtable[hash] = entries[index].next;
    } else {
      while (entries[k].next != index) {
        k = entries[k].next;
        do_assert(0 <= k && k < hm_int(entries.size()));
      }
      entries[k].next = entries[index].next;
    }

    hm_int back_idx = entries.size()-1;

    if (index != back_idx)
    {
      hm_int back_hash = do_hash(entries[back_idx].udata);

      k = hashtable[back_hash];
      if (k == back_idx) {
        hashtable[back_hash] = index;
      } else {
        while (entries[k].next != back_idx) {
          k = entries[k].next;
          do_assert(0 <= k && k < hm_int(entries.size()));
        }
        entries[k].next = index;
      }

      entries[index] = std::move(entries[back_idx]);
    }

    entries.pop_back();

    if (entries.empty())
      hashtable.clear();

    return 1;
  }

  /**
  * @brief Return hash and index for a key.
  *
  * @param [in]  key   The key value.
  * @param [out] hash  The hash, i.e. the hashtable index.
  *
  * @return The enties index.
  */
  hm_int do_lookup(const K &key, hm_int &hash) const
  {
    if (hashtable.empty())
      return -1;

    if (entries.size() * hashtable_size_trigger > hashtable.size()) {
      ((unordered_set*)this)->do_rehash();
      hash = do_hash(key);
    }

    hm_int index = hashtable[hash];

    while (index >= 0 && !ops.cmp(entries[index].udata, key)) {
      index = entries[index].next;
      do_assert(-1 <= index && index < hm_int(entries.size()));
    }

    return index;
  }

  /**
  * @brief Insert a pair consisting of a key and a default (empty) value.
  *
  * New values are always added to the end of the entries array.
  * Their index is inserted into the hashtable at the hash index.
  *
  * @param [in]  key   The key.
  * @param [out] hash  The hash index.
  *
  * @return Insert index.
  */
  hm_int do_insert(const K &value, hm_int &hash)
  {
    if (hashtable.empty()) {
      entries.push_back(entry_t(value, -1));
      do_rehash();
      hash = do_hash(value);
    }
    else {
      entries.push_back(entry_t(value, hashtable[hash]));
      hashtable[hash] = entries.size() - 1;
    }
    return entries.size() - 1;
  }

  public:
  class const_iterator : public std::iterator<std::forward_iterator_tag, K>
  {
    friend class unordered_set;

    protected:
    const unordered_set *ptr;
    hm_int index;
    const_iterator(const unordered_set *iptr, hm_int iindex) : ptr(iptr), index(iindex) { }

    public:
    const_iterator() { }
    const_iterator & operator++()  { index--; return *this; }  // pre-increment
    const_iterator operator++(int) { index--; return *this; }  // post-increment
    const_iterator & operator--()  { index++; return *this; }  // pre-decrement
    const_iterator operator--(int) { index++; return *this; }  // post-decrement
    bool operator==(const const_iterator &other) const { return index == other.index; }
    bool operator!=(const const_iterator &other) const { return index != other.index; }
    const K &operator*() const { return ptr->entries[index].udata; }
    const K *operator->() const { return &ptr->entries[index].udata; }
  };

  class iterator : public std::iterator<std::forward_iterator_tag, K>
  {
    friend class unordered_set;

    protected:
    unordered_set *ptr;
    hm_int index;
    iterator(unordered_set *iptr, hm_int iindex) : ptr(iptr), index(iindex) { }

    public:
    iterator() { }
    iterator & operator++()  { index--; return *this; }  // pre-increment
    iterator operator++(int) { index--; return *this; }  // post-increment
    iterator & operator--()  { index++; return *this; }  // pre-decrement
    iterator operator--(int) { index++; return *this; }  // post-decrement
    bool operator==(const iterator &other) const { return index == other.index; }
    bool operator!=(const iterator &other) const { return index != other.index; }
    K &operator*() { return ptr->entries[index].udata; }
    K *operator->() { return &ptr->entries[index].udata; }
    const K &operator*() const { return ptr->entries[index].udata; }
    const K *operator->() const { return &ptr->entries[index].udata; }
    operator const_iterator() const { return const_iterator(ptr, index); }
  };

  /// Empty constructor.
  unordered_set()
  { }

  /// Construct from another set.
  unordered_set(const unordered_set &other)
  {
    entries = other.entries;
    do_rehash();
  }

  /// Construct from another set.
  unordered_set(unordered_set &&other)
  {
    swap(other);
  }

  unordered_set &operator=(const unordered_set &other) {
    entries = other.entries;
    do_rehash();
    return *this;
  }

  unordered_set &operator=(unordered_set &&other) {
    clear();
    swap(other);
    return *this;
  }

  unordered_set(const std::initializer_list<K> &list)
  {
    for (auto &it : list)
      insert(it);
  }

  template<class InputIterator>
  unordered_set(InputIterator first, InputIterator last)
  {
    insert(first, last);
  }

  template<class InputIterator>
  void insert(InputIterator first, InputIterator last)
  {
    for (; first != last; ++first)
      insert(*first);
  }

  std::pair<iterator, bool> insert(const K &value)
  {
    hm_int hash = do_hash(value);
    hm_int i = do_lookup(value, hash);
    if (i >= 0)
      return std::pair<iterator, bool>(iterator(this, i), false);
    i = do_insert(value, hash);
    return std::pair<iterator, bool>(iterator(this, i), true);
  }

  hm_int erase(const K &key)
  {
    hm_int hash = do_hash(key);
    hm_int index = do_lookup(key, hash);
    return do_erase(index, hash);
  }

  iterator erase(iterator it)
  {
    hm_int hash = do_hash(*it);
    do_erase(it.index, hash);
    return ++it;
  }

  hm_int count(const K &key) const
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    return i < 0 ? 0 : 1;
  }

  hm_int count(const K &key, const_iterator it) const
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    return i < 0 || i > it.index ? 0 : 1;
  }

  iterator find(const K &key)
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    if (i < 0)
      return end();
    return iterator(this, i);
  }

  const_iterator find(const K &key) const
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    if (i < 0)
      return end();
    return const_iterator(this, i);
  }

  bool operator[](const K &key)
  {
    hm_int hash = do_hash(key);
    hm_int i = do_lookup(key, hash);
    return i >= 0;
  }

  template<typename Compare = std::less<K>>
  void sort(Compare comp = Compare())
  {
    std::sort(entries.begin(), entries.end(), [comp](const entry_t &a, const entry_t &b){ return comp(b.udata, a.udata); });
    do_rehash();
  }

  K pop()
  {
    iterator it = begin();
    K ret = *it;
    erase(it);
    return ret;
  }

  void swap(unordered_set &other)
  {
    hashtable.swap(other.hashtable);
    entries.swap(other.entries);
  }

  bool operator==(const unordered_set &other) const {
    if (size() != other.size())
      return false;
    for (auto &it : entries)
      if (!other.count(it.udata))
        return false;
    return true;
  }

  bool operator!=(const unordered_set &other) const {
    return !operator==(other);
  }

  void reserve(size_t n) { entries.reserve(n); }
  size_t size() const { return entries.size(); }
  bool empty() const { return entries.empty(); }
  void clear() { hashtable.clear(); entries.clear(); }

  iterator begin() { return iterator(this, hm_int(entries.size())-1); }
  iterator end() { return iterator(nullptr, -1); }

  const_iterator begin() const { return const_iterator(this, hm_int(entries.size())-1); }
  const_iterator end() const { return const_iterator(nullptr, -1); }
};

template<typename K, hm_int offset, typename OPS>
class idict
{
  unordered_set<K, OPS> database;

public:
  typedef typename unordered_set<K, OPS>::const_iterator const_iterator;

  hm_int operator()(const K &key)
  {
    hm_int hash = database.do_hash(key);
    hm_int i = database.do_lookup(key, hash);
    if (i < 0)
      i = database.do_insert(key, hash);
    return i + offset;
  }

  hm_int at(const K &key) const
  {
    hm_int hash = database.do_hash(key);
    hm_int i = database.do_lookup(key, hash);
    if (i < 0)
      throw std::out_of_range("idict::at()");
    return i + offset;
  }

  hm_int at(const K &key, hm_int defval) const
  {
    hm_int hash = database.do_hash(key);
    hm_int i = database.do_lookup(key, hash);
    if (i < 0)
      return defval;
    return i + offset;
  }

  hm_int count(const K &key) const
  {
    hm_int hash = database.do_hash(key);
    hm_int i = database.do_lookup(key, hash);
    return i < 0 ? 0 : 1;
  }

  void expect(const K &key, hm_int i)
  {
    hm_int j = (*this)(key);
    if (i != j)
      throw std::out_of_range("idict::expect()");
  }

  const K &operator[](hm_int index) const
  {
    return database.entries.at(index - offset).udata;
  }

  void swap(idict &other)
  {
    database.swap(other.database);
  }

  void reserve(size_t n) { database.reserve(n); }
  size_t size() const { return database.size(); }
  bool empty() const { return database.empty(); }
  void clear() { database.clear(); }

  const_iterator begin() const { return database.begin(); }
  const_iterator end() const { return database.end(); }
};

template<typename K, typename OPS>
class mfp
{
  mutable idict<K, 0, OPS> database;
  mutable std::vector<hm_int> parents;

public:
  typedef typename idict<K, 0, OPS>::const_iterator const_iterator;

  hm_int operator()(const K &key) const
  {
    hm_int i = database(key);
    parents.resize(database.size(), -1);
    return i;
  }

  const K &operator[](hm_int index) const
  {
    return database[index];
  }

  hm_int ifind(hm_int i) const
  {
    hm_int p = i, k = i;

    while (parents[p] != -1)
      p = parents[p];

    while (k != p) {
      hm_int next_k = parents[k];
      parents[k] = p;
      k = next_k;
    }

    return p;
  }

  void imerge(hm_int i, hm_int j)
  {
    i = ifind(i);
    j = ifind(j);

    if (i != j)
      parents[i] = j;
  }

  void ipromote(hm_int i)
  {
    hm_int k = i;

    while (k != -1) {
      hm_int next_k = parents[k];
      parents[k] = i;
      k = next_k;
    }

    parents[i] = -1;
  }

  hm_int lookup(const K &a) const
  {
    return ifind((*this)(a));
  }

  const K &find(const K &a) const
  {
    hm_int i = database.at(a, -1);
    if (i < 0)
      return a;
    return (*this)[ifind(i)];
  }

  void merge(const K &a, const K &b)
  {
    imerge((*this)(a), (*this)(b));
  }

  void promote(const K &a)
  {
    hm_int i = database.at(a, -1);
    if (i >= 0)
      ipromote(i);
  }

  void swap(mfp &other)
  {
    database.swap(other.database);
    parents.swap(other.parents);
  }

  void reserve(size_t n) { database.reserve(n); }
  size_t size() const { return database.size(); }
  bool empty() const { return database.empty(); }
  void clear() { database.clear(); parents.clear(); }

  const_iterator begin() const { return database.begin(); }
  const_iterator end() const { return database.end(); }
};

} /* namespace hashmap */

#endif
