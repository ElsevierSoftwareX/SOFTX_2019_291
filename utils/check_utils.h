/**
* @file check_utils.h
* @brief Sanity check functions
* @author Aurel Neic
* @version 
* @date 2018-12-13
*/


#ifndef _CHECK_UTILS_H
#define _CHECK_UTILS_H

#include <unistd.h>


enum mt_check_type {
  chk_nonzero,
  chk_fexists,
  dont_check
};

template<class U, class V>
inline void check_same_size(const mt_vector<U> & v1, const mt_vector<V> & v2, const char* caller)
{
  if(v1.size() != v2.size()) {
    fprintf(stderr, "%s error: vectors not of same size! Aborting!\n", caller);
    exit(1);
  }
}

template<class V>
inline void check_size(const mt_vector<V> & v, const size_t n, const char* caller)
{
  if(v.size() != n) {
    fprintf(stderr, "%s error: vector not of size %ld! Aborting!\n", caller, (long int) n);
    exit(1);
  }
}

inline void check_size(const size_t n1, const size_t n2, const char* caller)
{
  if(n1 != n2) {
    fprintf(stderr, "%s error: size missmatch %ld != %ld! Aborting!\n",
            caller, (long int) n1, (long int) n2);
    exit(1);
  }
}

inline void check_nonzero(const size_t n, const char* caller)
{
  if(n == 0) {
    fprintf(stderr, "%s error: Nonzero size check failed! Aborting!\n", caller);
    exit(1);
  }
}

inline void check_file_exists(const std::string file, const char* caller)
{
  bool exists = access(file.c_str(), F_OK) == 0;

  if(!exists) {
    fprintf(stderr, "%s error: File %s does not exist! Aborting!\n", caller, file.c_str());
    exit(1);
  }
}

inline void check_elem_type(elem_t ct, elem_t rt, const char* caller)
{
  if(ct != rt) {
    fprintf(stderr, "%s error: Unsupported element type! \n", caller);
    exit(1);
  }
}

inline void check_condition(bool cond, const char* cond_str, const char* caller)
{
  if(!cond) {
    fprintf(stderr, "%s error: Condition (%s) not fullfilled! \n", caller, cond_str);
    exit(1);
  }
}

#endif
