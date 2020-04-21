/**
* @file lookup_table.hpp
* @brief Minimal lookup table class implementation.
* @author Aurel Neic
* @version 
* @date 2017-09-12
*/


#ifndef _LOOKUP_TABLE
#define _LOOKUP_TABLE

#include <iostream>
#include <map>
#include <vector>

typedef long int la_int;

template<class V>
class lookup_table
{
  private:
  std::vector<V> _yvals;
  V   _x_min;
  V   _x_max;
  la_int _num_int;

  la_int _delta;

  inline la_int getIdx(V x) const
  {
    la_int i = la_int((x - _x_min) * _num_int) / _delta;
    if(i < 0 || i > (_num_int-1) ) {
      fprintf(stderr, "lookup table error: requested value out of bounds!"
                      "%g not in [%g, %g].\n", x, _x_min, _x_max);
      exit(1);
    }
    return i;
  }

  public:
  lookup_table(V xmin, V xmax, std::vector<V> & yvals) :
    _yvals(yvals),
    _x_min(xmin),
    _x_max(xmax),
    _num_int(yvals.size())
  {
    _delta = _x_max - _x_min;
  }

  inline V operator() (V x) const
  {
    const la_int i = getIdx(x);
    return _yvals[i];
  }

  inline V linear_interp(V x) const
  {
    const la_int i = getIdx(x);
    const V d = V(_delta) / V(_num_int);
    const V xs = _x_min + i*d, xe = xs + d, ys = _yvals[i], ye = _yvals[i+1];

    return ys + (ye - ys) * (x - xs) / (xe - xs);
  }
};


template<class V>
class dyn_lookup_table
{
  private:
  std::map<V,V> xypairs;
  typename std::map<V,V>::iterator start;
  typename std::map<V,V>::iterator end;

  public:
  /// construct interpolation from x-y pair vector
  dyn_lookup_table(const mt_vector<V> & xval, const mt_vector<V> & yval)
  {
    assert(xval.size() > 1 && xval.size() == yval.size());

    for(size_t i=0; i<xval.size(); i++)
      xypairs[xval[i]] = yval[i];

    start = xypairs.begin();
    end   = xypairs.end(); end--;
  }
  /// construct interpolation from tracefile
  dyn_lookup_table(const std::string tracefile)
  {
    FILE* fd = fopen(tracefile.c_str(), MT_FOPEN_READ);
    int numvals = 0, numread = 0;

    const int bufsize = 2048;
    char buffer[bufsize];
    char* ptr;
    float x, y;
    int line = 0;

    if(fd) {
      ptr = fgets(buffer, bufsize, fd);
      if(ptr != NULL) numread = sscanf(buffer, "%d", &numvals);
      if(numread != 1) {
        fprintf(stderr, "Bad read in file %s, line %d! Aborting!\n", tracefile.c_str(), line);
        fclose(fd);
        exit(1);
      }
      line++;

      for(int i=0; i<numvals; i++) {
        ptr = fgets(buffer, bufsize, fd);
        if(ptr != NULL) numread = sscanf(buffer, "%f %f", &x, &y);
        if(numread == 2) {
          xypairs[x] = y;
        }
        else {
          fprintf(stderr, "Bad read in file %s, line %d! Aborting!\n", tracefile.c_str(), line);
          fclose(fd);
          exit(1);
        }
        line++;
      }
      start = xypairs.begin();
      end   = xypairs.end(); end--;
    }
  }

  /// evaluate interpolation function
  V operator()(const V x) const
  {
    if(x <= start->first) return start->second;
    if(x >= end->first)   return end->second;

    typename std::map<V,V>::const_iterator up  = xypairs.upper_bound(x);
    typename std::map<V,V>::const_iterator low = up; --low;

    V xs = low->first , xe = up->first ;
    V ys = low->second, ye = up->second;

    V fctr = (x - xs) / (xe - xs);
    return ys + (ye - ys)*fctr;
  }

  void normalize()
  {
    V maxy = fabs(start->second);

    // compute maximum absolute value
    for(const std::pair<V,V> & xy : xypairs) {
      V val = fabs(xy.second);
      if(maxy < val) maxy = val;
    }

    // divide through maximum abs val
    for(typename std::map<V,V>::iterator it = xypairs.begin(); it != xypairs.end(); ++it)
      it->second /= maxy;
  }

  V start_x() {
    return start->first;
  }
  V end_x() {
    return end->first;
  }
};



#endif


