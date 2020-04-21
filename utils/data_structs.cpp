/**
* @file data_structs.cpp
* @brief Utility functions for geometric operations.
* @author Elias Karabelas
* @version
* @date 2019-07-05
*/
#include <cmath>
#include "data_structs.h"

float inv_sqrtf(const float number)
{
  float y;
#ifdef FAST_INV_SQRT
  uint32_t i;
  const float x2 = number * 0.5F;
  y  = number;
  i  = * ( uint32_t * ) &y;                       // evil floating point bit level hacking
  i  = 0x5f375a86 - ( i >> 1 );               // what the fuck?
  y  = * ( float * ) &i;
  y  = y * ( 1.5F - ( x2 * y * y ) );   // 1st iteration
  y  = y * ( 1.5F - ( x2 * y * y ) );   // 2nd iteration, this can be removed
#else
  y = 1.0F / std::sqrt(number);
#endif
  return y;
}

double inv_sqrtd(const double number)
{
  double y;
#ifdef FAST_INV_SQRT
  uint64_t i;
  const double x2 = number * 0.5;
  y = number;
  i = *(uint64_t *) &y;
  i = 0x5fe6eb50c7b537a9 - (i >> 1);
  y = *(double *) &i;
  y = y * (1.5 - (x2 * y * y));
  y = y * (1.5 - (x2 * y * y));
#else
  y = 1.0 / std::sqrt(number);
#endif
  return y;
}

