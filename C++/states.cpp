#include <stddef.h>
#include "states.h"

states states::operator+(const states& s)
{
  states out; 

  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < hiv_idx; m++)
	for(size_t u = 0; u < art_idx; u++)
	  out.X[g][a][m][u] = X[g][a][m][u] + s.X[g][a][m][u];

  out.hiv_idx = hiv_idx;
  out.art_idx = art_idx;

  return out;
}

states states::operator*(const double& c)
{
  states out; 

  for(size_t g = 0; g < NG; g++)
    for(size_t a = 0; a < AG; a++)
      for(size_t m = 0; m < hiv_idx; m++)
	for(size_t u = 0; u < art_idx; u++)
	  out.X[g][a][m][u] = X[g][a][m][u] * c;

  out.hiv_idx = hiv_idx;
  out.art_idx = art_idx;

  return out;
}

states & states::operator+=(const states& s)
{
  *this = *this + s;
  return *this;
}
