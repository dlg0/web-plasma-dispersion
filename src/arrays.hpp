#ifndef ARRAYS_HPP_
#define ARRAYS_HPP_

#include "boost/multi_array.hpp"
#include <complex>

typedef std::complex<double> dcomp;
typedef boost::multi_array<dcomp,2> array2d_dcomp;
typedef array2d_dcomp::index index2d_dcomp;

#endif
