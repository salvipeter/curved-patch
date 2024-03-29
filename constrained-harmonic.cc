#include "constrained-harmonic.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "curved-domain.hh"

ConstrainedHarmonic::ConstrainedHarmonic(size_t levels) : Harmonic(levels) {
}

ConstrainedHarmonic::~ConstrainedHarmonic() {
}

Point2D
ConstrainedHarmonic::mapToRibbon(size_t i, const Point2D &uv) const {
  Point2D sd = Harmonic::mapToRibbon(     i , uv);
  double s_1 = Harmonic::mapToRibbon(prev(i), uv)[0];
  double s1  = Harmonic::mapToRibbon(next(i), uv)[0];

  // As in Surface::blendSideSingular
  std::vector<double> blends, weights = { sd[1], 1.0 - sd[0], 1.0 - sd[1], sd[0] };
  size_t small = 0;
  for (const auto &w : weights)
    if (w < epsilon)
      ++small;
  if (small > 0) {
    double val = 1.0 / small;
    for (const auto &w : weights)
      blends.push_back(w < epsilon ? val : 0.0);
  } else {
    double denominator = 0.0;
    for (const auto &w : weights) {
      blends.push_back(std::pow(w, -2));
      denominator += blends.back();
    }
    std::transform(blends.begin(), blends.end(), blends.begin(),
                   [denominator](double x) { return x / denominator; });
  }

  sd[1] = sd[1] * (blends[0] + blends[2]) + s1 * blends[1] + (1.0 - s_1) * blends[3];
  return sd;
}
