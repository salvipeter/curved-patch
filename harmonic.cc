#include "harmonic.hh"

#include <algorithm>
#include <cmath>

#include "curved-domain.hh"

ConstrainedHarmonic::ConstrainedHarmonic(size_t levels) : levels_(levels) {
  size_ = std::pow(2, levels_);
}

ConstrainedHarmonic::~ConstrainedHarmonic() {
}

Point2D
ConstrainedHarmonic::harmonicMap(size_t i, const Point2D &uv) const {
  double x = uv[0] * size_, y = uv[1] * size_, value;
  int u = std::round(x), v = std::round(y);
  auto bc = [&](size_t j) {
              const auto &m = maps_[j];
              value = m[v*size_+u].value * (1.0 - y + v) * (1.0 - x + u);
              value += m[(v+1)*size_+u].value * (y - v) * (1.0 - x + u);
              value += m[v*size_+u+1].value * (1.0 - y + v) * (x - u);
              value += m[(v+1)*size_+u+1].value * (y - v) * (x - u);
              return value;
            };
  Point2D sd;
  double bi = bc(i), bi_1 = bc(prev(i));
  double denom = bi + bi_1;
  if (denom < epsilon)
    sd[0] = 0.0;                // kutykurutty
  else
    sd[0] = bi / denom;
  sd[1] = 1.0 - bi - bi_1;
  return sd;
}

Point2D
ConstrainedHarmonic::mapToRibbon(size_t i, const Point2D &uv) const {
  Point2D sd = harmonicMap(     i , uv);
  double s_1 = harmonicMap(prev(i), uv)[0];
  double s1  = harmonicMap(next(i), uv)[0];

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

namespace {

  void solve(HarmonicMap grid, size_t level) {
    size_t n = (size_t)std::pow(2, level);
    if (level > 3) {
      // Generate a coarser grid and solve that first to get good starting values
      size_t level1 = level - 1, n1 = (size_t)std::pow(2, level1);
      HarmonicMap grid1(n1 * n1);
      for (size_t i = 0; i < n1; ++i)
        for (size_t j = 0; j < n1; ++j) {
          grid1[j*n1+i].value = 0.0;
          int count = 0;
          if (grid[2*j*n+2*i].boundary) {
            ++count;
            grid1[j*n1+i].value += grid[2*j*n+2*i].value;
          }
          if (grid[2*j*n+2*i+1].boundary) {
            ++count;
            grid1[j*n1+i].value += grid[2*j*n+2*i+1].value;
          }
          if (grid[(2*j+1)*n+2*i].boundary) {
            ++count;
            grid1[j*n1+i].value += grid[(2*j+1)*n+2*i].value;
          }
          if (grid[(2*j+1)*n+2*i+1].boundary) {
            ++count;
            grid1[j*n1+i].value += grid[(2*j+1)*n+2*i+1].value;
          }
          if (count > 0) {
            grid1[j*n1+i].boundary = true;
            grid1[j*n1+i].value /= (double)count;
          }
          else
            grid1[j*n1+i].boundary = false;
        }
      solve(grid1, level1);
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
          if (!grid[j*n+i].boundary)
            grid[j*n+i].value = grid1[(j/2)*n1+i/2].value;
    }

    double change;
    do {
      change = 0.0;
      size_t count = 0, index = n + 1;
      for (size_t j = 1, n_1 = n - 1; j < n_1; ++j) {
        for (size_t i = 1; i < n_1; ++i, ++index)
          if (!grid[index].boundary) {
            double value = 0.0;
            value += grid[index-n].value;
            value += grid[index-1].value;
            value += grid[index+n].value;
            value += grid[index+1].value;
            value /= 4.0;
            change += std::abs(grid[index].value - value);
            grid[index].value = value;
            ++count;
          }
        index += 2;
      }
      change /= (double)count;
    } while (change > 1.0e-5);  // kutykurutty
  }

}

void
ConstrainedHarmonic::update() {
  const size_t resolution = 100;
  const auto &curves = dynamic_cast<CurvedDomain *>(domain_.get())->boundaries();
  n_ = curves.size();
  maps_.clear();
  for (size_t i = 0; i < n_; ++i) {
    HarmonicMap m(size_ * size_);
    for (auto &g : m) {
      g.boundary = false;
      g.value = 0.0;
    }
    for (size_t j = 0; j < n_; ++j) {
      const auto &c = curves[j];
      Point3D from, to = c.eval(0.0);
      if (j == next(i))
        to[2] = 1.0;
      for (size_t k = 1; k <= resolution; ++k) {
        from = to;
        double u = (double)k / resolution;
        to = c.eval(u);
        if (j == i)
          to[2] = u;
        else if (j == next(i))
          to[2] = 1.0 - u;
        // Line drawing:
        int x0 = from[0] * size_, y0 = from[1] * size_;
        int x1 = to[0] * size_, y1 = to[1] * size_;
        double v0 = from[2], v1 = to[2];
        int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
        int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
        int err = (dx > dy ? dx : -dy) / 2, e2;
        if (err == 0) {
          m[y0*size_+x0].boundary = true;
          m[y0*size_+x0].value = v0;
          m[y1*size_+x1].boundary = true;
          m[y1*size_+x1].value = v1;
          continue;
        }
        while (true) {
          double ratio;
          if (err > 0)
            ratio = (double)std::abs(x1 - x0) / (double)dx;
          else
            ratio = (double)std::abs(y1 - y0) / (double)dy;
          m[y0*size_+x0].boundary = true;
          m[y0*size_+x0].value = v0 * ratio + v1 * (1.0 - ratio);
          if (x0 == x1 && y0 == y1) break;
          e2 = err;
          if (e2 > -dx) { err -= dy; x0 += sx; }
          if (e2 <  dy) { err += dx; y0 += sy; }
        }
      }
    }
    solve(m, levels_);
    maps_.push_back(m);
  }
}
