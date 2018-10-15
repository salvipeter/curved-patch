#include <cstdlib>
#include <fstream>
#include <iostream>

#include "curved-gc.hh"

CurveVector readLOP(std::string filename) {
  size_t n, deg, nk, nc;
  DoubleVector knots;
  PointVector cpts;
  CurveVector result;
  std::ifstream f(filename);
  f >> n;
  result.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    f >> deg;
    f >> nk;
    knots.resize(nk);
    for (size_t j = 0; j < nk; ++j)
      f >> knots[j];
    f >> nc;
    cpts.resize(nc);
    for (size_t j = 0; j < nc; ++j)
      f >> cpts[j][0] >> cpts[j][1] >> cpts[j][2];
    result.push_back(std::make_shared<BSCurve>(deg, knots, cpts));
  }
  if (!f)
    return CurveVector();
  return result;
}

int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    std::cerr << "Usage: " << argv[0] << " infile.lop outfile.obj [resolution]" << std::endl;
    return 1;
  }
  CurveVector cv = readLOP(argv[1]);
  if (cv.empty()) {
    std::cerr << "Cannot read file: " << argv[1] << std::endl;
    return 2;
  }
  size_t resolution = 30;
  if (argc == 4)
    resolution = std::atoi(argv[3]);

  auto surf = std::make_shared<CurvedGC>();
  surf->setCurves(cv);
  surf->setupLoop();
  surf->update();
  surf->eval(resolution).writeOBJ(argv[2]);

  return 0;
}
