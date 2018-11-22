#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <domain.hh>
#include <parameterization.hh>
#include <surface-corner-based.hh>
#include <surface-generalized-coons.hh>

#include "curved-cb.hh"
#include "curved-gc.hh"
#include "perpendicular-cb.hh"

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

void fixMesh(TriMesh &mesh, const CurveVector &cv, size_t resolution) {
  size_t index = 0;
  for (const auto &c : cv)
    for (size_t i = 0; i < resolution; ++i) {
      double u = (double)i / resolution;
      mesh[index++] = c->eval(u);
    }
}

void
domainEval(const std::shared_ptr<Surface> &surf, size_t resolution, std::string filename) {
  auto mesh = surf->domain()->meshTopology(resolution);
  auto uvs = surf->domain()->parameters(resolution);
  std::vector<Point2DVector> points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) {
                   auto sds = surf->parameterization()->mapToRibbons(uv);
                   sds.push_back(uv);
                   return sds;
                 });
  std::ofstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }
  for (const auto &p : points) {
    f << 'v';
    for (auto coord : p)
      f << ' ' << coord[0] << ' ' << coord[1];
    f << std::endl;
  }
  for (const auto &t : mesh.triangles())
    f << "f " << t[0] + 1 << ' ' << t[1] + 1 << ' ' << t[2] + 1 << std::endl;
  f.close();
}

void surfaceTest(std::string name, std::shared_ptr<Surface> &&surf, const CurveVector &cv,
                 std::string filename, size_t resolution, bool fix_mesh = false) {
  std::cout << name << ":" << std::endl;
  std::chrono::steady_clock::time_point begin, end;

  begin = std::chrono::steady_clock::now();
  surf->setCurves(cv);
  surf->setupLoop();
  surf->update();
  end = std::chrono::steady_clock::now();
  std::cout << "  Setup time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;

  if (name == "CCB") {
    begin = std::chrono::steady_clock::now();
    domainEval(surf, resolution, filename + "-domain.obj");
    end = std::chrono::steady_clock::now();
    std::cout << "  Domain output time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << "ms" << std::endl;
  }
  
  begin = std::chrono::steady_clock::now();
  auto mesh = surf->eval(resolution);
  end = std::chrono::steady_clock::now();
  std::cout << "  Evaluation time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;

  if (fix_mesh)
    fixMesh(mesh, cv, resolution); // computes exact boundaries
  mesh.writeOBJ(filename + "-" + name + ".obj");
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0]
              << " basename [resolution]" << std::endl;
    return 1;
  }
  std::string fname(argv[1]);

  CurveVector cv = readLOP(fname + ".lop");
  if (cv.empty()) {
    std::cerr << "Cannot read file: " << argv[1] << std::endl;
    return 2;
  }
  
  size_t resolution = 30;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  // surfaceTest("CGC", std::make_shared<CurvedGC>(), cv, fname, resolution, true);
  surfaceTest("CCB", std::make_shared<CurvedCB>(), cv, fname, resolution, true);
  // surfaceTest("GC", std::make_shared<Transfinite::SurfaceGeneralizedCoons>(),
  //             cv, fname, resolution);
  // surfaceTest("CB", std::make_shared<Transfinite::SurfaceCornerBased>(),
  //             cv, fname, resolution);
  surfaceTest("PCB", std::make_shared<PerpCB>(), cv, fname, resolution);

  return 0;
}
