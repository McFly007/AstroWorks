#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <fastbvh>
#include <vector>

namespace py = pybind11;

std::vector<Tri>
make_tris(std::vector<float> const & V, std::vector<int> const & F)
{
  std::vector<Tri> tris;
  tris.reserve(F.size()/3);
  for (int i = 0; i < F.size(); i += 3) {
    int i0 = F[i], i1 = F[i + 1], i2 = F[i + 2];
    Vector3 v0 {V[3*i0], V[3*i0 + 1], V[3*i0 + 2]};
    Vector3 v1 {V[3*i1], V[3*i1 + 1], V[3*i1 + 2]};
    Vector3 v2 {V[3*i2], V[3*i2 + 1], V[3*i2 + 2]};
    Vector3 n = normalize((v1 - v0)^(v2 - v0));
    tris.emplace_back(v0, v1, v2, n, i/3);
  }
  return tris;
}

std::vector<Object *>
make_objects(std::vector<Tri> const & tris)
{
  std::vector<Object *> objects;
  objects.reserve(tris.size());
  for (int i = 0; i < tris.size(); ++i) {
    objects.push_back((Object *) &tris[i]);
  }
  return objects;
}

struct BVH_wrapper {
  BVH_wrapper(std::vector<float> const & V, std::vector<int> const & F):
    _tris {make_tris(V, F)},
    _objects {make_objects(_tris)},
    _bvh {&_objects}
  {}

  int intersect(Ray const & ray, IntersectionInfo * isect, bool occ) const {
    auto hit = _bvh.getIntersection(ray, isect, occ);
    return hit ? static_cast<Tri const *>(isect->object)->index : -1;
  }

private:
  std::vector<Tri> _tris;
  std::vector<Object *> _objects;
  BVH _bvh;
};

PYBIND11_MODULE(pybvh, m) {
  py::class_<BVH_wrapper>(m, "BVH")
    .def(py::init<std::vector<float> const &, std::vector<int> const &>())
    .def("intersect", [&] (BVH_wrapper const & bvh,
                           std::array<float, 3> const & p,
                           std::array<float, 3> const & n,
                           bool occlusion) {
                        Ray ray {Vector3 {p[0], p[1], p[2]},
                                 Vector3 {n[0], n[1], n[2]}};
                        IntersectionInfo isect;
                        return bvh.intersect(ray, &isect, occlusion);
                      });
}
