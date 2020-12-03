#include "simple_vis.h"

#define TINYOBJLOADER_IMPLEMENTATION

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include "./fastbvh/include/fastbvh"
#include "./tinyobjloader/tiny_obj_loader.h"

static std::vector<Object *> load_tris(char const * path);
static void free_tris(std::vector<Object *> & tris);

struct BVH_wrapper {
	BVH_wrapper(char const * path): _objs {load_tris(path)}, _bvh {&_objs} {}
	~BVH_wrapper() { free_tris(_objs); }
	bool query_occluded(Vector3 r, Vector3 d) const {
		IntersectionInfo unused;
		return !_bvh.getIntersection(Ray(r, d), &unused, true);
	}
private:
	std::vector<Object *> _objs;
	BVH _bvh;
};


static std::unordered_map<std::string, BVH_wrapper> _bvh_wrappers;

void sv_load_obj(char const * obj_key, char const * path) {
	_bvh_wrappers.emplace(std::string(obj_key), path);
}

void sv_free_obj(char const * obj_key) {
	_bvh_wrappers.erase(std::string(obj_key));
}

bool sv_query_occluded(char const * obj_key, float const * r, float const * d) {
	return _bvh_wrappers.at(obj_key).query_occluded(
		Vector3(r[0], r[1], r[2]),
		Vector3(d[0], d[1], d[2]));
}

static std::vector<Object *> load_tris(char const * path) {
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;

  std::string err;
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, path);
  if (!err.empty()) {
    std::cerr << err << std::endl;
  }
  if (!ret) {
    std::cerr << "Failed to load " << path << std::endl;
    std::exit(EXIT_FAILURE);
  }

  tinyobj::shape_t & shape = shapes[0];
  auto const & F = shape.mesh.indices;
  auto const & V = attrib.vertices;

  std::vector<Object *> objects;
  Vector3 vec[3];
  int vi[3];
  for (size_t i = 0; i < F.size(); i += 3) {
	for (int j = 0; j < 3; ++j) {
		vi[j] = F[i + j].vertex_index;
		for (int k = 0; k < 3; ++k) {
			vec[j][k] = V[3*vi[j] + k];
		}
	}
	// std::cerr << vec[0] << std::endl;
	
    Vector3 n = normalize((vec[1] - vec[0])^(vec[2] - vec[0]));
    objects.push_back(new Tri(vec[0], vec[1], vec[2], n, i/3));
  }
  return objects;
}

static void free_tris(std::vector<Object *> & tris) {
	for (Object * obj: tris) {
		delete static_cast<Tri *>(obj);
	}
}

void sv_load_obj_(char const *obj_key, char const *path,
						  int len_obj_key, int len_path)
{
	// Unused
	(void) len_obj_key;
	(void) len_path;

	sv_load_obj(obj_key, path);
}


void sv_free_obj_(char const *obj_key, int len_obj_key)
{
	// Unused
	(void) len_obj_key;

	sv_free_obj(obj_key);
}

bool sv_query_occluded_(char const *obj_key, float const * r, float const * d,
				   int len_obj_key)
{
	// Unused
	(void) len_obj_key;

	return sv_query_occluded(obj_key, r, d);
}
