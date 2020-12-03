#ifndef simple_vis_h
#define simple_vis_h

#include <stdbool.h>

// NOTE: functions below with an extra underscore at the end will be
// called from Fortran.

/**
 * Load the Wavefront OBJ file stored at `path' and associate it with
 * the key `obj_key' for later use. This function also builds the
 * acceleration structure used for visibility testing. After it is
 * called, the functions `sv_free_obj' and `sv_query_vis' can be
 * called using the key `obj_key'.
 *
 * Note: this function loads the first mesh in the OBJ file
 * (i.e. shape_index == 0).
 */
extern "C"
void sv_load_obj(char const *obj_key, char const *path);

extern "C"
void sv_load_obj_(char const *obj_key, char const *path,
				  int len_obj_key, int len_path);

/**
 * Free the Wavefront OBJ file associated to the key `obj_key'.
 */
extern "C"
void sv_free_obj(char const *obj_key);

extern "C"
void sv_free_obj_(char const *obj_key, int len_obj_key);

/**
 * For the Wavefront OBJ file associated with they key `obj_key', test
 * for occlusion. The float pointers are both assumed to point to
 * contiguous memory containing 3-vectors.
 *
 * For occlusion testing, a ray originating from the vector `r' in the
 * direction of `d' is cast: if it escapes into free space, then
 * `sv_query_occluded' returns 1, otherwise it returns 0.
 */
extern "C"
bool sv_query_occluded(char const *obj_key, float const *r, float const *d);

extern "C"
bool sv_query_occluded_(char const *obj_key, float const *r, float const *d,
						int len_obj_key);

#endif /* simple_vis_h */
