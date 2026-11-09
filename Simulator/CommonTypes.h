#ifndef __H__LIB__COMMON__
#define __H__LIB__COMMON__

#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <ranges>



typedef int32_t								i32;
typedef uint32_t							u32;
typedef double								f64;
typedef std::string							Str;
typedef std::vector<Str>					vec_Str;
typedef std::vector<i32>					vec_i32;
typedef std::vector<bool>					vec_bool;
typedef std::vector<f64>					vec_f64;
typedef std::vector<vec_i32>				vec2D_i32;
typedef std::vector<vec_f64>				vec2D_f64;
typedef std::vector<vec2D_f64>				vec3D_f64;
typedef std::map<Str, i32>					map_Str_i32;
typedef std::map<i32, Str>					map_i32_Str;
typedef std::map<i32, f64>					map_i32_f64;
typedef std::map<Str, Str>					map_Str_Str;
typedef std::set<Str>						set_Str;
typedef std::set<i32>						set_i32;
typedef std::ranges::iota_view<i32, i32>	irange;



#endif