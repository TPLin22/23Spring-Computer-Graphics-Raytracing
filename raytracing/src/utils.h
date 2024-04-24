#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <iostream>
#include <vector>

using std::shared_ptr;
using std::make_shared;
using std::sqrt;
using std::vector;

// Constants
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926;

inline double degrees_to_radians(double degrees)
{
	return degrees * pi / 180.0;
}

inline double clamp(double x, double min, double max)
{
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

inline double random_double()
{
	// Returns a random real in [0,1).
	return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max)
{
	// Returns a random real in [min,max).
	return min + (max - min) * random_double();
}

inline int random_int(int min, int max)
{
	// Returns a random integer in [min,max].
	return static_cast<int>(random_double(min, max + 1));
}

//using std::sqrt;
//using std::fabs;

class vec3
{
public:
	vec3() : e{ 0,0,0 } {}
	vec3(double e0, double e1, double e2) : e{ e0, e1, e2 } {}

	double x() const { return e[0]; }
	double y() const { return e[1]; }
	double z() const { return e[2]; }

	vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
	double operator[](int i) const { return e[i]; }
	double& operator[](int i) { return e[i]; }

	vec3& operator+=(const vec3& v)
	{
		e[0] += v.e[0];
		e[1] += v.e[1];
		e[2] += v.e[2];
		return *this;
	}

	vec3& operator*=(const double t)
	{
		e[0] *= t;
		e[1] *= t;
		e[2] *= t;
		return *this;
	}

	vec3& operator/=(const double t)
	{
		return *this *= 1 / t;
	}

	double length() const {
		return sqrt(length_squared());
	}

	double length_squared() const {
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
	}

	bool near_zero() const
	{
		// Return true if the vector is close to zero in all dimensions.
		const auto s = 1e-8;
		return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
	}

	inline static vec3 random()
	{
		return vec3(random_double(), random_double(), random_double());
	}

	inline static vec3 random(double min, double max) {
		return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
	}

public:
	double e[3];
};

#define point3 vec3
#define color vec3
// 重载
inline std::ostream& operator<<(std::ostream& out, const vec3& v)
{
	return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3& u, const vec3& v)
{
	return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3& u, const vec3& v)
{
	return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3& u, const vec3& v)
{
	return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3& v)
{
	return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline vec3 operator*(const vec3& v, double t)
{
	return t * v;
}

inline vec3 operator/(vec3 v, double t)
{
	return (1 / t) * v;
}

inline double dot(const vec3& u, const vec3& v) {
	return u.e[0] * v.e[0]
		+ u.e[1] * v.e[1]
		+ u.e[2] * v.e[2];
}

inline vec3 cross(const vec3& u, const vec3& v) {
	return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
		u.e[2] * v.e[0] - u.e[0] * v.e[2],
		u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v) {
	return v / v.length();
}

inline vec3 random_in_unit_disk() {
	while (true)
	{
		auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
		if (p.length_squared() >= 1) continue;
		return p;
	}
}

inline vec3 random_in_unit_sphere() {
	while (true) {
		auto p = vec3::random(-1, 1);
		if (p.length_squared() >= 1) continue;
		return p;
	}
}

inline vec3 random_unit_vector() {
	return unit_vector(random_in_unit_sphere());
}

inline vec3 random_in_hemisphere(const vec3& normal) {
	vec3 in_unit_sphere = random_in_unit_sphere();
	if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
		return in_unit_sphere;
	else
		return -in_unit_sphere;
}

// 反射
inline vec3 reflect(const vec3& v, const vec3& n)
{
	return v - 2 * dot(v, n) * n;
}

// 折射
inline vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat)
{
	auto cos_theta = fmin(dot(-uv, n), 1.0);
	vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
	vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
	return r_out_perp + r_out_parallel;
}

class mat3x3
{
public:
	mat3x3() {}
	mat3x3(const vec3& v0, const vec3& v1, const vec3& v2)
	{
		m[0] = v0;
		m[1] = v1;
		m[2] = v2;
	}

	vec3& operator[](int i) { return m[i]; }
	const vec3& operator[](int i) const { return m[i]; }

	vec3 m[3];

	static mat3x3 identity()
	{
		return mat3x3(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1));
	}

	vec3 operator*(const vec3& v) const
	{
		return vec3(dot(m[0], v), dot(m[1], v), dot(m[2], v));
	}
};


// rotate
inline vec3 rotate(const vec3& v, const vec3& axis, double angle)
{
	double rad = angle * pi / 180.0;
	double s = sin(rad);
	double c = cos(rad);
	double oc = 1.0 - c;
	double x = axis.x();
	double y = axis.y();
	double z = axis.z();
	double xs = x * s;
	double ys = y * s;
	double zs = z * s;
	double xy = x * y;
	double xz = x * z;
	double yz = y * z;
	double xoc = x * oc;
	double yoc = y * oc;
	double zoc = z * oc;

	mat3x3 m;
	m[0][0] = x * xoc + c;
	m[0][1] = xy * xoc + zs;
	m[0][2] = xz * xoc - ys;
	m[1][0] = xy * xoc - zs;
	m[1][1] = y * yoc + c;
	m[1][2] = yz * yoc + xs;
	m[2][0] = xz * xoc + ys;
	m[2][1] = yz * yoc - xs;
	m[2][2] = z * zoc + c;

	return m * v;
}

class ray {
public:
	ray() {}
	ray(const point3& origin, const vec3& direction)
		: orig(origin), dir(direction), tm(0)
	{}

	ray(const point3& origin, const vec3& direction, double time)
		: orig(origin), dir(direction), tm(time)
	{}


	point3 origin() const { return orig; }
	vec3 direction() const { return dir; }
	double time() const { return tm; }

	point3 at(double t) const
	{
		return orig + t * dir;
	}

public:
	point3 orig;
	vec3 dir;
	double tm;
};

class camera {
public:
	camera() : camera(point3(0, 0, -1), point3(0, 0, 0), vec3(0, 1, 0), 40, 1, 0, 10) {}

	camera(
		point3 lookfrom,
		point3 lookat,
		vec3   vup,
		double vfov, // vertical field-of-view in degrees
		double aspect_ratio,
		double aperture,
		double focus_dist,
		double _time0 = 0,
		double _time1 = 0
	) {
		auto theta = degrees_to_radians(vfov);
		auto h = tan(theta / 2);
		auto viewport_height = 2.0 * h;
		auto viewport_width = aspect_ratio * viewport_height;

		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);

		origin = lookfrom;
		horizontal = focus_dist * viewport_width * u;
		vertical = focus_dist * viewport_height * v;
		lower_left_corner = origin - horizontal / 2 - vertical / 2 - focus_dist * w;

		lens_radius = aperture / 2;
		time0 = _time0;
		time1 = _time1;
	}

	ray get_ray(double s, double t) const
	{
		vec3 rd = lens_radius * random_in_unit_disk();
		vec3 offset = u * rd.x() + v * rd.y();
		return ray(
			origin + offset,
			lower_left_corner + s * horizontal + t * vertical - origin - offset,
			random_double(time0, time1)
		);
	}

private:
	point3 origin;
	point3 lower_left_corner;
	vec3 horizontal;
	vec3 vertical;
	vec3 u, v, w;
	double lens_radius;
	double time0, time1;
};

#endif
