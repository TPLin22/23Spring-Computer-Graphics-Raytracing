#ifndef SPHERE_H
#define SPHERE_H

#include "utils.h"
#include <vector>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
class material;

struct hit_record
{
	// ����
	point3 p;
	// ����
	vec3 normal;
	// ����
	shared_ptr<material> mat_ptr;
	// ray ʱ��
	double t;
	//���滹������
	bool front_face;

	// ���ݹ��߷���normal,�۽�˵����������
	inline void set_face_normal(const ray& r, const vec3& outward_normal)
	{
		// dot С��0, ˵���Ƕ۽ǣ���ray�����෴������
		front_face = dot(r.direction(), outward_normal) < 0;
		normal = front_face ? outward_normal : -outward_normal;
	}
};

class hittable
{
public:
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

// �����е�����
class Scene : public hittable {
public:
	Scene()
	{}
	Scene(shared_ptr<hittable> object)
	{
		add(object);
	}

	void clear()
	{
		objects.clear();
	}
	void add(shared_ptr<hittable> object) { objects.push_back(object); }

	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

public:
	std::vector<shared_ptr<hittable>> objects;
};

bool Scene::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
	hit_record temp_rec;
	bool hit_anything = false;
	auto closest_so_far = t_max;

	for (const auto& object : objects)
	{
		if (object->hit(r, t_min, closest_so_far, temp_rec))
		{
			hit_anything = true;
			closest_so_far = temp_rec.t;
			rec = temp_rec;
		}
	}

	return hit_anything;
};

class material
{
public:
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const = 0;
};

// ������Ϊ�Ǵֲڲ��ʣ�ֱ�ӷ��ر�����ɫ
class lambertian : public material {
public:
	lambertian(const color& a) : albedo(a)
	{}

	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override
	{
		auto scatter_direction = rec.normal + random_unit_vector();

		// Catch degenerate scatter direction
		if (scatter_direction.near_zero())
		{
			scatter_direction = rec.normal;
		}

		scattered = ray(rec.p, scatter_direction);
		attenuation = albedo;
		return true;
	}

public:
	color albedo;
};

// �������ʣ�ȫ����
class metal : public material
{
public:
	metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1)
	{}

	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override
	{
		vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
		scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
		attenuation = albedo;
		return (dot(scattered.direction(), rec.normal) > 0);
	}

public:
	color albedo;
	double fuzz;
};

// ����ʣ��������ʣ��з���������
class dielectric : public material
{
public:
	dielectric(double index_of_refraction) : ir(index_of_refraction) {}

	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override
	{
		attenuation = color(1.0, 1.0, 1.0);
		double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

		vec3 unit_direction = unit_vector(r_in.direction());
		double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

		//���
		bool cannot_refract = refraction_ratio * sin_theta > 1.0;
		vec3 direction;

		//����
		if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
		{
			direction = reflect(unit_direction, rec.normal);
		}
		else // ����
		{
			direction = refract(unit_direction, rec.normal, refraction_ratio);
		}

		scattered = ray(rec.p, direction);
		return true;
	}

public:
	double ir; // Index of Refraction

private:
	static double reflectance(double cosine, double ref_idx) {
		// Use Schlick's approximation for reflectance.
		auto r0 = (1 - ref_idx) / (1 + ref_idx);
		r0 = r0 * r0;
		return r0 + (1 - r0) * pow((1 - cosine), 5);
	}
};

class lightsource : public material
{
public:
	lightsource(const color& a) : albedo(a)
	{}

	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override
	{
		attenuation = albedo;
		scattered = ray(rec.p, vec3(1e12, 1e12, 1e12));
		return true;
	}
private:
	color albedo;
};

class glittering : public material
{
public:
	glittering(const color& a) : albedo(a)
	{}

	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override
	{
		attenuation = albedo;
		scattered = ray(rec.p, vec3(1e12, 1e12, 1e12));
		return true;
	}
private:
	color albedo;
};

class sphere : public hittable
{
public:
	sphere() {}

	sphere(point3 cen, double r, shared_ptr<material> m)
		: center(cen), radius(r), mat_ptr(m) {};

	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

public:
	point3 center;
	double radius;
	shared_ptr<material> mat_ptr;
};


bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
	// �����������Ľ���
	vec3 oc = r.origin() - center;
	auto a = r.direction().length_squared();
	auto half_b = dot(oc, r.direction());
	auto c = oc.length_squared() - radius * radius;

	auto discriminant = half_b * half_b - a * c;
	if (discriminant < 0)
	{
		return false;
	}
	auto sqrtd = sqrt(discriminant);

	// Find the nearest root that lies in the acceptable range.
	auto root = (-half_b - sqrtd) / a;
	if (root < t_min || t_max < root)
	{
		root = (-half_b + sqrtd) / a;
		if (root < t_min || t_max < root)
		{
			return false;
		}
	}

	rec.t = root;
	rec.p = r.at(rec.t);
	vec3 outward_normal = (rec.p - center) / radius;
	rec.set_face_normal(r, outward_normal);
	rec.mat_ptr = mat_ptr;

	return true;
}

class plane : public hittable
{
public:
	plane() {}

	plane(point3 cen, vec3 dir, shared_ptr<material> m)
		: center(cen), normal(dir), mat_ptr(m) {};

	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

public:
	point3 center;
	vec3 normal;
	shared_ptr<material> mat_ptr;
};


bool plane::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
	// ���������ƽ��Ľ���
	double t = dot((center - r.origin()), normal) / dot(r.direction(), normal);
	if (t < t_min || t > t_max)
	{
		return false;
	}

	rec.t = t;
	rec.p = r.at(rec.t);
	rec.set_face_normal(r, normal);
	rec.mat_ptr = mat_ptr;

	return true;
}

color ray_color(const ray& r, const hittable& scene, color background, int depth)
{
	// ��¼�������Ϣ
	hit_record rec;

	// �����˷������޴���������0
	if (depth <= 0)
	{
		return color(1, 1, 1);
	}

	// ���߻�������Ļ�,�ݹ������material ����color��
	if (scene.hit(r, 0.001, infinity, rec))
	{
		// ÿ�η������ray
		ray scattered;
		// ������ɫ��
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
		{
			if (scattered.direction().length_squared() > 1e36)
			{
				return attenuation;
			}
			auto temp = ray_color(scattered, scene, background, depth - 1);
			return attenuation * temp;
		}
		return color(0, 0, 0);
	}

	// û�л��г����е����壬���ر���
	return background;
}

#endif
