#ifndef BVH_H
#define BVH_H
#include "hitsphere.h"
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

template <typename T>
T min(const T a, const T b)
{
	return a < b ? a : b;
}

template <typename T>
T max(const T& a, const T& b)
{
	return a > b ? a : b;
}

class aabb
{
public:
	aabb() { _min = vec3(1e12, 1e12, 1e12);  _max = vec3(-1e12, -1e12, -1e12); }
	aabb(const vec3& a, const vec3& b) { _min = a; _max = b; };

	vec3 min() const { return _min; }
	vec3 max() const { return _max; }

	bool hit(const ray& r, double tmin, double tmax) const
	{
		for (int a = 0; a < 3; a++)
		{
			double invD = 1.0f / r.direction()[a];
			double t0 = (min()[a] - r.origin()[a]) * invD;
			double t1 = (max()[a] - r.origin()[a]) * invD;
			if (invD < 0.0f)
				std::swap(t0, t1);
			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;
			if (tmax <= tmin)
				return false;
		}
		return true;
	}
private:
	vec3 _min;
	vec3 _max;

};

aabb surrounding_box(aabb box0, aabb box1)
{
	vec3 small(fmin(box0.min().x(), box1.min().x()),
		fmin(box0.min().y(), box1.min().y()),
		fmin(box0.min().z(), box1.min().z()));
	vec3 big(fmax(box0.max().x(), box1.max().x()),
		fmax(box0.max().y(), box1.max().y()),
		fmax(box0.max().z(), box1.max().z()));
	return aabb(small, big);
}

class bvh_node : public hittable
{
public:
	bvh_node() {};
	bvh_node(vector<point3>& vertices, vector<vec3> &normals, vector<int>::iterator begin, vector<int>::iterator end, shared_ptr<material> m);
	shared_ptr<bvh_node> left;
	shared_ptr<bvh_node> right;
	aabb box;
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
};

bool bvh_node::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
	if (!box.hit(r, t_min, t_max))
	{
		return false;
	}
	bool hit_left = false, hit_right = false;
	if (left != nullptr)
	{
		hit_left = left->hit(r, t_min, t_max, rec);
	}
	if (right != nullptr)
	{
		hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);
	}

	return hit_left || hit_right;
}

class triangle : public bvh_node
{
public:
	triangle() {}

	triangle(point3 a, point3 b, point3 c, vec3 na, vec3 nb, vec3 nc, shared_ptr<material> m)
		: A(a), B(b), C(c), mat_ptr(m), normalA(na), normalB(nb), normalC(nc)
	{
		left = nullptr;
		right = nullptr;
		box = aabb(vec3(min(min(A.x(), B.x()), C.x()), min(min(A.y(), B.y()), C.y()), min(min(A.z(), B.z()), C.z())),
			vec3(max(max(A.x(), B.x()), C.x()), max(max(A.y(), B.y()), C.y()), max(max(A.z(), B.z()), C.z())));

	};

	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
private:
	point3 A;
	point3 B;
	point3 C;
	vec3 normalA;
	vec3 normalB;
	vec3 normalC;
	shared_ptr<material> mat_ptr;
};

bool triangle::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
	vec3 e1 = B - A;
	vec3 e2 = C - A;
	vec3 s = r.origin() - A;
	vec3 s1 = cross(r.direction(), e2);
	vec3 s2 = cross(s, e1);
	double t = dot(s2, e2) / dot(s1, e1);
	double b1 = dot(s1, s) / dot(s1, e1);
	double b2 = dot(s2, r.direction()) / dot(s1, e1);
	if (t > t_min && t < t_max && b1 >= 0 && b2 >= 0 && b1 + b2 <= 1)
	{
		rec.t = t;
		rec.p = r.at(t);
		// bay centric coordinates
		rec.normal = unit_vector(b1 * normalB + b2 * normalC + (1 - b1 - b2) * normalA);
	    
		// or use the face normal
		//rec.normal = unit_vector(cross(e1, e2));
		rec.mat_ptr = mat_ptr;
		return true;
	}
	return false;
}

bvh_node::bvh_node(vector<point3>& vertices, vector<vec3> &normals, vector<int>::iterator begin, vector<int>::iterator end, shared_ptr<material> m)
{
	auto diff = end - begin;
	if (diff < 3)
	{
		left = nullptr;
		right = nullptr;
		box = aabb();
	}
	else if (end - begin == 3)
	{
		left = make_shared<triangle>(vertices[*(begin)], vertices[*(begin + 1)], vertices[*(begin + 2)], normals[*(begin)], normals[*(begin + 1)], normals[*(begin + 2)], m);
		right = nullptr;
		box = aabb(vec3(min(min(vertices[*(begin)].x(), vertices[*(begin + 1)].x()), vertices[*(begin + 2)].x()),
			min(min(vertices[*(begin)].y(), vertices[*(begin + 1)].y()), vertices[*(begin + 2)].y()),
			min(min(vertices[*(begin)].z(), vertices[*(begin + 1)].z()), vertices[*(begin + 2)].z())),
			vec3(max(max(vertices[*(begin)].x(), vertices[*(begin + 1)].x()), vertices[*(begin + 2)].x()),
				max(max(vertices[*(begin)].y(), vertices[*(begin + 1)].y()), vertices[*(begin + 2)].y()),
				max(max(vertices[*(begin)].z(), vertices[*(begin + 1)].z()), vertices[*(begin + 2)].z())));
	}
	else if (diff > 3 && diff < 6)
	{
		left = make_shared<triangle>(vertices[*(begin)], vertices[*(begin + 1)], vertices[*(begin + 2)], normals[*(begin)], normals[*(begin + 1)], normals[*(begin + 2)], m);
		right = nullptr;
		box = aabb(vec3(min(min(vertices[*(begin)].x(), vertices[*(begin + 1)].x()), vertices[*(begin + 2)].x()),
			min(min(vertices[*(begin)].y(), vertices[*(begin + 1)].y()), vertices[*(begin + 2)].y()),
			min(min(vertices[*(begin)].z(), vertices[*(begin + 1)].z()), vertices[*(begin + 2)].z())),
			vec3(max(max(vertices[*(begin)].x(), vertices[*(begin + 1)].x()), vertices[*(begin + 2)].x()),
				max(max(vertices[*(begin)].y(), vertices[*(begin + 1)].y()), vertices[*(begin + 2)].y()),
				max(max(vertices[*(begin)].z(), vertices[*(begin + 1)].z()), vertices[*(begin + 2)].z())));
	}
	else if (diff == 6)
	{
		left = make_shared<triangle>(vertices[*(begin)], vertices[*(begin + 1)], vertices[*(begin + 2)], normals[*(begin)], normals[*(begin + 1)], normals[*(begin + 2)], m);
		right = make_shared<triangle>(vertices[*(begin + 3)], vertices[*(begin + 4)], vertices[*(begin + 5)], normals[*(begin + 3)], normals[*(begin + 4)], normals[*(begin + 5)], m);
		box = aabb(vec3(min(min(vertices[*(begin)].x(), vertices[*(begin + 1)].x()), vertices[*(begin + 2)].x()),
			min(min(vertices[*(begin)].y(), vertices[*(begin + 1)].y()), vertices[*(begin + 2)].y()),
			min(min(vertices[*(begin)].z(), vertices[*(begin + 1)].z()), vertices[*(begin + 2)].z())),
			vec3(max(max(vertices[*(begin)].x(), vertices[*(begin + 1)].x()), vertices[*(begin + 2)].x()),
				max(max(vertices[*(begin)].y(), vertices[*(begin + 1)].y()), vertices[*(begin + 2)].y()),
				max(max(vertices[*(begin)].z(), vertices[*(begin + 1)].z()), vertices[*(begin + 2)].z())));
		box = surrounding_box(box, aabb(vec3(min(min(vertices[*(begin + 3)].x(), vertices[*(begin + 4)].x()), vertices[*(begin + 5)].x()),
			min(min(vertices[*(begin + 3)].y(), vertices[*(begin + 4)].y()), vertices[*(begin + 5)].y()),
			min(min(vertices[*(begin + 3)].z(), vertices[*(begin + 4)].z()), vertices[*(begin + 5)].z())),
			vec3(max(max(vertices[*(begin + 3)].x(), vertices[*(begin + 4)].x()), vertices[*(begin + 5)].x()),
				max(max(vertices[*(begin + 3)].y(), vertices[*(begin + 4)].y()), vertices[*(begin + 5)].y()),
				max(max(vertices[*(begin + 3)].z(), vertices[*(begin + 4)].z()), vertices[*(begin + 5)].z()))));
	}
	else
	{
		vector<int>::iterator mid = begin + (end - begin) / 2;
		auto count = mid - begin;
		if (count % 3 != 0)
		{
			mid = mid + 3 - (count % 3);
		}
		left = make_shared<bvh_node>(vertices, normals, begin, mid, m);
		right = make_shared<bvh_node>(vertices, normals, mid, end, m);

		//left = make_shared<bvh_node>(vertices, begin, begin + 3, m);
		//right = make_shared<bvh_node>(vertices, begin + 3, end, m);
		box = surrounding_box(left->box, right->box);
	}
}

class triangle_mesh : public hittable
{
public:
	triangle_mesh() {}
	triangle_mesh(const char* filename, shared_ptr<material> m, point3 pos, double scalar, vec3 axis, double degree);
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
private:
	std::vector<point3> vertices;
	std::vector<int> indices;
	std::vector<vec3> normals;
	aabb box;
	shared_ptr<material> mat_ptr;
	shared_ptr<bvh_node> root;
};

triangle_mesh::triangle_mesh(const char* filename, shared_ptr<material> m, point3 pos, double scalar, vec3 axis = vec3(1, 0, 0), double degree = 90)
{
	mat_ptr = m;
	point3 min;
	point3 max;
	max = point3(-1e12, -1e12, -1e12);
	min = point3(1e12, 1e12, 1e12);
	// use assimp
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(filename, aiProcess_Triangulate);
	if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
	{
		std::cout << "ERROR::ASSIMP::" << importer.GetErrorString() << std::endl;
		return;
	}

	// add vertices and indices
	for (unsigned int i = 0; i < scene->mNumMeshes; i++)
	{
		aiMesh* mesh = scene->mMeshes[i];
		for (unsigned int j = 0; j < mesh->mNumVertices; j++)
		{
			aiVector3D vertex = mesh->mVertices[j];
			point3 vertex_rotated = rotate(scalar * vec3(vertex.x, vertex.y, vertex.z), axis, degree) + pos;
			vertices.push_back(vertex_rotated);
			// normal
			aiVector3D normal = mesh->mNormals[j];
			normals.push_back(rotate(vec3(normal.x, normal.y, normal.z), axis, degree));
			//normals.push_back(vec3(normal.x, normal.y, normal.z));
		}
		for (unsigned int j = 0; j < mesh->mNumFaces; j++)
		{
			aiFace face = mesh->mFaces[j];
			indices.push_back(face.mIndices[0]);
			indices.push_back(face.mIndices[1]);
			indices.push_back(face.mIndices[2]);
		}
	}
	// normal on vertices

	root = make_shared<bvh_node>(vertices, normals, indices.begin(), indices.end(), m);
	box = root->box;
}

bool triangle_mesh::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
	double t = 1e12;
	// bound
	if (!box.hit(r, t_min, t_max))
	{
		return false;
	}

	return root->hit(r, t_min, t_max, rec);
}

#endif