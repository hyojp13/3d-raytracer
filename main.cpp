#include <iostream>
#include "FreeImage.h"
#include "readfile.h"
#include "variables.h"

using namespace std;

unsigned int height;
unsigned int width;
int maxdepth;
char* output_file;
int maxverts;
vec3 camera_eye;
vec3 center;
vec3 camera_up;
float camera_fov;
vec3 ambient;
vec3 diffuse;
vec3 specular;
vec3 emission;
float shininess;
vec3 attenuation;
int numused = 0;
float lightposn[30];
float lightcolor[30];
bool is_point[10];
vector<vec3> vertices;
const float EPSILON = 0.00001f;

struct triangle;
struct sphere;
std::vector<triangle> triangles;
std::vector<sphere> spheres;
int depth_counter;
int counter = 0;

vec3 rayCast(const int& x1, const int& y1) {
	float x = x1 + 0.5;
	float y = y1 + 0.5;

	vec3 w = glm::normalize(camera_eye - center);
	vec3 u = glm::normalize(glm::cross(camera_up, w));
	vec3 v = glm::normalize(glm::cross(w, u));	//makes it flipped
	if (counter == 1) {
		cout << glm::to_string(w);
		cout << glm::to_string(camera_up) << glm::to_string(w) << endl;

		cout << glm::to_string(u);

		cout << glm::to_string(v);
	}
	counter++;

	float camera_fovy = glm::radians(camera_fov);

	float alpha = glm::tan(camera_fovy / 2) * width / height * (2 * x - width) / width;
	float beta = glm::tan(camera_fovy / 2) * (height - y * 2) / height;

	vec3 ray = glm::normalize(alpha * u + beta * v - w);
	return ray;
}

bool tri_inter(const vec3& dir, const vec3& cent, const int& index, float& dist, vec3& normal, vec3& point) {
	vec3 A = triangles[index].vertices[0];
	vec3 B = triangles[index].vertices[1];
	vec3 C = triangles[index].vertices[2];
	/*vec3 ba = b - a;
	vec3 ca = c - a;
	vec3 normal = cross(ba, ca);
	float dr = glm::dot(normal, ray);
			// TODO: Handle single sided?
	if (dr == 0) return false;
	float d = dot(normal, a);
	float t = -(dot(normal, camera_eye) + d) / dr;	//camera_eye, not center
	
	if (t < 0) return false;
	vec3 p = camera_eye + (ray * t);
	vec3 pa = p - a;
	float v = dot(normal, cross(ba, pa));
	
	if (v < 0) return false;
	vec3 pb = p - b;
	vec3 cb = c - b;
	float w = dot(normal, cross(cb, pb));
	
	if (w < 0) return false;
	vec3 pc = p - c;
	vec3 ac = a - c;
	float u = dot(normal, cross(ac, pc));
	
	if (u < 0) return false;
	float nlen2 = dot(normal, normal);
	u /= nlen2;
	v /= nlen2;
	w /= nlen2;
	vec3 pt = camera_eye + t * ray;
	
	return true;*/
	vec3 n = glm::normalize(glm::cross(B - A, C - A));
	float denom = glm::dot(n, dir);

	if (denom == 0) return false;

	float t = (glm::dot(A, n) - glm::dot(cent, n)) / denom;
	vec3 P = cent + dir * t;

	float ABC = glm::dot(n, glm::cross(B - A, C - A));
	float PBC = glm::dot(n, glm::cross(B - P, C - P));
	float PCA = glm::dot(n, glm::cross(C - P, A - P));

	float alpha = PBC / ABC;
	float beta = PCA / ABC;
	float gamma = 1.0f - alpha - beta;
	//cout << alpha << ' ' << beta << ' ' << gamma << endl;
	if (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f) {
		dist = t;
		//normal = vec3(glm::transpose(triangles[index].inv_transform) * vec4(n, 1.0f));
		point = vec3(triangles[index].transform * vec4(alpha * A + beta * B + gamma * C, 1.0f));
		normal = glm::normalize(glm::transpose(glm::inverse(glm::mat3(triangles[index].transform))) * glm::normalize(n));
		//point = alpha * A + beta * B + gamma * C;
		return true;
	}
	return false;
}

bool sphere_inter(const vec3& dir, const vec3& cam_cent, const int& index, float& dist, vec3& normal, vec3& point) {

	vec3 cent = spheres[index].center;
	float radius = spheres[index].radius;

	float a = glm::dot(dir, dir);
	float b = 2 * glm::dot(dir, (cam_cent - cent));
	float c = glm::dot(cam_cent - cent, cam_cent - cent) - radius * radius;

	float disc = b * b - 4 * a * c;
	if (disc < 0) return false;
	
	float sol1 = (-1 * b + sqrt(disc)) / (2 * a);
	float sol2 = (-1 * b - sqrt(disc)) / (2 * a);
	if (sol1 < 0 && sol2 < 0) return false;

	if (sol1 > 0 && sol2 > 0) {
		dist = sol1 < sol2 ? sol1 : sol2;
	}
	else if (sol1 > 0 && sol2 < 0) {
		dist = sol1;
	}
	else if (sol1 < 0 && sol2 > 0) {
		dist = sol2;
	}
	vec3 p = dir * dist + cam_cent;
	point = vec3(spheres[index].transform * vec4(dir * dist + cam_cent, 1.0f));
	normal = glm::normalize(glm::transpose(glm::inverse(glm::mat3(spheres[index].transform))) * glm::normalize(p - cent));
	
	return true;
}

float check_inters(const vec3& ray, const vec3& ray_center, int& index, int& object, vec3& normal, vec3& point, bool p = false) {
	float shortest_dist = 100000000.0f;
	float dist = 0.0f;
	index = 0;
	object = 0; // 0 is none, 1 is triangle, 2 is sphere
	vec3 _normal, _point;
	for (int m = 0; m < triangles.size(); ++m) {
		vec4 cent4 = triangles[m].inv_transform * vec4(ray_center, 1.0f);
		vec3 cent(cent4.x, cent4.y, cent4.z);
		vec3 dir = vec3(triangles[m].inv_transform * vec4(ray, 0.0f));

		if (tri_inter(dir, cent, m, dist, _normal, _point)) {
			//cout << i % width << ' ' << (i - i % width) / width << endl;
			if (shortest_dist > dist && dist > 0) {
				normal = _normal;
				point = _point;
				shortest_dist = dist;
				index = m;
				object = 1;
				//if (p) cout << m << endl;
			}
		}
	}
	for (int m = 0; m < spheres.size(); ++m) {
		vec4 cent4 = spheres[m].inv_transform * vec4(ray_center, 1.0f);
		vec3 cent(cent4.x, cent4.y, cent4.z);
		vec4 dir_ = spheres[m].inv_transform * vec4(ray, 0.0f);
		vec3 dir = vec3(dir_.x, dir_.y, dir_.z);

		if (sphere_inter(dir, cent, m, dist, _normal, _point)) {
			if (shortest_dist > dist && dist > 0) {
				normal = _normal;
				point = _point;
				shortest_dist = dist;
				index = m;
				object = 2;
			}
		}
	}
	return dist;
}

float check_light_inter(vec3& ray, const vec3& ray_center, int& index) {
	// check if hits light before hitting object, if so, use light source not the object
	float shortest_dist = 1000000000.0f;
	float dist;
	for (int m = 0; m < numused; ++m) {
		vec3 loc = vec3(lightposn[m * 3], lightposn[m * 3 + 1], lightposn[m * 3 + 2]);
		if (glm::normalize(ray) == glm::normalize(loc - ray_center)) {
			dist = glm::distance(ray_center, loc);
			if (dist < shortest_dist) {
				shortest_dist = dist;
				index = m;
			}
		}
	}

	return shortest_dist;
}

vec3 findColor(const int& index, const vec3& point, const int& object, const vec3& normal, const vec3& init_point=camera_eye, const bool& color_light=true) {
	vec3 lighting(0, 0, 0);
	if (object != 0) {
		if (object == 1) {
			lighting = triangles[index].ambient + triangles[index].emission;
			for (int m = 0; m < numused; ++m) {
				vec3 loc = vec3(lightposn[m * 3], lightposn[m * 3 + 1], lightposn[m * 3 + 2]);	// fix for directional
				vec3 dir;
				if (is_point[m]) dir = glm::normalize(loc - point);
				else dir = glm::normalize(loc);

				vec3 _normal, _point;
				int _index = 0;
				int hits = 0;
				float dist = check_inters(dir, point + EPSILON * dir, _index, hits, _normal, _point);
				bool is_shadow = false;
				if (hits != 0 && dist > 0) {
					if (is_point[m] && distance(point, loc) > dist) {
						is_shadow = true;
					}
					else if (!is_point[m] && distance(point, loc) > 10000000) {
						is_shadow = true;
					}
				}

				vec3 color = vec3(lightcolor[m * 3], lightcolor[m * 3 + 1], lightcolor[m * 3 + 2]);
				float d = 0;
				float atten;
				if (is_point[m]) {
					d = glm::distance(point, loc);
					atten = 1.0f / (attenuation[0] + attenuation[1] * d + attenuation[2] * d * d);
				}
				else {
					atten = 1.0f;
				}

				float t = fmax(glm::dot(dir, normal), 0.0);
				vec3 factor1 = triangles[index].diffuse * t;

				vec3 half = glm::normalize(dir + glm::normalize(init_point - point));
				float y = fmax(glm::dot(half, normal), 0.0f);
				vec3 factor2 = triangles[index].specular * pow(y, triangles[index].shininess);

				// recursion
				dir = glm::normalize(point - init_point);
				vec3 ref_ray = glm::normalize(dir - 2 * glm::dot(dir, normal) * normal);
				vec3 recursion_lighting(0, 0, 0);
				if (depth_counter < maxdepth) {
					++depth_counter;
					float recur_dist = check_inters(ref_ray, point + EPSILON * ref_ray, _index, hits, _normal, _point);
					
					if (hits != 0 && recur_dist > 0) {

						recursion_lighting = findColor(_index, _point, hits, _normal, point, false);
					}
				}
				vec3 factor3 = triangles[index].specular * recursion_lighting;
				//if (!is_point[m]) factor3 = vec3(0, 0, 0);
				
				if (is_shadow) color = vec3(0, 0, 0);
				if (color_light) lighting += color / atten * (factor1 + factor2) + factor3;
				else return color / atten * (factor1 + factor2) + factor3 + lighting;
			}
		}
		else {
			lighting = spheres[index].ambient + spheres[index].emission;
			for (int m = 0; m < numused; ++m) {
				vec3 loc = vec3(lightposn[m * 3], lightposn[m * 3 + 1], lightposn[m * 3 + 2]);
				vec3 dir;
				if (is_point[m]) dir = glm::normalize(loc - point);
				else dir = glm::normalize(loc);
				vec3 _normal, _point;
				int hits = 0;
				int _index = 0;
				float dist = check_inters(dir, point + EPSILON * dir, _index, hits, _normal, _point);
				bool is_shadow = false;
				if (hits != 0 && dist > 0) {
					if (is_point[m] && distance(point, loc) > dist) {
						is_shadow = true;
					}
					else if (!is_point[m] && distance(point, loc) > 10000000) {
						is_shadow = true;
					}
				}

				vec3 color = vec3(lightcolor[m * 3], lightcolor[m * 3 + 1], lightcolor[m * 3 + 2]);
				float d = 0;
				float atten;
				if (is_point[m]) {
					d = glm::distance(point, loc);
					atten = 1.0f / (attenuation[0] + attenuation[1] * d + attenuation[2] * d * d);
				}
				else {
					atten = 1.0f;
				}
				//float atten = 1.0f;
				float t = fmax(glm::dot(dir, normal), 0.0);
				vec3 factor1 = spheres[index].diffuse * t;

				vec3 half = glm::normalize(dir + glm::normalize(init_point - point));
				float y = fmax(glm::dot(half, normal), 0.0f);
				vec3 factor2 = spheres[index].specular * pow(y, spheres[index].shininess);

				// recursion
				dir = glm::normalize(point - init_point);
				vec3 ref_ray = glm::normalize(dir - 2 * glm::dot(dir, normal) * normal);
				vec3 recursion_lighting(0, 0, 0);
				if (depth_counter < maxdepth) {
					++depth_counter;
					float recur_dist = check_inters(ref_ray, point + EPSILON * ref_ray, _index, hits, _normal, _point);
					if (hits != 0 && recur_dist > 0) {

						recursion_lighting = findColor(_index, _point, hits, _normal, point);
					}
				}
				vec3 factor3 = spheres[index].specular * recursion_lighting;
				//if (!is_point[m]) factor3 = vec3(0, 0, 0);	// ????
				if (is_shadow) color = vec3(0, 0, 0);

				if (color_light) lighting += color / atten * (factor1 + factor2) + factor3;
				else return color / atten * (factor1 + factor2) + factor3 + lighting;
			}
		}
	}

	return lighting;
}

int main() {
	readfile("scene7.test");
	//cout << triangles.size() << endl;
	unsigned char* pixels = (unsigned char*)malloc(static_cast<size_t>(3 * width * height));
	memset(pixels, 0, 3 * static_cast<size_t>(width* height));

	//cout << glm::to_string(triangles[0].transform);
	//BYTE pixels[3 * width * height] = { 0 };
	FreeImage_Initialise();
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			if (i % 10 == 0 && j == 0) cout << i << endl;
			vec3 ray = rayCast(j, i);
			//if (i == 0) std::cout << ray.x << ' ' << ray.y << ' ' << ray.z << endl;
			int index, object;
			vec3 normal, point;

			check_inters(ray, camera_eye, index, object, normal, point);
			depth_counter = 0;
			vec3 lighting = findColor(index, point, object, normal);

			unsigned int loc = (i * width + j) * 3;
			pixels[loc] = (int)(lighting[2] * 255);
			pixels[loc + 1] = (int)(lighting[1] * 255);
			pixels[loc + 2] = (int)(lighting[0] * 255);
		}
	}
	FIBITMAP* img = FreeImage_ConvertFromRawBits(pixels, width, height, width * 3, 24, 0, 0, 0, true);
	FreeImage_Save(FIF_PNG, img, output_file, 0);
	FreeImage_DeInitialise();
	
	std::cout << output_file;
}