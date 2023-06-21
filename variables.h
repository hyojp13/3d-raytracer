#pragma once

#include <vector>
#include "glm.0.9.7.1/build/native/include/glm/ext.hpp"
#include "glm.0.9.7.1/build/native/include/glm/glm.hpp"
#include "glm.0.9.7.1/build/native/include/glm/gtc/matrix_transform.hpp"

using namespace glm;

#ifdef MAINPROGRAM 
#define EXTERN 
#else 
#define EXTERN extern 
#endif

EXTERN unsigned int height, width;

EXTERN int maxdepth;
EXTERN char* output_file;
EXTERN int maxverts;
EXTERN vec3 camera_eye;
EXTERN vec3 center;
EXTERN vec3 camera_up;
EXTERN float camera_fov;
EXTERN vec3 ambient;
EXTERN vec3 diffuse;
EXTERN float shininess;
EXTERN vec3 specular;
EXTERN vec3 emission;
EXTERN vec3 attenuation;
EXTERN int numused;

/*EXTERN struct light {
	vec3 dir;
	vec3 color;
	mat4 transform;
	bool is_point;
};*/ //may need transformations, not sure if neccessary

EXTERN float lightposn[3 * 10];
EXTERN float lightcolor[3 * 10];
EXTERN bool is_point[10];
//EXTERN float lightransf[3 * 10];

EXTERN std::vector<vec3> vertices;

EXTERN struct triangle {
	std::vector<vec3> vertices;
	vec3 ambient;
	vec3 diffuse;
	vec3 specular;
	vec3 emission;
	float shininess;
	mat4 transform;
	mat4 inv_transform;

};

EXTERN struct sphere {
	vec3 center;
	float radius;
	vec3 ambient;
	vec3 diffuse;
	vec3 specular;
	vec3 emission;
	float shininess;
	mat4 transform;
	mat4 inv_transform;
};

EXTERN std::vector<triangle> triangles;
EXTERN std::vector<sphere> spheres;
