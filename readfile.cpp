#include "readfile.h"
#include "Transform.h"

using namespace std;
using namespace glm;

void rightmultiply(const mat4& M, stack<mat4>& transfstack)
{
	mat4& T = transfstack.top();
	T = M * T;
	transfstack.push(T);
}


void readfile(const string& filename) {
	stack<mat4> transfstack;
	transfstack.push(mat4(1.0));
	ifstream file{ filename };
	string line_in;
	ambient = vec3(0.2, 0.2, 0.2);
	maxdepth = 5;
	attenuation = vec3(1, 0, 0);
	output_file = (char*)"1.png";
	while (getline(file, line_in)) {
		stringstream line(line_in);
		string word;
		line >> word;
		if (word[0] == '#' || word[0] == ' ') continue;
		else if (word == "size") {
			int x, y;
			line >> x >> y;
			//width = 640;
			//height = 480;
			width = 1920;
			height = 1080;
			//width = x;
			//height = y;
		}
		else if (word == "output") {
			line >> word;
			output_file = _strdup(word.c_str());
		}
		else if (word == "maxdepth") {
			int d;
			line >> d;
			maxdepth = d;
		}
		else if (word == "camera") {
			float x, y, z;
			line >> x >> y >> z;
			camera_eye = vec3(x, y, z);
			line >> x >> y >> z;
			center = vec3(x, y, z);
			line >> x >> y >> z;
			camera_up = vec3(x, y, z);
			//camera_up = Transform::upvector(camera_up, camera_eye);
			line >> x;
			camera_fov = x;
		}
		else if (word == "tri") {
			float a, b, c;
			line >> a >> b >> c;
			triangle tri = triangle();
			vector<vec3> verts{ vertices[a], vertices[b], vertices[c] };
			tri.vertices = verts;
			tri.ambient = ambient;
			tri.diffuse = diffuse;
			tri.emission = emission;
			tri.shininess = shininess;
			tri.specular = specular;
			tri.transform = transfstack.top();
			tri.inv_transform = glm::inverse(transfstack.top());
			triangles.push_back(tri);
			//cout << transfstack.size();
		}
		else if (word == "sphere") {
			float x, y, z, r;
			line >> x >> y >> z >> r;
			sphere sph = sphere();
			sph.center = vec3(x, y, z);;
			sph.radius = r;
			sph.ambient = ambient;
			sph.diffuse = diffuse;
			sph.emission = emission;
			sph.shininess = shininess;
			sph.specular = specular;
			sph.transform = transfstack.top();
			sph.inv_transform = glm::inverse(transfstack.top());
			spheres.push_back(sph);
		}
		else if (word == "directional") {
			float a, b, c, d, e, f;
			line >> a >> b >> c >> d >> e >> f;

			lightposn[3 * numused] = a;
			lightposn[3 * numused + 1] = b;
			lightposn[3 * numused + 2] = c;

			lightcolor[3 * numused] = d;
			lightcolor[3 * numused + 1] = e;
			lightcolor[3 * numused + 2] = f;

			is_point[numused] = false;
			++numused;
		}
		else if (word == "point") {
			float x, y, z, r, g, b;
			line >> x >> y >> z >> r >> g >> b;

			lightposn[3 * numused] = x;
			lightposn[3 * numused + 1] = y;
			lightposn[3 * numused + 2] = z;

			lightcolor[3 * numused] = r;
			lightcolor[3 * numused + 1] = g;
			lightcolor[3 * numused + 2] = b;

			is_point[numused] = true;
			++numused;
		}
		else if (word == "attenuation") {
			float x, y, z;
			line >> x >> y >> z;
			attenuation = vec3(x, y, z);
		}
		else if (word == "maxverts") {
			line >> maxverts;
		}
		else if (word == "translate") {
			float x, y, z;
			line >> x >> y >> z;
			transfstack.top() = glm::translate(transfstack.top(), vec3(x, y, z));

			//transfstack.top() *= Transform::translate(x, y, z);
		}
		else if (word == "scale") {
			float x, y, z;
			line >> x >> y >> z;
			transfstack.top() = glm::scale(transfstack.top(), vec3(x, y, z));
			//transfstack.top() *= Transform::scale(x, y, z);
		}
		else if (word == "rotate") {
			float x, y, z, angle;
			line >> x >> y >> z >> angle;
			transfstack.top() *= glm::mat4(Transform::rotate(angle, vec3(x, y, z)));
			
			//vec3 axis(x, y, z);
			//normalize(axis);
			//transfstack.top() *= mat4(Transform::rotate(angle, axis));
			//mat3 rotmat = Transform::rotate(angle, axis);
			//mat4 rotate = mat4(rotmat);
			//rightmultiply(rotate, transfstack);
		}
		else if (word == "ambient") {
			float x, y, z;
			line >> x >> y >> z;
			ambient = vec3(x, y, z);
		}
		else if (word == "emission") {
			float x, y, z;
			line >> x >> y >> z;
			emission = vec3(x, y, z);
		}
		else if (word == "diffuse") {
			float x, y, z;
			line >> x >> y >> z;
			diffuse = vec3(x, y, z);
		}
		else if (word == "shininess") {
			float x;
			line >> x;
			shininess = x;
		}
		else if (word == "specular") {
			float x, y, z;
			line >> x >> y >> z;
			specular = vec3(x, y, z);
		}
		else if (word == "popTransform") {
			if (transfstack.size() <= 1) {
				cerr << "Stack has no elements.  Cannot Pop\n";
			}
			else {
				transfstack.pop();
			}
		}
		else if (word == "pushTransform") {
			transfstack.push(transfstack.top());
		}
		else if (word == "vertex") {
			float x, y, z;
			line >> x >> y >> z;
			vertices.push_back(vec3(x, y, z));
		}
		else{
		cout << word << endl;
		}
	}
}