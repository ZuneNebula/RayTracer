// The JSON library allows you to reference JSON arrays like C++ vectors and JSON objects like C++ maps.

#include "raytracer.h"

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "json.hpp"

using json = nlohmann::json;

const char *PATH = "scenes/";

#define M_PI 3.14159265358979323846264338327950288
double fov = 60;
double epsilon = 0.0001;
colour3 background_colour(0, 0, 0);

json scene;

/****************************************************************************/
glm::vec3 vector_to_vec3(const std::vector<float> &v) {
	return glm::vec3(v[0], v[1], v[2]);
}

class Object {
	char* type;
public:
	Object(char* type)
	{
		this->type = type;
	}
};

class Material {
public:
	colour3 ambient;
	colour3 diffuse;
	colour3 specular;
	float shininess;
	colour3 reflective;
	colour3 transmissive;
	float refraction;

	Material()
	{
		this->ambient = { 0,0,0 };
		this->diffuse = { 0,0,0 };
		this->specular = { 0,0,0 };
		this->reflective = { 0,0,0 };
		this->transmissive = { 0,0,0 };
		this->shininess = 0;
		this->refraction = 0;
	}
	Material(json M)
	{
		if (M.find("ambient") != M.end())
		{
			this->ambient = vector_to_vec3(M["ambient"]);
		}
		if (M.find("diffuse") != M.end())
		{
			this->diffuse = vector_to_vec3(M["diffuse"]);
		}
		if (M.find("specular") != M.end())
		{
			this->specular = vector_to_vec3(M["specular"]);
		}
		if (M.find("reflective") != M.end())
		{
			this->reflective = vector_to_vec3(M["reflective"]);
		}
		if (M.find("transmissive") != M.end())
		{
			this->transmissive = vector_to_vec3(M["transmissive"]);
		}
		if (M.find("shininess") != M.end())
		{
			this->shininess = M["shininess"];
		}
		if (M.find("refraction") != M.end())
		{
			this->refraction = M["refraction"];
		}

	}
};

class Sphere :Object{
public:
	float r;
	point3 pos;
	Material m;
	Sphere(float r, point3 pos, Material &m) : Object("Sphere")
	{
		this->r = r;
		this->pos = pos;
		this->m = m;
	}
};

class Plane :Object{
public:
	point3 pos;
	point3 n;
	Material m;
	Plane(point3 pos, point3 n, Material &m) : Object("Plane")
	{
		this->pos = pos;
		this->n = n;
		this->m = m;
	}
};

class Mesh :Object{
public:
	std::vector<std::vector<std::vector<float>>> triangles;
	Material m;
	Mesh(std::vector<std::vector<std::vector<float>>> triangles, Material &m) : Object("Mesh")
	{
		this->triangles = triangles;
		this->m = m;
	}
};

struct Objects {
	std::vector <Sphere> spheres;
	std::vector <Plane> planes;
	std::vector <Mesh> meshes;
};
Objects OBJs;

class Light {
	char *type;
public:
	Light(char *type)
	{
		this->type = type;
	}
};

class Ambient :Light {
public:
	colour3 Ia;
	Ambient(colour3 Ia) : Light("ambient")
	{
		this->Ia = Ia;
	}
};

class Directional :Light {
public:
	colour3 Id;
	point3 Dir;
	Directional(colour3 Ia, point3 Dir) : Light("ambient")
	{
		this->Id = Ia;
		this->Dir = Dir;
	}
};

class Point :Light {
public:
	colour3 Ip;
	point3 pos;
	Point(colour3 Ip, point3 pos) : Light("ambient")
	{
		this->Ip = Ip;
		this->pos = pos;
	}
};

class Spot :Light {
public:
	colour3 Is;
	point3 pos;
	point3 Dir;
	float cutoff;
	Spot(colour3 Is, point3 pos, point3 Dir, float cutoff) : Light("ambient")
	{
		this->Is = Is;
		this->pos = pos;
		this->Dir = Dir;
		this->cutoff = cutoff;
	}
};

struct Lights {
	std::vector <Ambient> ambience;
	std::vector <Directional> directionals;
	std::vector <Point> points;
	std::vector <Spot> spotlights;
	
};

Lights LIT;

// here are some potentially useful utility functions

json find(json &j, const std::string key, const std::string value) {
	json::iterator it;
	for (it = j.begin(); it != j.end(); ++it) {
		if (it->find(key) != it->end()) {
			if ((*it)[key] == value) {
				return *it;
			}
		}
	}
	return json();
}

colour3 Ambi(colour3 Ia, point3 Ka) //returns ambient component
{
	colour3 ambient = Ia * Ka;
	ambient = glm::clamp(ambient, { 0,0,0 }, { 1,1,1 });
	return ambient;
}
colour3 Diff(colour3 Id, point3 Kd, point3 N, point3 L) //returns diffuse component
{
	
	colour3 diffuse = Id * Kd*glm::dot(N, L);
	diffuse = glm::clamp(diffuse, { 0,0,0 }, { 1,1,1 });
	return diffuse;
}
colour3 Spec(colour3 Is, point3 Ks, point3 R, point3 V, float alpha) //returns specular component
{
	
	colour3 specular = Is * Ks*pow(glm::dot(R, V), alpha);
	specular = glm::clamp(specular, { 0,0,0 }, { 1,1,1 });
	return specular;
}

bool isCloser(point3 A, point3 B)
{
	float dist1 = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
	float dist2 = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);

	if (dist1 < dist2)
	{
		return true;
	}
	else
		return false;
}

void printVector(point3 A)
{
	std::cout << A[0] << "," << A[1] << "," << A[2] << std::endl;
}

bool hitSphere(point3 d, point3 e, point3 c, float r, double &Tnear, double&Tfar, point3 &hit) // sphere intersection func
{
	double disc = pow(glm::dot(d, e - c), 2) - glm::dot(d, d)*(glm::dot(e - c, e - c) - pow(r, 2));

	if (disc < 0)
	{
		return false;
	}

	double rest = glm::dot(-d, e - c) / glm::dot(d, d);
	double sqrtdisc = sqrt(disc);

	if (disc > 0)
	{
		double t1 = rest + sqrtdisc / glm::dot(d, d);
		double t2 = rest - sqrtdisc / glm::dot(d, d);

		if (t1 < t2)
		{
			Tnear = t1;
			Tfar = t2;
		}
		else
		{
			Tnear = t2;
			Tfar = t1;
		}
			
	}
	else {
		Tnear = rest;
		Tfar = rest;
	}

	hit = e + float(Tnear)*d;
	return true;

}

bool hitPlane(point3 d, point3 e, point3 n, point3 a, double &Tnear, double &Tfar, point3 &hitpos) // plane intersection func
{
	double denom = glm::dot(n, d);

	if (denom >= 0)
	{
		return false;
	}
	Tnear = glm::dot(n, a - e) / denom;
	Tfar = Tnear + epsilon;
	hitpos = e + float(Tnear)*d;
	return true;

}

bool hitTriangle(point3 d, point3 e, point3 a, point3 b, point3 c, point3 &n, double &Tnear, double &Tfar, point3 &hitpos) //triangle intersection func
{
	n = glm::normalize(glm::cross(b - a, c - a));
	if (hitPlane(d, e, n, a, Tnear, Tfar, hitpos))
	{
		double A = glm::dot(glm::cross(b - a, hitpos - a), n);
		double B = glm::dot(glm::cross(c - b, hitpos - b), n);
		double C = glm::dot(glm::cross(a - c, hitpos - c), n);

		if (A > 0 && B > 0 && C > 0 || A < 0 && B < 0 && C < 0)
		{
			return true;
		}
		else
			return false;
	}
	else
		return false;

}
/****************************************************************************/
void choose_scene(char const *fn) {
	if (fn == NULL) {
		std::cout << "Using default input file " << PATH << "c.json\n";
		fn = "f";
	}

	std::cout << "Loading scene " << fn << std::endl;

	std::string fname = PATH + std::string(fn) + ".json";
	std::fstream in(fname);
	if (!in.is_open()) {
		std::cout << "Unable to open scene file " << fname << std::endl;
		exit(EXIT_FAILURE);
	}

	in >> scene;

	json camera = scene["camera"];
	// these are optional parameters (otherwise they default to the values initialized earlier)
	if (camera.find("field") != camera.end()) {
		fov = camera["field"];
		std::cout << "Setting fov to " << fov << " degrees.\n";
	}
	if (camera.find("background") != camera.end()) {
		background_colour = vector_to_vec3(camera["background"]);
		std::cout << "Setting background colour to " << glm::to_string(background_colour) << std::endl;
	}
	json &objects = scene["objects"];

	for (json::iterator it = objects.begin(); it != objects.end(); ++it) { //loads from json objects to struct OBJs
		json &object = *it;
		if (object["type"] == "sphere") {
			Material m = Material(object["material"]);
			OBJs.spheres.push_back(Sphere(object["radius"], vector_to_vec3(object["position"]), m));
		}
		else if (object["type"] == "plane") {
			Material m = Material(object["material"]);
			OBJs.planes.push_back(Plane(vector_to_vec3(object["position"]),vector_to_vec3(object["normal"]), m));
		}
		else if (object["type"] == "mesh") {
			Material m = Material(object["material"]);
			OBJs.meshes.push_back(Mesh(object["triangles"], m));
		}
	}

	json &lights = scene["lights"];


	for (json::iterator jt = lights.begin(); jt != lights.end(); ++jt) { //loads from json objects to struct OBJs
		json &light = *jt;
		if (light["type"] == "ambient") {
			LIT.ambience.push_back(Ambient(vector_to_vec3(light["color"])));
		}
		else if (light["type"] == "directional") {
			LIT.directionals.push_back(Directional(vector_to_vec3(light["color"]),vector_to_vec3(light["direction"])));
		}
		else if (light["type"] == "point") {
			LIT.points.push_back(Point(vector_to_vec3(light["color"]), vector_to_vec3(light["position"])));
		}
		else if (light["type"] == "spot") {
			LIT.spotlights.push_back(Spot(vector_to_vec3(light["color"]), vector_to_vec3(light["position"]),vector_to_vec3(light["direction"]),light["cutoff"]));
		}
	}

}

bool shadowTrace(const point3 &e, const point3 &d, double T) {
	// NOTE 1: This is a demo, not ray tracing code! You will need to replace all of this with your own code...
  // NOTE 2: You can work with JSON objects directly (like this sample code), read the JSON objects into your own data structures once
	//and render from those (probably in choose_scene), or hard-code the objects in your own data structures and choose them by name in choose_scene;
	//e.g. choose_scene('e') would pick the same scene as the one in "e.json". Your choice.
  // If you want to use this JSON library, https://github.com/nlohmann/json for more information.
	//The code below gives examples of everything you should need: getting named values, iterating over arrays, and converting types.

	// traverse the objects
	double Tval = 1000;//stores closest hit
	double Tfar = 1000;
	double collide = Tval;
	point3 HITPOS;
	Material Mat;
	point3 Normal;

	// every object in the scene will have a "type"
	for (unsigned int i = 0; i < OBJs.spheres.size(); i++) {
		// This is NOT ray-sphere intersection
		// Every sphere will have a position and a radius

		point3 c = OBJs.spheres[i].pos;
		float r = OBJs.spheres[i].r;
		double t1;
		double tfar;
		point3 hitpos;
		bool didhit;

		didhit = hitSphere(d, e, c, r, t1, tfar, hitpos);

		if (didhit) {
			// Every object will have a material
			Material material = OBJs.spheres[i].m;
			point3 N = glm::normalize(hitpos - c);
			if (t1 < Tval)
			{
				Tval = t1;
				HITPOS = hitpos;
				Mat = material;
				Normal = N;
				Tfar = tfar;
			}

			// This is NOT correct: it finds the first hit, not the closest
		}
	}
	for (unsigned int i = 0; i < OBJs.planes.size(); i++) {

		point3 a = OBJs.planes[i].pos;
		point3 n = OBJs.planes[i].n;
		double t2;
		double tfar;
		point3 hitpos;
		bool didhit;

		didhit = hitPlane(d, e, n, a, t2, tfar, hitpos);

		if (didhit) {
			// Every object will have a material
			Material material = OBJs.planes[i].m;
			point3 N = glm::normalize(n);
			if (t2 < Tval)
			{
				Tval = t2;
				HITPOS = hitpos;
				Mat = material;
				Normal = N;
				Tfar = tfar;
			}

			// This is NOT correct: it finds the first hit, not the closest
		}

	}
	for (unsigned int i = 0; i < OBJs.meshes.size(); i++) {
		double t3;
		double tfar;

		for (unsigned int j = 0; j < OBJs.meshes[i].triangles.size(); j++) //triangles
		{
			std::vector<std::vector<float>> triangle = OBJs.meshes[i].triangles[j];

			std::vector<float> A = triangle[0];
			std::vector<float> B = triangle[1];
			std::vector<float> C = triangle[2];

			point3 a = vector_to_vec3(A);
			point3 b = vector_to_vec3(B);
			point3 c = vector_to_vec3(C);
			point3 n;
			point3 hitpos;

			bool didhit;

			didhit = hitTriangle(d, e, a, b, c, n, t3, tfar, hitpos);
			if (didhit) {
				// Every object will have a material
				Material material = OBJs.meshes[i].m;
				point3 N = glm::normalize(n);
				if (t3 < Tval)
				{
					Tval = t3;
					HITPOS = hitpos;
					Mat = material;
					Normal = N;
					Tfar = tfar;
				}

				// This is NOT correct: it finds the first hit, not the closest
			}
		}
	}
	if (Tval != collide && Tval > epsilon && Tval < T)
	{
		//colour = lighting(Tfar, e, s - e, HITPOS, Mat, Normal, refDepth);
		return true;
	}
	return false;
}

colour3 lighting(double t, point3 e, point3 d, point3 hitpos, Material material, point3 N, int refDepth) //lighting function
{
	//traverse the lights
	bool doesReflect = false;
	bool doesTransmit = false;
	float alpha =material.shininess;

	if (material.reflective != point3{0, 0, 0})
	{
		doesReflect = true;
	}
	if (material.transmissive != point3{ 0,0,0 })
	{
		doesTransmit = true;
	}

	point3 Ka,Kd,Ks,Kr,Kt;
	Ka = material.ambient;
	Kd = material.diffuse;
	Ks = material.specular;
	Kr = material.reflective;
	Kt = material.transmissive;

	colour3 Kt2 = point3(1 - Kt[0], 1 - Kt[1], 1 - Kt[2]);
	colour3 Kr2 = point3(1 - Kr[0], 1 - Kr[1], 1 - Kr[2]);

	colour3 color1, color2, color3, color4,color5;
	colour3 color;


	for(int i=0;i<LIT.ambience.size();i++)
	{
		colour3 Ia = LIT.ambience[i].Ia;
		color1 += Ambi(Ia, Ka);
	}
	for (int i = 0; i < LIT.directionals.size(); i++) //diffuse and specular
	{
		colour3 Id = LIT.directionals[i].Id;
		point3 Dir = LIT.directionals[i].Dir;
		point3 L = glm::normalize(-Dir);
		double T = 100;
			
		if (!shadowTrace(hitpos, L, T))
		{
			color2 = color2 + Diff(Id, Kd, N, L);
		}
			
		point3 V = glm::normalize(e-hitpos);
		point3 R = glm::normalize(2 * glm::dot(N, L)*N - L);
			
		if (!shadowTrace(hitpos, L, T))
		{
			color3 = color3 + Spec(Id, Ks, R, V, alpha);
		}
			
	}
	for (int i = 0; i < LIT.points.size(); i++)
	{
		colour3 Id = LIT.points[i].Ip;
		point3 Pos = LIT.points[i].pos;
		point3 Ray = Pos - hitpos;
		//double T = sqrt(pow(Pos[0] - hitpos[0], 2) + pow(Pos[1] - hitpos[1], 2) + pow(Pos[2] - hitpos[2], 2));
		//float r = sqrt(Ray[0] * Ray[0] + Ray[1] * Ray[1] + Ray[2] * Ray[2]);
		//Id = Id / (r*r);
		point3 L = glm::normalize(Pos - hitpos);
		double T = glm::length(Pos - hitpos);

		if (!shadowTrace(hitpos, L, T))
		{
			color2 = color2 + Diff(Id, Kd, N, L);
		}

		point3 V = glm::normalize(e-hitpos);
		point3 R = glm::normalize(2 * glm::dot(N, L)*N - L);
			
		if (!shadowTrace(hitpos, L, T))
		{
			color3 = color3 + Spec(Id, Ks, R, V, alpha);
		}
			
	}
	for (int i = 0; i < LIT.spotlights.size(); i++)
	{
		float cutoff = LIT.spotlights[i].cutoff;

		colour3 Id = LIT.spotlights[i].Is;
		point3 Pos = LIT.spotlights[i].pos;
		point3 Ray = glm::normalize(hitpos-Pos);
		point3 Dir = glm::normalize(LIT.spotlights[i].Dir);
		float angle = acos(glm::dot(Dir, Ray)) * 180 / M_PI;
		//double T = sqrt(pow(Pos[0] - hitpos[0], 2) + pow(Pos[1] - hitpos[1], 2) + pow(Pos[2] - hitpos[2], 2));
		double T = glm::length(Pos - hitpos);
		if (angle < cutoff)
		{
			point3 L = glm::normalize(Pos - hitpos);
				
			if (!shadowTrace(hitpos, L, T))
			{
				color2 = color2 + Diff(Id, Kd, N, L);
			}
				
			point3 V = glm::normalize(e-hitpos);
			point3 R = glm::normalize(2 * glm::dot(N, L)*N - L);
				
			if (!shadowTrace(hitpos, L, T))
			{
				color3 = color3 + Spec(Id, Ks, R, V, alpha);
			}
		}

	}


	if (doesReflect && refDepth > 1)
	{
		refDepth--;
		point3 V = glm::normalize(e - hitpos);
		point3 R = glm::normalize(2 * glm::dot(N, V)*N - V);
		trace(hitpos + point3(epsilon, epsilon, epsilon)*R, hitpos + point3(epsilon, epsilon, epsilon)*R + R, color4, refDepth, false);
		color4 = glm::clamp(Kr * color4, { 0,0,0 }, { 1,1,1 });
		color1 = Kr2 * color1;
		color2 = Kr2 * color2;
		color3 = Kr2 * color3;
	}

	if (!doesTransmit)
	{
		color = glm::clamp(color1 + color2 + color3 + color4, { 0,0,0 }, { 1,1,1 }); //sum of ambient/diffuse/specular/reflective
		return color;
	}
	if(doesTransmit && refDepth >1)
	{
		refDepth--;
		point3 Ray = glm::normalize(hitpos - e);
		point3 F = e + float(t)*Ray;
		//std::cout << t << std::endl;
		trace(F, F + Ray, color5, refDepth, false);
		color5 = glm::clamp(Kt * color5, { 0,0,0 }, { 1,1,1 });
		color = glm::clamp(Kt2*(color1 + color2 + color3 + color4) + color5, { 0,0,0 }, { 1,1,1 }); //transparency
		return color;
	}
	else
	{
		color5 = glm::clamp(Kt * color5, { 0,0,0 }, { 1,1,1 });
		color = glm::clamp(Kt2*(color1 + color2 + color3 + color4) + color5, { 0,0,0 }, { 1,1,1 }); //transparency
		return color;
	}
	return color;
}
bool trace(const point3 &e, const point3 &s, colour3 &colour, int refDepth, bool pick) {
	// NOTE 1: This is a demo, not ray tracing code! You will need to replace all of this with your own code...
  // NOTE 2: You can work with JSON objects directly (like this sample code), read the JSON objects into your own data structures once
	//and render from those (probably in choose_scene), or hard-code the objects in your own data structures and choose them by name in choose_scene;
	//e.g. choose_scene('e') would pick the same scene as the one in "e.json". Your choice.
  // If you want to use this JSON library, https://github.com/nlohmann/json for more information.
	//The code below gives examples of everything you should need: getting named values, iterating over arrays, and converting types.

	// traverse the objects
	double Tval = 1000;//stores closest hit
	double Tfar = 1000;
	double collide = Tval;
	point3 HITPOS;
	Material Mat;
	point3 Normal;

	// every object in the scene will have a "type"
	for (unsigned int i = 0; i < OBJs.spheres.size(); i++) {
		// This is NOT ray-sphere intersection
		// Every sphere will have a position and a radius

		point3 c = OBJs.spheres[i].pos;
		point3 d = s - e;
		float r = OBJs.spheres[i].r;
		double t1;
		double tfar;
		point3 hitpos;
		bool didhit;

		didhit = hitSphere(d, e, c, r, t1, tfar, hitpos);

		if (didhit) {
			// Every object will have a material
			Material material = OBJs.spheres[i].m;
			point3 N = glm::normalize(hitpos - c);
			if (t1 < Tval)
			{
				Tval = t1;
				HITPOS = hitpos;
				Mat = material;
				Normal = N;
				Tfar = tfar;
			}

			// This is NOT correct: it finds the first hit, not the closest
		}
	}
	for (unsigned int i = 0; i < OBJs.planes.size(); i++) {

		point3 a = OBJs.planes[i].pos;
		point3 d = s - e;
		point3 n = glm::normalize(OBJs.planes[i].n);
		double t2;
		double tfar;
		point3 hitpos;
		bool didhit;

		didhit = hitPlane(d, e, n, a, t2, tfar, hitpos);

		if (didhit) {
			// Every object will have a material
			Material material = OBJs.planes[i].m;
			point3 N = glm::normalize(n);
			if (t2 < Tval)
			{
				Tval = t2;
				HITPOS = hitpos;
				Mat = material;
				Normal = N;
				Tfar = tfar;
			}

			// This is NOT correct: it finds the first hit, not the closest
		}

	}
	for (unsigned int i = 0; i < OBJs.meshes.size(); i++) {
		double t3;
		double tfar;

		for (unsigned int j = 0; j < OBJs.meshes[i].triangles.size(); j++) //triangles
		{
			std::vector<std::vector<float>> triangle = OBJs.meshes[i].triangles[j];

			std::vector<float> A = triangle[0];
			std::vector<float> B = triangle[1];
			std::vector<float> C = triangle[2];

			point3 a = vector_to_vec3(A);
			point3 b = vector_to_vec3(B);
			point3 c = vector_to_vec3(C);
			point3 d = s - e;
			point3 n;
			point3 hitpos;

			bool didhit;

			didhit = hitTriangle(d, e, a, b, c, n, t3, tfar, hitpos);
			if (didhit) {
				// Every object will have a material
				Material material = OBJs.meshes[i].m;
				point3 N = glm::normalize(n);
				if (t3 < Tval)
				{
					Tval = t3;
					HITPOS = hitpos;
					Mat = material;
					Normal = N;
					Tfar = tfar;
				}

				// This is NOT correct: it finds the first hit, not the closest
			}
		}
	}
	if (Tval != collide)
	{
		colour = lighting(Tfar, e, s - e, HITPOS, Mat, Normal, refDepth);
		return true;
	}
	return false;
}
/*bool shadowTrace(const point3 &e, const point3 &d, double T) {
	// NOTE 1: This is a demo, not ray tracing code! You will need to replace all of this with your own code...
  // NOTE 2: You can work with JSON objects directly (like this sample code), read the JSON objects into your own data structures once
	//and render from those (probably in choose_scene), or hard-code the objects in your own data structures and choose them by name in choose_scene;
	//e.g. choose_scene('e') would pick the same scene as the one in "e.json". Your choice.
  // If you want to use this JSON library, https://github.com/nlohmann/json for more information.
	//The code below gives examples of everything you should need: getting named values, iterating over arrays, and converting types.

	// traverse the objects
	double Tval = 100;//stores closest hit
	double Tfar = 100;
	double collide = Tval;
	point3 HITPOS;
	json Mat;
	point3 Normal;

	json &objects = scene["objects"];
	for (json::iterator it = objects.begin(); it != objects.end(); ++it) {
		json &object = *it;

		// every object in the scene will have a "type"
		if (object["type"] == "sphere") {
			// This is NOT ray-sphere intersection
			// Every sphere will have a position and a radius
			std::vector<float> pos = object["position"];

			point3 c = vector_to_vec3(pos);
			float r = float(object["radius"]);
			double t1;
			double tfar;
			point3 hitpos;
			bool didhit;

			didhit = hitSphere(d, e, c, r, t1, tfar, hitpos);

			if (didhit) {
				// Every object will have a material
				json &material = object["material"];
				point3 N = glm::normalize(hitpos - c);
				if (t1 < Tval)
				{
					Tval = t1;
					HITPOS = hitpos;
					Mat = material;
					Normal = N;
					Tfar = tfar;
				}

				// This is NOT correct: it finds the first hit, not the closest
			}
		}
		else if (object["type"] == "plane") {
			std::vector<float> pos = object["position"];
			std::vector<float> normal = object["normal"];

			point3 a = vector_to_vec3(pos);
			point3 n = vector_to_vec3(normal);
			double t2;
			double tfar;
			point3 hitpos;
			bool didhit;

			didhit = hitPlane(d, e, n, a, t2, tfar, hitpos);

			if (didhit) {
				// Every object will have a material
				json &material = object["material"];
				point3 N = glm::normalize(n);
				if (t2 < Tval)
				{
					Tval = t2;
					HITPOS = hitpos;
					Mat = material;
					Normal = N;
					Tfar = tfar;
				}

				// This is NOT correct: it finds the first hit, not the closest
			}

		}
		else if (object["type"] == "mesh")
		{
			double t3;
			double tfar;
			json &triangles = object["triangles"];
			for (json::iterator jt = triangles.begin(); jt != triangles.end(); ++jt) //triangles
			{
				json& triangle = *jt;

				std::vector<float> A = triangle[0];
				std::vector<float> B = triangle[1];
				std::vector<float> C = triangle[2];

				point3 a = vector_to_vec3(A);
				point3 b = vector_to_vec3(B);
				point3 c = vector_to_vec3(C);
				point3 n;
				point3 hitpos;

				bool didhit;

				didhit = hitTriangle(d, e, a, b, c, n, t3, tfar, hitpos);
				if (didhit) {
					// Every object will have a material
					json &material = object["material"];
					point3 N = glm::normalize(n);
					if (t3 < Tval)
					{
						Tval = t3;
						HITPOS = hitpos;
						Mat = material;
						Normal = N;
						Tfar = tfar;
					}

					// This is NOT correct: it finds the first hit, not the closest
				}
			}
		}
	}
	if (Tval != collide && Tval > epsilon && Tval < T)
	{
		//colour = lighting(Tval, e, HITPOS, Mat, Normal);
		return true;
	}
	return false;
}*/
/*colour3 lighting(double t, point3 e, point3 d, point3 hitpos, json material, point3 N, int refDepth) //lighting function
{
	//traverse the lights
	json &lights = scene["lights"];
	std::vector<float> ambient = material["ambient"];
	std::vector<float> diffuse = { 0,0,0 };
	std::vector<float> specular = { 0,0,0 };
	std::vector<float> reflective = { 0,0,0 };
	std::vector<float> transmissive = { 0,0,0 };
	bool doesReflect = false;
	bool doesTransmit = false;
	float alpha = 0;

	if (material.find("diffuse") != material.end())
	{
		std::vector<float> temp = material["diffuse"];
		diffuse = temp;
	}
	if (material.find("specular") != material.end())
	{
		std::vector<float> temp = material["specular"];
		specular = temp;
		alpha = material["shininess"];
	}
	if (material.find("reflective") != material.end())
	{
		std::vector<float> temp = material["reflective"];
		reflective = temp;
		doesReflect = true;
	}
	if (material.find("transmissive") != material.end())
	{
		std::vector<float> temp = material["transmissive"];
		transmissive = temp;
		doesTransmit = true;
	}

	point3 Ka, Kd, Ks, Kr, Kt;
	Ka = vector_to_vec3(ambient);
	Kd = vector_to_vec3(diffuse);
	Ks = vector_to_vec3(specular);
	Kr = vector_to_vec3(reflective);
	Kt = vector_to_vec3(transmissive);
	colour3 Kt2 = point3(1 - Kt[0], 1 - Kt[1], 1 - Kt[2]);
	colour3 Kr2 = point3(1 - Kr[0], 1 - Kr[1], 1 - Kr[2]);

	colour3 color1, color2, color3, color4, color5;
	colour3 color;

	for (json::iterator it = lights.begin(); it != lights.end(); ++it)
	{
		json &light = *it;

		if (light["type"] == "ambient")
		{
			std::vector<float> color = light["color"];
			colour3 Ia = vector_to_vec3(color);
			color1 = Ambi(Ia, Ka);

		}
		else if (light["type"] == "directional") //diffuse and specular
		{
			std::vector <float> color = light["color"];
			std::vector <float> direction = light["direction"];
			colour3 Id = vector_to_vec3(color);
			point3 Dir = vector_to_vec3(direction);
			point3 L = glm::normalize(-Dir);
			double T = 100;

			if (!shadowTrace(hitpos, L, T))
			{
				color2 = color2 + Diff(Id, Kd, N, L);
			}

			point3 V = glm::normalize(e - hitpos);
			point3 R = glm::normalize(2 * glm::dot(N, L)*N - L);

			if (!shadowTrace(hitpos, L, T))
			{
				color3 = color3 + Spec(Id, Ks, R, V, alpha);
			}

		}
		else if (light["type"] == "point")
		{
			std::vector <float> color = light["color"];
			std::vector <float> position = light["position"];
			colour3 Id = vector_to_vec3(color);
			point3 Pos = vector_to_vec3(position);
			point3 Ray = Pos - hitpos;
			//double T = sqrt(pow(Pos[0] - hitpos[0], 2) + pow(Pos[1] - hitpos[1], 2) + pow(Pos[2] - hitpos[2], 2));
			//float r = sqrt(Ray[0] * Ray[0] + Ray[1] * Ray[1] + Ray[2] * Ray[2]);
			//Id = Id / (r*r);
			point3 L = glm::normalize(Pos - hitpos);
			double T = glm::length(Pos - hitpos);

			if (!shadowTrace(hitpos, L, T))
			{
				color2 = color2 + Diff(Id, Kd, N, L);
			}

			point3 V = glm::normalize(e - hitpos);
			point3 R = glm::normalize(2 * glm::dot(N, L)*N - L);

			if (!shadowTrace(hitpos, L, T))
			{
				color3 = color3 + Spec(Id, Ks, R, V, alpha);
			}

		}
		else if (light["type"] == "spot")
		{
			std::vector <float> color = light["color"];
			std::vector <float> direction = light["direction"];
			std::vector <float> position = light["position"];
			float cutoff = light["cutoff"];

			colour3 Id = vector_to_vec3(color);
			point3 Pos = vector_to_vec3(position);
			point3 Ray = glm::normalize(hitpos - Pos);
			point3 Dir = glm::normalize(vector_to_vec3(direction));
			float angle = acos(glm::dot(Dir, Ray)) * 180 / M_PI;
			//double T = sqrt(pow(Pos[0] - hitpos[0], 2) + pow(Pos[1] - hitpos[1], 2) + pow(Pos[2] - hitpos[2], 2));
			double T = glm::length(Pos - hitpos);
			if (angle < cutoff)
			{
				point3 L = glm::normalize(Pos - hitpos);

				if (!shadowTrace(hitpos, L, T))
				{
					color2 = color2 + Diff(Id, Kd, N, L);
				}

				point3 V = glm::normalize(e - hitpos);
				point3 R = glm::normalize(2 * glm::dot(N, L)*N - L);

				if (!shadowTrace(hitpos, L, T))
				{
					color3 = color3 + Spec(Id, Ks, R, V, alpha);
				}
			}

		}

	}

	if (doesReflect && refDepth > 1)
	{
		refDepth--;
		point3 V = glm::normalize(e - hitpos);
		point3 R = glm::normalize(2 * glm::dot(N, V)*N - V);
		trace(hitpos + point3(epsilon, epsilon, epsilon)*R, hitpos + point3(epsilon, epsilon, epsilon)*R + R, color4, refDepth, false);
		color4 = glm::clamp(Kr * color4, { 0,0,0 }, { 1,1,1 });
		color1 = Kr2 * color1;
		color2 = Kr2 * color2;
		color3 = Kr2 * color3;
	}

	if (!doesTransmit)
	{
		color = glm::clamp(color1 + color2 + color3 + color4, { 0,0,0 }, { 1,1,1 }); //sum of ambient/diffuse/specular/reflective
		return color;
	}
	else if (doesTransmit && refDepth > 1)
	{
		refDepth--;
		point3 Ray = glm::normalize(hitpos - e);
		point3 F = e + float(t)*Ray;
		//std::cout << t << std::endl;
		trace(F, F + Ray, color5, refDepth, false);
		color5 = glm::clamp(Kt * color5, { 0,0,0 }, { 1,1,1 });
		color = glm::clamp(Kt2*(color1 + color2 + color3 + color4) + color5, { 0,0,0 }, { 1,1,1 }); //transparency
		return color;
	}
	return color;
}*/
/*bool trace(const point3 &e, const point3 &s, colour3 &colour, int refDepth, bool pick) {
	// NOTE 1: This is a demo, not ray tracing code! You will need to replace all of this with your own code...
  // NOTE 2: You can work with JSON objects directly (like this sample code), read the JSON objects into your own data structures once
	//and render from those (probably in choose_scene), or hard-code the objects in your own data structures and choose them by name in choose_scene;
	//e.g. choose_scene('e') would pick the same scene as the one in "e.json". Your choice.
  // If you want to use this JSON library, https://github.com/nlohmann/json for more information.
	//The code below gives examples of everything you should need: getting named values, iterating over arrays, and converting types.
	
	// traverse the objects
	double Tval = 1000;//stores closest hit
	double Tfar = 1000;
	double collide = Tval;
	point3 HITPOS;
	json Mat;
	point3 Normal;

	json &objects = scene["objects"];
	for (json::iterator it = objects.begin(); it != objects.end(); ++it) {
		json &object = *it;

		// every object in the scene will have a "type"
		if (object["type"] == "sphere") {
			// This is NOT ray-sphere intersection
			// Every sphere will have a position and a radius
			std::vector<float> pos = object["position"];

			point3 c = vector_to_vec3(pos);
			point3 d = s - e;
			float r = float(object["radius"]);
			double t1;
			double tfar;
			point3 hitpos;
			bool didhit;

			didhit = hitSphere(d, e, c, r, t1, tfar, hitpos);

			if (didhit) {
				// Every object will have a material
				json &material = object["material"];
				point3 N = glm::normalize(hitpos - c);
				if (t1 < Tval)
				{
					Tval = t1;
					HITPOS = hitpos;
					Mat = material;
					Normal = N;
					Tfar = tfar;
				}
				
				// This is NOT correct: it finds the first hit, not the closest
			}
		}
		else if (object["type"] == "plane") {
			std::vector<float> pos = object["position"];
			std::vector<float> normal = object["normal"];

			point3 a = vector_to_vec3(pos);
			point3 d = s - e;
			point3 n = vector_to_vec3(normal);
			double t2;
			double tfar;
			point3 hitpos;
			bool didhit;

			didhit = hitPlane(d, e, n, a, t2, tfar, hitpos);

			if (didhit) {
				// Every object will have a material
				json &material = object["material"];
				point3 N = glm::normalize(n);
				if (t2 < Tval)
				{
					Tval = t2;
					HITPOS = hitpos;
					Mat = material;
					Normal = N;
					Tfar = tfar;
				}

				// This is NOT correct: it finds the first hit, not the closest
			}

		}
		else if (object["type"] == "mesh")
		{
			double t3;
			double tfar;
			json &triangles = object["triangles"];
			for (json::iterator jt = triangles.begin(); jt != triangles.end(); ++jt) //triangles
			{
				json& triangle = *jt;
				
				std::vector<float> A = triangle[0];
				std::vector<float> B = triangle[1];
				std::vector<float> C = triangle[2];

				point3 a = vector_to_vec3(A);
				point3 b = vector_to_vec3(B);
				point3 c = vector_to_vec3(C);
				point3 d = s - e;
				point3 n;
				point3 hitpos;
				
				bool didhit;

				didhit = hitTriangle(d, e, a, b, c, n, t3, tfar, hitpos);
				if (didhit) {
					// Every object will have a material
					json &material = object["material"];
					point3 N = glm::normalize(n);
					if(t3 < Tval)
					{
						Tval = t3;
						HITPOS = hitpos;
						Mat = material;
						Normal = N;
						Tfar = tfar;
					}

					// This is NOT correct: it finds the first hit, not the closest
				}
			}
		}
	}
	if (Tval != collide)
	{
		colour = lighting(Tfar, e, s-e, HITPOS, Mat, Normal,refDepth);
		return true;
	}
	return false;
}*/