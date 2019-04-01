#include <string.h>
#include "math.h"
#include <vector>
#include <list>
#include <iostream>
#include <cfloat>


#include "stb_image.h"

class Vector {
public:
	Vector(double x = 0, double y = 0, double z = 0) :
		x(x), y(y), z(z) {};
	double x, y, z;

	double norm2() {
		return x * x + y * y + z * z;
	}
	;
	void normalize() {
		double n = sqrt(norm2());
		x = x / n;
		y = y / n;
		z = z / n;
	}
	Vector& operator+=(const Vector& a) {
		x+=a.x;y+=a.y;z+=a.z;
		return *this;
	}
};

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(const Vector a, const Vector& b);
Vector operator*(double a, const Vector& b);
Vector operator*(const Vector& a, double b);
Vector operator/(const Vector& a, double b);
double dot(const Vector& a, const Vector& b);
Vector cross(const Vector& a, const Vector b);
Vector rotate(const Vector& a, const double& theta);

class Ray {
public:
	Ray(const Vector &o = Vector(), const Vector &d = Vector()) : origin(o), direction(d) {};
	Vector origin, direction;
	void reflechi(const Vector P, const Vector n) {
		direction = direction - 2 * dot(direction, n)* n;
		origin = P + 0.01*direction;
	}
	void refracte(const Vector P, const Vector n, const double n_air, const  double n_sphere, const double alea) {
		double k0 = (n_air - n_sphere)/(n_air + n_sphere); k0 = k0*k0;
		double D ;
		double R ;


		double A = dot(direction, n);
		double B = 1 - (n_air / n_sphere)*(n_air / n_sphere)*(1 - A*A);
		double C = 1 - (n_sphere / n_air)*(n_sphere / n_air)*(1 - A*A);

		if (A< 0 && B > 0) {
			D = 1 + dot(n,direction); D =D*D*D*D*D;
			R = k0 + (1-k0)*D;
			if (alea > R ){
				direction = (n_air / n_sphere)*direction - ((n_air / n_sphere)*A + sqrt(B))*n;
				origin = P + 0.01*direction;
			}
			else reflechi(P, n);
		}

		else if (A < 0) {
			reflechi(P, n);
		}
		else if (C > 0) {
			Vector direction_refracte = (n_sphere / n_air)*direction + (-(n_sphere / n_air)*A + sqrt(C))*n;
			D = 1 - dot(n,direction_refracte); D =D*D*D*D*D;
			R = 0.1;//k0 + (1-k0)*D;
			if (alea > R ){
				direction = direction_refracte;
				origin = P + 0.01*direction;
			}
			else reflechi(P, (-1)*n);

		}
		else {

			reflechi(P, (-1)*n);
		}
	}
};

enum MaterialType {
	DIFFUS,SPEC,TRANS,LAMPE
};

struct Diffus{};
struct Spec{};
struct Trans{};
struct Lampe{};

class Material {
public:
	Material() {};
	Material(const Diffus Type, const Vector &Albedo = Vector()) :  type(DIFFUS), albedo(Albedo) {};
	Material(const Spec Type):  type(SPEC){};
	Material(const Trans Type, const double &N_refract = 1):  type(TRANS), n_refract(N_refract){};
	Material(const Lampe Type, const double &Intensite = 1):  type(LAMPE), intensite(Intensite){};

	MaterialType type;
	Vector albedo;
	double n_refract;
	double intensite;
};


class Object {
public:
	Object(const Vector& o, double r, const Material& M) : origin(o), rayon(r), material(M) {};
	Object(const Vector A,const Vector B,const Vector C,const Material M) : A(A),B(B),C(C),material(M){};
	Object(const std::string object_directory, const Material M):   material(M) {};
	virtual bool intersection(const Ray &r, double &tInterSphere, Vector &N, Vector&albedo) = 0;
	Vector origin;
	double rayon;
	Vector A, B, C;
	Vector N;
	Material material;

};

class Sphere : public Object{
public:
	Sphere(const Vector& o, double r, const Material& M) : Object(o,r,M){};

	virtual bool intersection(const Ray &r, double &tInterSphere,  Vector &N, Vector&albedo) {
		double a = 1;
		double b = 2 * dot(r.direction, r.origin - origin);
		double c = (r.origin - origin).norm2() - (rayon)*(rayon);
		double delta = b * b - 4 * a*c;
		if (delta < 0) return false;
		double t1 = (-b - sqrt(delta)) / (2 * a);
		double t2 = (-b + sqrt(delta)) / (2 * a);

		if (t2 < 0) return false;

		if (t1 > 0) {
			tInterSphere = t1;
		}
		if (t1 < 0) {
			tInterSphere = t2;
		}
		N = r.origin + tInterSphere*r.direction - origin;
		//N.normalize();
		N = N/rayon;
		albedo = material.albedo;
		return true;
	}
};


class Triangle : public Object{
public:
	Triangle(const Vector A,const Vector B,const Vector C,const Material m) : Object(A,B,C,m){};
	virtual bool intersection(const Ray &r, double &tInterTriangle, Vector &N, Vector&albedo){

		N = cross(A - B,C - A); N.normalize();

		tInterTriangle = dot(C - r.origin,N)/dot(r.direction,N);
		if (tInterTriangle < 0) return false;
		Vector P = r.origin + tInterTriangle*r.direction;
		Vector AB = B - A; double norm2AB = AB.norm2();
		Vector AC = C - A; double norm2AC = AC.norm2();
		Vector AP = P - A;
		double dotABAC= dot(AB,AC);
		double dotAPAB= dot(AP,AB);
		double dotAPAC= dot(AP,AC);

		double detM = (norm2AB*norm2AC - dotABAC*dotABAC);
		double detBeta = (dotAPAB*norm2AC - dotABAC*dotAPAC);
		double detgamma = (dotAPAC*norm2AB - dotAPAB*dotABAC);

		double beta = detBeta /detM;
		double gamma = detgamma/detM;
		double alpha = 1 - gamma - beta;

		if (beta >= 0 && beta <=1 && gamma >= 0 && gamma <=1 && alpha >= 0 && alpha <=1){
			albedo = material.albedo;
			return true;
		}
		return false;

	}
	bool intersection_triangle(const Ray &r, double &tInterTriangle,double&alpha, double&beta,double&gamma, Vector&albedo ){
		Vector N = cross(A - B,C - A); N.normalize();

		tInterTriangle = dot(C - r.origin,N)/dot(r.direction,N);
		if (tInterTriangle < 0) return false;
		Vector P = r.origin + tInterTriangle*r.direction;
		Vector AB = B - A; double norm2AB = AB.norm2();
		Vector AC = C - A; double norm2AC = AC.norm2();
		Vector AP = P - A;
		double dotABAC= dot(AB,AC);
		double dotAPAB= dot(AP,AB);
		double dotAPAC= dot(AP,AC);

		double detM = (norm2AB*norm2AC - dotABAC*dotABAC);
		double detBeta = (dotAPAB*norm2AC - dotABAC*dotAPAC);
		double detgamma = (dotAPAC*norm2AB - dotAPAB*dotABAC);

		beta = detBeta /detM;
		gamma = detgamma/detM;
		alpha = 1 - gamma - beta;

		if (beta >= 0 && beta <=1 && gamma >= 0 && gamma <=1 && alpha >= 0 && alpha <=1){
			albedo = material.albedo;
			return true;
		}
		return false;
	}
};


class BBox {
public:
	BBox() {};
	BBox(const Vector &Top ,const Vector &Bottom) : top(Top),bottom(Bottom){}
	bool intersection(const Ray &r){
		double t1_x = (bottom.x - r.origin.x)/r.direction.x;
		double t2_x = (top.x - r.origin.x)/r.direction.x;
		if(t1_x > t2_x){double temp = t1_x;t1_x = t2_x;t2_x = temp;}


		double t1_y = (bottom.y - r.origin.y)/r.direction.y;
		double t2_y = (top.y - r.origin.y)/r.direction.y;
		if(t1_y > t2_y){double temp = t1_y;t1_y = t2_y;t2_y = temp;}


		double t1_z = (bottom.z - r.origin.z)/r.direction.z;
		double t2_z = (top.z - r.origin.z)/r.direction.z;
		if(t1_z > t2_z){double temp = t1_z;t1_z = t2_z;t2_z = temp;}

		double t_inter_min = std::max(std::max(t1_x,t1_y),t1_z);
		double t_inter_max = std::min(std::min(t2_x,t2_y),t2_z);

		if(t_inter_max > t_inter_min && t_inter_max > 0 ) return true;

		return false;

	}
	Vector top;
	Vector bottom;
};
class Boxs_tree {
public:
	int i0,i1;
	BBox box;
	Boxs_tree *left_branch ;
	Boxs_tree *right_branch ;
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};


class Geometry : public Object {
private:
	std::vector<stbi_uc*> textures;
	std::vector<int> textures_x;
	std::vector<int> textures_y;

	Boxs_tree boxs_tree;
public:
	~Geometry() {}
	Geometry(const std::string object_directory, const Material &M): Object(object_directory,M) {
		readOBJ((object_directory + "BeautifulGirl.obj").c_str());

		for (unsigned int i = 0; i < vertices.size(); i++) {
			vertices[i] = Vector(vertices[i].x, vertices[i].z, vertices[i].y) * 20 + Vector(0,-10,0);
		}

		for (unsigned int i = 0; i < normals.size(); i++) {
			normals[i] = Vector(normals[i].x, normals[i].z, normals[i].y);
		}


		add_textures(object_directory);
		build_tree(&boxs_tree,0, indices.size());


	}
	BBox build_box(int i0,int i1){
		BBox bbox;
		bbox.top.x = vertices[indices[i0].vtxi].x ;
		bbox.top.y = vertices[indices[i0].vtxi].y ;
		bbox.top.z = vertices[indices[i0].vtxi].z ;


		bbox.bottom =bbox.top;
		for(int i = i0 ; i < i1; i++){
			bbox.top.x = std::max(bbox.top.x,vertices[indices[i].vtxi].x);
			bbox.top.x = std::max(bbox.top.x,vertices[indices[i].vtxj].x);
			bbox.top.x = std::max(bbox.top.x,vertices[indices[i].vtxk].x);

			bbox.top.y = std::max(bbox.top.y,vertices[indices[i].vtxi].y);
			bbox.top.y = std::max(bbox.top.y,vertices[indices[i].vtxj].y);
			bbox.top.y = std::max(bbox.top.y,vertices[indices[i].vtxk].y);

			bbox.top.z = std::max(bbox.top.z,vertices[indices[i].vtxi].z);
			bbox.top.z = std::max(bbox.top.z,vertices[indices[i].vtxj].z);
			bbox.top.z = std::max(bbox.top.z,vertices[indices[i].vtxk].z);


			bbox.bottom.x = std::min(bbox.bottom.x,vertices[indices[i].vtxi].x);
			bbox.bottom.x = std::min(bbox.bottom.x,vertices[indices[i].vtxj].x);
			bbox.bottom.x = std::min(bbox.bottom.x,vertices[indices[i].vtxk].x);

			bbox.bottom.y = std::min(bbox.bottom.y,vertices[indices[i].vtxi].y);
			bbox.bottom.y = std::min(bbox.bottom.y,vertices[indices[i].vtxj].y);
			bbox.bottom.y = std::min(bbox.bottom.y,vertices[indices[i].vtxk].y);

			bbox.bottom.z = std::min(bbox.bottom.z,vertices[indices[i].vtxi].z);
			bbox.bottom.z = std::min(bbox.bottom.z,vertices[indices[i].vtxj].z);
			bbox.bottom.z = std::min(bbox.bottom.z,vertices[indices[i].vtxk].z);
		}
		return bbox;
	}
	void build_tree(Boxs_tree* boxs_tree ,int i0,int i1){
		BBox box = build_box(i0,i1);
		boxs_tree->i0 = i0;
		boxs_tree->i1 = i1;
		boxs_tree->box = box;


		//trouver la plus grande largeur
		int argmax_largeur = 0;
		double max_largeur = 0;

		std::vector<double> coord_top{box.top.x , box.top.y , box.top.z};
		std::vector<double> coord_bottom{box.bottom.x, box.bottom.y , box.bottom.z};

		for(int i = 0; i<3;i++){
			double largeur_i = coord_top[i] - coord_bottom[i];
			if(largeur_i > max_largeur){
				argmax_largeur = i;
				max_largeur = largeur_i;
			}
		}
		double middle = (coord_top[argmax_largeur] + coord_bottom[argmax_largeur])/2;
		// diviser la liste
		int ind_pivot = i0 - 1;
		TriangleIndices temp;
		for (int ind = i0; ind < i1;ind++){
			std::vector<double> A{vertices[indices[ind].vtxi].x,vertices[indices[ind].vtxi].y,vertices[indices[ind].vtxi].z};
			std::vector<double> B{vertices[indices[ind].vtxj].x,vertices[indices[ind].vtxj].y,vertices[indices[ind].vtxj].z};
			std::vector<double> C{vertices[indices[ind].vtxk].x,vertices[indices[ind].vtxk].y,vertices[indices[ind].vtxk].z};
			if(A[argmax_largeur] + B[argmax_largeur] + C[argmax_largeur] <= 3*middle){
				// swap avec le pivot
				ind_pivot++;
				std::swap(indices[ind_pivot],indices[ind]); // @suppress("Invalid arguments")

			}
		}
		if (ind_pivot < i0 || ind_pivot >= i1 - 1 || i1 == i0 + 1)return;

		// creer les deux branches, gauche et droite
		boxs_tree->left_branch = new Boxs_tree();
		build_tree(boxs_tree->left_branch,i0,ind_pivot+1);

		boxs_tree->right_branch = new Boxs_tree();
		build_tree(boxs_tree->right_branch,ind_pivot+1,i1);
	}
	void add_textures(std::string object_directory){
		std::vector<const char *> files_names;
		files_names.push_back("visage.bmp");
		files_names.push_back("cheveux.bmp");
		files_names.push_back("corps.bmp");
		files_names.push_back("pantalon.bmp");
		files_names.push_back("accessoires.bmp");
		files_names.push_back("mains.bmp");

		for (int i =  0; i <6; i++){
			int x,y, channels_in_file;
			stbi_uc* a = stbi_load((object_directory + files_names[i]).c_str(), &x, &y, &channels_in_file,3);
			textures.push_back(a);
			textures_x.push_back(x);
			textures_y.push_back(y);
		}
	}

	virtual bool intersection(const Ray &r, double &tInterGeometry,  Vector &N, Vector&albedo){

		if(!boxs_tree.box.intersection(r)) return false;

		tInterGeometry = DBL_MAX;
		Vector temp_albedo;

		bool geometryHasInter = false;
		std::list<const Boxs_tree*> l;
		l.push_back(&boxs_tree);

		while(!l.empty()){

			const Boxs_tree* current_branch = l.back();
			l.pop_back();
			if( current_branch->left_branch && current_branch->left_branch->box.intersection(r)){
				l.push_back(current_branch->left_branch);
			}
			if(current_branch->right_branch && current_branch->right_branch->box.intersection(r)){
				l.push_back(current_branch->right_branch);
			}
			if( !current_branch->left_branch){
				double tInterTriangle;
				for ( int ind = current_branch->i0; ind < current_branch->i1; ind++){

					Vector A = vertices[indices[ind].vtxi];
					Vector B = vertices[indices[ind].vtxj];
					Vector C = vertices[indices[ind].vtxk];

					Triangle T(A,B,C,Material(Diffus(),Vector(1,1,1)));
					double alpha,beta,gamma;
					bool hasInter = T.intersection_triangle(r, tInterTriangle,alpha,beta,gamma, temp_albedo);

					if (hasInter  && tInterTriangle < tInterGeometry) {
						Vector NA =normals[indices[ind].ni];
						Vector NB =normals[indices[ind].nj];
						Vector NC =normals[indices[ind].nk];

						N = alpha*NA + beta*NB + gamma*NC; N.normalize();
						tInterGeometry = tInterTriangle;

						int group = indices[ind].group;
						int H = textures_x[group];
						int W = textures_y[group];

						int x = (uvs[indices[ind].uvi].x * alpha + uvs[indices[ind].uvj].x * beta + uvs[indices[ind].uvk].x * gamma)*(W - 1);
						int y = H - (uvs[indices[ind].uvi].y * alpha + uvs[indices[ind].uvj].y * beta + uvs[indices[ind].uvk].y * gamma)*(H - 1);

						albedo.x = textures[group][(y * W + x) * 3 + 0]/255.;
						albedo.y = textures[group][(y * W + x) * 3 + 1]/255.;
						albedo.z = textures[group][(y * W + x) * 3 + 2]/255.;

						geometryHasInter = true;
					}
				}

			}
		}
		return geometryHasInter;
	}

	void readOBJ(const char* obj);

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
};


class Scene {
public:
	Scene(const double &Indice_refract, Sphere* lampe):indice_refract(Indice_refract),lampe(lampe) {
		objects.push_back(lampe);
	};
	double indice_refract;
	std::vector<Object*> objects;
	Sphere* lampe;

	void addSphere(Sphere *sphere) {
		objects.push_back(sphere);
	};

	void addGeometry(Geometry *geometry) {
		objects.push_back(geometry);
	};

	void addTriangle(Triangle *triangle) {
		objects.push_back(triangle);
	};

	bool intersection(const Ray &r, Vector &P, Vector &N, Material &M) {
		//*******trouver le premier contact*******
		Vector tempN;
		int indiceSphere;
		int i = 0;
		bool SceneHasInter = false;
		double tInterScene = DBL_MAX;
		Vector tempAlbedo;
		Vector albedo;
		for (auto & object : this -> objects) {
			double tInterSphere;
			bool SphereHasInter = object->intersection(r, tInterSphere,tempN,tempAlbedo);

			if (SphereHasInter) {
				if (tInterSphere < tInterScene) {
					N = tempN;
					albedo = tempAlbedo;
					tInterScene = tInterSphere;
					indiceSphere = i;
					SceneHasInter = true;
				}
			}
			i+=1;
		}

		if (!SceneHasInter) return false;


		M = objects[indiceSphere]->material;M.albedo = albedo;
		P = r.origin + tInterScene * r.direction;


		return true;
	}
	bool intersection_shadow(const Ray &r, double & tInterScene) {
		//*******trouver le premier contact*******
		bool SceneHasInter = false;
		tInterScene = DBL_MAX;
		for (auto & object : this -> objects) {
			double tInterSphere;
			Vector N;
			Vector albedo;
			bool SphereHasInter = object->intersection(r, tInterSphere,N,albedo);

			if (SphereHasInter) {
				if (tInterSphere < tInterScene) {

					tInterScene = tInterSphere;
					SceneHasInter = true;
				}
			}
		}
		if (!SceneHasInter) return false;
		return true;
	}
};

