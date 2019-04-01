// Raytracer.cpp : Defines the entry point for the console application.

#define _CRT_SECURE_NO_WARNINGS // for Visual Studio 2017 (maybe 2015 as well)

#include <iostream>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION

#include "Classes.h"
#include <algorithm>
#include <random>
#include <chrono>
#include <omp.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <thread>

std::default_random_engine engine;
std::uniform_real_distribution<double> distrib(0, 1);


//*******************************************************************************************
//*************************************global variables**************************************
//*******************************************************************************************
const bool anti_aliasing = true ;
const bool indirect_ligh = true;
const bool profondeur = true;

const int nb_ray = 10;
const int nb_rebounce = 4;


const int W = 512*2;
const int H = W;
const double fov = 60 * M_PI / 180;
const double depth = -H/(2*tan(fov*0.5));

const std::string images_directory = "";
const std::string girl_directory= "girl/";


//********************************************************************************************
//****************************************generateRay*****************************************
//********************************************************************************************
// generateRay prend en entrer un rayon et une distance focale et renvoie un rayon
// aleatoire autour de ce dernier en respectant la distance focale
// (si anti_aliasing et profondeur sont activés (voir lignes 26 et 28))
Ray generateRay( Ray principalRay, double distance_focale) {

	Vector origin = principalRay.origin;
	Vector direction = principalRay.direction;

	if (anti_aliasing){
		double x = distrib(engine), y = distrib(engine), R=sqrt(-2*log(x));
		double u =  R*cos(2*M_PI*y)*0.5;
		double v =  R*sin(2*M_PI*y)*0.5;

		direction.x +=  u - 0.5;
		direction.y +=  v - 0.5;

		direction.normalize();
	}
	if(profondeur) {
		double a = (distrib(engine) - .5)*.5;
		double b = (distrib(engine) - .5)*.5;
		if(!anti_aliasing) direction.normalize();
		Vector end = origin + distance_focale * direction;
		origin = origin + Vector(a,b,0);
		direction = (end - origin);direction.normalize();
	}
	if(!anti_aliasing && !profondeur) direction.normalize();
	return Ray(origin,direction);
};

//*********************************************************************************************
//***********************************randomDirection*******************************************
//*********************************************************************************************
// randomDirection prend en entrer un vecteur n et renvoie un vecteur
// aleatoire avec une probabilite cos(theta) autour de n
Vector randomDirection(const Vector & n){
	double r1 = distrib(engine);
	double r2 = distrib(engine);

	double A = std::sqrt(1 - r2);

	double x = cos(2 * M_PI * r1) * A;
	double y = sin(2 * M_PI * r1) * A;
	double z = std::sqrt(r2);

//	Vector alea(distrib(engine),distrib(engine),distrib(engine));
	Vector alea(r1,r2,1);// en choisisant les nombres aleatoires deja generes on optimise jusqu'a 20% le temps de calcul

	Vector tangentiel1 = cross(n,alea);tangentiel1.normalize();
	Vector tangentiel2 = cross(tangentiel1,n);

	return x * tangentiel1 + y * tangentiel2 + z * n;
}

//**********************************************************************************************
//*****************************************geColor**********************************************
//**********************************************************************************************
Vector getColor(Ray &r, int nb_rebound, Scene* scene_pointer) {
	if (nb_rebound == 0) {
		return Vector();
	}
	Vector P;
	Vector n;
	Material m;
	Vector couleur_pixel;

	bool sceneHasInter = scene_pointer->intersection(r, P, n, m);

	if (!sceneHasInter) return Vector();

	if (m.type == LAMPE){
		return 200000*Vector(1,1,1);
	}

	else if (m.type == SPEC) {
		// reflechir le rayon
		r.reflechi(P, n);

		return getColor(r, nb_rebound - 1,scene_pointer);
	}
	else if (m.type == TRANS) {
		double alea = distrib(engine);
		r.refracte(P, n, scene_pointer->indice_refract, m.n_refract,alea);

		return getColor(r, nb_rebound - 1,scene_pointer);
	}
	else /*if (m.type == DIFFUS)*/ {

		Vector albedo_materiau = m.albedo;
		Vector directionAlea;
		// *******contribution direct (source_position spherique)*********
		Vector Ox = P - scene_pointer->lampe->origin; Ox.normalize();
		Vector N_prime = randomDirection(Ox);
		Vector P_prime = N_prime*scene_pointer->lampe->rayon + scene_pointer->lampe->origin;
		directionAlea = P_prime - P;

		double dist_P_P_prime2 = directionAlea.norm2();
		directionAlea.normalize();
		Ray rayAlea(P + 0.001*directionAlea,directionAlea);
		double tInter2;
		// il y a toujours une itersection avec la lampe, c'est pas la peine de renvoyer la valeur scene.intersection2
		scene_pointer->intersection_shadow(rayAlea, tInter2);

		if (tInter2 * tInter2 > dist_P_P_prime2*0.99) {
			couleur_pixel = (scene_pointer->lampe->material.intensite*std::max(0.,dot(directionAlea,n))*dot(-1*directionAlea,N_prime))/
					(4*M_PI*dist_P_P_prime2*dot(N_prime,Ox)) * albedo_materiau;

		}
		//******contribution indirect*********
		if (!indirect_ligh){return couleur_pixel;}

		directionAlea = randomDirection(n);
		Ray rayonAlea(P + 0.001*directionAlea, directionAlea);

		couleur_pixel +=  getColor(rayonAlea, nb_rebound - 1,scene_pointer)*albedo_materiau;
		}

		return couleur_pixel;
	return Vector();
}
//******************************************************************************************
//************************************create scene *****************************************
//******************************************************************************************
Scene* create_scene(){
	double I = 2550000000;
	double lampe_rayon = 3;
	Vector source_position(-10, 20, 30);
	Sphere* lampe = new Sphere(source_position,lampe_rayon,Material(Lampe(),I));

	Scene* scene = new Scene(1.,lampe);

	//*************les materiaux utilisés****************
	Material blue_diffus(Diffus(),Vector(0, 0, 1));
	Material blanc_diffus(Diffus(),Vector(1, 1, 1));
	Material rouge_diffus(Diffus(),Vector(1, 0, 0));
	Material vert_diffus(Diffus(),Vector(0, 1, 0));
	Material gris_diffus(Diffus(),Vector(.5, .5, .5));
	Material jaune_diffus(Diffus(),Vector(1, 1, 0));

	Material transparent(Trans(),1.4);

	//*********ajouter les objets à la la scene************
	scene->addSphere(new Sphere(Vector(0, -1000, 0), 990, blue_diffus));//sphere bas
	scene->addSphere(new Sphere(Vector(0, 1000, 0), 960, blanc_diffus));//sphere haut
	scene->addSphere(new Sphere(Vector(0, 0, -1000), 940, vert_diffus));//sphere devant
	scene->addSphere(new Sphere(Vector(1000, 0, 0), 940, rouge_diffus));//sphere droite
	scene->addSphere(new Sphere(Vector(-1000, 0, 0), 940, jaune_diffus));//sphere gauche
	scene->addSphere(new Sphere(Vector(0, 0, 1000), 940, gris_diffus));//sphere derriere

	scene->addSphere(new Sphere(Vector(9, 4, 0), 4, Material(Spec())));
	scene->addSphere(new Sphere(Vector(-9, 4, 0), 4, transparent));
	scene->addGeometry(new Geometry(girl_directory,rouge_diffus));

	return scene;
}

//******************************************************************************************
//************************************generate image****************************************
//******************************************************************************************

void generate_image(int image_nb, double theta,double distance_focale,Scene* scene_pointer) {
	const Vector centreCamera(0, 0, 55);
	// centre de camera apres rotation
	Vector centreCameraRotated = rotate(centreCamera,-theta);

	// tracer l'image
	std::vector<unsigned char> image(W*H * 3, 0);
	#pragma omp parallel for schedule(dynamic,1) //num_threads(10)
	for (int i = 0; i < H; i++){
		for (int j = 0; j < W; j++) {
			// le rayon emis de l'ecran ( qui pase par le centre du pixel i,j)
			Vector direction(j - W / 2 + 0.5, i - H / 2 + 0.5, depth);
			// on a pas besoin de le normaliser puisqu'il le sera dans la fonction generateRay
//			direction.normalize();

			direction = rotate(direction,-theta);

			Ray centralRay(centreCameraRotated, direction);

			// trouver la couleur du pixel
			Vector couleur_pixel;

			for (int rayInd = 0; rayInd < nb_ray; rayInd++) {
				Ray r = generateRay(centralRay,distance_focale);
				couleur_pixel +=  getColor(r, nb_rebounce,scene_pointer);
			}

			couleur_pixel = couleur_pixel/nb_ray;

			// tracer le pixel 
			image[((H - i - 1) * W + j) * 3 + 0] = std::min(std::pow(couleur_pixel.x, 1 / (2.2)), 255.);
			image[((H - i - 1) * W + j) * 3 + 1] = std::min(std::pow(couleur_pixel.y, 1 / (2.2)), 255.);
			image[((H - i - 1) * W + j) * 3 + 2] = std::min(std::pow(couleur_pixel.z, 1 / (2.2)), 255.);
		}
	}

	std::string image_name = images_directory+"image" + std::to_string(image_nb) + ".png";

	stbi_write_png(image_name.c_str(), W, H, 3, &image[0], 0); // @suppress("Invalid arguments")
};


//******************************************************************************************
//*****************************************main*********************************************
//******************************************************************************************
int main(){
	Scene* scene_pointer = create_scene();
	auto start = std::chrono::high_resolution_clock::now();


	//************** rotate the scene *****************************
	const int nb_images = 1;
	Vector intial_lampe_postion = scene_pointer->lampe->origin;
	for(int i = 0; i<nb_images; i++){
		double theta = 2*M_PI*i/nb_images;
		// rotate light source
		scene_pointer->lampe->origin = rotate(intial_lampe_postion,-theta);

		generate_image(i,theta,55,scene_pointer);
	}

	//************** change focale distance *****************************
//	for(int distance_focale = 30; distance_focale<120; distance_focale += 3){
//
//		generate_image((distance_focale - 30)/3,0,distance_focale,scene_pointer);
//	}


	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
	std::cout << "Elapsed time: " << elapsed.count()/60 << " min\n";
	//637.006 s pour une image
}

