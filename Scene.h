#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include <iostream>
#include <vector>

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<std::vector<Vec3 *>> transformed_vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName, int osType);
	void transformVertices();
	void projection(Camera * camera);
	void perspective_divide(Camera * camera);
	void barycentric_alg(Mesh * mesh_ptr, Camera * camera, const vector<Vec4 *> & tr_vertices);
	void midpoint_alg(Mesh * mesh_ptr, Camera * camera, const vector<Vec4 *> & tr_vertices);
	void forwardRenderingPipeline(Camera *camera);
	void viewport_trans(Camera * camera);
	void draw_line_midpoint(int x0, int y0, double z0,
                             int x1, int y1, double z1,
                             Camera *camera, 
                             Color c0, Color c1);
	void clip(Camera * camera);
	bool clip_line(Camera * camera, double w, Vec4 * v0, Vec4 * v1);
	bool visible(double den, double num, double & te, double & tl);
	void bfc(Camera *camera);
};
#endif
