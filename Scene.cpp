#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		

		for (int y = 0; y < camera->verRes; y++)
		{
			vector<Color> rowOfColors;
			for (int x = 0; x < camera->horRes; x++)
			{
				rowOfColors.push_back(backgroundColor);
			}
			image.push_back(rowOfColors);
		}

	}
	else
	{
		for (int j = 0; j < camera->verRes; j++)
		{
			for (int i = 0; i < camera->horRes; i++)
			{
				this->image[j][i].r = this->backgroundColor.r;
				this->image[j][i].g = this->backgroundColor.g;
				this->image[j][i].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}


/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int y = camera->verRes - 1; y >= 0; y--)
	{
		for (int x = 0; x < camera->horRes; x++)
		{
			fout << makeBetweenZeroAnd255(image[y][x].r) << " "
				<< makeBetweenZeroAnd255(image[y][x].g) << " "
				<< makeBetweenZeroAnd255(image[y][x].b) << " ";
		}
		fout << endl;
	}

	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}



void Scene::transformVertices() 
{

	for (int mesh_idx=0; mesh_idx < meshes.size(); mesh_idx++) 
	{
		Mesh * mesh_ptr = meshes[mesh_idx];

		vector<Vec3 *> mesh_vertices;
		for (Vec3 * vertex : vertices) {
			mesh_vertices.push_back(new Vec3(*vertex));
		}
		transformed_vertices.push_back(mesh_vertices);


		for (int face_idx=0; face_idx < mesh_ptr->numberOfTriangles; face_idx++) 
		{
			Triangle tri = (mesh_ptr->triangles)[face_idx];
			Vec3 vec1 = *vertices[tri.vertexIds[0] - 1];
			Vec3 vec2 = *vertices[tri.vertexIds[1] - 1];
			Vec3 vec3 = *vertices[tri.vertexIds[2] - 1];


			for (int trans=0; trans < mesh_ptr->numberOfTransformations; trans++) 
			{				
				int transformation_id = mesh_ptr->transformationIds[trans];

				switch (mesh_ptr->transformationTypes[trans])
				{
					case 't':
					{
						Translation translation = *translations[transformation_id-1];
						Vec3 t_values = Vec3(translation.tx, translation.ty, translation.tz);
						vec1 = addVec3(vec1, t_values);
						vec2 = addVec3(vec2, t_values);
						vec3 = addVec3(vec3, t_values);
						break;
					}
					case 's':
					{
						Scaling scaling = *scalings[transformation_id-1];
						Vec3 s_values = Vec3(scaling.sx, scaling.sy, scaling.sz);
						vec1 = multiplyVec3(vec1, s_values);
						vec2 = multiplyVec3(vec2, s_values);
						vec3 = multiplyVec3(vec3, s_values);
						break;
					}
					case 'r':
					{
						Rotation rotation = *rotations[transformation_id-1];
						double angle = rotation.angle * M_PI / 180.0;
						Vec3 u = Vec3(rotation.ux, rotation.uy, rotation.uz);
						u = normalizeVec3(u);

						vec1 = rotateVec3(vec1, u, angle);
						vec2 = rotateVec3(vec2, u, angle);
						vec3 = rotateVec3(vec3, u, angle);
						break;
					}
				}
			}

			transformed_vertices[mesh_idx][tri.vertexIds[0] - 1]->x = vec1.x;
			transformed_vertices[mesh_idx][tri.vertexIds[0] - 1]->y = vec1.y;
			transformed_vertices[mesh_idx][tri.vertexIds[0] - 1]->z = vec1.z;

			transformed_vertices[mesh_idx][tri.vertexIds[1] - 1]->x = vec2.x;
			transformed_vertices[mesh_idx][tri.vertexIds[1] - 1]->y = vec2.y;
			transformed_vertices[mesh_idx][tri.vertexIds[1] - 1]->z = vec2.z;

			transformed_vertices[mesh_idx][tri.vertexIds[2] - 1]->x = vec3.x;
			transformed_vertices[mesh_idx][tri.vertexIds[2] - 1]->y = vec3.y;
			transformed_vertices[mesh_idx][tri.vertexIds[2] - 1]->z = vec3.z;

		}
	}

}



void Scene::projection(Camera * camera) {

	double r = camera->right, l = camera->left;;
	double t = camera->top, b = camera->bottom;
	double f = camera->far, n = camera->near;

	Matrix4 M_project;

	if (camera->projectionType == 0) //orthogonal
	{
		M_project = Matrix4({
			{2/(r-l), 0, 0, -(r+l)/(r-l)},
			{0, 2/(t-b), 0, -(t+b)/(t-b)},
			{0, 0, -2/(f-n), -(f+n)/(f-n)},
			{0, 0, 0, 1}
		});
	}
	else //perspective
	{
		M_project = Matrix4({
			{2*n/(r-l), 0, (r+l)/(r-l), 0},
			{0, 2*n/(t-b), (t+b)/(t-b), 0},
			{0, 0, -(f+n)/(f-n), -2*f*n/(f-n)},
			{0, 0, -1, 0}
		});
	}


	int num_of_meshes = meshes.size();
	int num_of_vecs = vertices.size();

	for (int mesh_idx = 0; mesh_idx < num_of_meshes; mesh_idx++)
	{
		for (int vec_idx = 0; vec_idx < num_of_vecs; vec_idx++) 
		{
			Vec4 * vec_ptr = camera->cam_trans_vertices[mesh_idx][vec_idx];

			*vec_ptr = multiplyMatrixWithVec4(M_project, *vec_ptr);
		}
	}
}


void Scene::perspective_divide(Camera * camera) {
	int num_of_meshes = meshes.size();
	int num_of_vecs = vertices.size();

	for (int mesh_idx = 0; mesh_idx < num_of_meshes; mesh_idx++)
	{
		for (int vec_idx = 0; vec_idx < num_of_vecs; vec_idx++) 
		{
			Vec4 * vec_ptr = camera->cam_trans_vertices[mesh_idx][vec_idx];
			double w = vec_ptr->t; 

			if (w==0) w = EPSILON;

			vec_ptr->x = vec_ptr->x / w;
			vec_ptr->y = vec_ptr->y / w;
			vec_ptr->z = vec_ptr->z / w;
			vec_ptr->t = 1;
		}
	}

}

void Scene::viewport_trans(Camera * camera) {
	double nx = camera->horRes, ny = camera->verRes;

	Matrix4 M_vp = Matrix4({
		{nx/2, 0, 0, (nx-1)/2},
		{0, ny/2, 0, (ny-1)/2},
		{0, 0, (double)1/2, (double)1/2},
		{0, 0, 0, 0} //not important
	});
	
	int num_of_meshes = meshes.size();
	int num_of_vecs = vertices.size();

	for (int mesh_idx = 0; mesh_idx < num_of_meshes; mesh_idx++)
	{
		if (meshes[mesh_idx]->type == WIREFRAME_MESH) continue;
		for (int vec_idx = 0; vec_idx < num_of_vecs; vec_idx++) 
		{
			Vec4 * vec_ptr = camera->cam_trans_vertices[mesh_idx][vec_idx];

			*vec_ptr = multiplyMatrixWithVec4(M_vp, *vec_ptr);
		}
	}

	int num_of_pairs = camera->vertex_pairs.size();

	for (int pair_idx=0; pair_idx < num_of_pairs; pair_idx++)
	{
		Vec4 * vec_ptr0 = camera->vertex_pairs[pair_idx].first;
		*vec_ptr0 = multiplyMatrixWithVec4(M_vp, *vec_ptr0);

		Vec4 * vec_ptr1 = camera->vertex_pairs[pair_idx].second;
		*vec_ptr1 = multiplyMatrixWithVec4(M_vp, *vec_ptr1);
	}
}




void Scene::barycentric_alg(Mesh * mesh_ptr, Camera * camera, const vector<Vec4 *> & tr_vertices) {
	for (Triangle &triangle : mesh_ptr->triangles) 
	{
		if (triangle.backface == true) continue;

		Vec4 *v1 = tr_vertices[triangle.vertexIds[0] - 1];
        Vec4 *v2 = tr_vertices[triangle.vertexIds[1] - 1];
        Vec4 *v3 = tr_vertices[triangle.vertexIds[2] - 1];

        Color c1 = *(this->colorsOfVertices[v1->colorId - 1]);
        Color c2 = *(this->colorsOfVertices[v2->colorId - 1]);
        Color c3 = *(this->colorsOfVertices[v3->colorId - 1]);

        int x1 = (int)(v1->x + 0.5);
        int y1 = (int)(v1->y + 0.5);
        double z1 = v1->z;

        int x2 = (int)(v2->x + 0.5);
        int y2 = (int)(v2->y + 0.5);
        double z2 = v2->z;

        int x3 = (int)(v3->x + 0.5);
        int y3 = (int)(v3->y + 0.5);
        double z3 = v3->z;

		int y_min = min3(y1,y2,y3);
		int x_min = min3(x1,x2,x3);

		int y_max = max3(y1,y2,y3);
		int x_max = max3(x1,x2,x3);

		double f12 = line_eqn(*v1, *v2, v3->x, v3->y);
		double f23 = line_eqn(*v2, *v3, v1->x, v1->y);
		double f31 = line_eqn(*v3, *v1, v2->x, v2->y);

		for (int y=y_min-5; y < y_max+5; y++) 
		{
			if ((y < 0)||(y > (camera->verRes - 1))) continue;

			for (int x=x_min-5; x < x_max+5; x++) 
			{
				if ((x < 0)||(x > (camera->horRes - 1))) continue;

				double alpha = line_eqn(*v1, *v2, x, y) / f12;
				double beta = line_eqn(*v2, *v3, x, y) / f23;
				double omega = line_eqn(*v3, *v1, x, y) / f31;


				double z = alpha * v3->z + beta * v1->z + omega * v2->z;


				if (((alpha>=0)&&(beta>=0)&&(omega>=0)) && (z < camera->depth_buff[y][x])) { //check inside triangle and in front
					camera->depth_buff[y][x] = z;


					double r = alpha * c3.r + beta * c1.r + omega * c2.r;
					double g = alpha * c3.g + beta * c1.g + omega * c2.g;
					double b = alpha * c3.b + beta * c1.b + omega * c2.b;

					Color color = Color((int)(r+0.5), (int)(g+0.5), (int)(b+0.5));
					camera->color_buff[y][x] = color;
				}
			}

		}
	}
}

bool Scene::visible(double den, double num, double & te, double & tl) {
	if (den > 0)
	{
		double t = num/den;
		if (t > tl) return false;
		if (t > te) te = t;
	}
	else if (den < 0)
	{
		double t = num/den;
		if (t < te) return false;
		if (t < tl) tl = t;
	}
	else if (num > 0) return false;

	return true;
}

bool Scene::clip_line(Camera * camera, double w, Vec4 * v0, Vec4 * v1) { //returns true if object is fully or partially inside

	double min_val = -w, max_val = w;

	double dx = v1->x - v0->x;
	double dy = v1->y - v0->y;
	double dz = v1->z - v0->z;

	Vec4 * new_v0 = new Vec4(*v0);
	Vec4 * new_v1 = new Vec4(*v1);

	Color c0 = *colorsOfVertices[v0->colorId - 1];
	Color c1 = *colorsOfVertices[v1->colorId - 1];
	Color dcolor = c1 - c0;

	double te = 0, tl = 1;
	bool visible = false;
	if (this->visible(dx, (min_val - v0->x), te, tl))
	{
		if (this->visible(-dx, (v0->x - max_val), te, tl))
		{
			if (this->visible(dy, (min_val - v0->y), te, tl))
			{
				if (this->visible(-dy, (v0->y - max_val), te, tl))
				{
					if (this->visible(dz, (min_val - v0->z), te, tl))
					{
						if (this->visible(-dz, (v0->z - max_val), te, tl)) 
						{
							visible = true;

							if (tl < 1) 
							{
								new_v1->x = v0->x + dx*tl;
								new_v1->y = v0->y + dy*tl;
								new_v1->z = v0->z + dz*tl;		

								Color * new_c1 = new Color(c0 + dcolor*tl);
								colorsOfVertices.push_back(new_c1);
								new_v1->colorId = colorsOfVertices.size();

							}
							if (te > 0)
							{
								new_v0->x = v0->x + dx*te;
								new_v0->y = v0->y + dy*te;
								new_v0->z = v0->z + dz*te;					

								Color * new_c0 = new Color(c0 + dcolor*te);
								colorsOfVertices.push_back(new_c0);
								new_v0->colorId = colorsOfVertices.size();	
							}
						}
					}
				}
			}
		}
	}

	if (visible) camera->vertex_pairs.push_back({new_v0, new_v1});

	return visible;
}

void Scene::clip(Camera * camera) {

	int num_of_meshes = meshes.size();

	for (int mesh_idx=0; mesh_idx < num_of_meshes; mesh_idx++) 
	{
		Mesh * mesh_ptr = meshes[mesh_idx];
		if (mesh_ptr->type == SOLID_MESH) continue;

		vector<Vec4 *> & mesh_tr_vertices = camera->cam_trans_vertices[mesh_idx];

		for (Triangle &triangle : mesh_ptr->triangles) 
		{

			if (triangle.backface == true) continue;

			Vec4 *v0 = mesh_tr_vertices[triangle.vertexIds[0] - 1];
			Vec4 *v1 = mesh_tr_vertices[triangle.vertexIds[1] - 1];
			Vec4 *v2 = mesh_tr_vertices[triangle.vertexIds[2] - 1];

			clip_line(camera, 1, v0, v1);
			clip_line(camera, 1, v1, v2);
			clip_line(camera, 1, v2, v0);

		}

	}


}

void Scene::draw_line_midpoint(int x0, int y0, double z0,
                             int x1, int y1, double z1,
                             Camera *camera, 
                             Color c0, Color c1) {
   
    bool steep = abs(y1 - y0) > abs(x1 - x0);

    if (steep) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        std::swap(z0, z1);
    }

    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
        std::swap(z0, z1);
		swap(c0, c1);
    }

    int dx = x1 - x0;
    int dy = abs(y1 - y0);
    int yStep = (y0 < y1) ? 1 : -1;
    int err = dx / 2;
    int y = y0;

    double dz = dx != 0 ? (z1 - z0) / dx : 0.0;
    double z = z0;

    double dr = dx != 0 ? (c1.r - c0.r) / dx : 0.0;
    double dg = dx != 0 ? (c1.g - c0.g) / dx : 0.0;
    double db = dx != 0 ? (c1.b - c0.b) / dx : 0.0;

    double r = c0.r;
    double g = c0.g;
    double b = c0.b;

    for (int x = x0; x <= x1; x++) {
        int plotX = steep ? y : x;
        int plotY = steep ? x : y;

        if (plotX >= 0 && plotX < camera->horRes && plotY >= 0 && plotY < camera->verRes) {
            
            if (z < camera->depth_buff[plotY][plotX]) { //depth check
                camera->depth_buff[plotY][plotX] = z;
				Color color = Color((int)(r+0.5), (int)(g+0.5), (int)(b+0.5));
                camera->color_buff[plotY][plotX] = color;
            }
        }

        err -= dy;
        if (err < 0) {
            y += yStep;
            err += dx;
        }

        r += dr; 
		g += dg; 
		b += db;
		z += dz;
    }
}



void Scene::bfc(Camera *camera) {

	int num_of_meshes = meshes.size();

	for (int mesh_idx=0; mesh_idx < num_of_meshes; mesh_idx++) 
	{
		Mesh * mesh_ptr = meshes[mesh_idx];

		vector<Vec4 *> mesh_tr_vertices = camera->cam_trans_vertices[mesh_idx];

		for (Triangle & triangle : mesh_ptr->triangles) 
		{
			Vec4 *v0 = mesh_tr_vertices[triangle.vertexIds[0] - 1];
			Vec4 *v1 = mesh_tr_vertices[triangle.vertexIds[1] - 1];
			Vec4 *v2 = mesh_tr_vertices[triangle.vertexIds[2] - 1];

			Vec3 vec3_v0(v0->x, v0->y, v0->z);
			Vec3 vec3_v1(v1->x, v1->y, v1->z);
			Vec3 vec3_v2(v2->x, v2->y, v2->z);

			Vec3 edge1 = normalizeVec3(subtractVec3(vec3_v0, vec3_v2));
			Vec3 edge2 = normalizeVec3(subtractVec3(vec3_v1, vec3_v2));

			Vec3 normal = crossProductVec3(edge1, edge2);
			normal = normalizeVec3(normal);

			Vec3 ray = normalizeVec3(vec3_v0); //v0 - camera_position

			if (dotProductVec3(normal, ray) < 0) triangle.backface = true;
			else triangle.backface = false;
		}

	}
}



/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	Matrix4 M_cam({
		{camera->u.x, camera->u.y, camera->u.z, -(camera->u.x*camera->position.x + camera->u.y*camera->position.y + camera->u.z*camera->position.z)},
		{camera->v.x, camera->v.y, camera->v.z, -(camera->v.x*camera->position.x + camera->v.y*camera->position.y + camera->v.z*camera->position.z)},
		{camera->w.x, camera->w.y, camera->w.z, -(camera->w.x*camera->position.x + camera->w.y*camera->position.y + camera->w.z*camera->position.z)},
		{0,0,0,1}
	});

	int num_of_meshes = meshes.size();
	int num_of_vecs = vertices.size();

	for (int mesh_idx = 0; mesh_idx < num_of_meshes; mesh_idx++)
	{
		vector<Vec4 *> mesh_vertices;

		for (int vec_idx = 0; vec_idx < num_of_vecs; vec_idx++)
		{
			Vec3 current_vec = *transformed_vertices[mesh_idx][vec_idx];
			Vec4 current_vec4(current_vec.x, current_vec.y, current_vec.z, 1, current_vec.colorId);

			Vec4 res = multiplyMatrixWithVec4(M_cam, current_vec4);

			mesh_vertices.push_back(new Vec4(res));
		}

		camera->cam_trans_vertices.push_back(mesh_vertices);
	}

	projection(camera);

	if (camera->projectionType == PERSPECTIVE_PROJECTION) perspective_divide(camera);

	if (cullingEnabled) bfc(camera);

	clip(camera);

	viewport_trans(camera);

	//init depth and color buff
	camera->depth_buff = new double*[camera->verRes];
	camera->color_buff = new Color*[camera->verRes];

	for (int j = 0; j < camera->verRes; ++j) 
	{
        camera->depth_buff[j] = new double[camera->horRes];
		camera->color_buff[j] = new Color[camera->horRes];

        for (int i = 0; i < camera->horRes; ++i) 
		{
            camera->depth_buff[j][i] = numeric_limits<double>::max();
			camera->color_buff[j][i] = Color(backgroundColor);
        }
    }

	int num_of_pairs = camera->vertex_pairs.size();

	for (int pair_idx=0; pair_idx < num_of_pairs; pair_idx++)
	{
		Vec4 * v0 = camera->vertex_pairs[pair_idx].first;
		Vec4 * v1 = camera->vertex_pairs[pair_idx].second;

		Color c0 = *(colorsOfVertices[v0->colorId - 1]);
		Color c1 = *(colorsOfVertices[v1->colorId - 1]);

        int x0 = (int)(v0->x + 0.5);
        int y0 = (int)(v0->y + 0.5);
        double z0 = v0->z;

		int x1 = (int)(v1->x + 0.5);
        int y1 = (int)(v1->y + 0.5);
        double z1 = v1->z;

		draw_line_midpoint(x0, y0, z0, x1, y1, z1, camera, c0, c1);
		
	}

	for (int mesh_idx=0; mesh_idx < num_of_meshes; mesh_idx++) 
	{
		vector<Vec4 *> mesh_tr_vertices = camera->cam_trans_vertices[mesh_idx];
		Mesh * mesh_ptr = meshes[mesh_idx];

		if (mesh_ptr->type == SOLID_MESH) barycentric_alg(mesh_ptr, camera, mesh_tr_vertices); 
	}


	//write to image
	for (int j = 0; j < camera->verRes; ++j) 
	{
        for (int i = 0; i < camera->horRes; ++i) 
		{
			image[j][i] = camera->color_buff[j][i];
		}
	}


}
