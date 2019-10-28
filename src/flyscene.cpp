#include "flyscene.hpp"
#include <cmath>
#include <GLFW/glfw3.h>
#include <limits>
#include <tucano/shapes/box.hpp>

void Flyscene::initialize(int width, int height) {
  // initiliaze the Phong Shading effect for the Opengl Previewer
  phong.initialize();

  // set the camera's projection matrix
  flycamera.setPerspectiveMatrix(60.0, width / (float)height, 0.1f, 100.0f);
  flycamera.setViewport(Eigen::Vector2f((float)width, (float)height));

  // load the OBJ file and materials
  Tucano::MeshImporter::loadObjFile(mesh, materials,
                                    "resources/models/cube.obj");

  // normalize the model (scale to unit cube and center at origin)
  mesh.normalizeModelMatrix();

  // pass all the materials to the Phong Shader
  for (int i = 0; i < materials.size(); ++i)
    phong.addMaterial(materials[i]);

  // set the color and size of the sphere to represent the light sources
  // same sphere is used for all sources
  lightrep.setColor(Eigen::Vector4f(1.0, 1.0, 0.0, 1.0));
  lightrep.setSize(0.15);

  // create a first ray-tracing light source at some random position
  lights.push_back(Eigen::Vector3f(-1.0, 1.0, 1.0));

  // scale the camera representation (frustum) for the ray debug
  camerarep.shapeMatrix()->scale(0.2);

  // craete a first debug ray pointing at the center of the screen
  createDebugRay(Eigen::Vector2f(width / 2.0, height / 2.0));

  glEnable(GL_DEPTH_TEST);
}

void Flyscene::paintGL(void) {

  // update the camera view matrix with the last mouse interactions
  flycamera.updateViewMatrix();
  Eigen::Vector4f viewport = flycamera.getViewport();

  // clear the screen and set background color
  glClearColor(0.9, 0.9, 0.9, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // position the scene light at the last ray-tracing light source
  scene_light.resetViewMatrix();
  scene_light.viewMatrix()->translate(-lights.back());

  ////glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // render the scene using OpenGL and one light source
  phong.render(mesh, flycamera, scene_light);

  // generate bounding box (not rendered)
  generateBoundingBox();

  // render the ray and camera representation for ray debug
  ray.render(flycamera, scene_light);
  camerarep.render(flycamera, scene_light);

  // render ray tracing light sources as yellow spheres
  for (int i = 0; i < lights.size(); ++i) {
    lightrep.resetModelMatrix();
    lightrep.modelMatrix()->translate(lights[i]);
    lightrep.render(flycamera, scene_light);
  }

  // render coordinate system at lower right corner
  flycamera.renderAtCorner();
}

void Flyscene::simulate(GLFWwindow *window) {
  // Update the camera.
  // NOTE(mickvangelderen): GLFW 3.2 has a problem on ubuntu where some key
  // events are repeated: https://github.com/glfw/glfw/issues/747. Sucks.
  float dx = (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS ? 1.0 : 0.0) -
             (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS ? 1.0 : 0.0);
  float dy = (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS
                  ? 1.0
                  : 0.0) -
             (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS
                  ? 1.0
                  : 0.0);
  float dz = (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS ? 1.0 : 0.0) -
             (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS ? 1.0 : 0.0);
  flycamera.translate(dx, dy, dz);
}

void Flyscene::createDebugRay(const Eigen::Vector2f &mouse_pos) {
  ray.resetModelMatrix();
  // from pixel position to world coordinates
  Eigen::Vector3f screen_pos = flycamera.screenToWorld(mouse_pos);

  // direction from camera center to click position
  Eigen::Vector3f dir = (screen_pos - flycamera.getCenter()).normalized();

  // set cylinder length to ray colission distance 
  float rayLength = 10;
  ray.setSize(0.005, rayLength); 

  // set cylinder length to collision distance
  traceRay(screen_pos, dir);
  
  // position and orient the cylinder representing the ray
  ray.setOriginOrientation(flycamera.getCenter(), dir);

  // place the camera representation (frustum) on current camera location, 
  camerarep.resetModelMatrix();
  camerarep.setModelMatrix(flycamera.getViewMatrix().inverse());
}

void Flyscene::raytraceScene(int width, int height) {

  // if no width or height passed, use dimensions of current viewport
  Eigen::Vector2i image_size(width, height);
  if (width == 0 || height == 0) {
    image_size = flycamera.getViewportSize();
  }

  // create 2d vector to hold pixel colors and resize to match image size
  vector<vector<Eigen::Vector3f>> pixel_data;
  pixel_data.resize(image_size[1]);
  for (int i = 0; i < image_size[1]; ++i)
    pixel_data[i].resize(image_size[0]);

  // origin of the ray is always the camera center
  Eigen::Vector3f origin = flycamera.getCenter();
  Eigen::Vector3f screen_coords;
	float progress = image_size[1];

	int barWidth = 45;

	std::cout << "Ray tracing scene..." << std::endl;

  std::cout << "hey: " << origin.transpose() << std::endl;

  // for every pixel shoot a ray from the origin through the pixel coords
  for (int j = 0; j < image_size[1]; ++j) {
	for (int i = 0; i < image_size[0]; ++i) {
		// create a ray from the camera passing through the pixel (i,j)
		screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
		// launch raytracing for the given ray and write result to pixel data
		pixel_data[i][j] = traceRay(origin, screen_coords);

		// Kowalski... Status!
		std::cout << "[";
		int pos = barWidth * (j / progress);
		for (int i = 0; i < barWidth; ++i) {
			if (i <= pos) std::cout << "#";
			else std::cout << "-";
		}
		std::cout << "] " << int((j/ progress) * 100.0) << " %\r";
		std::cout.flush();
    }
  } 
  // write the ray tracing result to a PPM image
  Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
}

bool Flyscene::triangleIntersect(float &t, const Eigen::Vector3f origin, const Eigen::Vector3f dir, const Eigen::Vector3f v0, const Eigen::Vector3f v1, const Eigen::Vector3f v2) {
  // Create matrices and the normalized direction
  Eigen::Matrix<float, 3, 3> mat;
  Eigen::Vector3f rv0 = v0 - origin;

	// Calculate a, b, and t
  mat << v0[0] - v1[0], v0[0] - v2[0], dir[0]
       , v0[1] - v1[1], v0[1] - v2[1], dir[1]
       , v0[2] - v1[2], v0[2] - v2[2], dir[2];

  Eigen::Vector3f solution = (mat.inverse() * rv0); // a, b, and t are stored in first to last indeces of the 'solution' matrix respectively
  t = solution[2];

	// check if the point of intersection lies within the triangle
  if(solution[0] < 0 || solution[1] < 0 || solution[0] + solution[1] - 1 > 0) return false;

	// check if the ray is not pointing backwards
	if(solution[2] < 0) return false;

  return true;
}
Eigen::Vector3f Flyscene::pointShading(const Tucano::Material::Mtl material, const Eigen::Vector3f p, const Eigen::Vector3f n, const Eigen::Vector3f dir, const Eigen::Vector3f light_position) {
	// Normalize normal
	Eigen::Vector3f normal = n.normalized();

	// Assume all lights have the same intesity, and color
	Eigen::Vector3f light_value = Eigen::Vector3f(1.0, 1.0, 1.0);
	
	// Getting all the needed vectors read
	Eigen::Vector3f light_dir = (light_position - p).normalized();
	Eigen::Vector3f reflection = 2 * (light_dir.dot(normal)) * normal - light_dir;
	Eigen::Vector3f eye_dir = -dir;

	// Ambient term
	Eigen::Vector3f ambient = light_value.cwiseProduct(material.getAmbient());

	// Diffuse term
	float costh = (light_dir).dot(normal);
	Eigen::Vector3f diffuse = light_value.cwiseProduct(std::max(0.0f, costh) * material.getDiffuse());

	// Specular term
	float pre = reflection.dot(eye_dir);
	float cosph = std::max(0.0f, pre);
	float spec_term = std::pow(cosph, material.getShininess());
	Eigen::Vector3f specular = light_value.cwiseProduct(spec_term * material.getSpecular());

	// Add it all up
	Eigen::Vector3f color = ambient + diffuse + specular;
	return color;
}

Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f &origin,
                                   Eigen::Vector3f &dest) {
  float t, tmin;
  t = tmin = INFINITY;

	// Face that's closest to the camera
	Tucano::Face min_face;
 
	// Compute ray direction
  Eigen::Vector3f dir = (dest - origin).normalized();
	
	// Get modelMatrix for the given mesh
  Eigen::Affine3f modelMatrix = mesh.getShapeModelMatrix();
 
	// Loop over all the faces in the scene
	for (int i = 0; i<mesh.getNumberOfFaces(); ++i){
		Tucano::Face face = mesh.getFace(i);    

		 // Convert mesh coordinates to world coordinates
     Eigen::Vector3f v1 = (modelMatrix * mesh.getVertex(face.vertex_ids[0])).head<3>();
     Eigen::Vector3f v2 = (modelMatrix * mesh.getVertex(face.vertex_ids[1])).head<3>();
     Eigen::Vector3f v3 = (modelMatrix * mesh.getVertex(face.vertex_ids[2])).head<3>();

		 // Update tmin if the ray has an intersection with the given face, and t is smaller than tmin
		 // or: if we found a face closer the ray
     if(triangleIntersect(t, origin, dir, v1, v2, v3) && t < tmin) {
       	tmin = t;
			 	min_face = face;
			 } 
  }
	// Return the shaded pixel if tmin updated, else return a nice sky color
	if (tmin != INFINITY) {
		// Calculate point of intersection
		// TODO Multiple light sources?
		Eigen::Vector3f p = origin + tmin * dir;
		return pointShading(materials[min_face.material_id], p, min_face.normal, dir, lights[0]);
	} else {
		// Show background color
		return Eigen::Vector3f(85.0/255.0, 191.0/255.0, 225.0/255.0);
	}
}

void Flyscene::generateBoundingBox() {
	// Get max values of the mesh
	Eigen::Vector3f maxVector3 = mesh.getBoundingMax();
	Eigen::Vector3f minVector3 = mesh.getBoundingMin();

	Eigen::Affine3f modelMatrix = mesh.getShapeModelMatrix();

	// Scale to mesh
	Eigen::Vector4f minVector = modelMatrix * Eigen::Vector4f(minVector3[0], minVector3[1], minVector3[2], 1.0f);
	Eigen::Vector4f maxVector = modelMatrix * Eigen::Vector4f(maxVector3[0], maxVector3[1], maxVector3[2], 1.0f);

	// Initialize bounding box width, height, depth
	float w = maxVector[0] - minVector[0];
	float h = maxVector[1] - minVector[1];
	float d = maxVector[2] - minVector[2];
	
	// Generate a box using the size parameters 
	Tucano::Shapes::Box bBox = Tucano::Shapes::Box(w, h, d);
}

