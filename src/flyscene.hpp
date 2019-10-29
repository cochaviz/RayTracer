#ifndef __FLYSCENE__
#define __FLYSCENE__

// Must be included before glfw.
#include <GL/glew.h>

#include <GLFW/glfw3.h>

#include <tucano/effects/phongmaterialshader.hpp>
#include <tucano/mesh.hpp>
#include <tucano/shapes/camerarep.hpp>
#include <tucano/shapes/cylinder.hpp>
#include <tucano/shapes/sphere.hpp>
#include <tucano/utils/flycamera.hpp>
#include <tucano/utils/imageIO.hpp>
#include <tucano/utils/mtlIO.hpp>
#include <tucano/utils/objimporter.hpp>

class Flyscene {

public:

  Flyscene(void) {}

  /**
   * @brief Initializes the shader effect
   * @param width Window width in pixels
   * @param height Window height in pixels
   */
  void initialize(int width, int height);

  /**
   * Repaints screen buffer.
   **/
  virtual void paintGL();

  /**
   * Perform a single simulation step.
   **/
  virtual void simulate(GLFWwindow *window);

  /**
   * Returns the pointer to the flycamera instance
   * @return pointer to flycamera
   **/
  Tucano::Flycamera *getCamera(void) { return &flycamera; }

  /**
   * @brief Add a new light source
   */
  void addLight(void) { lights.push_back(flycamera.getCenter()); }

  /**
   * @brief Create a debug ray at the current camera location and passing
   * through pixel that mouse is over
   * @param mouse_pos Mouse cursor position in pixels
   */
  void createDebugRay(const Eigen::Vector2f &mouse_pos);

  /**
   * @brief raytrace your scene from current camera position   
   */
  void raytraceScene(int width = 0, int height = 0);

  /**
   * @brief trace a single ray from the camera passing through dest
   * @param origin Ray origin
   * @param dest Other point on the ray, usually screen coordinates
   * @return a RGB color
   */
  Eigen::Vector3f traceRay(Eigen::Vector3f &origin, Eigen::Vector3f &dest);
  
  Eigen::Vector3f traceRayRecursive(Eigen::Vector3f &origin, Eigen::Vector3f &dest, int depth);

  /**
   * @briefs Check if the current ray intersects the given triangle, and returns the length of the ray (stored in variable t)
   * @param the length of the ray (valid if function returns true)
   * @param origin Ray origin
   * @param dest Other point on the ray, usually screen coordinates
   * @param v0 first vector of the triangle
   * @param v1 second vector of the triangle
   * @param v2 third vector of the triangle
   * @return boolean check; whether or not the ray intersects with the triangle; validity of t
   */
  bool triangleIntersect(float &t, const Eigen::Vector3f origin, const Eigen::Vector3f dest, const Eigen::Vector3f v0, const Eigen::Vector3f v1, const Eigen::Vector3f v2);

  Eigen::Vector3f pointShading(float& t, const int iterations, const Tucano::Material::Mtl material, const Eigen::Vector3f p, const Eigen::Vector3f n, const Eigen::Vector3f dir, const Eigen::Vector3f light_position);
  
  	/**
	* @brief Loops over all the faces and checks whether the ray (cast to the light(s)) intersects with any face
	* @param t the length of the ray (valid if function returns true)
	* @param p Ray origin - the closest intersection point from traceRay
	* @param dest The location of the light source
	* @param face The face containg the point of intersection p
	* @return whether a point is in shadow
	*/
	bool isInShadow(float& t, const Eigen::Vector3f p, const Eigen::Vector3f dest, const Eigen::Vector3f n);
  
	Eigen::Vector3f random_unit_vector();
  
	void generateBoundingBox();

private:
	// A simple phong shader for rendering meshes
	Tucano::Effects::PhongMaterial phong;

	// A fly through camera
	Tucano::Flycamera flycamera;

	// the size of the image generated by ray tracing
	Eigen::Vector2i raytracing_image_size;

	// A camera representation for animating path (false means that we do not
	// render front face)
	Tucano::Shapes::CameraRep camerarep = Tucano::Shapes::CameraRep(false);

	// a frustum to represent the camera in the scene
	Tucano::Shapes::Sphere lightrep;

	// light sources for ray tracing
	vector<Eigen::Vector3f> lights;

	// Scene light represented as a camera
	Tucano::Camera scene_light;

	/// A very thin cylinder to draw a debug ray
	Tucano::Shapes::Cylinder ray = Tucano::Shapes::Cylinder(0.1, 1.0, 16, 64);

	// Scene meshes
	Tucano::Mesh mesh;

	/// MTL materials
	vector<Tucano::Material::Mtl> materials;
};

#endif // FLYSCENE
