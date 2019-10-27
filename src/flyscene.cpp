#include "flyscene.hpp"
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

  // the debug ray is a cylinder, set the radius and length of the cylinder
  // ray.setSize(0.005, 10.0);

  // craete a first debug ray pointing at the center of the screen
  createDebugRay(Eigen::Vector2f(width / 2.0, height / 2.0));

  glEnable(GL_DEPTH_TEST);


  /**
  * @brief Loop over all vertices and get the min and max vertex.
  *
  */
  min = Eigen::Vector3f(INFINITY, INFINITY, INFINITY);
  max = Eigen::Vector3f(-INFINITY, -INFINITY, -INFINITY);
  for (int i = 0; i < mesh.getNumberOfVertices(); i++) {
	  if (mesh.getVertex(i)[0] < min[0]) {
		  std::cout << "OLD MIN: " << min[0] << "  ::  NEW MIN: " << mesh.getVertex(i)[0] << std::endl;
		  min[0] = mesh.getVertex(i)[0];
	  }
	  else if (mesh.getVertex(i)[0] > max[0]) {
		  std::cout << "OLD MAX: " << max[0] << "  ::  NEW MAX: " << mesh.getVertex(i)[0] << std::endl;
		  max[0] = mesh.getVertex(i)[0];
	  }

	  if (mesh.getVertex(i)[1] < min[1]) {
		  std::cout << "OLD MIN: " << min[1] << "  ::  NEW MIN: " << mesh.getVertex(i)[1] << std::endl;
		  min[1] = mesh.getVertex(i)[1];
	  }
	  else if (mesh.getVertex(i)[1] > max[1]) {
		  std::cout << "OLD MAX: " << max[1] << "  ::  NEW MAX: " << mesh.getVertex(i)[1] << std::endl;
		  max[1] = mesh.getVertex(i)[1];
	  }

	  if (mesh.getVertex(i)[2] < min[2]) {
		  std::cout << "OLD MIN: " << min[2] << "  ::  NEW MIN: " << mesh.getVertex(i)[2] << std::endl;
		  min[2] = mesh.getVertex(i)[2];
	  }
	  else if (mesh.getVertex(i)[2] > max[2]) {
		  std::cout << "OLD MAX: " << max[2] << "  ::  NEW MAX: " << mesh.getVertex(i)[2] << std::endl;
		  max[2] = mesh.getVertex(i)[2];
	  }
  }
  std::cout << "Min vector: " << min.transpose() << std::endl;
  std::cout << "Max vector: " << max.transpose() << std::endl;


  for (int i = 0; i<mesh.getNumberOfFaces(); ++i){
    Tucano::Face face = mesh.getFace(i);    
    for (int j =0; j<face.vertex_ids.size(); ++j){
      std::cout<<"vid "<<j<<" "<<face.vertex_ids[j]<<std::endl;
      std::cout<<"vertex "<<mesh.getVertex(face.vertex_ids[j]).transpose()<<std::endl;
      std::cout<<"normal "<<mesh.getNormal(face.vertex_ids[j]).transpose()<<std::endl;
    }
    std::cout<<"mat id "<<face.material_id<<std::endl<<std::endl;
    std::cout<<"face   normal "<<face.normal.transpose() << std::endl << std::endl;
  }
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

  drawBoundingBox();

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
  std::cout << "Intersection: " << boundingIntersection(screen_pos, dir) << std::endl;
  
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

  std::cout << "hey: " << origin.transpose() << std::endl;

  // for every pixel shoot a ray from the origin through the pixel coords
  for (int j = 0; j < image_size[1]; ++j) {
    for (int i = 0; i < image_size[0]; ++i) {
      // create a ray from the camera passing through the pixel (i,j)
      screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
      // launch raytracing for the given ray and write result to pixel data
      pixel_data[i][j] = traceRay(origin, screen_coords);
    }
  } 
  // write the ray tracing result to a PPM image
  Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
}

bool Flyscene::triangleIntersect(float &t, const Eigen::Vector3f origin, const Eigen::Vector3f dest, const Eigen::Vector3f v0, const Eigen::Vector3f v1, const Eigen::Vector3f v2) {
  Eigen::Matrix<float, 3, 3> mat;
  Eigen::Vector3f dir = (dest - origin).normalized();
  Eigen::Vector3f rv0 = v0 - origin;

  mat << v0[0] - v1[0], v0[0] - v2[0], dir[0]
       , v0[1] - v1[1], v0[1] - v2[1], dir[1]
       , v0[2] - v1[2], v0[2] - v2[2], dir[2];

  Eigen::Vector3f solution = (mat.inverse() * rv0);
  t = solution[2];

  if(solution[0] < 0) return false;
  if(solution[1] < 0) return false;
  if(solution[0] + solution[1] - 1 > 0) return false;

  return true;
}

/*
Ray can be described as U = Origin (O) + t*Dest (D)
For each plane (x,y,z) check for which t the ray intersects its lower and upper bounds.
	so U = target coordinate, for instance lower bound x = 2
Then you get 2 = O + t*D
Solve for t.

The minimal t of the 2 per plane gives the one hit first.
The max of the 3 min t's gives the plane that intersects the box.
*/
bool Flyscene::boundingIntersection (Eigen::Vector3f& origin,
												Eigen::Vector3f& destination) {

	Eigen::Vector3f dest = (destination - origin).normalized();
	float tmin = 0;
	float tylb = 0;
	float tzlb = 0;
	float tmax = 0;
	float tyub = 0;
	float tzub = 0;
	//Calculate the t for intersecting the lower bound and upper bound planes of 3 dimensions
	//txlb and txub are named tmin and max for efficiency later

	//Set x plane upper and lower bounds
	if (dest[0] > 0) {
		tmin = (min[0] - origin[0]) / dest[0];
		tmax = (max[0] - origin[0]) / dest[0];
	}
	else if (dest[0] == 0) {
		//Handle edge case of dest = -0
		/*Eigen::Vector3f invdest = dest.inverse();
		tmin = (min[0] - origin[0]) / invdest[0];
		tmax = (max[0] - origin[0]) / invdest[0];*/
	}
	else {
		tmin = (max[0] - origin[0]) / dest[0];
		tmax = (min[0] - origin[0]) / dest[0];
	}

	//Set y plane upper and lower bounds
	if (dest[1] > 0) {
		tylb = (min[1] - origin[1]) / dest[1];
		tyub = (max[1] - origin[1]) / dest[1];
	}
	else if (dest[1] == 0) {
		////Handle edge case of dest = -0
		//Eigen::Vector3f invdest = dest.inverse();
		//tmin = (min[1] - origin[1]) / invdest[1];
		//tmax = (max[1] - origin[1]) / invdest[1];
	}
	else {
		tylb = (max[1] - origin[1]) / dest[1];
		tyub = (min[1] - origin[1]) / dest[1];
	}

	//If a lower bound t is greater than an upper bound t of another dimension we do not hit the box
	if ((tmin > tyub) || (tylb > tmax)) return false;

	//test for new highest lower bound and lowest higher bound
	if (tylb > tmin) tmin = tylb;
	if (tyub < tmax) tmax = tyub;


	//Set z plane upper and lower bounds
	if (dest[1] > 0) {
		tzlb = (min[2] - origin[2]) / dest[2];
		tzub = (max[2] - origin[2]) / dest[2];
	}
	else if (dest[2] == 0) {
		////Handle edge case of dest = -0
		//Eigen::Vector3f invdest = dest.inverse();
		//tmin = (min[2] - origin[2]) / invdest[2];
		//tmax = (max[2] - origin[2]) / invdest[2];
	}
	else {
		tzlb = (max[2] - origin[2]) / dest[2];
		tzub = (min[2] - origin[2]) / dest[2];
	}

	if ((tmin > tzub) || (tylb > tzub)) return false;

	//once more for z test for new highest low and lowest high
	if (tzlb > tmin) tmin = tzlb;
	if (tzub < tmax) tmax = tzub;

	std::cout << "Lower t: " << tmin << " : " << tylb << " : " << tzlb << std::endl;
	std::cout << "Upper t: " << tmax << " : " << tyub << " : " << tzub << std::endl;

	return true;
}

Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f &origin,
                                   Eigen::Vector3f &dest) {
  Eigen::Affine3f modelMatrix = mesh.getShapeModelMatrix();

  float t, tmin;
  t = tmin = INFINITY;
   
  for (int i = 0; i<mesh.getNumberOfFaces(); ++i){
     Tucano::Face face = mesh.getFace(i);    

     Eigen::Vector3f v1 = (modelMatrix * mesh.getVertex(face.vertex_ids[0])).head<3>();
     Eigen::Vector3f v2 = (modelMatrix * mesh.getVertex(face.vertex_ids[1])).head<3>();
     Eigen::Vector3f v3 = (modelMatrix * mesh.getVertex(face.vertex_ids[2])).head<3>();

     if(triangleIntersect(t, origin, dest, v1, v2, v3))
       tmin = (t < tmin) ? t : tmin;
     
    //std::cout<<"mat id "<<face.material_id<<std::endl<<std::endl;
    //std::cout<<"face   normal "<<face.normal.transpose() << std::endl << std::endl;
  }
  return (tmin != INFINITY) ? Eigen::Vector3f(1.0, 1.0, 1.0) : Eigen::Vector3f(0.0, 0.0, 0.0);
}

void Flyscene::drawBoundingBox() {
	Eigen::Vector3f maxVector = mesh.getBoundingMax();
	Eigen::Vector3f minVector = mesh.getBoundingMin();

	float w = max[0] - min[0];
	float h = max[1] - min[1];
	float d = max[2] - min[2];
	Tucano::Shapes::Box bBox = Tucano::Shapes::Box(w, h, d);
	bBox.setColor(Eigen::Vector4f(1, 0, 0, 1));

	Eigen::Affine3f mMatrix = bBox.getModelMatrix();
	mMatrix.translate(Eigen::Vector3f(1.0, 0, 0));
	bBox.setModelMatrix(mMatrix);
	bBox.normalizeModelMatrix();

	Eigen::Affine3f sMatrix = mesh.getShapeModelMatrix();
	sMatrix = sMatrix.linear();
	sMatrix = sMatrix.inverse();
	

	bBox.setModelMatrix(sMatrix);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	bBox.render(flycamera, scene_light);
	bBox.normalizeModelMatrix();

	//Temp
	//std::cout << "Scalars: " << w << " : " << h << " : " << d << std::endl;

	Tucano::Shapes::Box testBox = Tucano::Shapes::Box(1, 1, 1);
	testBox.setColor(Eigen::Vector4f(1, 1, 0, 1));

	mMatrix = testBox.getModelMatrix();
	mMatrix.translate(Eigen::Vector3f(-1.0, 0, 0));
	
	testBox.setModelMatrix(mMatrix);
	testBox.render(flycamera, scene_light);
	// end temp

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}
