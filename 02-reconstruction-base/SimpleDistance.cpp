#include "SimpleDistance.h"
#include <iostream>


/* Initialize everything to be able to compute the implicit distance of [Hoppe92] 
   at arbitrary points that are close enough to the point cloud.
 */

void SimpleDistance::init(const PointCloud *pointCloud, float samplingRadius)
{
	cloud = pointCloud;
	knn.setPoints(&pointCloud->getPoints());
	radius = samplingRadius;
}


/* This operator returns a boolean that if true signals that the value parameter
   has been modified to contain the value of the implicit function of [Hoppe92]
   at point P.
 */

bool SimpleDistance::operator()(const glm::vec3 &P, float &value) const
{
	std::vector<size_t> neighbors;
	std::vector<float> sqDists;
	knn.getKNearestNeighbors(P,1,neighbors,sqDists);
	//std::cout << "index: " << index << std::endl;

	std::vector<glm::vec3> points = cloud->getPoints();
	std::vector<glm::vec3> normals = cloud->getNormals();

	glm::vec3 pi = points.at(neighbors.at(0));
	//std::cout << "point index ok" << std::endl;
	glm::vec3 ni = normals.at(neighbors.at(0));
	//std::cout << "normal index ok" << std::endl;

	//projection of p onto tangent plane at pi
	glm::vec3 z = pi - ((P - pi)*ni)*ni;

	if(glm::distance(z,pi)<=this->radius){
		value = glm::dot((P - pi),ni);
		return true;
	}
	return false;
}






