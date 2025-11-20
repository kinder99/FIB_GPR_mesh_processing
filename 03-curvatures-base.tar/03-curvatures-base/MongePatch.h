#ifndef _MONGE_PATCH_INCLUDE
#define _MONGE_PATCH_INCLUDE


#include <vector>
#include <eigen/Core>
#include <eigen/Dense>
#include "glm/glm.hpp"


using namespace std;


class MongePatch
{

public:
	void init(const glm::vec3 &P, const glm::vec3 &normal, const vector<glm::vec3> &closest);
	
	void principalCurvatures(float &kmin, float &kmax) const;

	Eigen::VectorXf computeQ(float u, float v, float w) const;

	std::vector<glm::vec3> getNeighbors(){return closest_t;}
	void setNeighbors(std::vector<glm::vec3> vec){closest_t = vec;}

private:
	std::vector<glm::vec3> closest_t;
	
};


#endif // _MONGE_PATCH_INCLUDE


