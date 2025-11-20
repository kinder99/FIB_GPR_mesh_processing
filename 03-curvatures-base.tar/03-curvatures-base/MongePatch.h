#ifndef _MONGE_PATCH_INCLUDE
#define _MONGE_PATCH_INCLUDE


#include <vector>
#include "glm/glm.hpp"


using namespace std;


class MongePatch
{

public:
	void init(const glm::vec3 &P, const glm::vec3 &normal, const vector<glm::vec3> &closest);
	
	void principalCurvatures(float &kmin, float &kmax) const;

	void setNeighbors(std::vector<glm::vec3> vec){closest_t = vec;}

private:
	std::vector<glm::vec3> closest_t;
	
};


#endif // _MONGE_PATCH_INCLUDE


