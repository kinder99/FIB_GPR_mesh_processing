#include <iostream>
#include "MongePatch.h"


// Given a point P, its normal, and its closest neighbors (including itself) 
// compute a quadratic Monge patch that approximates the neighborhood of P.
// The resulting patch will be used to compute the principal curvatures of the 
// surface at point P.

void MongePatch::init(const glm::vec3 &P, const glm::vec3 &normal, const vector<glm::vec3> &closest)
{
	//create coordinate system from normal
	glm::vec3 w = -normal;
	glm::vec3 u = glm::cross(P,w);
	glm::vec3 v = glm::cross(w,u);

	//transform all neighbors to created system using dot products
	for(glm::vec3 pi : closest){
		float ui = glm::dot(u, (pi - P));
		float vi = glm::dot(v, (pi - P));
		float wi = glm::dot(w, (pi - P));

		glm::vec3 pi_t = glm::vec3(ui,vi,wi);
		closest_t.push_back(pi_t);
	}
}

// Return the values of the two principal curvatures for this patch

void MongePatch::principalCurvatures(float &kmin, float &kmax) const
{
	kmin = 0.f;
	kmax = 0.f;
}


