#include <iostream>
#include "MongePatch.h"


// Given a point P, its normal, and its closest neighbors (including itself) 
// compute a quadratic Monge patch that approximates the neighborhood of P.
// The resulting patch will be used to compute the principal curvatures of the 
// surface at point P.
void MongePatch::init(const glm::vec3 &P, const glm::vec3 &normal, const vector<glm::vec3> &closest)
{
	//create coordinate system from normal
	glm::vec3 w = -glm::normalize(normal);
	glm::vec3 u = glm::normalize(glm::cross(glm::vec3(1.0,0.0,0.0),w));
	glm::vec3 v = glm::normalize(glm::cross(w,u));

	//transform all neighbors to created system using dot products
	closest_t.resize(closest.size());
	for(int i = 0; i < closest.size(); i++){
		double ui = glm::dot(u, (closest[i] - P));
		double vi = glm::dot(v, (closest[i] - P));
		double wi = glm::dot(w, (closest[i] - P));

		glm::vec3 pi_t = glm::vec3(ui,vi,wi);
		closest_t[i] = pi_t;
	}
	
	//As = b; with A = sum(qi*qi^T) and b = sum(wi,qi)
	Eigen::VectorXd Q(6);
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6,6);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(6);
	for(int i = 0; i < closest.size(); i++){
		Q(0) = closest_t[i].x * closest_t[i].x;
		Q(1) = closest_t[i].x * closest_t[i].y;
		Q(2) = closest_t[i].y * closest_t[i].y;
		Q(3) = closest_t[i].x;
		Q(4) = closest_t[i].y;
		Q(5) = 1;

		A += Q * Q.transpose();
		b += closest_t[i].z * Q;
	}
	
	Eigen::VectorXd s = A.ldlt().solve(b);

	//assign values of H according to computed s
	H(0,0) = (double)(2*s(0));
	H(0,1) = (double)s(1);
	H(1,0) = (double)s(1);
	H(1,1) = (double)2*s(2);
}

// Return the values of the two principal curvatures for this patch
void MongePatch::principalCurvatures(float &kmin, float &kmax) const
{
	//solve for eigenvalues of H
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
	if(solver.info() != Eigen::Success){
		std::cout << "eigensolver failed" << std::endl;
		kmin = 0.f;
		kmax = 0.f;
		return;
	}

	Eigen::VectorXd evalues = solver.eigenvalues();
	kmin = evalues(0);
	kmax = evalues(1);
}

