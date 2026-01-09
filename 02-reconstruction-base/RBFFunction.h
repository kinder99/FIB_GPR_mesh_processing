#ifndef _RBF_FUNCTION_INCLUDE
#define _RBF_FUNCTION_INCLUDE

#include "ImplicitFunction.h"
#include "PointCloud.h"
#include "NearestNeighbors.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseQR>

class RBFFunction : public ImplicitFunction
{

public:
	void init(const PointCloud *pointCloud, float standardDeviation, float supportRadius);

	bool operator()(const glm::vec3 &P, float &value) const;
	
private:
	double computeGaussian(glm::vec3 pos1, glm::vec3 pos2, float c) const;
	double computeGaussian(float r) const;
	Eigen::VectorXd m_c;
	const PointCloud* m_pointCloud;
	const vector<glm::vec3> *m_positions;
	const vector<glm::vec3> *m_normals;
	std::vector<glm::vec3> m_positionsMatrix;
	NearestNeighbors m_neighbors;
	float m_std;
	float m_radius;
	float m_c_squared;
};

#endif // _RBF_FUNCTION_INCLUDE