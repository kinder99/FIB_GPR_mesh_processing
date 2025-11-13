#include "NormalEstimator.h"
#include "NearestNeighbors.h"

// This method has to compute a normal per point in the 'points' vector and put it in the 
// 'normals' vector. The 'normals' vector already has the same size as the 'points' vector. 
// There is no need to push_back new elements, only to overwrite them ('normals[i] = ...;')
// In order to compute the normals you should use PCA. The provided 'NearestNeighbors' class
// wraps the nanoflann library that computes K-nearest neighbors effciently. 

void NormalEstimator::computePointCloudNormals(const vector<glm::vec3> &points, vector<glm::vec3> &normals)
{
	//initialize useful stuff
	NearestNeighbors nn = NearestNeighbors();
	nn.setPoints(&points);
	glm::vec3 sum = glm::vec3(0,0,0);
	std::vector<size_t> neighbors;
	std::vector<float> dists;
	int c = 0;

	for(glm::vec3 p : points){
		sum = glm::vec3(0,0,0);
		neighbors.clear();
		dists.clear();
		c++;
		nn.getKNearestNeighbors(p,10,neighbors,dists); //get neighbors of current point
		for (int i = 0; i<neighbors.size(); i++){
			sum += points.at(neighbors.at(i));
		}
		auto const count = static_cast<float>(neighbors.size());
		glm::vec3 p_centroid = sum/count; //compute centroid of current point subset

		Eigen::Matrix3f covariance_matrix = Eigen::Matrix3f::Zero();
		for (int i = 0; i<neighbors.size(); i++){ //for all neighbors of current point
			glm::vec3 point = points.at(neighbors.at(i));
			point = point - p_centroid; //translate all points by the centroid coordinates 
			Eigen::Vector3f p_eigen = Eigen::Vector3f(point.x,point.y,point.z);
			covariance_matrix += p_eigen * p_eigen.transpose(); //compute covariance matrix
		}

		//compute spectral decomposition
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(covariance_matrix);
		Eigen::Matrix3f eigenvecs = eigensolver.eigenvectors();

		//define normal
		glm::vec3 temp = glm::vec3(eigenvecs.col(0)[0], eigenvecs.col(0)[1], eigenvecs.col(0)[2]);
		temp = glm::normalize(temp);

		//check for orientation
		if(temp.z < 0){temp = -temp;};

		//assignment
		normals[c] = temp;
	}	
}