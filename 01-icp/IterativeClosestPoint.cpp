#include <iostream>
#include <algorithm>
#include "IterativeClosestPoint.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <cmath>
#include <glm/gtc/matrix_transform.hpp>


void IterativeClosestPoint::setClouds(PointCloud *pointCloud1, PointCloud *pointCloud2)
{
	cloud1 = pointCloud1;
	cloud2 = pointCloud2;
	knn.setPoints(&(cloud1->getPoints()));
}

// This method should mark the border points in cloud 1. It also changes their color (for example to red).
// You will need to add an attribute to this class that stores this property for all points in cloud 1. 
// maxDeltaAlpha is supposed to be between 90 and 180 degrees apparently
// compute border points and store them locally, plus change colors of border points using PointCloud::getColors()
void IterativeClosestPoint::markBorderPoints()
{
	// TODO
	std::cout << "marking border points" << std::endl;
	std::vector<size_t> neighbors;
	std::vector<float> dists;
	std::vector<glm::vec3>& points = (cloud1->getPoints());
	//initialize values of border vector
	for (int i = 0; i < points.size(); i++){
		borderPoints.push_back(false);
	}
	int counter = 0;
	for(glm::vec3 point : points){
		neighbors.clear();
		dists.clear();
		//get nearest neighbors of current point	
		knn.getKNearestNeighbors(point,10,neighbors,dists);

		//PCA
		glm::vec3 sum = glm::vec3(0,0,0);
		for (int i = 0; i<neighbors.size(); i++){
			sum += points.at(neighbors.at(i));
		}
		auto const count = static_cast<float>(neighbors.size());
		glm::vec3 p_centroid = sum/count; //compute centroid of current point subset
		Eigen::Matrix3f covariance_matrix = Eigen::Matrix3f::Zero();
		for (int i = 0; i<neighbors.size(); i++){ //for all neighbors of current point
			glm::vec3 p = points.at(neighbors.at(i));
			p = p - p_centroid; //translate all points by the centroid coordinates 
			Eigen::Vector3f p_eigen = Eigen::Vector3f(p.x,p.y,p.z);
			covariance_matrix += p_eigen * p_eigen.transpose(); //compute covariance matrix
		}
		//compute spectral decomposition
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(covariance_matrix);
		Eigen::Matrix3f eigenvecs = eigensolver.eigenvectors();

		//transform the neighbors to the PCA computed local frame, then discard z
		std::vector<float> neighbor_angles;
		for (int i = 1; i < neighbors.size(); i++){
			glm::vec3 current_neighbor = points.at(neighbors.at(i));
			Eigen::Vector3f substract = Eigen::Vector3f(current_neighbor.x - point.x, current_neighbor.y - point.y, current_neighbor.z - point.z);
			current_neighbor.x = eigenvecs.col(2).dot(substract);
			current_neighbor.y = eigenvecs.col(1).dot(substract);
			current_neighbor.z = 0;//discard Z component to project on xy plane
			float angle = atan2(current_neighbor.y, current_neighbor.x) * (180/M_PI); //get angle value in degrees
			neighbor_angles.push_back(angle);
		}

		neighbors.erase(neighbors.begin() + 0); //erase current point from neighbors

		//sort angles then look for maxDeltaAlpha (the largest)
		std::sort(neighbor_angles.begin(), neighbor_angles.end());
		float maxDeltaAlpha = 0;
		float difference = 360 - (neighbor_angles.back() - neighbor_angles.front());
		difference > maxDeltaAlpha ? maxDeltaAlpha = difference : maxDeltaAlpha = maxDeltaAlpha;
		for(int i = 0; i < neighbor_angles.size()-1; i++){
			difference = neighbor_angles.at(i+1) - neighbor_angles.at(i);
			difference > maxDeltaAlpha ? maxDeltaAlpha = difference : maxDeltaAlpha = maxDeltaAlpha; // accumulate maxDeltaAlpha difference value
		}
		maxDeltaAlpha > 90 ? borderPoints.at(counter) = true : borderPoints.at(counter) = false; //check value of maxDeltaAlpha angle to set border point status
		counter++;
	}
	for(int i = 0; i < points.size(); i++){
		if(borderPoints.at(i) == true){
			(cloud1->getColors()).at(i) = glm::vec4(1,0,0,1); //set color as red if point is a border point
		}
	}
}



// This method should compute the closest point in cloud 1 for all non border points in cloud 2. 
// This correspondence will be useful to compute the ICP step matrix that will get cloud 2 closer to cloud 1.
// Store the correspondence in this class as the following method is going to need it.
// As it is evident in its signature this method also returns the correspondence. The application draws this if available.
vector<int> *IterativeClosestPoint::computeCorrespondence()
{
    // reset the correspondance vector
    correspondance.resize(this->cloud2->size(), -1);
    this->numberOfCorrespondance = 0;

    std::vector<size_t> neighbors;
    std::vector<float> distances;
    for(int i = 0 ; i < this->cloud2->size() ; ++i) {
        glm::vec3 q = this->cloud2->getPoints()[i];
        // find the nearest neighbor
        this->knn.getKNearestNeighbors(q, 1, neighbors, distances);
        if(this->borderPoints[neighbors[0]]) continue;

        correspondance[i] = (int)neighbors[0];
        this->numberOfCorrespondance++;
    }
	return &(this->correspondance);
}



// This method should compute the closest point in cloud 1 for all non border points in cloud 2. 
// This correspondence will be useful to compute the ICP step matrix that will get cloud 2 closer to cloud 1.
// Store the correspondence in this class as the following method is going to need it.
// As it is evident in its signature this method also returns the correspondence. The application draws this if available.
// store correspondances into vector of int, if border point assign value -1
// vector<int> *IterativeClosestPoint::computeCorrespondence()
// {
// 	correspondance.resize(this->cloud2->size(), -1);
// 	this->numberOfCorrespondance = 0;

// 	std::vector<size_t> neighbors;
// 	std::vector<float> dists;

// 	for(int i = 0; i<this->cloud2->size();++i){
// 		glm::vec3 q = this->cloud2->getPoints()[i];
// 		this->knn.getKNearestNeighbors(q, 1, neighbors, dists);
// 		if(this->borderPoints[neighbors[0]]){
// 			continue;
// 		}
// 		correspondance[i] = (int)neighbors[0];
// 		this->numberOfCorrespondance;
// 	}
// 	return &(this->correspondance);
// }


// This method should compute the rotation and translation of an ICP step from the correspondence
// information between clouds 1 and 2. Both should be encoded in the returned 4x4 matrix.
// To do this use the SVD algorithm in Eigen.
glm::mat4 IterativeClosestPoint::computeICPStep()
{
    // 1. Compute centroids
    glm::vec3 pCentroid = glm::vec3(0.f, 0.f, 0.f);
    glm::vec3 qCentroid = glm::vec3(0.f, 0.f, 0.f);
    for(int i = 0 ; i < this->cloud2->size() ; ++i) {
        if(this->correspondance[i] == -1) continue;

        pCentroid += this->cloud1->getPoints()[this->correspondance[i]];
        qCentroid += this->cloud2->getPoints()[i];
    }
    pCentroid = pCentroid / (float)this->numberOfCorrespondance;
    qCentroid = qCentroid / (float)this->numberOfCorrespondance;

    Eigen::Matrix3f covarianceMatrix = Eigen::Matrix3f::Zero();
    for(int i = 0 ; i < this->cloud2->size() ; ++i) {
        if(this->correspondance[i] == -1) continue;

        // 1.2. shift points
        glm::vec3 shiftedP = this->cloud1->getPoints()[this->correspondance[i]] - pCentroid;
        glm::vec3 shiftedQ = this->cloud2->getPoints()[i] - qCentroid;

        // 2. compute covariance matrix
        Eigen::Vector3f eigenPi(shiftedP.x, shiftedP.y, shiftedP.z);
        Eigen::Vector3f eigenQi(shiftedQ.x, shiftedQ.y, shiftedQ.z);
        covarianceMatrix += eigenQi * eigenPi.transpose();
    }

    // 3. Compute SVD decomposition
    Eigen::JacobiSVD<Eigen::Matrix3f> solver;
    solver.compute(covarianceMatrix, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Matrix3f U(solver.matrixU());
    Eigen::Matrix3f Sigma(solver.singularValues().asDiagonal());
    Eigen::Matrix3f V(solver.matrixV());
    Eigen::Matrix3f S =  U * Sigma * V.transpose();

    // 4. Compute the optimal rotation
    Eigen::Matrix3f R = V * U.transpose();

    // 5. If we have a reflexion (det S = -1), we change the sign of the last element of R
    if(S.determinant() < 0.f) { // float equal comparison is bad
        R.col(2).row(2) *= -1;
    }

    // 6. Compute translation vector
    Eigen::Vector3f eigenPCentroid(pCentroid.x, pCentroid.y, pCentroid.z);
    Eigen::Vector3f eigenQCentroid(qCentroid.x, qCentroid.y, qCentroid.z);

    Eigen::Vector3f t = eigenPCentroid - R * eigenQCentroid;

    // glm is column major
    glm::mat4 res{
        R(0,0), R(1,0), R(2,0), 0.f,
        R(0,1), R(1,1), R(2,1), 0.f,
        R(0,2), R(1,2), R(2,2), 0.f,
        t(0), t(1), t(2), 1.0f
    };

    /*
    for (int row = 0; row < 4; ++row)
    {
        for (int col = 0; col < 4; ++col)
        {
            std::cout << res[col][row] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    */

	return res;
}

// This method should perform the whole ICP algorithm with as many steps as needed.
// It should stop when maxSteps are performed, when the Frobenius norm of the transformation matrix of
// a step is smaller than a small threshold, or when the correspondence does not change from the 
// previous step.

vector<int> *IterativeClosestPoint::computeFullICP(unsigned int maxSteps)
{
	int step = 0;

    // Update correspondances
    (void)computeCorrespondence();

    std::vector<int> previousCorrespondance = this->correspondance;
    while(step < maxSteps) {
        step++;



        // Compute ICP Step
        glm::mat4 ICPStep(computeICPStep());

        // Apply transformation
        this->cloud2->transform(ICPStep);

        // early stop: if correspondance didn't change a bit
        (void)computeCorrespondence();
        if(previousCorrespondance == this->correspondance) break;
        previousCorrespondance = this->correspondance;
    }

    std::cout << "Number of steps: " << step << std::endl;

	return &(this->correspondance);
}