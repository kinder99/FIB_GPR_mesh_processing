#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "LaplacianSmoothing.h"


void LaplacianSmoothing::setMesh(TriangleMesh *newMesh)
{
	mesh = newMesh;
}

glm::vec3 LaplacianSmoothing::compute1DLaplacian(std::vector<glm::vec3>& vertices, const std::vector<unsigned int> &neighbors, const glm::vec3& Pi)
{
	float factor = 1.0f/neighbors.size();
	glm::vec3 sum = glm::vec3(0);
	for(int i = 0; i < neighbors.size(); i++){
		glm::vec3 Pj = vertices.at(neighbors.at(i));
		sum += Pj - Pi;
	}
	return sum*factor;
}

/* This method should apply nIterations iterations of the laplacian vector multiplied by lambda 
   to each of the vertices. */

void LaplacianSmoothing::iterativeLaplacian(int nIterations, float lambda)
{
	std::vector<glm::vec3> vertices = mesh->getVertices();
	std::vector<glm::vec3> smooth = vertices;
	//for all iterations
	for (int curIteration = 0; curIteration < nIterations; curIteration++)
	{
		//iterate over all the vertices
		for(int i = 0; i < vertices.size(); i++){
			//compute neighbors
			vector<unsigned int> neighbors;
			mesh->getNeighbors(i,neighbors);

			//update vertex positions
			smooth.at(i) = vertices.at(i) + lambda * compute1DLaplacian(vertices, neighbors, vertices.at(i));
		}
		//assign new positions
		vertices = smooth;
	}
	//set smoothed vertices into mesh
	for(int i = 0; i < smooth.size(); i++){
		mesh->getVertices().at(i) = smooth.at(i);
	}
}

/* This method should apply nIterations iterations of the bilaplacian operator using lambda 
   as a scaling factor. */

void LaplacianSmoothing::iterativeBilaplacian(int nIterations, float lambda)
{
}

/* This method should apply nIterations iterations of Taubin's operator using lambda 
   as a scaling factor, and computing the corresponding nu value. */

void LaplacianSmoothing::iterativeLambdaNu(int nIterations, float lambda)
{
}

/* This method should compute new vertices positions by making the laplacian zero, while
   maintaing the vertices marked as constraints fixed. */

void LaplacianSmoothing::globalLaplacian(const vector<bool> &constraints)
{
}

/* This method has to optimize the vertices' positions in the least squares sense, 
   so that the laplacian is close to zero and the vertices remain close to their 
   original locations. The constraintWeight parameter is used to control how close 
   the vertices have to be to their original positions. */

void LaplacianSmoothing::globalBilaplacian(const vector<bool> &constraints, float constraintWeight)
{
}








