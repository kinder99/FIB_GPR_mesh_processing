#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "LaplacianSmoothing.h"
#include "progressbar.hpp"

void LaplacianSmoothing::setMesh(TriangleMesh *newMesh)
{
	mesh = newMesh;
}

glm::vec3 LaplacianSmoothing::computeLaplacian(std::vector<glm::vec3>& vertices, const std::vector<unsigned int> &neighbors, const glm::vec3& Pi)
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
	std::vector<glm::vec3> meshVertices = mesh->getVertices();
	std::vector<glm::vec3> smooth = meshVertices;
	//for all iterations
	for (int curIteration = 0; curIteration < nIterations; curIteration++)
	{
		//iterate over all the vertices
		for(int i = 0; i < meshVertices.size(); i++){
			//compute neighbors
			std::vector<unsigned int> neighbors;
			mesh->getNeighbors(i,neighbors);

			//update vertex positions
			smooth.at(i) = meshVertices.at(i) + lambda * computeLaplacian(meshVertices, neighbors, meshVertices.at(i));
		}
		//assign new positions
		meshVertices = smooth;
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
	std::vector<glm::vec3> meshVertices = mesh->getVertices();
	std::vector<glm::vec3> smooth = meshVertices;
	//for all iterations
	for (int curIteration = 0; curIteration < nIterations; curIteration++)
	{
		//iterate over all the vertices
		for(int i = 0; i < meshVertices.size(); i++){
			//compute neighbors
			std::vector<unsigned int> neighbors;
			mesh->getNeighbors(i,neighbors);

			//update vertex positions
			glm::vec3 temp;
			temp = meshVertices.at(i) + lambda * computeLaplacian(meshVertices, neighbors, meshVertices.at(i));
			//then correct using opposite sign
			smooth.at(i) = temp - lambda * computeLaplacian(meshVertices, neighbors, temp);
		}
		//assign new positions
		meshVertices = smooth;
	}
	//set smoothed vertices into mesh
	for(int i = 0; i < smooth.size(); i++){
		mesh->getVertices().at(i) = smooth.at(i);
	}
}

/* This method should apply nIterations iterations of Taubin's operator using lambda 
   as a scaling factor, and computing the corresponding nu value. */

void LaplacianSmoothing::iterativeLambdaNu(int nIterations, float lambda)
{
	std::vector<glm::vec3> meshVertices = mesh->getVertices();
	std::vector<glm::vec3> smooth = meshVertices;
	//compute nu using given formula
	float nu = lambda/(0.1*lambda-1);
	//for all iterations
	for (int curIteration = 0; curIteration < nIterations; curIteration++)
	{
		//iterate over all the vertices
		for(int i = 0; i < meshVertices.size(); i++){
			//compute neighbors
			std::vector<unsigned int> neighbors;
			mesh->getNeighbors(i,neighbors);

			//update vertex positions
			glm::vec3 temp;
			temp = meshVertices.at(i) + lambda * computeLaplacian(meshVertices, neighbors, meshVertices.at(i));
			//then correct using opposite sign
			smooth.at(i) = temp + nu * computeLaplacian(meshVertices, neighbors, temp);
		}
		//assign new positions
		meshVertices = smooth;
	}
	//set smoothed vertices into mesh
	for(int i = 0; i < smooth.size(); i++){
		mesh->getVertices().at(i) = smooth.at(i);
	}
}

/* This method should compute new vertices positions by making the laplacian zero, while
   maintaing the vertices marked as constraints fixed. */
void LaplacianSmoothing::globalLaplacian(const vector<bool> &constraints)
{
	//FYI this takes quite a bit of time depending on the mesh loaded
	std::vector<glm::vec3> &meshVertices = mesh->getVertices();
	size_t size = meshVertices.size();
	
	//get matrix containing all vertices as rows
	Eigen::MatrixXf P(size,3);
	for(int i = 0; i < size; i++){
		for(int j = 0; j < 3; j++){
			P(i,j) = meshVertices[i][j];
		}
	}

	Eigen::SparseMatrix<float> C(size, size);
	Eigen::SparseMatrix<float> M(size, size);
	Eigen::SparseMatrix<float> L(size, size);

	//fill C and M matrices
	std::vector<unsigned int> neighbors;
	progressbar bar(size); //progress bar by Luigi Perdolti
	for(int i = 0; i < size; i++){
		mesh->getNeighbors(i, neighbors);
		float neighborNum = (float)neighbors.size();

		C.coeffRef(i,i) = -neighborNum;
		M.coeffRef(i,i) = 1/neighborNum;

		for(int j = 0; j < neighborNum; j++){
			C.coeffRef(i, neighbors[j]) = 1;
		}

		//update progress
		bar.update();
	}

	// //compute Laplacian matrix L
	L = M * C;

	//modify vertices according to constraint condition and computed Laplacian, with formula P' = P + LP
	P += L * P;
	for(int i = 0; i < size; i++){
		if(constraints[i]){
			for(int j = 0; j < 3; j++){
				meshVertices[i][j] = P(i,j);
			}
		}
	}
}

/* This method has to optimize the vertices' positions in the least squares sense, 
   so that the laplacian is close to zero and the vertices remain close to their 
   original locations. The constraintWeight parameter is used to control how close 
   the vertices have to be to their original positions. */
void LaplacianSmoothing::globalBilaplacian(const vector<bool> &constraints, float constraintWeight)
{
	std::vector<glm::vec3> &meshVertices = mesh->getVertices();
	size_t size = meshVertices.size();
	
	//get matrix containing all vertices as rows
	Eigen::MatrixXf P(size,3);
	for(int i = 0; i < size; i++){
		for(int j = 0; j < 3; j++){
			P(i,j) = meshVertices[i][j];
		}
	}

	Eigen::SparseMatrix<float> C(size, size);
	Eigen::SparseMatrix<float> M(size, size);
	Eigen::SparseMatrix<float> L(size, size);

	//fill C and M matrices
	std::vector<unsigned int> neighbors;
	progressbar bar(size);
	for(int i = 0; i < size; i++){
		mesh->getNeighbors(i, neighbors);
		float neighborNum = (float)neighbors.size();

		C.coeffRef(i,i) = -neighborNum;
		M.coeffRef(i,i) = 1/neighborNum;

		for(int j = 0; j < neighborNum; j++){
			C.coeffRef(i, neighbors[j]) = 1;
		}

		//update progress bar
		bar.update();
	}

	//compute Laplacian matrix L
	L = M * C;

	//modify vertices according to constraint condition and weight, and computed Laplacian 
	P += constraintWeight * L * P;
	for(int i = 0; i < size; i++){
		if(constraints[i]){
			for(int j = 0; j < 3; j++){
				meshVertices[i][j] = P(i,j);
			}
		}
	}
}








