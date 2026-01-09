#include <iostream>
#include "RBFFunction.h"
#include <math.h>

const float ARTIFICIAL_DISTANCE = 0.01;
const float LAMBDA = 1.00;

//helper function to compute the Gaussian RBF, two definitions depending on the use case
double RBFFunction::computeGaussian(glm::vec3 pos1, glm::vec3 pos2, float c) const
{
	double r = glm::distance(pos1, pos2);
	//check for the radius condition 
	if(r >= 3*c)
		return 0;
	return std::exp(-(r*r)/(double)m_c_squared);
}

double RBFFunction::computeGaussian(float r) const
{
	return std::exp(-(r*r)/(double)m_c_squared);
}

/* Initialize everything to be able to compute the implicit distance to the reconstructed
   point cloud at arbitrary points that are close enough to the point cloud. As should be
   obvious by the name of the class, the distance has to be computed using RBFs.
 */

void RBFFunction::init(const PointCloud *pointCloud, float standardDeviation, float supportRadius)
{
	//initialize all relevant attributes of the class
	m_pointCloud = pointCloud;
	m_positions = &m_pointCloud->getPoints();
	m_normals = &m_pointCloud->getNormals();
	m_std = standardDeviation;
	m_radius = supportRadius;
	m_c.resize(3*m_positions->size());
	m_c_squared = 2*m_std*m_std;
	Eigen::VectorXd v(3*m_positions->size());
	m_positionsMatrix.resize(3*m_positions->size());
	
	//compute the v vector by adding artificial points according to the definition
	for(int i = 0; i <m_positions->size(); ++i)
	{
		glm::vec3 vec = m_positions->at(i);
		m_positionsMatrix[i] = vec;
		v[i] = 0;
		//set the directions
		m_positionsMatrix[i+m_positions->size()] = vec + glm::normalize(m_normals->at(i))*ARTIFICIAL_DISTANCE;
		m_positionsMatrix[i+m_positions->size()*2] = vec - glm::normalize(m_normals->at(i))*ARTIFICIAL_DISTANCE;
		//set the amplitudes
		v[i+m_positions->size()] = ARTIFICIAL_DISTANCE;
		v[i+m_positions->size()*2] = -ARTIFICIAL_DISTANCE;
	}

	//initialize relevant stuff to prepare for A matrix computation
	m_neighbors.setPoints(&m_positionsMatrix);
	Eigen::SparseMatrix<double> A(m_positionsMatrix.size(), m_positionsMatrix.size());
	//we use a vector of triplets which will be computed individually, then used to compute A
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(m_positionsMatrix.size());

    //compute A matrix using neighbors in the support radius for all points
	for(int i=0; i<m_positionsMatrix.size();++i)
	{
		//get neighbors in radius
		std::vector<std::pair<size_t, float>> current_neighborsIndex;
		m_neighbors.getNeighborsInRadius(m_positionsMatrix[i], m_radius, current_neighborsIndex);

		//for all neighbors
		for(int j=0; j<current_neighborsIndex.size();++j)
		{
			int indexNeigh = current_neighborsIndex[j].first;
			//compute rbf value for the current point
            double RBF = computeGaussian(m_positionsMatrix[i], m_positionsMatrix[j], standardDeviation);
            if(i != indexNeigh && RBF != 0.0)
            {
                double rbfValue = RBF;
                tripletList.emplace_back(i,indexNeigh,rbfValue);
            }
		}
        tripletList.emplace_back(i,i,1.0f + LAMBDA);
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	A.makeCompressed();

	//solve A and check status to know if the computation failed or not
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>solver;
	solver.compute(A);
	if(solver.info() != Eigen::Success)
	{
		std::cout << "decomposition failed" << std::endl;
		return;
	}
	m_c.setZero();
	m_c = solver.solve(v);
	if(solver.info() != Eigen::Success)
	{
		std::cout << "solving failed" << std::endl;
		return;
	}
	std::cout << "RBF initialization complete" << std::endl;
}


/* This operator returns a boolean that if true signals that the value parameter
   has been modified to contain the value of the RBF implicit distance at point P.
 */

bool RBFFunction::operator()(const glm::vec3 &P, float &value) const
{
	std::vector<std::pair<size_t, float>> current_neighborsIndex;
	m_neighbors.getNeighborsInRadius(P, m_radius, current_neighborsIndex);
	
	value = 0;
	for(int i=0; i < current_neighborsIndex.size();++i)
	{
		int indexNeigh = current_neighborsIndex[i].first;
		if(indexNeigh < m_positionsMatrix.size())
		{
            double rbf = computeGaussian(P, m_positionsMatrix[indexNeigh], m_std);
            value += rbf*m_c(indexNeigh);
		}
	}
	return true;
}

