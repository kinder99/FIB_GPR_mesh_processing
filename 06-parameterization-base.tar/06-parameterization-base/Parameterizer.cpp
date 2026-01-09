#include <iostream>
#include <cstdlib>
#include <set>
#include <map>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "Parameterizer.h"
#include "timing.h"
#include "progressbar.hpp"


// This method should compute new texture coordinates for the input mesh
// using the harmonic coordinates approach with uniform weights.
// The 'TriangleMesh' class has a method 'getTexCoords' that may be used 
// to access and update its texture coordinates per vertex.
void Parameterizer::harmonicCoordinates(TriangleMesh *mesh)
{
	// First we need to determine the border edges of the input mesh
	std::set<std::pair<int, int>> halfEdges;
	std::vector<unsigned int> &triangles = mesh->getTriangles();
	size_t trianglesNum = triangles.size() / 3;
	size_t verticesNum = mesh->getVertices().size();
	progressbar bar(trianglesNum);

	//look through all triangles
	for (int i = 0; i < trianglesNum; i++)
	{
		//look for half edges in the current triangle, three edges per triangle
		for (int j = 0; j < 3; j++)
		{
			int first = triangles[i * 3 + j]; //we get the j edge in the i triangle
			int second = triangles[i * 3 + (j+1)%3]; //and the following one, taking care not to overflow
			std::pair<int, int> halfEdge(first, second);
			std::pair<int, int> opposed(second, first);

			//if the opposite half edge is in the set, remove it
			if (halfEdges.find(opposed) != halfEdges.end()){
				halfEdges.erase(opposed);
			} 
			else { //else we add the current hald edge to the set
				halfEdges.insert(halfEdge);
			}
		}
		//send progress to standard output
		bar.update();
	}
	
	// Then, these edges need to be arranged into a single cycle, the border polyline

	// Each of the vertices on the border polyline will receive a texture coordinate
	// on the border of the parameter space ([0,0] -> [1,0] -> [1,1] -> [0,1]), using 
	// the chord length approach
	
	
	// We build an equation system to compute the harmonic coordinates of the 
	// interior vertices
	// TODO
	
	// Finally, we solve the system and assign the computed texture coordinates
	// to their corresponding vertices on the input mesh
	// TODO
}
