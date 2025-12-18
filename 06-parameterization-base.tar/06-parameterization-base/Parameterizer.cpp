#include <iostream>
#include <cstdlib>
#include <set>
#include <map>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "Parameterizer.h"
#include "timing.h"


using namespace std;


// This method should compute new texture coordinates for the input mesh
// using the harmonic coordinates approach with uniform weights.
// The 'TriangleMesh' class has a method 'getTexCoords' that may be used 
// to access and update its texture coordinates per vertex.

void Parameterizer::harmonicCoordinates(TriangleMesh *mesh)
{
	// First we need to determine the border edges of the input mesh
	// note that this is probably not working but at least i have smth lmao
	set<unsigned int, unsigned int> halfEdges;
	vector<unsigned int> triangles = mesh->getTriangles();
	for(unsigned int id : triangles){
		//since i'm not on cpp 20 or higher, can't use contains method :(
		if(auto search = halfEdges.find((id, mesh->next(id))); search != halfEdges.end()){
			halfEdges.erase((id, mesh->next(id)));
		}
		else{
			halfEdges.insert((id, mesh->next(id)));
		}
	}
	
	// Then, these edges need to be arranged into a single cycle, the border polyline
	// TODO
	for(set<unsigned int, unsigned int>::iterator it = halfEdges.begin(); it != halfEdges.end(); it++){
		//it*
	}

	// Each of the vertices on the border polyline will receive a texture coordinate
	// on the border of the parameter space ([0,0] -> [1,0] -> [1,1] -> [0,1]), using 
	// the chord length approach
	// TODO
	
	// We build an equation system to compute the harmonic coordinates of the 
	// interior vertices
	// TODO
	
	// Finally, we solve the system and assign the computed texture coordinates
	// to their corresponding vertices on the input mesh
	// TODO
}






