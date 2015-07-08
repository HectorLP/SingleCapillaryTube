#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>

#include "DynamicSingleTube.h"

SingleCapillaryTube::SingleCapillaryTube()
{
	loadTubeGeometry();
	loadFluidProperties();
	loadPressureValues();
	loadPositionAndTime();
}

void SingleCapillaryTube::loadTubeGeometry()
{
	using std::ifstream;
	ifstream tubeGeometryFile;
	tubeGeometryFile.open("tubegeometry.dat");
	while (!tubeGeometryFile.is_open())
	{
		std::cout << "Error happened when trying to open the file on tube" << 
					"geometry.\n";
		system("pause");
		exit(1);
	}
	tubeGeometryFile >> geometry.radius >> geometry.length >> 
						geometry.angleToHorizontal;
	while (geometry.radius <= 0 || geometry.length <= 0 || \
			geometry.angleToHorizontal <= 0)
	{
		
	}
}