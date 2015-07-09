#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>

#include "DynamicSingleTube.h"
#include "ScantMethod.h"

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
	tubeGeometryFile >> Geometry.radius >> Geometry.length >> 
						Geometry.angleToHorizontal;
	while (Geometry.radius <= 0 || Geometry.length <= 0 || \
			Geometry.angleToHorizontal <= 0)
	{
		std::cout << "The values on tube geometry should be larger than 0.\n";
		system("pause");
		exit(1);
	}
	tubeGeometryFile.close();
}

void SingleCapillaryTube::loadFluidProperties()
{
	using std::ifstream;
	ifstream fluidPropertiesFile;
	fluidPropertiesFile.open("fluidproperties.dat");
	while (!fluidPropertiesFile.is_open())
	{
		std::cerr << "Error happened when trying to open the file on fluids.\n";
		system("pause");
		exit(1);
	}
	fluidPropertiesFile >> Fluids.densityWetting >> Fluids.densityNonWetting >> \
						Fluids.dynamicViscosityWetting >> \
						Fluids.dynamicViscosityNonWetting >> Fluids.contactAngle;
	while (Fluids.densityWetting <= 0 || Fluids.densityNonWetting <= 0 || \
			Fluids.dynamicViscosityWetting <=0 || Fluids.dynamicViscosityNonWetting <= 0 || \
			Fluids.contactAngle < 0)
	{
		std::cerr << "The values on fluid properties should be larger than 0.\n";
		system("pause");
		exit(1);
	}
	fluidPropertiesFile.close();
}

void SingleCapillaryTube::loadPressureValues()
{
	using std::ifstream;
	ifstream pressureValuesFile;
	pressureValuesFile.open("pressurevalue.dat");
	while (!pressureValuesFile.is_open())
	{
		std::cerr << "Error happened when trying to open the file on pressure boundary.\n";
		system("pause");
		exit(1);
	}
	pressureValuesFile >> leftPressure >> rightPressure;
	pressureValuesFile.close();
}

void SingleCapillaryTube::loadPositionAndTime()
{
	using std::ifstream;
	ifstream positionAndTimeFile;
	positionAndTimeFile.open("positiontime.dat");
	while (!positionAndTimeFile.is_open())
	{
		std::cerr << "Error happened when trying to open the file on pressure boundary.\n";
		system("pause");
		exit(1);
	}
	pressureValuesFile >> leftPressure >> rightPressure;
	pressureValuesFile.close();
}

double SingleCapillaryTube::calCoefficientA(const TubeGeometry &TG, const FluidProperties &FP)
{
	
}

double SingleCapillaryTube::