#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>

#include "DynamicSingleTube.h"
#include "ScantMethod.h"

#define PI 3.14159265
#define G 9.80

SingleCapillaryTube::SingleCapillaryTube()
{
	loadTubeGeometry();
	loadFluidProperties();
	loadPressureValues();
	loadPositionAndTime();
}

double calLocationInterface(const TubeGeometry &TG, const FluidProperties &FP)
{
	
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
	pressureValuesFile >> initialLocation >> initialTime >> timeStep >> timeEndPoint;
	pressureValuesFile.close();
}

double SingleCapillaryTube::calLocationFunctionWithAngle(const TubeGeometry &TG, \
							const FluidProperties &FP, double tempLocation, \
							double initialLocationValue)
{
	auto coefficientA = calCoefficientA(TG, FP);
	auto coefficientB = calCoefficientB(TG, FP);
	auto coefficientC = calCoefficientC(TG, FP);
	auto coefficientD = calCoefficientD(TG, FP);
	tempValue = coefficientA * (tempLocation - initialLocationValue) + \
				(coefficientB - ())
}

double SingleCapillaryTube::calCoefficientA(const TubeGeometry &TG, const FluidProperties &FP)
{
	double coefficientA; 
	coefficientA = 8.0 * (FP.dynamicViscosityWetting - FP.dynamicViscosityNonWetting) / \
				(TG.radius * TG.radius);
	return coefficientA;
}

double SingleCapillaryTube::calCoefficientB(const TubeGeometry &TG, const FluidProperties &FP)
{
	double coefficientB;
	coefficientB = 8.0 * (FP.dynamicViscosityNonWetting * TG.length) / \
					pow(TG.radius, 2.);
	return coefficientB;
}

double SingleCapillaryTube::calCoefficientC(const TubeGeometry &TG, const FluidProperties &FP)
{
	double coefficientC;
	coefficientC = (FP.densityWetting - FP.densityNonWetting) * G * \
					sin(TG.angleToHorizontal * PI / 180.);
	return coefficientC;
}

double SingleCapillaryTube::calCoefficientD(const TubeGeometry &TG, const FluidProperties &FP)
{
	double coefficientD;
	double tempCapillary;
	tempCapillary = calCapillaryPressure(TG, FP);
	coefficientD = FP.densityNonWetting * TG.length * G * \
				sin(TG.angleToHorizontal * PI / 180.) + \
				(rightPressure - leftPressure) - tempCapillary;
	return coefficientD;
}

double SingleCapillaryTube::calCapillaryPressure(const TubeGeometry &TG, \
							const FluidProperties &FP)
{
	//double capillaryPressure;
	auto capillaryPressure = 2.0 * FP.surfaceTension * cos(FP.contactAngle * PI / \
							180.) / TG.radius;
	return capillaryPressure;
}
