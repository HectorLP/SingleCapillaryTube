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

double SingleCapillaryTube::calLocationInterface()
{
	long numTimeSteps;
	numTimeSteps = long ((timeEndPoint - initialTime) / timeStep);
	double *interfaceLocation;
	interfaceLocation = new double [numTimeSteps + 1];
	interfaceLocation[0] = 0.0;
	double timePoint = 0.0;
	
	double tempValue0, tempValue1;
	
	double tempL0, tempL1, tempL;
	tempL0 = interfaceLocation[0];
	tempL1 = Geometry.length / 25.0;
	
	double valueForL0, valueForL1, valueL;
	if (typeFunction == "havingAngle")
	{
		for (int i = 1; i < (numTimeSteps + 1); ++i)
		{
			timePoint += timeStep;
			if (i > 1)
			{
				tempL0 = interfaceLocation[i - 1];
				tempL1 = interfaceLocation[i - 1] + (interfaceLocation[i -1] - \
						interfaceLocation[i - 2]) / 2.;
			}
			valueForL0 = calLocationFunctionWithAngle(Geometry, Fluids, 
							tempL0, initialLocation, timePoint, initialTime);
			valueForL1 = calLocationFunctionWithAngle(Geometry, Fluids,
							tempL1, initialLocation, timePoint, initialTime);
			interfaceLocation[i] = useScantMethod(tempL1, tempL0, valueForL1,
												valueForL0, 1.0e-7,  timePoint);
		}
	}
	else if (typeFunction == "noAngle")
	{
		
	}
}

double SingleCapillaryTube::calLocationFunctionWithAngle(const TubeGeometry &TG, 
							const FluidProperties &FP, double tempLocation, 
							double initialLocationValue, double tempTime, 
							const double initialTimePoint) 
														)
{
	auto coefficientA = calCoefficientA(TG, FP);
	auto coefficientB = calCoefficientB(TG, FP);
	auto coefficientC = calCoefficientC(TG, FP);
	auto coefficientD = calCoefficientD(TG, FP);
	auto tempValue = coefficientA * (tempLocation - initialLocationValue) + \
				(coefficientB - (coefficientA * coefficientD / coefficientC)) * \
				log(fabs(coefficientC * tempLocation + coefficientD)/ \
					fabs(coefficientC * tempLocation + coefficientD)) + \
					+ coefficientC * (tempTime - initialTimePoint);
	return    tempValue;
}

double SingleCapillaryTube::calLocationFunctionWithoutAngle(const TubeGeometry &TG,
								const FluidProperties &FP, double tempLocation,
								double initialLocationValue, double tempTime,
								const double initalTimePoint)
{
	auto coefficientB = calCoefficientB(TG, FP);
	auto coefficientA = calCoefficientA(TG, FP);
	auto coefficientD = calCoefficientD(TG, FP);
	auto tempValue = -coefficientB * (tempLocation - initialLocationValue) - \
					1./2. * coefficientA * (pow(tempLocation, 2.) - \
					pow(initialLocationValue, 2.)) - coefficientD * \
					(tempTime - initialTimePoint);
	return tempValue;
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

double SingleCapillaryTube::useScantMethod(double x1, double x0, double f0, 
										   double f1, const double specificErr,
										   const t)
{
	double tempL1, tempL0, tempL;
	double valueForL0, valueForL1, valueL;
	tempL1 = x1;
	tempL0 = x0;
	valueForL0 = f0;
	valueForL1 = f1;
	while (fabs(tempL0 - tempL1) > specificErr)
	{
		tempL = tempL1 - valueForL1 * (tempL1 - tempL0) / (valueForL1 - valueForL0);
		valueL = calLocationFunctionWithAngle(Geometry, Fluids,
							tempL, initialLocation, t, initialTime);
		tempL0 = tempL1;
		tempL1 = tempL;
		valueL0 = ValueL1;
		valueL1 = valueL;
	}
	return tempL1;
}