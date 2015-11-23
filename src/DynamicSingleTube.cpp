#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <functional>
#include <utility>

#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint.hpp>

#include "DynamicSingleTube.h"
#include "tolerance.h"
//#include "tolerance.h"
//#include "ScantMethod.h"

#define PI 3.14159265
#define G 9.80

using namespace boost::numeric::odeint;

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
			Geometry.angleToHorizontal < 0)
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
						Fluids.dynamicViscosityNonWetting >> \
						Fluids.surfaceTension >> Fluids.contactAngle;
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
	positionAndTimeFile >> initialLocation >> initialTime >> timeStep >> timeEndPoint >> typeFunction;
	positionAndTimeFile.close();
}

void SingleCapillaryTube::calLocationInterfaceScant()
{
	long numTimeSteps;
	numTimeSteps = long ((timeEndPoint - initialTime) / timeStep);
	std::vector<double> interfaceLocation;
	//interfaceLocation = new double [numTimeSteps + 1];
	interfaceLocation.push_back(0.0);
	double timePoint = 0.0;

	double tempValue0, tempValue1;

	double tempL0, tempL1, tempL;
	tempL0 = interfaceLocation[0];
	tempL1 = Geometry.length / 25.0;
	double valueForL0, valueForL1, valueL;
	int i = 1;
	if (typeFunction == "havingAngle")
	{
		while (i < (numTimeSteps + 1) && fabs(tempL0) < Geometry.length)
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
// 			interfaceLocation[i] = useScantMethodWithAngle(tempL1, tempL0, valueForL1,
// 												valueForL0, 1.0e-7,  timePoint);
			double tempScant = useScantMethodWithAngle(tempL1, tempL0, valueForL1,
												valueForL0, 1.0e-7,  timePoint);
			interfaceLocation.push_back(tempScant);
			i += 1;
		}
	}
	else if (typeFunction == "noAngle")
	{
		while (i < (numTimeSteps + 1) && fabs(tempL0) < Geometry.length)
		{
			timePoint += timeStep;
			if (i > 1)
			{
				tempL0 = interfaceLocation[i - 1];
				tempL1 = interfaceLocation[i - 1] + (interfaceLocation[i -1] - \
						interfaceLocation[i - 2]) / 2.;
			}
			valueForL0 = calLocationFunctionWithoutAngle(Geometry, Fluids,
							tempL0, initialLocation, timePoint, initialTime);
			valueForL1 = calLocationFunctionWithoutAngle(Geometry, Fluids,
							tempL1, initialLocation, timePoint, initialTime);
// 			interfaceLocation[i] = useScantMethodWithoutAngle(tempL1, tempL0, valueForL1,
// 												valueForL0, 1.0e-7,  timePoint);
			double tempScant = useScantMethodWithoutAngle(tempL1, tempL0, valueForL1,
												valueForL0, 1.0e-7,  timePoint);
			interfaceLocation.push_back(tempScant);
			i += 1;
		}
	}
	outputInterfaceLocation(interfaceLocation);
}

double SingleCapillaryTube::calLocationFunctionWithAngle(const TubeGeometry &TG,
							const FluidProperties &FP, double tempLocation,
							double initialLocationValue, double tempTime,
							const double initialTimePoint)
{
	auto coefficientA = calCoefficientA();
	auto coefficientB = calCoefficientB();
	auto coefficientC = calCoefficientC();
	auto coefficientD = calCoefficientD();
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
								const double initialTimePoint)
{
	auto coefficientB = calCoefficientB();
	auto coefficientA = calCoefficientA();
	auto coefficientD = calCoefficientD();
	auto tempValue = -coefficientB * (tempLocation - initialLocationValue) - \
					1./2. * coefficientA * (pow(tempLocation, 2.) - \
					pow(initialLocationValue, 2.)) - coefficientD * \
					(tempTime - initialTimePoint);
	return tempValue;
}

void SingleCapillaryTube::calLocationInterfaceBisect()
{
	//TODO use boost::math::tools::bisect(F, min, max, iteration_times)\
	to solve the interface location
	using namespace std::placeholders;
	std::vector<double> interfaceLocation;
	long numTimeSteps;
	numTimeSteps = long ((timeEndPoint - initialTime) / timeStep);
	//tolerance Tol = 1.0e-7; 
	interfaceLocation.push_back(0.0);
	double timePoint = 0.0;
	
	double tempValue0, tempValue1;
	tolerance Tol = 1.0e-7;
	double tempL0, tempL1, tempL;
	tempL0 = interfaceLocation[0];
	tempL1 = Geometry.length / 25.0;
	double valueForL0, valueForL1, valueL;
	int i = 1;	
	if (typeFunction == "havingAngle")
	{
		while (i < (numTimeSteps + 1) && fabs(tempL) < Geometry.length)
		{
			timePoint += timeStep;
			if (i > 1)
			{
				tempL0 = interfaceLocation[i];
				tempL1 = interfaceLocation[i - 1] + (interfaceLocation[i -1] - \
				interfaceLocation[i - 2]) / 2.;
			}
			typedef std::pair<double, double> resultEachStep;
			resultEachStep stepResult;
			stepResult = boost::math::tools::bisect(std::bind(&SingleCapillaryTube::calLocationFunctionWithAngle,this,\
											Geometry, Fluids, _1, initialLocation, \
											timePoint, initialTime), 0.0, \
											-Geometry.length, Tol);
			if (stepResult.first == stepResult.second)
			{
				tempL = stepResult.first;
			}
			interfaceLocation.push_back(tempL);
			i += 1;
		}
	}
	else if (typeFunction == "noAngle")
	{
		while (i < (numTimeSteps + 1) && fabs(tempL) < Geometry.length)
		{
			timePoint += timeStep;
			if (i > 1)
			{
				tempL0 = interfaceLocation[i];
				tempL1 = interfaceLocation[i - 1] + (interfaceLocation[i -1] - \
				interfaceLocation[i - 2]) / 2.;
			}
			typedef std::pair<double, double> resultEachStep;
			resultEachStep stepResult;
			stepResult = boost::math::tools::bisect(std::bind(&SingleCapillaryTube::calLocationFunctionWithoutAngle,this,\
						Geometry, Fluids, _1, initialLocation, timePoint, initialTime), \
						0.0, -Geometry.length, Tol);
			if (stepResult.first == stepResult.second)
			{
				tempL = stepResult.first;
			}
			interfaceLocation.push_back(tempL);
			i += 1;
		}
	}
	outputInterfaceLocation(interfaceLocation);
}

void SingleCapillaryTube::calThreeCoefficients()
{
	CoefficientB = calCoefficientB();
	CoefficientA = calCoefficientA();
	CoefficientD = calCoefficientD();
	std::cout << "The coefficient value is " << CoefficientA << CoefficientB << CoefficientD << "\n";
	std::cin.ignore();
}

void SingleCapillaryTube::calInterfaceLocationSolveOde()
{
	double finalLocation;
	finalLocation = useOdeSolver();
	std::cout << "The final position is " << finalLocation << std::endl;
}

void SingleCapillaryTube::rhsPart(const double location, double& dldt, const double t)
{
	calThreeCoefficients();
	dldt = - (CoefficientD) / (CoefficientA * location + CoefficientB);
}



inline void observeProcess(const double &x, const double t)
{
	std::cout << t << "\t" << x << std::endl;
}

double SingleCapillaryTube::useOdeSolver()
{
	typedef runge_kutta_dopri5<double, double, double, double, vector_space_algebra> stepperType;
	namespace pl = std::placeholders;
	double tempLocation;
	tempLocation = initialLocation;
	int numTimeSteps;
	numTimeSteps = (timeEndPoint - initialTime) / timeStep;
	std:: cout << "Time step is " << timeStep << std::endl;
	std::cin.ignore();
	double tempInitialTime = initialTime;
	int tempStep = 0;
	while (tempStep < numTimeSteps and fabs(tempLocation) < Geometry.length)
	{
		tempStep += 1;
		double t1 = initialTime + tempStep * timeStep;
		tempInitialTime += (tempStep - 1) * timeStep;
		integrate_adaptive(make_controlled(1e-12, 1e-12, stepperType() ), std::bind(&SingleCapillaryTube::rhsPart,*this, pl::_1 , pl::_2 , pl::_3), \
						tempLocation, tempInitialTime, t1, timeStep, observeProcess);
	}
	return tempLocation;
}


double SingleCapillaryTube::calCoefficientA()
{
	double coefficientA;
	coefficientA = 8.0 * (Fluids.dynamicViscosityWetting - Fluids.dynamicViscosityNonWetting) / \
				(Geometry.radius * Geometry.radius);
	return coefficientA;
}

double SingleCapillaryTube::calCoefficientB()
{
	double coefficientB;
	coefficientB = 8.0 * (Fluids.dynamicViscosityNonWetting * Geometry.length) / \
					pow(Geometry.radius, 2.);
	return coefficientB;
}

double SingleCapillaryTube::calCoefficientC()
{
	double coefficientC;
	coefficientC = (Fluids.densityWetting - Fluids.densityNonWetting) * G * \
					sin(Geometry.angleToHorizontal * PI / 180.);
	return coefficientC;
}

double SingleCapillaryTube::calCoefficientD()
{
	double coefficientD;
	double tempCapillary;
	tempCapillary = calCapillaryPressure();
	coefficientD = Fluids.densityNonWetting * Geometry.length * G * \
				sin(Geometry.angleToHorizontal * PI / 180.) + \
				(rightPressure - leftPressure) - tempCapillary;
	return coefficientD;
}

double SingleCapillaryTube::calCapillaryPressure()
{
	//double capillaryPressure;
	auto capillaryPressure = 2.0 * Fluids.surfaceTension * cos(Fluids.contactAngle * PI / \
							180.) /Geometry.radius;
	return capillaryPressure;
}

double SingleCapillaryTube::useScantMethodWithAngle(double x1, double x0, double f1,
										   double f0, const double specificErr,
										   double t)
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
		valueForL0 = valueForL1;
		valueForL1 = valueL;
	}
	return tempL1;
}

double SingleCapillaryTube::useScantMethodWithoutAngle(double x1, double x0, double f1,
										   double f0, const double specificErr,
										   double t)
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
		valueL = calLocationFunctionWithoutAngle(Geometry, Fluids,
							tempL, initialLocation, t, initialTime);
		tempL0 = tempL1;
		tempL1 = tempL;
		valueForL0 = valueForL1;
		valueForL1 = valueL;
	}
	return tempL1;
}

void SingleCapillaryTube::outputInterfaceLocation(const std::vector<double> &interfaceLocation)
{
	using std::ofstream;
	ofstream outputInterfaceLocationFile;
	outputInterfaceLocationFile.open("location.dat");
	//int tempSize;
	//tempSize = sizeof(interfaceLocation) / sizeof(int);
	for (auto locationValue : interfaceLocation)
	{
		outputInterfaceLocationFile << locationValue << std::endl;
	}
	outputInterfaceLocationFile.close();
}
