#ifndef DYNAMICSINGLETUBE_H_
#define DYNAMICSINGLETUBE_H_

#include <string>
#include <vector>

struct TubeGeometry
{
	double radius;
	double length;
	double angleToHorizontal;
};

struct FluidProperties
{
	double densityWetting;
	double densityNonWetting;
	double dynamicViscosityWetting;
	double dynamicViscosityNonWetting;
	double surfaceTension;
	double contactAngle;
};

class SingleCapillaryTube
{
private:
	TubeGeometry Geometry;
	FluidProperties Fluids;
	
	double leftPressure;
	double rightPressure;
	
	double initialLocation;
	double initialTime;
	
	double timeStep;
	double timeEndPoint;
	
	double CoefficientA;
	double CoefficientB;
	double CoefficientD;

	std::string typeFunction;
	
	void loadTubeGeometry();
	void loadFluidProperties();
	void loadPressureValues();
	void loadPositionAndTime();
	
	double calCoefficientA();
	double calCoefficientB();
	double calCoefficientC();
	double calCoefficientD();
	void calThreeCoefficients();
	
	double calCapillaryPressure();
	
	double calLocationFunctionWithAngle(const TubeGeometry &TG, const FluidProperties &FP,
					double tempLocation, double initialLocation,
					double tempTime, const double initialT);
	double calLocationFunctionWithoutAngle(const TubeGeometry &TG, const FluidProperties &FP, 
					double tempLocation, double initialLocation,
					double tempTime, const double initialT);
	
	double useScantMethodWithAngle(double x1, double x2, double f1, double f0, 
					const double specificErr, double timePoint);
	double useScantMethodWithoutAngle(double x1, double x2, double f1, double f0,
					const double specificErr, double timePoint);
	void rhsPart(const double location, double& dldt, const double t = 0.);
	double useOdeSolver();
	void outputInterfaceLocation(const std::vector<double> &interfaceLocation);
public:
	SingleCapillaryTube();
	~SingleCapillaryTube() {}
	void calLocationInterfaceScant();	//use Scant method to solve interface location
	void calLocationInterfaceBisect();	//use bisect method from boost library.
	void calLocationInterfaceBrent();
	void calInterfaceLocationSolveOde();	//use the ode solver from boost library
 	//double useScantMethod();
	//friend double useFixedMethod();
};
#endif 
