#ifndef DYNAMICSINGLETUBE_H_
#define DYNAMICSINGLETUBE_H_

#include <string>

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

	std::string typeFunction;
	
	void loadTubeGeometry();
	void loadFluidProperties();
	void loadPressureValues();
	void loadPositionAndTime();
	
	double calCoefficientA(const TubeGeometry &TG, const FluidProperties &FP);
	double calCoefficientB(const TubeGeometry &TG, const FluidProperties &FP);
	double calCoefficientC(const TubeGeometry &TG, const FluidProperties &FP);
	double calCoefficientD(const TubeGeometry &TG, const FluidProperties &FP);
	
	double calCapillaryPressure(const TubeGeometry &TG, const FluidProperties &FP);
	
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
public:
	SingleCapillaryTube();
	~SingleCapillaryTube() {}
	double calLocationInterface();
	//double useScantMethod();
	//friend double useFixedMethod();
};
#endif 
