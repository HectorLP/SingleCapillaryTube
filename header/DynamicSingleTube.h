#ifndef DYNAMICSINGLETUBE_H_
#define DYNAMICSINGLETUBE_H_

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

	const char *typeFunction;
	
	double leftPressure;
	double rightPressure;
	
	double initalLocation;
	double initialTime;
	
	double timeStep;
	double timeEndPoint;
	
	void loadTubeGeometry();
	void loadFluidProperties();
	void loadPressureValues();
	void loadPositionAndTime();
	
	double calCoefficientA(const TubeGeometry &TG, const FluidProperties &FP);
	double calCoefficientB(const TubeGeometry &TG, const FluidProperties &FP);
	double calCoefficientC(const TubeGeometry &TG, const FluidProperties &FP);
	double calCoefficientD(const TubeGeometry &TG, const FluidProperties &FP);
	
	double calCapillaryPressure(const TubeGeometry &TG, const FluidProperties &FP);
	
	double calLocationFunctionWithAngle(const TubeGeometry &TG, const FluidProperties &FP);
	double calLocationFunctionWithoutAngle(const TubeGeometry &TG, const FluidProperties &FP);
	
	double useScantMethod();
public:
	SingleCapillaryTube();
	~SingleCapillaryTube() {}
	double calLocationInterface();
	//double useScantMethod();
	//friend double useFixedMethod();
};
#endif 
