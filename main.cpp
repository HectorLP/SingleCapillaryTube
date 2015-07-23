
#include <iostream>
#include <cmath>
#include <iomanip>

#include "DynamicSingleTube.h"

int main()
{
	std::cout << "Start to calculate the location of interface.\n";
	SingleCapillaryTube singleTube;
	//singleTube.calLocationInterfaceScant();
	singleTube.calLocationInterfaceBisect();
	std::cout << "finish the calculation.\n";
	return 0;
}