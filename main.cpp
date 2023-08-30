/*
 * main.cpp
 *
 *  Created on: Aug 28, 2023
 *      Author: misterchief53
 */

#include <cmath>
#include <iostream>
#include "./include/pso.h"

double rosenborok(const std::vector<double>& position, double a, double b) {
	return std::pow(a - position[0], 2.0) + b * std::pow(position[1] - std::pow(position[0], 2.0), 2.0);
}
/*
double rosenborok(double x, double y, double a, double b) {
	return std::pow(a - x, 2.0) + b * std::pow(y - std::pow(x, 2.0), 2.0);
}
*/
/*
double rosenborok(double x, double y) {
	return std::pow(1 - x, 2) + 100 * std::pow(y - std::pow(x, 2), 2);
}
*/

void computeRosenborok() {
	double res;
	//a usually set to 1, b usually set to 100,
	//global minima is where x and y = a and a^2
	std::vector<double> position {1.0,1,0};
	res = rosenborok(position, 1, 100);
	std::cout << "Rosenborok sample: " << res << std::endl;
	pso optimizer = pso(2.0, -100.0, 100.0, 10000.0, 10.0, 1.0, 2.0, 2.0, rosenborok, 1.0, 100.0);
	optimizer.print_particles();
	optimizer.compute();
	optimizer.print_particles();
}

#include <iostream>

int main(){
	std::cout << "hello world\n";
	computeRosenborok();
}

