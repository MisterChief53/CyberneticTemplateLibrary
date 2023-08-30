/*
 * pso.h
 *
 *  Created on: Aug 28, 2023
 *      Author: misterchief53
 *  Notes: Can only use this with numbers that
 *  can be converted from "real numbers", due to
 *  the random distribution.
 */

#ifndef INCLUDE_PSO_H_
#define INCLUDE_PSO_H_

#pragma once

#include <cstddef>
#include <vector>
#include <random>
#include <iostream>
#include <tuple>

template <typename T, typename Func, typename... Args>
class particle
{
public:
	std::vector<T> position; // Vector since it can be on N dimensions.
	std::vector<T> velocity; // Assuming the same here
	T cost;
	std::vector<T> best_position;
	T best_cost;
	T& global_best_cost;
	std::vector<T>& global_best_position;

	particle(const std::vector<T>& position, T& global_best_cost,
			std::vector<T>& global_best_position,Func func, Args... args)
		: position(position),
		  velocity(position.size(),0), //I rather specify values whenever I can.
		  cost(func(position, args...)),
		  best_position(position),
		  best_cost(cost),
		  global_best_cost(global_best_cost), //Initialize ref member
		  global_best_position(global_best_position)
	{
		/*
		 * For the cost above
		 * Call the func function with position and optional args
		 * Remember that the args have to be of the same type that the function
		 * accepts, as well as the positions!
		 */

		if(best_cost < global_best_cost){
			this->global_best_cost = best_cost;
			this->global_best_position = best_position;
		}else{
			best_cost = global_best_cost;
			best_position = global_best_position;
		}
	}
};

template <typename T, typename Func, typename... Args>
class pso
{
public:
	size_t nVar; //Unknown variables
	//std::vector<double> varSize; //Matrix size of decision variables
	T lowerBound;
	T upperBound;
	size_t maxIterations;
	size_t nPop;
	T w; //Innertial coefficient, usually 1
	T wDamp = 0.99; //Damping ration of innertia weight
	T c1; //Personal acceleration coefficient, usually 2
	T c2; //Global acceleration coefficient, ususally 2
	std::vector<particle <T, Func, Args...>> particles;
	std::random_device rd;
	std::mt19937 gen;
	std::uniform_real_distribution<T> dist;
	T global_best_cost;
	std::vector<T> global_best_position;
	Func func;
	std::tuple<Args...> argsTuple;

	pso(const size_t& nVar, const T& lowerBound, const T& upperBound,
			const size_t& maxIterations, const size_t& nPop, const T& w,
			const T& c1, const T& c2, Func func, Args... args)
			: argsTuple(std::make_tuple(args...)),
			  gen(std::random_device()()){
		this->nVar = nVar;
		this->lowerBound = lowerBound;
		this->upperBound = upperBound;
		this->maxIterations = maxIterations;
		this->nPop = nPop;
		this->w = w;
		this->c1 = c1;
		this->c2 = c2;
		dist = std::uniform_real_distribution<T>(lowerBound, upperBound);
		global_best_cost = std::numeric_limits<T>::max();

		this->func = func;

		/*
		 * Create population array.
		 * We use auto to avoid dealing with nested types and other complexity.
		 * without it, we have to figure out complex types in case that T is something
		 * more than a simple primitive.
		 */
		for(size_t i = 0; i < nPop; i++){
			std::vector<T> position(nVar,0);
			for(size_t j = 0; j < nVar; j++){
				position[j] = dist(gen); //Supposedly more performant that emplacing back
			}
			particles.emplace_back(position, global_best_cost, global_best_position,
					func, args...);
		}
	}

	void compute(){
		std::cout << "started compute func\n";
		for(size_t i = 0; i < maxIterations; i++){
			for(size_t j = 0; j < nPop; j++){
				// TODO: for paralelization, we should probably offload this computation
				// to each particle.

				// We use vectors since std::array needs size known at compile time
				std::vector<T> random1(nVar,0), random2(nVar,0), velocityTerm(nVar,0),
					c1Term(nVar,0), c2Term(nVar,0), c1Subtraction(nVar,0), c2Subtraction(nVar,0);
				for(size_t k = 0; k < nVar; k++){
					//This divides the random number by the magnitude of the upper and lower bounds.
					//Not doing this resulted in this being extremely erratic.
					random1[k] = dist(gen) / abs(upperBound - lowerBound);
					random2[k] = dist(gen) / abs(upperBound - lowerBound);

					velocityTerm[k] = w * particles[j].velocity[k];
					c1Term[k] = c1*random1[k];
					c2Term[k] = c2*random2[k];
					c1Subtraction[k] = particles[j].best_position[k] - particles[j].position[k];
					c2Subtraction[k] = global_best_position[k] - particles[j].position[k];

					// So this should translate to calculating each term... hopefully.
					particles[j].velocity[k] = velocityTerm[k] + c1Term[k] * c1Subtraction[k] + c2Term[k] * c2Subtraction[k];

					// update the position, element wise as well.
					particles[j].position[k] += particles[j].velocity[k];
				}


				auto argsWithPosition = std::tuple_cat(std::make_tuple(particles[j].position), argsTuple);
				//particles[i].cost = func(particles[i].position, argsTuple);
				particles[j].cost = std::apply(func, argsWithPosition);
				//std::cout << "particle's " << j+1 << " cost: " << particles[j].cost << '\n';
				//std::cout << "position: " << particles[j].position[0] << '\n';
				if(particles[j].cost < particles[j].best_cost){
					//std::cout << "new local best cost found\n";
					particles[j].best_position = particles[j].position;
					particles[j].best_cost = particles[j].cost;

					if(particles[j].best_cost < global_best_cost){
						//std::cout << "new global best" << particles[j].best_cost << '\n';
						global_best_cost = particles[j].best_cost;
						global_best_position = particles[j].best_position;
					}
				}
			}

			std::cout << "Epoch: " << i+1 << " Best Cost: " << global_best_cost << "\n";

			// Damping innertia
			w *= wDamp;
		}
	}

	void print_particles(){
		for(const auto& x : particles){
			std::cout << "position: ";
			for(const auto& c: x.position){
				std::cout << c << " ";
			}
			std::cout << "\n";
			std::cout << "cost: " << x.cost << "\n";
		}
		std::cout << std::endl;
		std::cout << "global best position: ";
		for(const auto& x : global_best_position){
			std::cout << x << " ";
		}
		std::cout << std::endl;
		std::cout << "global best cost: " << global_best_cost << '\n';
		std::cout << std::endl;
	}
};

#endif /* INCLUDE_PSO_H_ */
