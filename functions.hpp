#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <random>
using namespace std;

const int N = 60700; //Januaray 20th 2020 to September 20th 2021
const float h = 0.01;

const string saveLocation = "";
const int numberOfFiles = 3;

void rungeKutta1(vector <float> &S, vector <float> &I, vector <float> &R, float Params[4]);
void rungeKutta2(vector <float> &S, vector <float> &I, vector <float> &R, float Params[4]);
void rungeKutta4(vector <float> &S, vector <float> &I, vector <float> &R, float Params[4]);

void rungekutta1(vector <float> &S, vector <float> &I, vector <float> &R, float Params[6], vector <float> &E, vector <float> &D);
void rungekutta2(vector <float> &S, vector <float> &I, vector <float> &R, float Params[6], vector <float> &E, vector <float> &D);
void rungeKutta4(vector <float> &S, vector <float> &I, vector <float> &R, float Params[6], vector <float> &E, vector <float> &D, vector<float> &Q);

void resetArrays(vector <float> &S, vector <float> &I, vector <float> &R);
void resetArrays(vector <float> &S, vector <float> &I, vector <float> &R, vector <float> &E, vector <float> &D, vector <float> &Q);

void write(string location,string fileName,vector <float> S,vector <float> I,vector <float> R);
void write(string location,string fileName,vector <float> S,vector <float> I,vector <float> R, vector <float> E, vector <float> &D, vector <float> &Q);

void generateSERData(string fileList[numberOfFiles], float Params[4]);
void generateSEIRData(string fileList[numberOfFiles], float Params[7]);


#endif /* functions_hpp */
