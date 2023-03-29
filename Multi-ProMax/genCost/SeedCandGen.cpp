#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <math.h>		//exp function

#include <fstream>		// std::file

#include "head.h"
#include "graph.h"

#define e exp(1)

using namespace std;

class Argument{
public:
	string dataset;
	string model;	//used for spread estimation
	int type;		//0 for rand, 1 for outdegree
	double epsilon; //control the accurarcy of spread accuracy
	int num;	//the percentage of the seed candidate set
	string index;	//index of the seedfile
	string name;
};

#include "SeedGen.h"

int main(int argc, char* argv[])
{
	if (argc < 5)
	{
		cout << "Usgae: ./SeedGen dataset model type epsilon" << endl; //// ../../dataset/hep IC 1 0.1
		exit(0);
	}

	Argument arg;
	//arg.type = 1;
	arg.epsilon = 0.1;
	//arg.num = 100;
	
	arg.dataset = argv[1];
	arg.model = argv[2];
	arg.type = atoi(argv[3]); ////gamma
	arg.epsilon = atof(argv[4]);
	//arg.scale = atoi(argv[4]);
	//arg.num = atof(argv[5]);
	//arg.index = argv[6];
	//arg.name = argv[7];

	ASSERT(arg.dataset != "");
	ASSERT(arg.model == "IC" || arg.model == "LT");
	string seedfile = arg.name;
	arg.dataset = arg.dataset + "/";	
	string graph_file;
    graph_file = arg.dataset + "graph.txt";
//	if (arg.model == "IC")
//		graph_file = arg.dataset + "graph_ic.inf";
//	else if (arg.model == "LT")
//		graph_file = arg.dataset + "graph_lt.inf";
//	else
//		ASSERT(false);


	SeedGen sg(arg.dataset, graph_file);
	if (arg.model == "IC")
		sg.setInfuModel(SeedGen::IC);
	else if (arg.model == "LT")
		sg.setInfuModel(SeedGen::LT);
	else
		ASSERT(false);

	sg.Gen(arg, seedfile);


	
	return 0;
}
