#define HEAD_INFO

#include "sfmt/SFMT.h"
#include "head.h"
#include "string.h"

class Argument{
public:	
    string dataset;
    double epsilon;
    string model;   
	string type;			//0 for rand, 1 for outdegree
	double scale;			// scale of cost
	//unsigned int batch;	  //this might be useless.
    string seedfile;
    int time;  //this is the number of times needed to measure the result.
	double theta;
    int algo; // 0 for single company; 1 for one-by-one; 2 for Fill; 3 for Iter
    ////
    double gammaR, gammaP; // paramether of reward and penalty
    int batchPS;
    int batchIS;
    int h;

};

#include "graph.h"
#include "infgraph.h"
#include "pm.h"

//static unsigned int rr_num=0;
void OutputSeedSetToFile(vector<int> seed_set, string seedfile)
{
    ofstream of;
    of.open(seedfile);
    for (int seed: seed_set)
    {
        of << seed << " ";
    }
	of << endl;
    of.close();
}

void run_with_parameter(InfGraph &g, const Argument & arg)
{
    if (arg.h == 1) {
            cout << "\n********************************** single company (ROI) ***************************************" << endl;
            APM::ProfitMaximizeROI(g, arg);
            cout << "******************************** single company (ROI) end ***************************************" << endl;

    }
    else {
        if (arg.algo == 1) {
            cout << "\n************************* multiple companies (Fill method) *******************************" << endl;
            APM::MultiProfitMaxFill(g, arg);
            cout << "************************* multiple companies (Fill method) end *****************************" << endl;
        }

        //////////////////////////////////////  balance influence distribution   /////////////////////////////////
        if (arg.algo == 2) {
            cout << "\n************************* multiple companies (one by one) *******************************" << endl;
            APM::MultiProfitMaximize(g, arg);
            cout << "************************* multiple companies (one by one) end *****************************" << endl;
        }
        if (arg.algo == 3){
            cout << "\n************************** multiple companies (Iter method) ****************************" << endl;
            APM::MultiProfitMaxIter(g, arg);
            cout << "\n************************ multiple companies (Iter method) end **************************" << endl;
        }
    }
    //INFO(g.seedSet);
    //OutputSeedSetToFile(g.seedSet, arg.seedfile);
}
void Run(int argn, char **argv)
{
    Argument arg;

  	arg.scale = 0.2;
  	arg.type = "1";
    arg.algo=1;
  	arg.theta = 10000;
    arg.epsilon=0.2;
    arg.batchPS = 10;
    arg.batchIS = 5;
    arg.gammaR = 1; // reward parameter gamma
    arg.gammaP = 1; // penalty parameter gamma
    arg.h = 5; // number of companies

    for (int i = 0; i < argn; i++)
    {

        if (argv[i] == string("-help") || argv[i] == string("--help") || argn == 1)
        {
            cout << "./apm -dataset *** -epsilon ***  -model IC|LT -time *** -theta *** -h *** -scale *** -type *** -gamma *** -algo *** -batchPS *** -batchIS ***" << endl;
            return ;
        }
		
		if (argv[i] == string("-dataset"))
			arg.dataset = argv[i + 1];
		if (argv[i] == string("-epsilon"))
			arg.epsilon = atof(argv[i + 1]);		
		if (argv[i] == string("-model"))
			arg.model = argv[i + 1];
		if (argv[i] == string("-seedfile"))
			arg.seedfile = argv[i + 1];
        if (argv[i] == string("-scale"))
            arg.scale = atof(argv[i + 1]);
		if (argv[i] == string("-type"))
			arg.type = argv[i + 1];
        if (argv[i] == string("-h"))
            arg.h = atoi(argv[i + 1]);
		if (argv[i] == string("-time"))
			arg.time = atoi(argv[i + 1]);
		if (argv[i] == string("-algo"))
			arg.algo = atoi(argv[i + 1]);
		if (argv[i] == string("-theta"))
			arg.theta = atof(argv[i + 1]);
        if (argv[i] == string("-gammaR"))
            arg.gammaR = atof(argv[i + 1]);
        if (argv[i] == string("-gammaP"))
            arg.gammaP = atof(argv[i + 1]);
        if (argv[i] == string("-batchPS"))
            arg.batchPS = atof(argv[i + 1]);
        if (argv[i] == string("-batchIS"))
            arg.batchIS = atof(argv[i + 1]);
    }
    ASSERT(arg.dataset != "");
    ASSERT(arg.model == "IC" || arg.model == "LT");

    arg.dataset = arg.dataset + "/";
    string graph_file = arg.dataset + "graph.txt";
    string cost_file = arg.dataset + "/cost_" + arg.type + ".txt";
    string budget_file = arg.dataset +"/budgets_" + to_string(arg.h) + ".txt";

    InfGraph g(arg.dataset, graph_file, cost_file, budget_file, arg.scale ,arg.gammaR, arg.gammaP, arg.h);

    if (arg.model == "IC")
        g.setInfuModel(InfGraph::IC);
    else if (arg.model == "LT")
        g.setInfuModel(InfGraph::LT);
    else
        ASSERT(false);	

    run_with_parameter(g, arg);
}


int main(int argn, char **argv)
{
    __head_version = "v1";
    OutputInfo info(argn, argv);
    
    Run( argn, argv );
    // disp_mem_usage();
    cout<<"Memory(MB) "<<getProcMemory()<<endl;
}

