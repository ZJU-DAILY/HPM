#include <numeric>		// std::itoa
#include <ctime>		// time()
#include <utility>		// std::pair
#include <algorithm>	// std::sort
#include <math.h>       /* pow */

class SeedGen : public Graph
{
private:	
	vector<bool> visit;
	vector<int> visit_mark;	
	//vector<int>seq;		
	sfmt_t sfmtSeed;

public:
	
	vector<int>SeedSet;
	vector<double>cst;

	SeedGen(string folder, string graph_file) :Graph(folder, graph_file)
	{
		srand(time(NULL));
		sfmt_init_gen_rand(&sfmtSeed, rand());

		visit = vector<bool>(max_n+1);
		visit_mark = vector<int>(max_n+1);

		//seq = vector<int>(n);
		cst = vector<double>(n, 0);
	}

	void Gen(Argument &arg, string seedfile)
	{
		//SeedGenOutdeg(size);
		//seedfile += "_" + to_string(size);
		//LoadSeedSet(seedfile, size);
                
		SeedSet.resize(n);
        std::copy(nodes.begin(), nodes.end(), SeedSet.begin());
        ////iota(SeedSet.begin(), SeedSet.end(), 0); //initialize from 0 to n;
		
		double lambda_prime=0;
		if(arg.type!=3){
		    int max_node=0;
		    double threshold = SpreadEst(arg, max_node);

		    //random_shuffle(SeedSet.begin(), SeedSet.end()); //shuffle this first
 		    lambda_prime=0.95*threshold/outdeg[max_node];

		    cout << "node: " << max_node << " spread: " << threshold <<" outdeg: " << outdeg[max_node] << " lambda_prime: " << lambda_prime << endl;
            //getchar();
        }
        lambda_prime=0.95;
		//arg.type = 0;		
		//CostGenNor(lambda_prime, arg.type);
		//Output(arg);

		//arg.type = 1;
		if(arg.type==0){
            arg.type = 2;
            CostGenNor(lambda_prime, arg.type);
            Output(arg);

		    arg.type = 4;
		    CostGenNor(lambda_prime, arg.type);
		    Output(arg);

		    arg.type = 5;
		    CostGenNor(lambda_prime, arg.type);
		    Output(arg);

            arg.type = 6;
            CostGenNor(lambda_prime, arg.type);
            Output(arg);

            arg.type = 7;
            CostGenNor(lambda_prime, arg.type);
            Output(arg);

            arg.type = 8;
            CostGenNor(lambda_prime, arg.type);
            Output(arg);

            arg.type = 9;
            CostGenNor(lambda_prime, arg.type);
            Output(arg);
		}
		else{
		    CostGenNor(lambda_prime, arg.type);
		    Output(arg);
		}
		//if (arg.type == 0)SeedGenRand(size);
		//else if (arg.type == 1)SeedGenOutdeg(size);
		//else ASSERT(false);
        //Output_graph(arg);
	}



	//how to generate cost, should the cost be aligned with its outdegree?
	void CostGenNor(double factor, int type)
	{
		//random generate cost
			/*
			double sum = 0;
			for (auto it : SeedSet)
			{
				cst[it] = sfmt_genrand_real1(&sfmtSeed);
				sum += cst[it];
			}
			double scale = threshold / sum;
			for (auto it : SeedSet)cst[it] *= scale;
			*/
		//generate based on outdegree

		if (type == 1)
		{
			//double degree = 0;
			//for (auto it : SeedSet)degree += outdeg[it];
			//double factor = threshold / degree;
			for (auto it : SeedSet){
				if(outdeg[it]==0){
					cst[it]=factor;
				}
				else{
					cst[it] = outdeg[it] * factor;
				}
			}
		}
		else if (type == 3) // uniform
		{
			//double cstsum=0;
			//for (auto it : SeedSet) cstsum += outdeg[it] * factor;
			
			//double avgcst=cstsum/SeedSet.size();
			for (auto it : SeedSet) cst[it] = 1;//avgcst;
                        
		}
		else if (type ==2) //degree^0
		{
			for (auto it : SeedSet)cst[it] = factor;
		}
		else if (type==4) // degree^0.2
		{
			for (auto it : SeedSet){
				if(outdeg[it]==0){
					cst[it]=factor;
				}
				else{
					cst[it] = pow(outdeg[it],0.2) * factor;
				}
			}
		}
		else if (type==5) // degree^0.4
		{
                        for (auto it : SeedSet)
			{
                                if(outdeg[it]==0){
                                     cst[it]=factor;
                                }
                                else{
				     cst[it] = pow(outdeg[it],0.4) * factor;
				}
			}
		}
		else if (type==6) // degree^0.6
		{
                        for (auto it : SeedSet){
                                if(outdeg[it]==0){
                                    cst[it]=factor;
                                }
                                else{
				    cst[it] = pow(outdeg[it],0.6) * factor;
				}
			}
		}
		else if (type==7) // degree^0.8
		{
                        for (auto it : SeedSet){
                                if(outdeg[it]==0){
                                    cst[it]=factor;
                                }
                                else{
				    cst[it] = pow(outdeg[it],0.8) * factor;
				}
			}
		}
		else if (type==8) // degree^1.2
		{
                        for (auto it : SeedSet){
                                if(outdeg[it]==0){
                                    cst[it]=factor;
                                }
                                else{
				    cst[it] = pow(outdeg[it],1.2) * factor;
				}
			}
		}
		else if (type==9) // degree^1.4
		{
                        for (auto it : SeedSet){
                                if(outdeg[it]==0){
                                    cst[it]=factor;
                                }
                                else{
				    cst[it] = pow(outdeg[it],1.4) * factor;
				}
			}
		}
		else
			ASSERT(false);				

	}


	void Output(Argument arg)
	{
		ofstream out(arg.dataset + "/new/cost_" + to_string(arg.type) + ".txt");
		cout << "writing " << arg.dataset + "new_cost_" + to_string(arg.type) + ".txt" << endl;
		for (auto node : SeedSet)
        {
            out << node << " " << cst[node] << endl;
            cout << node << " " << cst[node] << endl;
        }

		out.close();
	}

    void Output_graph(Argument arg)
    {
        ofstream out(arg.dataset + "graph_prob.txt");
        cout << "writing " << arg.dataset + "graph_prob.txt" << endl;
        set <int> :: iterator iter;
        for(iter = nodes.begin(); iter != nodes.end(); iter ++){
            auto node = *iter;
            if (gT[node].size()==0) continue;
            else{
                for (int j = 0; j < gT[node].size(); j++){
                    out << gT[node][j] << " " << node <<" "<< probT[node][j] << endl;
                }
            }
        }
        out.close();
    }




    double SpreadEst(const Argument & arg, int& max_node)
	{
        cout << "spread estimation ... " << endl;
		vector<int>spread(max_n, 0);
		//for (auto it : SeedSet)tag[it] = true;

		double epsilon = arg.epsilon, delta = 1. / n;
		double lambda = 4 * (e - 2)*log(2 / delta) / epsilon / epsilon;
		double tau = 1 + (1 + epsilon)*lambda;
		int cnt = 0, max_spread=0;
        cout << "tau: " << tau<< endl;
		while (max_spread < tau)
		{
		  BuildHypergraph(spread, cnt++);		  
		  auto max_spread_ptr=max_element(spread.begin(), spread.end());
		  max_spread=*max_spread_ptr;
		  max_node=max_spread_ptr-spread.begin();
		}
        cout <<"cnt: " << cnt << "\tmax_spread: " << max_spread << endl;
		return 1.*max_spread / cnt*n / (1 + epsilon);
	}

	void BuildHypergraph(vector<int>& spread, int hyperId)
	{
        //cout << "build rr...\t" ;
		////const auto uStart = sfmt_genrand_uint32(&sfmtSeed) % n;
        int uStart;
        //// con1: uStart has no information (uStart does not exist)
        int flag = 0;
        set<int>::iterator iter;
        while (flag == 0){
            uStart = sfmt_genrand_uint32(&sfmtSeed) % max_n;
            iter = nodes.find(uStart);
            if(iter != nodes.end()) flag = 1;
            else continue;
        }
        //cout << "current uStart: " <<uStart<<endl;
        //// con2: gT[uStart].size=0/t (uStart, x) /(x, uStart)
		unsigned int n_visit_mark = 0, curIdx = 0;
		visit_mark[n_visit_mark++] = uStart;
		visit[uStart] = true;		

		if (influModel == IC)
		{
			while (curIdx < n_visit_mark)
			{
				int i = visit_mark[curIdx++];
				for (int j = 0; j < (int)gT[i].size(); j++)
				{
					int v = gT[i][j];
					if (visit[v])continue;

					double randDouble = sfmt_genrand_real1(&sfmtSeed);
					if (randDouble > probT[i][j])continue;

					visit[v] = true;
					visit_mark[n_visit_mark++] = v;				
				}				
			}
		}
		else if (influModel == LT)
		{
			while (curIdx < n_visit_mark)
			{
				int expand = visit_mark[curIdx++];

				if (gT[expand].size() == 0)
					continue;

				int index = sfmt_genrand_uint32(&sfmtSeed) % gT[expand].size();
				int v = gT[expand][index];
				if (visit[v])continue;
				/*
				if (tag[v])
				{
					for (unsigned int i = 0; i < n_visit_mark; i++)visit[visit_mark[i]] = false;
					return true;
				}
				*/
				visit[v] = true;
				visit_mark[n_visit_mark++] = v;				
			}
		}
		else
			ASSERT(false);	

		for (unsigned int i = 0; i < n_visit_mark; i++)
		{
		  visit[visit_mark[i]] = false;
		  ++spread[visit_mark[i]];
		}				
	}
};
