#include <chrono>
#include <ctime>
#include <ratio>
#include <cmath>




#define e exp(1)
#define c 2*(exp(1)-2)

using namespace std::chrono;

class Math{
    public:
        static double log2(int n){
            return log(n) / log(2);
        }
        static double logcnk(int n, int k) {
            double ans = 0;
            for (int i = n - k + 1; i <= n; i++)
            {
                ans += log(i);
            }
            for (int i = 1; i <= k; i++)
            {
                ans -= log(i);
            }
            return ans;
        }
        static bool min(double a, double b)
        {
            return a < b;
        }
};

class APM
{
public:

	static void ProfitMaximizeROI(InfGraph &g, const Argument &arg) //// single-ROI greedy
    {
        double total_cost = 0;
        double total_time = 0;
        double total_profit = 0;
        unsigned int totalseednum = 0;
        double delta=1./g.n;
        double epsilon=arg.epsilon;
        double bound=(8+2*epsilon)*g.n*(log(6/delta)+g.n*log(2))/epsilon/epsilon;
        double theta;
        for (int i = 0; i < arg.time; i++)
        {
            double theta_done=0;
            theta=arg.theta;
            high_resolution_clock::time_point startTime = high_resolution_clock::now();
            g.init_hyper_graph(arg); //contains the index of node sequence
            unsigned int cnt=0;
            double estprofit;
            while(theta<bound*(1+g.eps1)/max(1.0,g.EstPro3))
            {
                cnt++;
                g.buildhypergraph((unsigned int)(theta+1-theta_done));
                theta_done += (theta+1-theta_done);  // the number of RR sets that are already generated
                if(g.build_seedset_ROI((unsigned int)(theta+1), cnt, epsilon, delta, estprofit))break;
                theta*=2;
            }
//            cout << " seed: ";
//            for (int i = 0; i < g.seedSet.size(); i ++){
//                cout << g.seedSet[i] << ", ";
//            }cout << endl;

            high_resolution_clock::time_point endTime = high_resolution_clock::now();
            duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
            //unsigned int singlespread = g.simulation();
            double singlecost = 0;
            for (auto it : g.seedSet)singlecost += g.Cost[it];
            //double singleprofit = estprofit - singlecost;
            total_cost += singlecost;
            total_time += (double)interval.count();
            total_profit += estprofit;//cout << "total_profit: " << total_profit << endl;
            totalseednum += g.seedSet.size();
        }
        cout << "Sample size " << theta << endl;
        cout << "RunningTime(s) " << total_time / arg.time << endl;
        cout << "TotalSeedNum " << totalseednum / arg.time << endl;
        cout << "AvgCost " << total_cost / arg.time << endl;
        cout << "AvgProfit " << total_profit / arg.time << endl;

    }

    ////************************************************Multi-companys**********************************************


    static void MultiProfitMaxFill(InfGraph &g, const Argument &arg) ////Fill-greedy
    {
        vector<double> adv_cost(g.nrCompanies,0.0);
        vector<double> adv_profit(g.nrCompanies,0.0);
        vector<double> adv_inf(g.nrCompanies,0.0);
        double total_profit = 0.0, total_cost = 0.0;
        unsigned int totalseednum = 0;
        double total_time = 0;
        double delta=1./g.n;
        double epsilon=arg.epsilon;
        //////////////////////////sum_tau
        vector<double> cost_temp;
        vector<double> tau_i;
        double T1 = 0.,T2 = 0., gamma1 = max(arg.gammaP,arg.gammaR), gamma2 = min(arg.gammaP,arg.gammaR);
        for(int i = 0; i < g.nrCompanies; i++)
        {
            company *adv = g.advList.at(i);
            cost_temp.assign(g.Cost.begin(), g.Cost.end());
            sort(cost_temp.begin(), cost_temp.end());//sort-increasing order
            int size=cost_temp.size();
            float sum=0.0;
            bool flag=true;
            for(int j=0;j<g.n;j++)
            {
                sum+=cost_temp[j];
                if(sum> adv->budget)
                {
                    tau_i.push_back(j); //
                    flag=false;
                    cost_temp.erase(cost_temp.begin(),cost_temp.begin()+j);
                    break;
                }
            }
            if(flag)
                tau_i.push_back(g.n);

            T1+=adv->bpi*gamma1;
            T2+=adv->bpi*gamma2;
        }
        double sum_tau = 0.0;
        for(int k = 0;k < g.nrCompanies;k++)
        {
            sum_tau += tau_i[k]*log(2*(double)g.n/tau_i[k]);
        }
        if (sum_tau >= g.n * log(2*g.n)) sum_tau=g.n * log(2*g.n);
        //////////////////////////sum_tau

        double bound=(8+2*epsilon)*g.n*(log(6/delta)+sum_tau)/epsilon/epsilon;
        double theta_done;
        double theta;
        for (int i = 0; i < arg.time; i++)
        {
            theta_done = 0;
            theta = arg.theta; //arg.theta
            high_resolution_clock::time_point startTime = high_resolution_clock::now();
            g.init_hyper_graph(arg); //contains the index of node sequence
            unsigned int cnt=0;
            double EstRev1, EstRev2, EstCost, EstPro1, EstPro2, eps1, eps2, EstRev3;
            vector<bool> edgeMarkR2;
            while(theta<bound*(1+eps1)/max(1.0,EstRev3))//
            {
                cnt++;
                cout << "\nround " << cnt << " RR set: " << theta << endl;
                EstRev1 = 0.; EstRev2 = 0.; EstCost = 0.; EstPro1 = 0.; EstPro2 = 0.; eps1=0.; eps2=0.; EstRev3=0;
                ////R1
                g.buildhypergraph((unsigned int)(theta+1));
                g.candidate_greedy((theta+1));
                g.allocate_Fill((theta+1), arg.gammaR, arg.gammaP);
                EstRev1 = g.estimateRevenue(g.advList, theta+1, arg.gammaR, arg.gammaP);//F_R(S*,R1)
                EstCost = g.estimateCost(g.advList);
                EstPro1 = EstRev1 - EstCost;
                cout << "==========> current total profit on R1: " << EstPro1 << endl;
                ////R2
                g.updateInfoR2(g.advList,arg.gammaP); // initial RR set covered
                g.buildhypergraph_R2((unsigned int)(theta+1),g.advList, edgeMarkR2);//update rr set covered and node influenced in R2
                //g.buildhypergraphR_R2((unsigned int)(theta+1));//update rr set covered and node influenced in R2
                //g.updateInfoR2(g.advList,arg.gammaP); // initial RR set covered
                EstRev2 = g.estimateRevenueR2(g.advList, theta+1, arg.gammaR, arg.gammaP);//F_R(S*,R2)
                EstCost = g.estimateCost(g.advList);
                EstPro2 = EstRev2 - EstCost;
                cout << "==========> current total profit on R2: " << EstPro2 << endl;
                double x = EstRev2 * theta / log(5.*cnt*cnt/delta)/(g.n*T1);
                double y = (EstRev2-(1.+eps1)*EstCost) * theta / log(5.*cnt*cnt/delta)/(g.n*T2);
                eps1=4./(sqrt(1.+8.*x)-3);
                //ASSERT(eps1>0);
                eps2=sqrt((2.*eps1+2.)/y);
                //ASSERT(eps2>0);
                double beta=EstPro1/EstPro2;
                EstRev3 = EstRev2- (1+eps1)*EstCost;
                cout << "x: " << x<< " y: " << y<< " eps1: " << eps1 << " eps2: " <<eps2 << " beta: " << beta << " EstRev3: " << EstRev3 << endl;
                cout << "condition 1: " << (beta-1.)/beta+eps1+eps2 << " epsilon: " << epsilon << endl;
                if( (beta-1.)/beta+eps1+eps2 <=epsilon && eps1+eps2<=epsilon) break;
                theta_done += (theta+1-theta_done);  // the number of RR sets that are already generated
                theta*=2;
            }
            total_profit += EstPro2;
            total_cost += EstCost;
            high_resolution_clock::time_point endTime = high_resolution_clock::now();
            duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
            total_time += (double)interval.count();
            cout << "SingleRuntime " << (double)interval.count() << endl << endl;

            for (int i=0; i<g.nrCompanies; i++){
                company *adv = g.advList.at(i);
                adv_cost[adv->companyID] += adv->totalSeedCosts;
                adv_inf[adv->companyID] += adv->totalInf;
                adv_profit[adv->companyID] += adv->totalProfit;
                totalseednum += adv->seedSet.size();
            }

        }
        cout << "Sample size " << theta << endl;
        cout << "RunningTime(s) " << total_time / arg.time << endl;
        cout << "TotalSeedNum " << totalseednum / arg.time << endl;
        cout << "AvgCost " << total_cost / arg.time << endl;
        cout << "AvgProfit " << total_profit / arg.time << endl;
    }


    //// ******************************************* Balance Multi-companys **********************************************

    static void MultiProfitMaximize(InfGraph &g, const Argument &arg) //// one-by-one
    {
        vector<double> adv_cost(g.nrCompanies,0.0);
        vector<double> adv_profit(g.nrCompanies,0.0);
        vector<double> adv_inf(g.nrCompanies,0.0);
        double total_profit = 0.0, total_cost = 0.0;
        unsigned int totalseednum = 0;
        double total_time = 0;
        double theta = arg.theta; //stable rr set as input number
        for (int i = 0; i < arg.time; i++) {
            high_resolution_clock::time_point startTime = high_resolution_clock::now();
            g.init_hyper_graph(arg); //contains the index of node sequence
            g.buildhypergraph((unsigned int) (theta+1));
            g.build_multi_seedsets_obo((unsigned int) (theta+1),arg.gammaR, arg.gammaP); // select seed node one-by-one
            double single_profit = g.estimateProfit(g.advList, theta, arg.gammaR, arg.gammaP);
            cout << "========================> total profit (one by one): " << single_profit << endl<<endl;
            total_profit += single_profit;
            high_resolution_clock::time_point endTime = high_resolution_clock::now();
            duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
            total_time += (double) interval.count();

            ///--------------------------------------------------------
            for (int i=0; i<g.nrCompanies; i++){
                company *adv = g.advList.at(i);
                adv_cost[adv->companyID] += adv->totalSeedCosts;
                adv_inf[adv->companyID] += adv->totalInf;
                adv_profit[adv->companyID] += adv->totalProfit;
                totalseednum += adv->seedSet.size();
            }
            ///--------------------------------------------------------
        }
        for (int i=0; i<g.nrCompanies; i++){company *adv = g.advList.at(i);total_cost += adv_cost[adv->companyID];}
        cout << "Sample size " << theta << endl;
        cout << "RunningTime(s) " << total_time / arg.time << endl;
        cout << "TotalSeedNum " << totalseednum / arg.time << endl;
        cout << "AvgCost " << total_cost / arg.time << endl;
        cout << "AvgProfit " << total_profit / arg.time << endl;
    }


    static void MultiProfitMaxIter(InfGraph &g, const Argument &arg) //// Iteration
    {
        cout <<"profit batch: " << arg.batchPS << "\tinfluence batch: " << arg.batchIS << endl;
        vector<double> adv_cost(g.nrCompanies,0.0);
        vector<double> adv_profit(g.nrCompanies,0.0);
        vector<double> adv_inf(g.nrCompanies,0.0);
        double total_profit = 0.0, total_cost = 0.0;
        unsigned int totalseednum = 0;
        double total_time = 0;
        double theta = arg.theta; //stable rr set as input number
        for (int i = 0; i < arg.time; i++) {
            high_resolution_clock::time_point startTime = high_resolution_clock::now();
            g.init_hyper_graph(arg); //contains the index of node sequence
            g.buildhypergraph((unsigned int) (theta+1));
            g.allocate_Iter((theta+1),arg.batchPS, arg.batchIS, arg.gammaR, arg.gammaP);
            double single_profit = g.estimateProfit(g.advList, theta, arg.gammaR, arg.gammaP);
            cout << "========================> total profit (Iter method): " << single_profit << endl<<endl;
            total_profit += single_profit;
            high_resolution_clock::time_point endTime = high_resolution_clock::now();
            duration<double> interval = duration_cast<duration<double>>(endTime - startTime);
            total_time += (double) interval.count();
            ///--------------------------------------------------------
            for (int i=0; i<g.nrCompanies; i++){
                company *adv = g.advList.at(i);
                adv_cost[adv->companyID] += adv->totalSeedCosts;
                adv_inf[adv->companyID] += adv->totalInf;
                adv_profit[adv->companyID] += adv->totalProfit;
                totalseednum += adv->seedSet.size();
            }
            ///--------------------------------------------------------
        }
        for (int i=0; i<g.nrCompanies; i++){company *adv = g.advList.at(i);total_cost += adv_cost[adv->companyID];}
        cout << "Sample size " << theta << endl;
        cout << "RunningTime(s) " << total_time / arg.time << endl;
        cout << "TotalSeedNum " << totalseednum / arg.time << endl;
        cout << "AvgCost " << total_cost / arg.time << endl;
        cout << "AvgProfit " << total_profit / arg.time << endl;
    }

};

