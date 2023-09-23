#include "iheap.h"
#include <queue>	//priority_queue
#include <utility>  // pair
#include <numeric>
#include <unordered_map>
#include <limits.h>
#include "company.h"
#include <algorithm>
#include <float.h>

#include <iostream>
#include <sstream>
#include <iomanip>


double round(double number, unsigned int bits) {
    stringstream ss;
    ss << fixed << setprecision(bits) << number;
    ss >> number;
    return number;
}

bool comparison(company* a,company* b){
    return a->bpi > b->bpi;
}

class InfGraph: public Graph
{
private:
    vector<bool> visit;
    vector<int> visit_mark;


public:

    vector<vector<int>> hyperG;
    vector<vector<int>> hyperGT;

    vector<vector<int>> hyperG_R2;
    vector<vector<int>> hyperGT_R2;

    vector<vector<int>>PO; //useless
    double EstPro1, EstPro2, EstPro3, eps1, eps2;

    vector<double>Cost;
    vector<int> seedSet;

    //// 0
    typedef std::vector<company*> companyList;
    companyList advList;
    companyList order_advList;

    //// 1 obo_greedy
    typedef struct{
        int advMaxId;
        double Maxweight;
    }adoptInfo;
    vector <adoptInfo> Adoption;
    //int timeCount; // record user influenced time

    /// 2 Fill
    vector<int> candidateSetT;
    vector<int> isSelect; // disjointness
    int candidateNode;
    set <int> usersExpand;
    //multimap<double, int> allocQueueM; // allocation priority, ordered by first value (double) in increase order
    //multimap<double,int> criterQueue;
    //typedef std::pair<int, double> infPair; // adv - node - profit marginal gain
    typedef struct{
        int adv;
        int node;
        double mpg; //return a.mpg<b.mpg;//从小到大排<，若要从大到小排则>
    } infPair;
    static bool cmp(const infPair &a, const infPair &b)
    {
        return a.mpg<b.mpg;//从小到大排<，若要从大到小排则>
    }
    vector<infPair> allocQueueM;
    vector<infPair> criterQueue;



    sfmt_t sfmtSeed;


    InfGraph(string folder, string graph_file, string cost_file, string budget_file, double scale, double gammaR, double gammaP, int h): Graph(folder, graph_file)
    {
        std::cout << "Initial InfGraph ... " << endl;
        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed , rand());
        //sfmt_init_gen_rand(&sfmtSeed, 95082); //By using a determined seed number, we could debug without randomness.

        visit = vector<bool> (max_n);
        visit_mark = vector<int> (max_n);

        Cost = vector<double>(max_n, 0);
        hyperG.resize(max_n, vector<int>());
        hyperG_R2.resize(max_n,vector<int>());
        nrCompanies = h;
        cout << "n: " <<n<<"\tm: "<<m<<"\tcompanies: "<< nrCompanies << endl;

        for(int i = 0; i < nrCompanies; i++)
        {

            company *aa = new company(i);
            advList.push_back(aa);

        }
        load_cost(cost_file, scale);
        load_budget(budget_file,gammaP);
    }

    ////load cost_file (cost_scale)
    void load_cost(string filename, double scale)
    {
        cout << "Read Incentive Costs ..." << endl;
        std::ifstream  infile(filename);
        if (!infile.is_open())
        {
            std::cout << "The file \"" + filename + "\" can NOT be opened\n";
            return;
        }
        for (unsigned int i = 0; i < n; i++){
            int node;
            double cst;
            infile >> node >> cst;
            //ASSERT( node < n );
            cout << "node: " << node << "\tcst: " << cst << endl;
            ASSERT(cst > 0);
            Cost[node] = scale*cst;

        }
        infile.close();
    }

    ////load budget_file (budgets)
    void load_budget(string filename, double gammaP)
    {
        cout << "Read Budget and Influence ..." << endl;
        cout << "Read Incentive Costs ..." << endl;
        std::ifstream  infile(filename);
        if (!infile.is_open())
        {
            std::cout << "The file \"" + filename + "\" can NOT be opened\n";
            return;
        }
        company *adv;
        int advIndex = 0;
        for (unsigned int i = 0; i < nrCompanies; i++){
            float budget, influence;
            infile >> budget >> influence;
            adv = advList.at(advIndex++);
            adv->budget = budget;
            adv->influence = influence;
            adv->bpi = (float) budget/influence;
            ////adv->bpi = adv->bpi * gammaP;
        }
        infile.close();
        //sort(advList.begin(),advList.end(),comparison); //// order companies based on decreasing order of Bi/Ii
        //cout <<" After sort advList based on Bi/Ii " << endl;
        for (unsigned int i = 0; i < nrCompanies; i++){
            adv = advList.at(i);
            cout << "companyID: " << adv->companyID <<"\tbudget: " << adv->budget << "\tinfluence: " << adv->influence << endl;
        }
    }

    char* map_file(const char* fname, size_t& length)
    {
        int fd = open(fname, O_RDONLY);
        if (fd == -1)
            handle_error("open");

        // obtain file size
        struct stat sb;
        if (fstat(fd, &sb) == -1)
            handle_error("fstat");

        length = sb.st_size;

        char* addr = static_cast<char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
        if (addr == MAP_FAILED)
            handle_error("mmap");

        // TODO close fd at some point in time, call munmap(...)
        close(fd);
        return addr;
    }


    void init_hyper_graph(const Argument & arg)
    {
        //fill(spread.begin(), spread.end(), 0);
        std::cout << "Initial hyper_graph ..." << endl;
        cout <<"gamma_reward: " << arg.gammaR << "\tgamma_penalty: " << arg.gammaP << endl;
        for (auto& hyper : hyperG)hyper.clear();
        for (auto& hyperT : hyperGT)vector<int>().swap(hyperT);
        hyperGT.clear();

        for (auto& hyper : hyperG_R2)hyper.clear();
        for (auto& hyperT : hyperGT_R2)vector<int>().swap(hyperT);
        hyperGT_R2.clear();

        seedSet.clear();
        candidateSetT.clear();
        allocQueueM.clear();
        isSelect.clear();
        usersExpand.clear();
        criterQueue.clear();

        EstPro1=1.0, EstPro2=0, EstPro3=0, eps1=0., eps2=0.;
        /*adopt_set = vector<int> (n, -1);
        max_weight = vector<vector<double>> (n);
        for (int i = 0; i < n; i++){
            max_weight[i] = vector<double> (nrCompanies, 0.0);
        }*/

        Adoption = vector<adoptInfo> (max_n);
        for (int i = 0; i < max_n; i ++){
            Adoption[i].advMaxId = -1; Adoption[i].Maxweight = 0.0;
        }

        for (int i = 0; i < nrCompanies; i ++){ //clear adv seedset
            company *adv;
            adv = advList.at(i);

            adv->seedSet.clear();
            adv->seedSet.shrink_to_fit();

            adv->num_coveredRR.clear();
            adv->num_coveredRR2.clear();
            adv->totalProfit = 0.0; adv->totalRevenue = 0.0; adv->totalSeedCosts = 0.0;
            adv->currentInf = 0.0; adv->totalInf = 0.0; adv->isFull = 0;
        }
    }


#include "discrete_rrset.h"

    void buildhypergraph(unsigned int R)
    {
        if (R > INT_MAX){
            cout << "Error:R too large" << endl;
            exit(1);
        }
        //ASSERT(hyperGT.size() == 0);
        unsigned int size = hyperGT.size();
        while (size < R)BuildHypergraphNode(size++, max_n, nodes);
    }
/*
    void buildhypergraphR_R2(unsigned int R)
    {
        if (R > INT_MAX){
            cout << "Error:R too large" << endl;
            exit(1);
        }
        //ASSERT(hyperGT.size() == 0);
        unsigned int size = hyperGT_R2.size();
        while (size < R)BuildHypergraphNode(size++ ,hyperG_R2, hyperGT_R2);
    }
*/
    void buildhypergraph_R2(unsigned int R, companyList advList, vector<bool> &edgeMarkR2)
    {
        if (R > INT_MAX){
            cout << "Error:R2 too large" << endl;
            exit(1);
        }
        //ASSERT(hyperGT_R2.size() == 0);
        unsigned int size = hyperGT_R2.size();
        while (size < R) { //size: hyperR2Id
            //BuildHypergraphNode(size++ ,hyperG_R2, hyperGT_R2);
            //BuildHypergraphNodeSingle(advList,hyperG_R2, hyperGT_R2);
            //cout << "current rr set id: " << size<< endl;
            unsigned int hyperId = size;
            edgeMarkR2.push_back(true);

            int uStart;
            int flag = 0;
            set<int>::iterator iter;
            while (flag == 0){
                uStart = sfmt_genrand_uint32(&sfmtSeed) % max_n;
                iter = nodes.find(uStart);
                if(iter != nodes.end()) flag = 1;
                else continue;
            }
            //cout << "current uStart: " <<uStart<<endl;

            visit.clear();visit_mark.clear();
            visit = vector<bool> (max_n);
            visit_mark = vector<int> (max_n);
            company *adv;
            // step1: check whether uStart is a sampled node
            bool flag1 = false;
            for (int i = 0; i < nrCompanies; i ++){
                if (flag1) break;
                adv = advList.at(i);
                for (int j = 0; j <adv->seedSet.size(); j++){
                    if (uStart == adv->seedSet[j]){
                        adv->num_coveredRR2.insert(pair<int, int>((int)hyperId, 1));
                        hyperG_R2[uStart].push_back(hyperId);
                        vector<int> temp_visit_mark; temp_visit_mark.push_back(uStart);
                        hyperGT_R2.push_back(temp_visit_mark);
                        //cout << "[uStart is seed] current seed is: " << uStart << "\tadv is: " << adv->companyID <<endl;
                        flag1 = true;
                        break;
                    }
                }
            }
            if (flag1) {size ++; continue;}
            // step2: uStart is not a sampled node, begin  to generate rr set
            unsigned int n_visit_mark = 0, curIdx = 0;
            visit_mark[n_visit_mark++] = uStart;
            visit[uStart] = true;
            hyperG_R2[uStart].push_back(hyperId);
            bool flag2=false;

            if (influModel == LT)
            {
                while (curIdx < n_visit_mark)
                {
                    int expand = visit_mark[curIdx++];

                    if (gT[expand].size() == 0)
                        continue;

                    int index = sfmt_genrand_uint32(&sfmtSeed) % gT[expand].size();
                    int v = gT[expand][index];
                    if (visit[v])continue;
                    for (int i = 0; i < nrCompanies; i ++){ // find seed of adv
                        if(flag2) break;
                        adv = advList.at(i);
                        for (int j = 0; j <adv->seedSet.size(); j++){
                            if (v == adv->seedSet[j]){
                                //cout << "[v is seed] current seed is: " << v << "\tadv is: " << adv->companyID <<endl;
                                candidateNode = v; // v is seed
                                flag2 = true;
                                break;
                            }
                        }
                    }
                    visit[v] = true;
                    visit_mark[n_visit_mark++] = v;
                    hyperG_R2[v].push_back(hyperId);
                    //++spread[v];
                    if(flag2) break;
                }
            }
            else
            ASSERT(false);

            if (!flag2){
                for (unsigned int i = 0; i < n_visit_mark; i++){
                    visit[visit_mark[i]] = false;
                }
                hyperGT_R2.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark));
            }
            else{
                int coveredRR2_count = 0;
                adv->num_coveredRR2.insert(pair<int, int>((int)hyperId, 0));
                for (unsigned int i = 0; i < n_visit_mark; i++){
                    // step3: assign node v before seed v company product
                    coveredRR2_count += estimateInfnode(candidateNode, visit_mark[i], adv, edgeMarkR2, 2);
                    // step4: reset visit[node]
                    visit[visit_mark[i]] = false;
                }
                adv->num_coveredRR2[hyperId] = coveredRR2_count;
                hyperGT_R2.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark));
                //if(coveredRR2_count>hyperGT_R2[hyperId].size())cout << "[adv " << adv->companyID << "]\tcoveredRRnum: " << coveredRR2_count << "\trr set size: " <<hyperGT_R2[hyperId].size() <<endl;
            }
            size ++;
        }

    }


////-------------------------------------------------------- basic function begin -----------------------------------------------------

    //// this function select company with the max weight adopted by v's active in-neighbor (nnbr/node(seed))
    int search_maxweight_comy (int v, int currentID){
        unordered_map<int,double> com_weight;
        int flag = 0;
        for (int i = 0; i < nrCompanies; i ++)
            com_weight[i] = 0.0;

        for (int j = 0; j < gT[v].size(); j ++){
            int innbr = gT[v][j];
            //cout << "[search max_weight adv]innbr: " << innbr <<" adopt " << Adoption[innbr].advMaxId<< "\ttime_tag: " << Adoption[innbr].time_tag << "\t current time: " << timeCount << endl;
            if (Adoption[innbr].advMaxId != -1 ) //node v's in-neighbor actived[is active] in last time_tag
            {
                //cout << "[yes]innbr "<< innbr << " adopt " << Adoption[innbr].advMaxId <<"  weight " << probT[v][innbr] << endl;
                com_weight[Adoption[innbr].advMaxId] += probT[v][j];
                flag = 1; // node v has active in-neighbor
            }
        }
        if (!flag)  return currentID;

        //for (int i = 0; i < nrCompanies; i ++) cout << "company " << i << " weight is: " << com_weight[i] << endl;
        std::pair<int,double> max = *max_element(com_weight.begin(), com_weight.end(), [](const std::pair<int, double> &p1, const std::pair<int, double> &p2){ return p1.second < p2.second;});
        int maxcomId = max.first;
        Adoption[v].advMaxId = maxcomId;
        Adoption[v].Maxweight = max.second;

        return maxcomId;
    }

    //// this function calculate the sum weight of current with v's active in-neighbor
    double cal_sum_weight (int v, int currComID){
        double sum_weight = 0.0;
        company *curr_comy;
        for (int k = 0; k < advList.size(); k ++) {
            if (advList.at(k)->companyID == currComID) {
                curr_comy = advList.at(k);
                break;
            }
        }
        for (int j = 0; j < gT[v].size(); j ++){
            int innbr = gT[v][j];
            //cout << "[calculate weight]innbr "<< innbr << " adopt " << Adoption[innbr].advMaxId <<"  time  " << Adoption[innbr].time_tag <<"  current time  " << timeCount<< endl;
            if (Adoption[innbr].advMaxId == curr_comy->companyID) //node v's active in-neighbor adopted current Company
                sum_weight += probT[v][j];
        }
        return sum_weight;
    }

    //// [R1] this function decrease num of influenced node in rr set that old company has covered
    void decrease_infRR (int node, int oldComID, vector<bool> &RRedgeMark){
        company *old_company;
        for (int k = 0; k < advList.size(); k ++){// step1: find the old_company
            if (advList.at(k)->companyID == oldComID){
                old_company = advList.at(k);
                break;
            }
        }
        for (auto idr : hyperG[node]){
            if (!RRedgeMark[idr])   continue;
            map<int, int>::iterator iter;
            iter = old_company->num_coveredRR.find(idr); // find the rr set covered by old company that contained node
            if(iter != old_company->num_coveredRR.end()){
                old_company->num_coveredRR[idr] = iter->second -1;
                if(old_company->num_coveredRR[idr] < 0) old_company->num_coveredRR[idr]=0;
                //cout <<"update old company"" [" <<idr<<"] number: " << old_company->num_coveredRR[idr] << endl;
            }

            else continue;
        }
    }

    //// [R1] this function increase num of new influenced node in maxComID's rr set (maxComID: most of node in-neightbor adopted)
    void increase_infRR (int node, int maxComID, vector<bool> &RRedgeMark){
        company *max_comy;
        for (int k = 0; k < advList.size(); k ++) {
            if (advList.at(k)->companyID == maxComID) {
                max_comy = advList.at(k);
                break;
            }
        }
        for (auto idr : hyperG[node]){
            if (!RRedgeMark[idr])   continue;
            map<int, int>::iterator iter;
            iter = max_comy->num_coveredRR.find(idr); // find the rr set covered by old company that contained node
            if(iter != max_comy->num_coveredRR.end())
                max_comy->num_coveredRR[idr] = iter->second +1;
            if(max_comy->num_coveredRR[idr] > hyperGT[idr].size()) max_comy->num_coveredRR[idr]=hyperGT[idr].size();
        }
    }

    //// [R2] this function decrease num of influenced node in rr set that old company has covered
    void decrease_infRR2 (int node, int oldComID, vector<bool> &RRedgeMark){
        company *old_company;
        for (int k = 0; k < advList.size(); k ++){
            if (advList.at(k)->companyID == oldComID){
                old_company = advList.at(k);
                break;
            }
        }
        for (auto idr : hyperG_R2[node]){
            if (!RRedgeMark[idr])   continue;
            map<int, int>::iterator iter;
            iter = old_company->num_coveredRR2.find(idr); // find the rr set covered by old company that contained node
            if(iter != old_company->num_coveredRR2.end()){
                old_company->num_coveredRR2[idr] = iter->second -1;
                if(old_company->num_coveredRR2[idr] < 0) old_company->num_coveredRR2[idr]=0;
                //cout <<"update old company"" [" <<idr<<"] number: " << old_company->num_coveredRR2[idr] << endl;
            }

            else continue;
        }
    }

    //// [R2] this function increase num of new influenced node in maxComID's rr set (R2) (maxComID: most of node in-neightbor adopted)
    void increase_infRR2 (int node, int maxComID, vector<bool> &RRedgeMark){
        company *max_comy;
        for (int k = 0; k < advList.size(); k ++) {
            if (advList.at(k)->companyID == maxComID) {
                max_comy = advList.at(k);
                break;
            }
        }
        for (auto idr : hyperG_R2[node]){
            if (!RRedgeMark[idr])   continue;
            map<int, int>::iterator iter;
            iter = max_comy->num_coveredRR2.find(idr); // find the rr set covered by old company that contained node
            if(iter != max_comy->num_coveredRR2.end())
                max_comy->num_coveredRR2[idr] = iter->second +1;
            if(max_comy->num_coveredRR2[idr] > hyperGT_R2[idr].size()) max_comy->num_coveredRR2[idr]=hyperGT_R2[idr].size();
        }
    }

////----------------------------------------------------------------------------------------------------------

    //// this function estimate influenced node in a rr set covered by a new selected seed
    // node: seed  it: user
    int estimateInfnode (int node, int it, company *adv, vector<bool> &edgeMark, int RR){
        int currentRR_adopt = 0;
        //// step2: v is activated by node and adopt company based the maximum probability of v's in-neighbors
        // Not switch: [it] is seed (including current company and other companies) and already adopted current company
        //// if [it] is already another company's seed, it cannot switch: isSelect[it] = 0
        if ((!isSelect[it]) && it == node) {currentRR_adopt++; return currentRR_adopt;}
        if ((!isSelect[it]) && it != node)  return currentRR_adopt;
        if (it == node){ // seed must adopt current company
            isSelect[it] = 0;
            if (Adoption[it].advMaxId == adv->companyID) {currentRR_adopt ++; }
            else if (Adoption[it].advMaxId == -1)
            {   Adoption[it].advMaxId = adv->companyID; Adoption[it].Maxweight = cal_sum_weight(it, adv->companyID);
                currentRR_adopt ++;}
            else{ // seed must switch to current company
                int oldID = Adoption[it].advMaxId;
                if (RR == 1) decrease_infRR (it, oldID, edgeMark); // update old_company's totalInf
                else    decrease_infRR2 (it, oldID, edgeMark);
                Adoption[it].advMaxId = adv->companyID;
                Adoption[it].Maxweight = cal_sum_weight(it, adv->companyID);
                currentRR_adopt ++;
                //cout << "seed " <<it <<" switch " << "old company "<< oldID << " ===> " << "new company "<< adv->companyID << endl;
            }
            return currentRR_adopt;
        }

        if (Adoption[it].advMaxId == adv->companyID) {currentRR_adopt ++; return currentRR_adopt;}

        // active + adopt (max influence probability company adopted in last time by it's in-neighbours)
        if ( Adoption[it].advMaxId == -1){
            int adoptId = search_maxweight_comy(it,adv->companyID);
            Adoption[it].advMaxId = adoptId;
            Adoption[it].Maxweight = cal_sum_weight(it, adoptId);
            //cout << "max weight: " << Adoption[it].Maxweight<<"\tmax weight company ID is: " << adoptId  <<"\t current company ID is: " << adv->companyID<< endl;

            if (adoptId != adv->companyID){
                if (RR == 1) increase_infRR(it, adoptId, edgeMark);
                else  increase_infRR2(it, adoptId, edgeMark);
            }
            else
                currentRR_adopt ++;
        }
        //// step3: switch
        else{
            int oldID = Adoption[it].advMaxId;
            company *old_company;
            for (int k = 0; k < advList.size(); k ++){
                if (advList.at(k)->companyID == oldID){
                    old_company = advList.at(k);
                    break;
                }
            }
            // old_company = advList.at(oldID);//不能这样访问，因为advList已经根据bpi sort过
            // switch condition 1: current company's sum weight of adopted by active neighbor of nbr > max_weight [v][oldID]
            // switch condition 2: old adopted company Bi/Ii < current company Bj/Ij
            bool condition1 = false, condition2 = false, condition3 = false;
            double current_weight = cal_sum_weight(it, adv->companyID);
            //cout <<"maybe switch node "<<it<< "  old weight: " << Adoption[it].Maxweight << "\tnew weight: " << current_weight << "  old bpi: " << old_company->bpi << "\tnew bpi: " << adv->bpi << endl;
            if (Adoption[it].Maxweight < current_weight && old_company->bpi < adv->bpi)  {condition1 = true;}
            if (Adoption[it].Maxweight < current_weight && old_company->bpi == adv->bpi) {condition2 = true;}
            if (Adoption[it].Maxweight == current_weight && old_company->bpi < adv->bpi) {condition3 = true;}
            if (condition1 || condition2 || condition3){
                Adoption[it].advMaxId = adv->companyID;
                Adoption[it].Maxweight = current_weight;
                currentRR_adopt ++;
                if (RR == 1) decrease_infRR (it, oldID, edgeMark); // update old_company's totalInf
                else    decrease_infRR2 (it, oldID, edgeMark);
                //cout << "switch node: " << it <<"\t"<<oldID << " ===> " << adv->companyID << endl;
            }
        }//end else
        return currentRR_adopt;
    }

    //// this function estimate (upper) bound of total profit (R1) for all companies based on their seed set
    double estimateProfit (companyList advList, unsigned int R, double gammaR, double gammaP){
        company *adv;
        double estiProfit = 0.0;
        int total_seed = 0;
        for (int i = 0; i < nrCompanies; i ++){
            adv = advList.at(i);
            // step1: estimate total influence
            double spreadRR = 0.0;
            map<int, int>::iterator iter;
            cout << "adv "<< adv->companyID << " seed"<<"["<<adv->seedSet.size()<<"]: " <<"{ " ;
            for (int seed : adv->seedSet){
                cout <<  seed << "," ;
                total_seed ++;
            }cout <<" }"<< endl;
            for(iter = adv->num_coveredRR.begin(); iter != adv->num_coveredRR.end(); iter++)
            {
                spreadRR += iter->second * 1.0 / hyperGT[iter->first].size(); //spreadRR==>spread[node]

                //cout << "adv "<<adv->companyID<<"  covered rrId is: " << iter->first << "  covered ratio is: " << spreadRR << "  covered num is: " << iter->second
                //    << "  rr size: " <<hyperGT[iter->first].size() << endl;

            }
            adv->totalInf = spreadRR /  R * 1.0 * n;
            cout << "adv total influence spread :" << adv->totalInf << " influence threshold: " << adv->influence << endl;
            // step2: estimate total profit
            if (adv->totalInf < adv->influence) //// host: penalty
                adv->totalRevenue = adv->budget * (1.0 + gammaP * ((adv->totalInf - adv->influence) / adv->influence));
            else //// host: reward
                adv->totalRevenue = adv->budget * (1.0 + gammaR * ((adv->totalInf - adv->influence) / adv->influence));

            ////20230106: in fact, this case not happens in obo and iter
            //if (adv->totalInf == 0) {adv->totalRevenue = 0.0; adv->totalProfit = 0.0;};

            adv->totalSeedCosts = 0.0;
            for (int seed : adv->seedSet){
                adv->totalSeedCosts += Cost[seed];
            }
            adv->totalProfit = adv->totalRevenue - adv->totalSeedCosts;
            cout << "adv total Profit :" << adv->totalProfit <<"\tadv total revenue :" << adv->totalRevenue << "\tadv total SeedCosts :" << adv->totalSeedCosts << endl << endl;
            estiProfit += adv->totalProfit;
        }
        cout << "total seed node :" << total_seed << endl;
        return estiProfit;
    }

    //// [R1] this function estimate (upper) bound of total revenue (R1) for all companies based on their seed set
    double estimateRevenue (companyList advList, unsigned int R, double gammaR, double gammaP){
        company *adv;
        double EstRevR1 = 0.0;
        int total_seed = 0;
        for (int i = 0; i < nrCompanies; i ++){
            adv = advList.at(i);
            // step1: estimate total influence
            double spreadRR = 0.0;
            map<int, int>::iterator iter;
            total_seed += adv->seedSet.size();
            for(iter = adv->num_coveredRR.begin(); iter != adv->num_coveredRR.end(); iter++)
            {
                spreadRR += iter->second * 1.0 / hyperGT[iter->first].size(); //spreadRR==>spread[node]
                //cout << "adv "<<adv->companyID<<"  covered rrId is: " << iter->first << "  covered ratio is: " << spreadRR << "  covered num is: " << iter->second
                //   << "  rr size: " <<hyperGT[iter->first].size() << endl;
            }
            adv->totalInf = spreadRR /  R * 1.0 * n;
            cout << "adv "<< adv->companyID << " total influence spread :" << adv->totalInf << "\tinfluence threshold: " << adv->influence;
            // step2: estimate total profit
            if (adv->totalInf < adv->influence) //// host: penalty
                adv->totalRevenue = adv->budget * (1.0 + gammaP * ((adv->totalInf - adv->influence) / adv->influence));
            else //// host: reward
                adv->totalRevenue = adv->budget * (1.0 + gammaR * ((adv->totalInf - adv->influence) / adv->influence));

            ////20230106: in fact, this case not happens in obo and iter
            //if (adv->totalInf == 0) {adv->totalRevenue = 0.0; adv->totalProfit = 0.0;};

            ////20230105 revenue=max{0, r}
            //if (adv->totalRevenue<0) adv->totalRevenue=0;

            cout << "\tadv total revenue :" << adv->totalRevenue << endl;
            EstRevR1 += adv->totalRevenue;
        }
        cout << "total seed node :" << total_seed << endl;
        return EstRevR1;
    }

    //// this function estimate (upper) bound of total seed cost (R1) for all companies based on their seed set
    double estimateCost (companyList advList){
        company *adv;
        double EstCost1 = 0.0;
        for (int i = 0; i < nrCompanies; i ++){
            adv = advList.at(i);

            adv->totalSeedCosts = 0.0;
            cout << "adv "<< adv->companyID << " seed"<<"["<<adv->seedSet.size()<<"]: " <<"{ " ;
            for (int seed : adv->seedSet){
                adv->totalSeedCosts += Cost[seed];
                cout <<  seed << "," ;
            }cout <<" }"<< endl;
            adv->totalProfit = adv->totalRevenue - adv->totalSeedCosts;
            cout << "adv total Profit :" << adv->totalProfit <<"\tadv total revenue :" << adv->totalRevenue << "\tadv total SeedCosts :" << adv->totalSeedCosts << endl;
            EstCost1 += adv->totalSeedCosts;
        }
        return EstCost1;
    }

    //// [R2] this function estimate (lower) bound of total revenue [R2] for all companies based on their seed set
    double estimateRevenueR2 (companyList advList, unsigned int R, double gammaR, double gammaP){
        company *adv;
        double EstRevR2 = 0.0;
        int total_seed = 0;
        for (int i = 0; i < nrCompanies; i ++){
            adv = advList.at(i);
            /// step1: estimate total influence
            double spreadRR = 0.0;
            map<int, int>::iterator iter;
            total_seed += adv->seedSet.size();
            for(iter = adv->num_coveredRR2.begin(); iter != adv->num_coveredRR2.end(); iter++)
            {
                auto rrid = iter->first;
                spreadRR += adv->num_coveredRR2[rrid] * 1.0 / hyperGT_R2[rrid].size(); //spreadRR==>spread[node]

                // cout << "adv "<<adv->companyID<<"  covered rrId is: " << iter->first << "  covered ratio is: " << spreadRR << "  covered num is: " << iter->second
                //     << "  rr size: " <<hyperGT_R2[iter->first].size() << endl;
            }
            adv->totalInf = spreadRR /  R * 1.0 * n;
            cout << "[R2]adv "<< adv->companyID << " total influence spread :" << adv->totalInf << "\tinfluence threshold: " << adv->influence;

            /// step2: estimate total profit
            if (adv->totalInf < adv->influence) //// host: penalty
                adv->totalRevenue = adv->budget * (1.0 + gammaP * ((adv->totalInf - adv->influence) / adv->influence));
            else //// host: reward
                adv->totalRevenue = adv->budget * (1.0 + gammaR * ((adv->totalInf - adv->influence) / adv->influence));

            ////20230106: in fact, this case not happens in obo and iter
            //if (adv->totalInf == 0) {adv->totalRevenue = 0.0; adv->totalProfit = 0.0;};

            cout << "\tadv total revenue :" << adv->totalRevenue << endl;
            EstRevR2 += adv->totalRevenue;

            /// step3: estimate total cost
            double singlecost =0.0;
            for (auto it : adv->seedSet)   singlecost += Cost[it];
            adv->totalSeedCosts=singlecost;

            /// step4: estimate total profit
            adv->totalProfit = adv->totalRevenue-adv->totalSeedCosts;
        }
        cout << "total seed node :" << total_seed << endl;
        return EstRevR2;
    }

    //// this function estimate current influence spread for company based on its current seed set
    bool estimateCurrInf (company *adv, unsigned int R){
        double currSpreadRR = 0.0;
        map<int, int>::iterator iter;
        /*cout << "current adv "<< adv->companyID << " seed: " <<"{ " ;
        for (int seed : adv->seedSet){
                cout <<  seed << "," ;
        }cout <<" }"<< endl;*/
        for(iter = adv->num_coveredRR.begin(); iter != adv->num_coveredRR.end(); iter++)
        {
            currSpreadRR += iter->second * 1.0 / hyperGT[iter->first].size(); //spreadRR==>spread[node]
            //cout << "adv "<<adv->companyID<<"  covered rrId is: " << iter->first << "  covered ratio is: " << spreadRR << "  covered num is: " << iter->second
            //    << "  rr size: " <<hyperGT[iter->first].size() << endl;
        }
        adv->currentInf = currSpreadRR /  R * 1.0 * n;
        //cout << "[adv " <<adv->companyID<<" current influence spread] :" << adv->currentInf << " influence threshold: " << adv->influence << endl;
        if (adv->currentInf > 1.0 * adv->influence) return true;
        else   return false;
    }

    //// this function estimate lower bound of total profit (R2) for all companies based selected seed set
    void updateInfoR2 (companyList advList,double gammaP){
        //step1: generate rr set (R): hyperG_R2, hyperGT_R2 in main()
        //step2: initial advertiser covered information
        Adoption = vector<adoptInfo> (max_n);
        for (int i = 0; i < max_n; i ++){
            Adoption[i].advMaxId = -1; Adoption[i].Maxweight = 0.0;
        }
        isSelect = std::vector<int>(max_n,1); // only disjointness

        company *adv;
        for (int i = 0; i < nrCompanies; i ++){ //clear adv unless seedset
            adv = advList.at(i);
            adv->num_coveredRR.clear();
            adv->totalProfit = 0.0; adv->totalRevenue = 0.0; adv->totalSeedCosts = 0.0;
            adv->currentInf = 0.0; adv->totalInf = 0.0; adv->isFull = 0;

            adv->bpi = (float) adv->budget/adv->influence*1.;
            adv->bpi = adv->bpi * gammaP;

            int currSeed;
            for (int j = 0; j <adv->seedSet.size(); j++){
                currSeed = adv->seedSet[j];
                isSelect[currSeed] = 0;
                Adoption[currSeed].advMaxId=adv->companyID;
            }
        }

    }


////----------------------------------------------------------------------------------------------------------

    //// this function returns the best only influence-sensitive node (maxVal node, minInf adv)
    bool selectBestISNode (infPair &bestIS, vector<int> &hyper_degree, unsigned int R){

        double mpg;
        int id = -1;
        double maxVal = 0.0;

        // select the company with minimal influence spread (balance)
        company *adv;
        int minInfId = -1;
        double currMinInf_ratio = 1.0;
        double currInf;
        for (int i = 0; i < nrCompanies; i ++){
            adv = advList.at(i);
            //estimate current total influence
            double currSpreadRR = 0.0;
            map<int, int>::iterator iter;
            for(iter = adv->num_coveredRR.begin(); iter != adv->num_coveredRR.end(); iter++)
            {
                currSpreadRR += iter->second * 1.0 / hyperGT[iter->first].size(); //spreadRR==>spread[node]
            }
            currInf = currSpreadRR /  R * 1.0 * n;
            if(currInf >= adv->influence)    continue;
            if(currInf/ (1.0 *adv->influence) < currMinInf_ratio){
                currMinInf_ratio = currInf/ (1.0 *adv->influence) ;
                minInfId = adv->companyID;
            }
            //cout << "current acompany: " << adv->companyID << " current influence: " << currInf << " influence threshold: " << adv->influence << endl;
        }
        if (minInfId == -1) // no company's current influence spread is smaller than threshold
            return false;

        for(int node = 0; node < (int) hyper_degree.size(); node++) { // "node" : node id
            if((isSelect[node] == 1) && (maxVal < hyper_degree[node])) {
                maxVal = hyper_degree[node]; //hyper_degree[] is updated
                id = node;
            }
        }

        adv = advList.at(minInfId);
        mpg = adv->bpi * (( double) n * (( double) maxVal / R)) / (1.0 * Cost[id]);
        this->candidateNode = id;
        bestIS.adv = adv->companyID; bestIS.node = id; bestIS.mpg = mpg;
        //bestIS = std::make_pair(id,mpg);
        //hyper_degree[id]=0;  // assign zero since either it will be allocated soon or will be out of the game due to attention bound
        return true;
    }

    //// this function returns the best only profit-sensitive node (marginal revenue gain/cost)
    bool selectBestPSNode (infPair &bestPS, vector<int> hyper_degree, vector<infPair> &criterQueue, unsigned int R){

        if(criterQueue.empty()) return false;
        //multimap<double,int>::iterator it = criterQueue.end(); //it 指向 multimap 最后一个元素的后一个位置
        //it--;
        sort (criterQueue.begin(), criterQueue.end(), cmp);
        infPair it = criterQueue[criterQueue.size()-1]; // it 指向 criterQueue 最后一个元素（最大mpg）

        while (true) {
            while(isSelect[it.node] == 0) {//cout << "criterQueue erase node: " << it.node << "\terase adv: " << it.adv<< endl;
                criterQueue.pop_back();
                if(criterQueue.empty()) return false;
                //it = criterQueue.end();
                //it--;
                it = criterQueue[criterQueue.size()-1];
            }
            //// it->second(user id) is not in usersExpand
            if(usersExpand.find(it.node) == usersExpand.end()) { //this node is not explored after assignBestNode (ie. hyper_degree unchanged)

                this->candidateNode = it.node;
                bestPS.adv = it.adv; bestPS.node = it.node; bestPS.mpg = it.mpg;
                //bestPS = std::make_pair(this->candidateNode,it->first); // should contain the ratio for CS
                ////hyper_degree[candidateNode]=0;  // assign zero since either it will be allocated soon or will be out of the game due to
                //criterQueue.erase(it);
                criterQueue.pop_back();
                return true;
                //break;
            }
            //// it->second(user id) in usersExamined (update hyper_degree[user id])
            else {
                int idTemp = it.node;
                company *adv;
                for(vector<infPair>::iterator it_temp=criterQueue.begin();it_temp!=criterQueue.end();it_temp++){
                    if(it_temp->node == idTemp){
                        adv = advList.at(it_temp->adv);
                        double mpg_new = (( double) n * (( double) hyper_degree[idTemp] / R)) * adv->bpi / (1.0 * Cost[idTemp]);//hyper_degree[] has updated
                        //cout << "[update margin profit gain] adv: "<<adv->companyID << "\tnode id: " << idTemp <<"\told gain: "<<it_temp->mpg <<"\t new gain: " <<mpg_new << endl;
                        it_temp->mpg = mpg_new;// update all elements in M'
                    }
                }
                //criterQueue.erase(it++);
                usersExpand.erase(idTemp);//已完成更新，则从待更新结点set中删除该结点
                //cout << "[update margin profit gain] ended! " << endl;
                sort (criterQueue.begin(), criterQueue.end(), cmp);
                it = criterQueue[criterQueue.size()-1];
                //criterQueue.insert(pair<double,int>(mpg_new ,idTemp));
                //it = criterQueue.end();
                //it--;
            }
        }

    }

    //// this function update hyperG/hyperGT and adv_information after selecting a new seed
    void assignBestNode (company *adv, vector<int> &hyper_degree, vector<bool> &edgeMark, unsigned int R) {
        // step1: update hyperG/hyperGT
        isSelect[candidateNode] = 0;
        hyper_degree[candidateNode] = 0;
        for (auto Idx : hyperG[candidateNode])  // Idx: rr set
        {//Remove from R all rr sets that are covered by node

            if (edgeMark[Idx])continue;
            adv->num_coveredRR.insert(pair<int, int>((int)Idx, 0));
            int coveredRR_count = 0;
            for (auto it : hyperGT[Idx]){ // it: node id
                --hyper_degree[it]; ////
                usersExpand.insert(it);
                if(hyper_degree[it] < 0) {
                    hyper_degree[it] = 0;
                }
                coveredRR_count += estimateInfnode(candidateNode, it, adv, edgeMark, 1);
            }
            edgeMark[Idx] = true; // best node covers it
            adv->num_coveredRR[Idx] = coveredRR_count;

        }
        // step2: update adv_information
        adv->seedSet.push_back(candidateNode);
    }



////-------------------------------------------------------- basic function end -----------------------------------------------------

////======================================================= main function  begin ====================================================

    // this function needs to be modified. We could directly minus cost from the spread.
    // And then select the node with the largest profit.
    //// ROI-Greedy
    bool build_seedset_ROI(unsigned int R, unsigned int cnt, double epsilon, double delta, double &estprofit)
    {
        EstPro1=1.0, EstPro2=0, EstPro3=0, eps1=0., eps2=0.;

        seedSet.clear();
        double cov=0, CstS=0;

        vector<bool>tag(max_n, false);
        vector<int>inactive(max_n, 0);
        iota(inactive.begin(), inactive.end(), 0); //initialize from 0 to n;
        vector<int>spread(max_n, 0);
        for (unsigned int k = 0; k < max_n; k++)spread[k] = (int)hyperG[k].size();
        long long numEdge = hyperGT.size();
        vector<bool> edgeMark(numEdge, false);

        if(true){ //algo>0 ROI greedy

            unordered_map<int,double> sc_map;
            for(unsigned int it=0; it<max_n; it++){
                ASSERT(spread[it] >= 0);
                sc_map[it]=spread[it] * 1.0 /  Cost[it];
            }

            while (!sc_map.empty())
            {
                std::pair<int,double> top = *max_element(sc_map.begin(), sc_map.end(), [](const std::pair<int, double> &p1, const std::pair<int, double> &p2){ return p1.second < p2.second;});
                int node = top.first;
                if (spread[node] * 1.0 / R*n <= Cost[node])break;
                ASSERT(node >= 0);
                seedSet.push_back(node);
                sc_map.erase(node);
                cov+=spread[node];
                CstS+=Cost[node];
                tag[node] = true;
                for (auto Idx : hyperG[node])  // Idx: rr set
                {//Remove from R all rr sets that are covered by node

                    if (edgeMark[Idx])continue;
                    for (auto it : hyperGT[Idx]){ // it: node id
                        --spread[it];
                        if(!tag[it]) sc_map[it] = spread[it] * 1.0 /  Cost[it];
                    }
                    edgeMark[Idx] = true;
                }
            }
        }

        EstPro1=cov/R*n-CstS;
        unsigned int i=0;
        while(i++<R)EstPro2+=BuildHypergraphNodeSingle(tag);
        double x = (EstPro2/R*n)*R/log(6.*cnt*cnt/delta)/n;
        eps1=4./(sqrt(1.+8.*x)-3);
        //ASSERT(eps1>0);
        eps2=sqrt((2.*eps1+2.)/x);
        //ASSERT(eps2>0);
        EstPro2=EstPro2/R*n-CstS;
        double t=EstPro1/EstPro2;
        EstPro3=EstPro2/R*n-(1+eps1)*CstS;
        //if(1-(1-eps2)/(t+t*eps1)<=epsilon && eps2<=epsilon)return true;
        estprofit = EstPro2;
        if( (t-1.)/t+eps1+eps2 <=epsilon && eps1+eps2<=epsilon) return true;
        return false;

    }



////******************************* multiple companies ****************************************

    //// this function is candidate-select method: select candidate seedset T for allocation
    //// note that candidateSetT is the maximal union seedset of all the companies
    void candidate_greedy(unsigned int R)
    {
        candidateSetT.clear();
        vector<bool>tag(max_n, false);
        isSelect = std::vector<int>(max_n,1); // only disjointness
        vector<int>spread(max_n, 0);
        for (unsigned int k = 0; k < max_n; k++)spread[k] = (int)hyperG[k].size(); /////!!!!!
        long long numEdge = hyperGT.size();
        vector<bool> edgeMark(numEdge, false);
        double max_bpi = 0.0; // the most relaxed conditions
        company *maxComy;
        for (int i = 0; i < nrCompanies; i ++){
            maxComy = advList.at(i);
            maxComy->bpi = maxComy->budget / maxComy->influence * 1.0;
            if (max_bpi <= maxComy->bpi)
                max_bpi = maxComy->bpi;
        }
        //ROI greedy
        unordered_map<int,double> sc_map;
        for(unsigned int it=0; it<max_n; it++){
            ASSERT(spread[it] >= 0);
            sc_map[it]=spread[it] * 1.0 /  Cost[it];
        }
        while (!sc_map.empty())
        {
            std::pair<int,double> top = *max_element(sc_map.begin(), sc_map.end(), [](const std::pair<int, double> &p1, const std::pair<int, double> &p2){ return p1.second < p2.second;});
            int node = top.first;
            //cout << "candidateT node: " << node << "\tmarginal profit gain: " << top.second <<endl;
            if (spread[node] * 1.0 / R*n <= Cost[node]/(1.0 * max_bpi))break; //// focus point
            ASSERT(node >= 0);
            candidateSetT.push_back(node);
            sc_map.erase(node);
            tag[node] = true;
            for (auto Idx : hyperG[node])  // Idx: rr set
            {
                //Remove from R all rr sets that are covered by node
                if (edgeMark[Idx])continue;
                for (auto it : hyperGT[Idx]){ // it: node id
                    --spread[it];
                    if(!tag[it]) sc_map[it] = spread[it] * 1.0 /  Cost[it];
                }
                edgeMark[Idx] = true;
            }
        }
        cout << "candidate seed set T ["<< candidateSetT.size()<<"]: {" ;
        for (int candidateSeed: candidateSetT ){
            cout << candidateSeed << ", ";
        }cout << "}"<< endl;
    }

    //// 1. this function is one-by-one method: select seedset for multi-companies based on ROI-Greedy
    void build_multi_seedsets_obo(unsigned int R, double gammaR, double gammaP)
    {
        cout << "one-by-one method select seedsets..." << endl;
        companyList advList_sort;
        advList_sort.assign(advList.begin(),advList.end());
        sort(advList_sort.begin(),advList_sort.end(),comparison); //// order companies based on decreasing order of Bi/Ii

        for (int i = 0; i <nrCompanies; i++){
            company *adv = advList.at(i);
            adv->bpi = gammaP * adv->budget / adv->influence *1.0;
        }
        vector<bool>tag(max_n, false); //seed tag
        vector<int>spread(max_n, 0);
        isSelect = std::vector<int>(max_n,1); // only disjointness
        for (unsigned int k = 0; k < max_n; k++)spread[k] = (int)hyperG[k].size();
        long long numEdge = hyperGT.size();
        vector<bool> edgeMark(numEdge, false);

        int flag = nrCompanies; // initial nrCompanies, set 0 when no user satisfy any company
        company *adv;

        //ROI greedy
        unordered_map<int,double> sc_map;
        for(unsigned int it=0; it<max_n; it++){
            ASSERT(spread[it] >= 0);
            sc_map[it]=spread[it] * 1.0 /  Cost[it];
        }
        vector<bool> overInf(nrCompanies, false);
        while (!sc_map.empty() && flag)
        {
            for (int i = 0; i < nrCompanies; i ++) //select seed for companies one by one
            {
                adv = advList_sort.at(i);
                if (adv->isFull) continue; // there is no node satisfy adv's profit > 0 condition
                std::pair<int,double> top = *max_element(sc_map.begin(), sc_map.end(), [](const std::pair<int, double> &p1, const std::pair<int, double> &p2){ return p1.second < p2.second;});
                int node = top.first;
                if (spread[node] * 1.0 / R*n <= Cost[node] / adv->bpi){ //if there is no user satisfies any companies (profit gain >= 0)
                    adv->isFull = 1;
                    flag --;
                    continue;
                }
                ASSERT(node >= 0);
                //cout <<"current seed: " <<node << endl;

                adv->seedSet.push_back(node);

                if (overInf[adv->companyID] == false){
                    if (estimateCurrInf(adv,R)){
                        adv->bpi = (adv->budget/adv->influence) *gammaR * 1.0; overInf[adv->companyID] = true;
                    }

                }

                sc_map.erase(node);

                tag[node] = true;
                for (auto Idx : hyperG[node])  // Idx: rr set covered by node
                {
                    //Remove from R all rr sets that are covered by node
                    if (edgeMark[Idx])continue;
                    adv->num_coveredRR.insert(pair<int, int>((int)Idx, 0));

                    int coveredRR_count = 0;

                    for (auto it : hyperGT[Idx]){ // it: node id
                        //// step1: nodes in rr set Idx spread--
                        --spread[it];
                        if(!tag[it]) sc_map[it] = spread[it] * 1.0 /  Cost[it];
                        //// step2: v is activated by node and adopt company based the maximum probability of v's in-neighbors
                        //cout << "isSelect["<<it<<"] =  " << isSelect[it] << "\tadopt company: " << Adoption[it].advMaxId << endl;
                        coveredRR_count += estimateInfnode(node, it, adv, edgeMark,1);


                    }//end for it

                    edgeMark[Idx] = true;

                    adv->num_coveredRR[Idx] = coveredRR_count;
                    //cout <<"current rr id: " << Idx << "\tcurrent rr covered num is: " << coveredRR_count << "\t current rr size is : " << hyperGT[Idx].size() << endl;
                }//end for Idx

            }//end for one batch seed select for nrCompanies (one by one)

        }// end for while (all companies seed selection ended)


    }

    //// 2. this function is allocation method: M = [T] * [h] metroid
    void allocate_Fill(unsigned int R, double gammaR, double gammaP)
    {
        allocQueueM.clear(); // (adv id, node id, max gain) (<int, int, double>)
        usersExpand.clear();
        criterQueue.clear();
        isSelect.clear();

        Adoption = vector<adoptInfo> (max_n);
        for (int i = 0; i < max_n; i ++){
            Adoption[i].advMaxId = -1; Adoption[i].Maxweight = 0.0;
        }

        for (int i = 0; i < nrCompanies; i ++){ //clear adv seedset
            company *adv;
            adv = advList.at(i);

            adv->seedSet.clear();
            adv->seedSet.shrink_to_fit();

            adv->num_coveredRR.clear();
            adv->totalProfit = 0.0; adv->totalRevenue = 0.0; adv->totalSeedCosts = 0.0;
            adv->currentInf = 0.0; adv->totalInf = 0.0; adv->isFull = 0;

            adv->bpi = (float) adv->budget/adv->influence;
            adv->bpi = adv->bpi * gammaP *1.0;
        }

        vector<bool> overInf (nrCompanies, false);
        isSelect = std::vector<int>(max_n,1); // only disjointness
        vector<int>spread(max_n, 0);
        for (unsigned int k = 0; k < max_n; k++)spread[k] = (int)hyperG[k].size();
        long long numEdge = hyperGT.size();
        vector<bool> edgeMark(numEdge, false);
        infPair best;  // (adv id, node id, marginal profit gain) (<int,int, double>)
        company *adv;
        int flag_isfull = nrCompanies;
//        int PorI = 0; // 0: profit max, 1: inf max

        for (int iterC = 0; iterC < nrCompanies; iterC ++){ // M' = [h] * [n]
            adv = advList.at(iterC);
            //for (unsigned int k = 0; k < max_n; k++){  ////M' = [V] * [h]
            for (int k: candidateSetT) { ////candidate set T M = [T] * [h]
                double mpg_temp = ((( double) n * (( double) spread[k] / R)) * adv->bpi) / (1.0 * Cost[k]);
                best.adv = adv->companyID; best.node = k; best.mpg = mpg_temp;
                if (mpg_temp > 1.0) // marginal profit gain based on company adv > 0
                    criterQueue.push_back(best);
                //cout << "initial criterQueue company id: " << adv->companyID << "\tnode: "<< k << "\t mpg: "<< mpg_temp << endl;
            }
        }
        cout << "allocQueueM initial begin!" <<"\t criterQueue.size: "<< criterQueue.size()<< endl;
        if (selectBestPSNode(best, spread, criterQueue, R)) // maximum profit marginal gain
            allocQueueM.push_back(best);//(adv id, node id, marginal profit gain)
        while (!allocQueueM.empty() && flag_isfull){ //stop allocation when no more company is available for allocation
            //multimap< double, int>::iterator iter = allocQueueM.end(); // the maximal marginal gain node
            //iter--;
            infPair iter = allocQueueM[allocQueueM.size()-1];
            adv = advList.at(iter.adv);
            //allocQueueM.erase(iter);
            allocQueueM.pop_back();
            if (adv->isFull){ // find the next best node
                if (!criterQueue.empty()){
                    if (selectBestPSNode(best, spread, criterQueue, R)) // maximum profit marginal gain
                        allocQueueM.push_back(best);//(adv id, node id, marginal profit gain)
                }
                continue;
            }
            if(isSelect[candidateNode] == 1) { // if the best candidate still has chance for assignment
                if(iter.mpg > 1.0) { // if the best profit gain > 0
                    //cout << "current seed: " << candidateNode << "\tadv: " << adv->companyID<<"\tprofit gain: " << iter.mpg << "\tspread: " << spread[candidateNode] << endl;
                    if (overInf[adv->companyID] == false){ //// if adv_currentinf> adv_influence, change gammaP to gammaR
                        if (estimateCurrInf(adv,R))
                        { adv->bpi = (adv->budget/adv->influence) *gammaR; overInf[adv->companyID] = true;}
                    }
                    assignBestNode(adv, spread, edgeMark,R); // inserts into advertiser's seed set and decreases the node's attention quota
                    if (!criterQueue.empty()){
                        if (selectBestPSNode(best, spread, criterQueue, R)) // maximum profit marginal gain
                            allocQueueM.push_back(best);//(adv id, node id, marginal profit gain)
                    }
                } // end of marginal profit>0 feasible part
                else{ // if not marginal profit>0 feasible, we stop the allocation for this advertiser so sktir et
                    cout << "adv " <<adv->companyID <<" [marginal profit <= 0] " << endl;
                    ////isSelect[candidateNode] = 0;
                    adv->isFull = 1;
                    flag_isfull--;
                    if (!criterQueue.empty()){
                        if (selectBestPSNode(best, spread, criterQueue, R)) // maximum profit marginal gain
                            allocQueueM.push_back(best);//(adv id, node id, marginal profit gain)
                    }
                }

            }// isSelect[candidateNode] > 0
            else {  //select another best node for this advertiser cause this seed is already allocated to another advertiser due to allocation priority
                if (!criterQueue.empty()){
                    if (selectBestPSNode(best, spread, criterQueue, R)) // maximum profit marginal gain
                        allocQueueM.push_back(best);//(adv id, node id, marginal profit gain)
                }
            }
        }
    }

    //// 3 this function is allocation method: M = [V] * [h] metroid
    void allocate_Iter(unsigned int R, int batchPS, int batchIS, double gammaR, double gammaP) {
        allocQueueM.clear(); // (max gain, adv id) (<double, int>)
        usersExpand.clear();
        criterQueue.clear();
        isSelect = std::vector<int>(max_n, 1); // only disjointness
        vector<int> spread(max_n, 0); for (unsigned int k = 0; k < max_n; k++){spread[k] = (int)hyperG[k].size();}
        long long numEdge = hyperGT.size();
        vector<bool> edgeMark(numEdge, false);
        infPair bestPS, bestIS, tempPS;  // (node id, marginal gain) (<int,double>)
        company *adv;
        vector<bool> overInf (nrCompanies, false);
        int flag_isfull = nrCompanies;
        for (int iterC = 0; iterC < nrCompanies; iterC++) { // initial PSNode selection
            adv = advList.at(iterC);
            //for (int candNode: candidateSetT) {
            for (unsigned int k = 0; k < max_n; k++){
                double mpg_temp = ((double) n * ((double) spread[k] / R)) * adv->bpi / (1.0 * Cost[k]);
                tempPS.adv = adv->companyID; tempPS.node = k; tempPS.mpg = mpg_temp;
                if (mpg_temp > 1.0)
                    criterQueue.push_back(tempPS);
                //cout << "initial criterQueue company id: " << adv->companyID << "\tcandNode: "<< k << "\t mpg: "<< mpg_temp << endl;
            }
        }
        cout << "criterQueue.size: " << criterQueue.size() <<endl;
        if (selectBestPSNode(bestPS, spread, criterQueue, R)) // maximum profit marginal gain
            allocQueueM.push_back(bestPS); //(adv id, node id, marginal profit gain)
        while (!allocQueueM.empty() && flag_isfull) {
            //// in profit-considering batches
            for (int i = 0; i < batchPS; i++) {
                if (allocQueueM.empty()) break;
                infPair iter = allocQueueM[allocQueueM.size()-1];
                adv = advList.at(iter.adv);
                allocQueueM.pop_back();
                if (adv->isFull){
                    if (!criterQueue.empty()){
                        if (selectBestPSNode(bestPS, spread, criterQueue, R)) // maximum profit marginal gain
                            allocQueueM.push_back(bestPS);//(adv id, node id, marginal profit gain)
                    }
                    continue;
                }
                if (isSelect[candidateNode] > 0) { // if the best candidate still has chance for assignment
                    if (iter.mpg > 1.0) { // if the best profit gain > 0
                        if (overInf[adv->companyID] == false){ //// if adv_currentinf> adv_influence, change gammaP to gammaR
                            if (estimateCurrInf(adv,R))
                            { adv->bpi = (adv->budget/adv->influence*1.0) *gammaR; overInf[adv->companyID] = true;}
                        }
                        //cout << "[profit batch] current seed: " << candidateNode << " current company: " << adv->companyID << "  profit gain: " << iter.mpg<< endl;
                        assignBestNode(adv, spread, edgeMark,R); // inserts into advertiser's seed set and decreases the node's attention quota
                        if (!criterQueue.empty()){
                            if (selectBestPSNode(bestPS, spread, criterQueue, R)) // maximum profit marginal gain
                                allocQueueM.push_back(bestPS);//(adv id, node id, marginal profit gain)
                        }
                    } // end of marginal profit>0 feasible part
                    else { // if not marginal profit>0 feasible, we stop the allocation for this advertiser so sktir et
                        cout << "adv " <<adv->companyID <<" [marginal profit <= 0] " << endl;
                        ////isSelect[candidateNode] = 0;
                        adv->isFull = 1;
                        flag_isfull--;
                        if (!criterQueue.empty()){
                            if (selectBestPSNode(bestPS, spread, criterQueue, R)) // maximum profit marginal gain
                                allocQueueM.push_back(bestPS);//(adv id, node id, marginal profit gain)
                        }
                    }

                }// isSelect[candidateNode] > 0
                else {  //select another best node for this advertiser cause this seed is already allocated to another advertiser due to allocation priority
                    if (!criterQueue.empty()) {
                        if (selectBestPSNode(bestPS,spread, criterQueue, R)) // maximum profit marginal gain
                            allocQueueM.push_back(bestPS);
                    }
                }
            }
            //// in influence-considering batch
            for (int j = 0; j < batchIS; j++) {

                if (allocQueueM.empty()) break;
                infPair iter = allocQueueM[allocQueueM.size()-1];
                adv = advList.at(iter.adv);
                allocQueueM.pop_back();
                if (adv->isFull){
                    if (!criterQueue.empty()){
                        if (selectBestISNode(bestIS,  spread, R)) // maximum profit marginal gain
                            allocQueueM.push_back(bestIS);//(adv id, node id, marginal profit gain)
                        else{
                            ////20221129
                            if (selectBestPSNode(bestPS,spread, criterQueue, R)) // maximum profit marginal gain
                                allocQueueM.push_back(bestPS);
                            batchIS = 0;
                            //break;
                        }
                    }
                    continue;
                }
                if (isSelect[candidateNode] == 1) { // if the best influence candidate still has chance for assignment
                    if (iter.mpg > 1.0) { // if the best profit gain > 0
                        if (overInf[adv->companyID] == false){ //// if adv_currentinf> adv_influence, change gammaP to gammaR
                            if (estimateCurrInf(adv,R))
                            { adv->bpi = (adv->budget / adv->influence *1.0) *gammaR; overInf[adv->companyID] = true;}
                        }
                        //cout << "[influence batch] current seed: " << candidateNode << " current company: " << adv->companyID << "  profit gain: " << iter.mpg << endl;
                        assignBestNode(adv, spread, edgeMark,R); // inserts into advertiser's seed set and decreases the node's attention quota
                        if (!criterQueue.empty()) {
                            if (selectBestISNode(bestIS, spread, R))//// exists bestIS company
                                allocQueueM.push_back(bestIS);
                            else{
                                ////20221129 IS->PS
                                if (selectBestPSNode(bestPS,spread, criterQueue, R)) // maximum profit marginal gain
                                    allocQueueM.push_back(bestPS);
                                batchIS = 0;
                                //break; ////20221129
                            }
                        }

                    } // end of marginal profit>0 feasible part
                    else {
                        isSelect[candidateNode] = 2; // may still be another company's best profit batch node
                        if (!criterQueue.empty()) {
                            if (selectBestISNode(bestIS, spread, R))//// exists bestIS company
                                allocQueueM.push_back(bestIS);
                            else{
                                if (!criterQueue.empty()) {
                                    if (selectBestPSNode(bestPS,spread, criterQueue, R)) // maximum profit marginal gain
                                        allocQueueM.push_back(bestPS);
                                }
                                batchIS = 0;
                                //break;
                            }
                        }

                    }
                }// isSelect[candidateNode] = 1
                else if (isSelect[candidateNode] == 0){//select another best node for this advertiser cause this seed is already allocated to another advertiser due to allocation priority
                    if (!criterQueue.empty()) {
                        if (selectBestISNode(bestIS, spread, R))//// exists bestIS company
                            allocQueueM.push_back(bestIS);
                        else{
                            if (selectBestPSNode(bestPS,spread, criterQueue, R)) // maximum profit marginal gain
                                allocQueueM.push_back(bestPS);
                            batchIS = 0;
                            //break;
                        }
                    }
                }
                else{ //isSelect[candidateNode] = 2
                    if (!criterQueue.empty()) {
                        if (selectBestISNode(bestIS, spread, R))//// exists bestIS company
                            allocQueueM.push_back(bestIS);
                        else{
                            if (selectBestPSNode(bestPS,spread, criterQueue, R)) // maximum profit marginal gain
                                allocQueueM.push_back(bestPS);
                            batchIS = 0;
                            //break;
                        }
                    }
                }
            }
        }
    }




////======================================================= main function end =====================================================


};
