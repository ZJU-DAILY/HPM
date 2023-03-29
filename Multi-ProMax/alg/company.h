

#ifndef ALG_COMPANY_H
#define ALG_COMPANY_H
using namespace std;

class company{

public:

    int companyID;
    std::vector<int> seedSet;
    map<int, int> num_coveredRR; //<rrId, covered user in rrId>
    map<int, int> num_coveredRR2;
    int isFull; //initial: 0

    double totalProfit;
    double totalRevenue;
    double totalSeedCosts;
    double currentInf, totalInf;

    double budget, influence;
    double bpi;

    company(int id)
    {
        this->companyID = id;
        this->bpi = 1.0 * this->budget / this->influence;

    }

    ~company() {}

};






#endif //ALG_COMPANY_H


