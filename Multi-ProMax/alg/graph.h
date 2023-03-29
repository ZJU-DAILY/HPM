#define HEAD_INFO
#include "sfmt/SFMT.c"
#include "sfmt/SFMT.h"
#include "head.h"
#include <iostream>

#include "MeasureM.h"
#include "memoryusage.h"

using namespace std;
typedef double (*pf)(int, int);
void handle_error(const char* msg);

class Graph
{
public:    	
	unsigned int n, m, nrCompanies;
    ////
    unsigned int max_n;

    vector<vector<int>> gT;
	vector<vector<double>> probT;
    set<int> nodes;


    enum InfluModel {IC, LT};
    InfluModel influModel;
    void setInfuModel(InfluModel p)
    {
        influModel = p;
        TRACE(influModel == IC);        
        TRACE(influModel == LT);
    }

    string folder;
    string graph_file;
    void readNM()
    {
        cout << "Read attribute ..." << endl;
        ifstream cin((folder + "attribute.txt").c_str());
        ASSERT(!cin == false);
        string s;
        while (cin >> s)
        {
            if (s.substr(0, 2) == "n=")
            {
                n = atoi(s.substr(2).c_str());
                continue;
            }
            if (s.substr(0, 2) == "m=")
            {
                m = atoi(s.substr(2).c_str());
                continue;
            }
//            if (s.substr(0, 2) == "r=")
//            {
//                nrCompanies = atoi(s.substr(2).c_str());
//                continue;
//            }
            ASSERT(false);

        }
        TRACE(n, m);
        cin.close();
    }

    void readGraph()
    {
        max_n = 0;
        cout << "Read graph ..." << endl;
        std::ifstream infile(graph_file);
        if (!infile.is_open())
        {
            std::cout << "The file \"" + graph_file + "\" can NOT be opened\n";
            return;
        }
        for (unsigned int i = 0; i < m; i++){
            unsigned int a, b;
            double p;
            infile >> a >> b >> p;
            cout << "a: " <<a<<"\tb: "<<b <<"\tp: " << p<< endl;
            //ASSERT( a < n );
            //ASSERT( b < n );
            if (a == b)
                continue;

            probT[b].push_back(p);
            gT[b].push_back(a);

            if (a>=max_n) max_n = a;
            if (b>=max_n) max_n = b;
            nodes.insert(a); nodes.insert(b);
        }
        max_n ++;
        cout << "max_n: " << max_n <<endl;
        infile.close();

    }

    void readGraph_()
    {
        cout << "read Graph..." <<endl;
        max_n = 0;
        size_t length;
        int fd = open((graph_file).c_str(), O_RDWR);
        if (fd == -1)
            handle_error("open");
        struct stat sb;
        int rc = fstat(fd, &sb);
        if (rc == -1)
            handle_error("fstat");

        length = sb.st_size;
        auto ptr = static_cast<char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));  //byte by byte
        auto f = ptr;

        int gap = 2 * sizeof(int) + sizeof(double);
        //ASSERT(fin != false);
        unsigned int readCnt = 0;
        for (unsigned int i = 0; i < m; i++)
        {
            readCnt ++;
            unsigned int a, b;
            double p;
            memcpy(&a, f, sizeof(int));
            memcpy(&b, f + sizeof(int), sizeof(int));
            memcpy(&p, f + 2 * sizeof(int), sizeof(double));
            f += gap;
            cout << "a: " <<a<<"\tb: "<<b <<"\tp: " << p<< endl;
            //ASSERT( a < n );
            //ASSERT( b < n );

            //probT[b].push_back((unsigned int)(p*UI_MAX));
            probT[b].push_back(p);
            gT[b].push_back(a);

            if (a>=max_n) max_n = a;
            if (b>=max_n) max_n = b;
            nodes.insert(a); nodes.insert(b);
        }

        ASSERT(readCnt == m);
        rc = munmap(ptr, length);
        close(fd);
        //fclose(fin);
        cout << "read graph end"<< endl;
    }


    Graph(string folder, string graph_file): folder(folder), graph_file(graph_file)
    {
		//UI_MAX = 4294967295U;		

		readNM();

		gT = vector<vector<int>>(n*2, vector<int>()); ////n*2
		probT = vector<vector<double>>(n*2, vector<double>());////n*2
		//probT = vector<vector<unsigned int>>(n, vector<unsigned int>());
        nodes.clear();

		readGraph();
    }
};

double sqr(double t)
{
    return t * t;
}

void handle_error(const char* msg) 
{
	perror(msg);
	exit(255);
}

