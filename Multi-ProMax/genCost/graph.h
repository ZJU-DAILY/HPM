#define HEAD_INFO
#include "sfmt/SFMT.c"
#include "sfmt/SFMT.h"
#include "head.h"

using namespace std;
typedef double (*pf)(int, int);
void handle_error(const char* msg);

class Graph
{
public:    	
	unsigned int n, m;
    unsigned int max_n=0;

    vector<vector<int>> gT;
	vector<vector<double>> probT; 
	vector<int>outdeg;
	//vector<vector<unsigned int>> probT;
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
            ASSERT(false);
        }
        cout << "n: " <<n<<"\tm: "<<m<<endl;
        TRACE(n, m );
        cin.close();
    }

/*    void readGraph() ////probT no use
    {
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

            ASSERT( a < n );
            ASSERT( b < n );
			
			//probT[b].push_back((unsigned int)(p*UI_MAX));
			probT[b].push_back(p);
			gT[b].push_back(a);		
			++outdeg[a];
        }        

        ASSERT(readCnt == m);
		rc = munmap(ptr, length);
		close(fd);
        //fclose(fin);
    }
*/
    void readGraph()
    {
        cout << "Read graph ..." << endl;
        std::ifstream infile(graph_file);
        if (!infile.is_open())
        {
            std::cout << "The file \"" + graph_file + "\" can NOT be opened\n";
            return;
        }
        for (unsigned int i = 0; i < m; i++){
            unsigned int a, b;
            infile >> a >> b;
            ASSERT( a < n );
            ASSERT( b < n );
            if (a == b)
                continue;
            gT[b].push_back(a);
            ++outdeg[a];
            if(a>=max_n) max_n=a;
            if(b>=max_n) max_n=b;
            nodes.insert(a); nodes.insert(b);
            cout << "a: " << a << "\tb: " <<b<< endl;
        }

        infile.close();

        //// probability on edge : Weighted cascade （1/d_v^in）
        set<int>::iterator Iter;
        cout << "max_n: " << max_n<< endl;
        cout << "nodes.size: "<<nodes.size()<<endl;
        //probT = vector<vector<double>>(max_n, vector<double>(max_n,0.0));
        probT = vector<vector<double>>(max_n+1, vector<double>());
        for(Iter = nodes.begin(); Iter != nodes.end(); Iter ++)
        {
            int Idx = *Iter;
            cout << "node " << Idx;
            double weight;
            if (gT[Idx].size() == 0) {weight = 0.; continue;}
            else weight = 1.0 / gT[Idx].size();
            cout << " probability weight: " << weight << endl;
            for (int i = 0; i < gT[Idx].size(); i ++)
            {
                int Idy = gT[Idx][i];
                probT[Idx].push_back(weight);
                //probT[Idx][Idy] = weight;  //// edge u1--u2(weight) according to code graphT[u2].push_back(u1)
            }
        }

    }

	Graph(string folder, string graph_file) : folder(folder), graph_file(graph_file)
    {
		//UI_MAX = 4294967295U;		

		readNM();

		gT = vector<vector<int>>(n*2, vector<int>());
        //probT = vector<vector<double>>(n, vector<double>(n,0.0));
		outdeg = vector<int>(n*2, 0);
		//probT = vector<vector<unsigned int>>(n, vector<unsigned int>());

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

