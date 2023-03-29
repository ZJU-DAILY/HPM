
//// * BFS starting from one node

bool BuildHypergraphNodeSingle(vector<bool>&tag) //while(i++<R)EstPro2+=BuildHypergraphNodeSingle(tag);
{
    auto uStart = sfmt_genrand_uint32(&sfmtSeed) % n;
    if(tag[uStart]) return true;

    unsigned int n_visit_mark = 0, curIdx = 0;
    visit_mark[n_visit_mark++] = uStart;
    visit[uStart] = true;
    //hyperG[uStart].push_back(hyperId);
    bool flag=false;

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
                if(tag[v])
                {
                    flag=true;
                    break;
                }
                visit[v] = true;
                visit_mark[n_visit_mark++] = v;
                //hyperG[v].push_back(hyperId);
                //++spread[v];
            }
            if(flag)break;
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
            if(tag[v])
            {
                flag=true;
                break;
            }
            visit[v] = true;
            visit_mark[n_visit_mark++] = v;
            //hyperG[v].push_back(hyperId);
            //++spread[v];
        }
    }
    else
        ASSERT(false);

    for (unsigned int i = 0; i < n_visit_mark; i++)visit[visit_mark[i]] = false;	//// reset visit vector

    //hyperGT.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark));

    return flag;
}

void BuildHypergraphNode(unsigned int hyperId,  unsigned int max_n, set<int> &nodes)
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
	hyperG[uStart].push_back(hyperId);
	//++spread[uStart];

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
				hyperG[v].push_back(hyperId);
				//++spread[v];
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

			visit[v] = true;
			visit_mark[n_visit_mark++] = v;						
			hyperG[v].push_back(hyperId);
			//++spread[v];
		}
	}

	else
		ASSERT(false);

 	for (unsigned int i = 0; i < n_visit_mark; i++)visit[visit_mark[i]] = false; //方便构造下个rr set

	hyperGT.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark));

	return;
}

/*
//Possible world simulation, useless
unsigned int simulation()
{
	unsigned int n_visit_mark = 0, curIdx = 0;
	vector<bool>active(n, false);
	for (auto it : seedSet)
	{
		visit_mark[n_visit_mark++] = it;
		active[it] = true;
	}

	while (curIdx < n_visit_mark)
	{
		int expand = visit_mark[curIdx++];
		for (auto v : PO[expand])
		{
			if (active[v])continue;
			visit_mark[n_visit_mark++] = v;
			active[v] = true;
		}
	}
	return n_visit_mark;
}*/
//new
/*
unsigned int simulation()
{
    unsigned int n_visit_mark = 0, curIdx = 0;
    vector<bool>active(n, false);
    for (auto it : seedSet)
    {
        visit_mark[n_visit_mark++] = it;
        active[it] = true;
    }

    while (curIdx < n_visit_mark)
    {
        int expand = visit_mark[curIdx++];
        for (auto v : gT[expand])
        {
            if (active[v])continue;
            visit_mark[n_visit_mark++] = v;
            active[v] = true;
        }
    }
    return n_visit_mark;
}
*/

/*
void realization(int curnode)
{	
	unsigned int n_visit_mark = 0, curIdx = 0;	
	
	active_set[curnode] = true;
	visit_mark[n_visit_mark++] = curnode;
	--NumcurNode;
	curNodeIdx[curNode[NumcurNode]] = curNodeIdx[curnode];
	curNode[curNodeIdx[curnode]] = curNode[NumcurNode];
	

	while (curIdx < n_visit_mark)
	{
		int expand = visit_mark[curIdx++];		
		for (auto v : PO[expand])
		{
			if (active_set[v])continue;
			visit_mark[n_visit_mark++] = v;
			active_set[v] = true;
			--NumcurNode;
			curNodeIdx[curNode[NumcurNode]] = curNodeIdx[v];
			curNode[curNodeIdx[v]] = curNode[NumcurNode];
		}
	}
	fill(spread.begin(), spread.end(), 0);
}
*/
