#include "expandedmarkovchain.h"

EMkChain::EMkChain()
{

}

EMkChain::EMkChain(const Graph &originalGraph):m_graph(originalGraph)
{
    R1 = m_graph.edgeNum(1);
    R2 = m_graph.edgeNum(2);

#ifdef _DEBUG
    qDebug() << "R(2) = " << R2;
#endif
//    srand(time(0));
}

EMkStat EMkChain::init(int d, int sid)
{
    m_estate.clear();

    QVector<bool> visited(m_graph.nodeNum(), false);
    visited[sid] = true;

    QQueue<int> q;
    q.push_back(sid);

    MkStat stat0;
    while(!q.empty())
    {
        int tmp = q.front();
        q.pop_front();

        stat0.push_back(tmp);
        if (stat0.size() == d)
        {
            break;
        }

        Node* tmpNode = m_graph.GetNode(tmp);
        for (int i = 0; i < tmpNode->degree(); i++)
        {
            if (!visited[tmpNode->childs[i]])
            {
                visited[tmpNode->childs[i]] = true;
                q.push_back(tmpNode->childs[i]);
            }
        }
    }
    // initialize current state
    if (stat0.size() == d) {
        m_estate.push_back(stat0);
        return m_estate;
    }
    else
    {
        qDebug() << "Not enough connected nodes to form a state!";
        assert(0);
    }
}


EMkStat EMkChain::NBInit(int k, int d, int sid)
{
    // first find d-node state
   init(d, sid);
   // then update to next d-node state, and it ensures not the same one as before
   update(k, d);

   // now we have e ij as our initial state
   return m_estate;
}


EMkStat EMkChain::NBUpdate(int k, int d)
{
    // current two states together combines the NB state
    NBMkStat eij = QPair<MkStat, MkStat>(m_estate[m_estate.size()-2], m_estate[m_estate.size()-1]);
    NBMkStat ejk = eij;


    QVector<MkStat> Nj = popNeighbors(eij.second, d);
    int dj = Nj.size();

    if (dj > 0)
    {
        if (dj == 1)
        {
    //        qDebug() << "d(j) = 1!";
            ejk = QPair<MkStat, MkStat>(eij.second, Nj[0]);
        }
        else
        {
            while (true)
            {
                int r = rand()%dj;
                if (Nj[r] != eij.first) // Might bug here because of order in the MkStat
                {
                    ejk = QPair<MkStat, MkStat>(eij.second, Nj[r]);
                    break;
                }
            }
        }

    }

    if (m_estate.size() == k-d+1)
    {
        m_estate.pop_front();
    }

    m_estate.push_back(ejk.second);
    return m_estate;
}

//bool ExpandedMarkovChain::statEqual(const MkStat &xi, const MkStat &xj)
//{
//    return (xi[0] == xj[0] && xi[1] == xj[1] ) || (xi[0] == xj[1] && xi[1] == xj[0] );
//}

QVector<MkStat> EMkChain::popNeighbors(MkStat xi, int d)
{
    QVector<MkStat> neighbors;

    switch (d) {
    case 1:
    {
        Node* nodeu = m_graph.GetNode(xi[0]);

        // Uniformly choose one of u's neighbor
        for (int i = 0; i < nodeu->degree(); i++)
        {
            neighbors.append(MkStat({nodeu->childs[i]}));
        }
        break;
    }
    case 2:
    {
        int u, v;
        u = xi[0];
        v = xi[1];

        Node *nodeU, *nodeV;
        nodeU = m_graph.GetNode(u);
        nodeV = m_graph.GetNode(v);

        int du, dv;
        du = nodeU->degree();
        dv = nodeV->degree();

        // push all u's neighbor
        for (int i = 0; i < du; i++)
        {
            int w = nodeU->childs[i];
            if (w != v)
            {
                neighbors.append(MkStat({u, w}));
            }
        }

        // push all v's neighbor
        for (int i = 0; i < dv; i++)
        {
            int w = nodeV->childs[i];
            if (w != u)
            {
                neighbors.append(MkStat({w, v}));
            }
        }

        break;
    }
    default:
        assert(false);
        break;
    }


    return neighbors;
}

int EMkChain::R(int d)
{
    if (d == 1)
    {
        return R1;
    }
    else if (d == 2)
    {
        return R2;
    }
}

int EMkChain::isoIndicator(const EMkStat& xl, int k, int i)
{
    QVector<int> V = xl.nodeSet();

    if (V.size() < k)
    {
        return 0;
    }
    else
    {
        Graph s = m_graph.inducedSubgraph(V);
        int type = graphType(s);

        if (type != i)
        {
            return 0;
        }
        else {
            return 1;
        }
    }
}

int EMkChain::sizeIndicator(const EMkStat &xl, int k)
{
    QVector<int> V = xl.nodeSet();
//    qDebug() << int(V.size() == k);
    return V.size() == k;//// might bug here
}

float EMkChain::stableDistr(const EMkStat &xl, int d)
{
    float deno = 2*R(d);

    if (xl.size() == 1)
    {
        return dgree(xl[0])/deno;
    }
    else if(xl.size() == 2)
    {
        return 1.f/deno;
    }
    else
    {
        float product = 1.f/deno;
        for (int i = 1; i < xl.size()-1; i++)
        {
            product *= (1.f / dgree(xl[i]));
        }
        return product;
    }
}

float EMkChain::sampleProb(const EMkStat &xl, int k, int d)
{
    int l = k-d+1;

    // subgraph induced by nodes in Xl
    QVector<int> V = xl.nodeSet();
    Graph s = m_graph.inducedSubgraph(V);

    QTime time;
    time.start();
    // find all d-node connected induced subgraph
    QVector<MkStat> S = s.connectedInducedGraph(d);////Slow step

//    qDebug() << "d-node connected induced subgraph " << time.elapsed() << "ms";

//    QVector<Graph> S;
//    for (int i = 0; i < dgraphs.size(); i++)
//    {
//        QVector<int> dgraph = dgraphs[i];
//        S.append(m_graph.inducedSubgraph(dgraph));
//    }


    float PXl = 0.f;

    time.start();
    // combination of l elements from S
    QVector<QVector<MkStat>> Sls = combine(S, l);     ////Slow step
//    qDebug() << "combine " << time.elapsed() << "ms";

    time.start();
    for (int i =0; i < Sls.size(); i++)
    {
        QVector<MkStat> Sl = Sls[i];
        // sort nodes in each Markov state
        for (int j = 0; j < Sl.size(); j++)
        {
            std::sort(Sl[j].begin(), Sl[j].end());
        }

        QSet<int> U;
        for (int j = 0; j < Sl.size(); j++)
        {
            for (int k = 0; k < Sl[j].size(); k++)
            {
                U.insert(Sl[j][k]);
            }
        }


        // if size of node set ∪ s∈S (l) V (s) equals to k
        if (U.size() == k)
        {
            // all permutations of Sl
            while(std::next_permutation(Sl.begin(), Sl.end()))
            {
                // corresponding state
                bool testflag = true;
                for (int j = 0; j <= l-2; j++)
                {
                    QVector<int> inter;
                    // if Sl[j] and Sl[j+1] shares d-1 nodes for i in [0, l-1)
                    std::set_intersection(Sl[j].begin(), Sl[j].end(), Sl[j+1].begin(), Sl[j+1].end(), std::back_inserter(inter));
                    if (inter.size() != d-1)
                    {
                        testflag = false;
                        break;
                    }
                }
                // corresponding state X 0 ← (S x 1 , · · · , S x l )
                if (testflag)
                {
                    EMkStat c_xl;
                    for (int j = 0; j < Sl.size(); j++)
                    {
                        c_xl.push_back(Sl[j]);
                    }

                    PXl += stableDistr(c_xl, d);
                }
            }
        }
    }
//    qDebug() << "loop over all combinations " << time.elapsed() << "ms";

    return PXl;
}

int EMkChain::coff(int k, int i, int d)
{
    // Hard code for d=1, 2, for 29 types of graphlets
    QVector<QVector<int>> hardCode = {
        {
            2, 6, 2, 0, 8, 4, 12, 24, 2, 0, 0, 2, 4, 0, 10, 4, 4, 8, 8, 12, 14, 12, 12, 20, 28, 36, 48, 72, 120
        },
        {
            2, 6, 2, 6, 8, 10, 24, 48, 4, 24, 10, 8, 32, 10, 12, 48, 48, 24, 36, 30, 108, 72, 12, 84, 68, 164, 152, 288, 480
        },
    };

    if (k == 3)
    {
        return hardCode[d-1][i];
    }
    else if (k == 4)
    {
        return hardCode[d-1][2+i];
    }
    else if (k == 5)
    {
        return hardCode[d-1][8+i];
    }
    qDebug() << "error";
    assert(0);
}

QVector<float> EMkChain::SRW(int n, int d, int k, bool norm)
{
//    qDebug() << QString("SRW(n: %1, d: %2, k: %3)").arg(n).arg(d).arg(k);

    EMkStat xl;
    // initialize state
    xl = init(d);

    // walk l-1 steps to get the initial state Xl
    int l = k - d + 1;
    for (int i = 0; i < l-1; i++)
    {
        xl = update(k, d);
    }

    int gk = gletTypeNum(k);
    QVector<float> _c(gk, 0.f);

    int t = 0;
    while(t < n)
    {
        // subgraph induced by nodes in Xl
        QVector<int> V = xl.nodeSet();
        if (V.size() == k)//collect k distinct nodes
        {
            Graph s = m_graph.inducedSubgraph(V);
            int id = graphType(s);
            if (id != -1)
            {
                _c[id] += 1.f/(coff(k, id, d)*stableDistr(xl, d));
            }
        }
        // uniformly choose neighbor and update
        xl = update(k, d);
//        xl.printStatus();

        t++;
    }

    if (norm)
    {
        math_utils::normalize(_c);
    }

    return _c;
}

QVector<float> EMkChain::SRW(int n)
{
    bool normInner = false;

    QVector<float> gfd3 = SRW(n, 1, 3, normInner);
    QVector<float> gfd4 = SRW(n, 2, 4, normInner);
    QVector<float> gfd5 = SRW(n, 2, 5, normInner);

    QVector<float> gfd;
    gfd.append(gfd3);
    gfd.append(gfd4);
    gfd.append(gfd5);

    if (!normInner)
    {
//        qDebug() << "normalize together!";
//        utils::printVector(gfd);
//        utils::normalize(gfd);
    }
    assert(gfd.size() == ALL_GRAPHLET);

    return gfd;
}

float EMkChain::NB_stableDistr(const EMkStat &xl, int d)
{
    float deno = 2*R(d);

    // we discard |Rd| here because it will be emiminated during the division
    if (xl.size() == 1)
    {
        return NB_dgree(xl[0])/deno;
    }
    else if(xl.size() == 2)
    {
        return 1.f/deno;
    }
    else
    {
        float product = 1.f/deno;
        for (int i = 1; i < xl.size()-1; i++)
        {
            product *= (1.f / NB_dgree(xl[i]));
        }
        return product;
    }
}

QVector<float> EMkChain::NB_SRW(int n, int d, int k, bool norm)
{
#ifdef _DEBUG
    qDebug() << QString("NB_SRW(n: %1, d: %2, k: %3)").arg(n).arg(d).arg(k);
#endif

    EMkStat xl;
    // initialize state
    xl = NBInit(k, d);

    // walk l-1 steps to get the initial state Xl
    int l = k - d + 1;
    for (int i = 0; i < l-1; i++)
    {
        xl = NBUpdate(k, d);
    }

//#ifdef _DEBUG
//        xl.printStatus();
//#endif

    int gk = gletTypeNum(k);
    QVector<float> _c(gk, 0.f);

    int t = 0;
    while(t < n)
    {
        // subgraph induced by nodes in Xl
        QVector<int> V = xl.nodeSet();
        if (V.size() == k)  //collect k distinct nodes
        {
            Graph s = m_graph.inducedSubgraph(V);
            int id = graphType(s);
            if (id != -1)
            {
                _c[id] += 1.f/(coff(k, id, d)*NB_stableDistr(xl, d));
//                _c[id] += 1.f/NB_stableDistr(xl, d);
            }
        }

        // uniformly choose neighbor and update
        xl = NBUpdate(k, d);

#ifdef _DEBUG
//        xl.printStatus();
#endif

        t++;
    }

//    for (int i = 0; i < _c.size(); i++)
//    {
//        _c[i] *= coff(k, i, d);
//    }

    if (norm)
    {
        math_utils::normalize(_c);
    }

    return _c;
}

QVector<float> EMkChain::NB_SRW(int n)
{
    bool normInner = false;
    QVector<float> gfd3 = NB_SRW(n, 1, 3, normInner);
    QVector<float> gfd4 = NB_SRW(n, 2, 4, normInner);
    QVector<float> gfd5 = NB_SRW(n, 2, 5, normInner);

    QVector<float> gfd;
    gfd.append(gfd3);
    gfd.append(gfd4);
    gfd.append(gfd5);

    if (!normInner)
    {
//        qDebug() << "normalize together!";
//        utils::printVector(gfd);
        math_utils::normalize(gfd);
    }
    assert(gfd.size() == ALL_GRAPHLET);

    return gfd;
}

QSet<GraphLet> EMkChain::NB_SRW_GLETS(int n, int d, int k)
{
#ifdef _DEBUG
    qDebug() << QString("NB_SRW(n: %1, d: %2, k: %3)").arg(n).arg(d).arg(k);
#endif

    QVector<QVector<int>> ccs = m_graph.CCs();
//    QVector<int> ccs = m_graph.CCs();

    QSet<GraphLet> glets;
    for (int cc = 0; cc < ccs.size(); cc++)  // for each connected component
    {
        EMkStat xl;
        // initialize state, select one node in current cc as starting node
//        xl = NBInit(k, d, ccs[cc]);
        xl = NBInit(k, d, ccs[cc][0]);

        // walk l-1 steps to get the initial state Xl
        int l = k - d + 1;
        for (int i = 0; i < l-1; i++)
        {
            xl = NBUpdate(k, d);
        }

    #ifdef _DEBUG
            xl.printStatus();
    #endif


        int t = 0;
        while(t < n)
        {
            // subgraph induced by nodes in Xl
            QVector<int> V = xl.nodeSet();
            // if collect k distinct nodes from state xl
            if (V.size() == k)
            {
                GraphLet s = m_graph.inducedGraphlet(V);
                glets.insert(s);
            }

            // uniformly choose neighbor and update
            xl = NBUpdate(k, d);

    #ifdef _DEBUG
    //        xl.printStatus();
    #endif

            t++;
        }
    }

    return glets;
}

QSet<GraphLet> EMkChain::NB_SRW_GLETS(int n)
{
    //相同节点集，induced graph肯定一样，则对应的graphlet一定是同一个
    QSet<GraphLet> glets3 = NB_SRW_GLETS(n, 1, 3);
    QSet<GraphLet> glets4 = NB_SRW_GLETS(n, 2, 4);
    QSet<GraphLet> glets5 = NB_SRW_GLETS(n, 2, 5);

    QSet<GraphLet> glets = glets3 + glets4 + glets5;
    return glets;
}

QVector<float> EMkChain::SRW_CSS(int n, int d, int k)
{
    qDebug() << QString("SRW_CSS(n: %1, d: %2, k: %3)").arg(n).arg(d).arg(k);

    EMkStat xl;

    // initialize state
    xl = init(d);
    qDebug() << xl.size();

    // walk l-1 steps to get the initial state Xl
    int l = k - d + 1;
    for (int i = 0; i < l-1; i++)
    {
        xl = update(k, d);
    }
    qDebug() << xl.size();


    int gk = gletTypeNum(k);
    QVector<float> numers(gk, 0.f);
    QVector<float> denos(gk, 0.f);

    int t = 0;
    while(t < n)
    {
        for (int i = 0; i < gk; i++)
        {
            float prob = sampleProb(xl, k, d);      // why zero?
            if (prob == 0.f)
            {
                qDebug() << "impossble!" ;
            }
            qDebug() << "sample probability: " << prob;
            // Accumulate
            numers[i] += isoIndicator(xl, k, i) / prob;
            denos[i] += sizeIndicator(xl, k) / prob;
        }

        // uniformly choose neighbor and update
        xl = update(k, d);
//        xl.printStatus();

        t++;
    }

    // the concentration of graphlets
    QVector<float> _c;
    for (int i = 0; i < numers.size(); i++)
    {
        _c.append(numers[i]/denos[i]);
    }
    return _c;
}

QVector<float> EMkChain::SRW_CSS(int n)
{
    QVector<float> gfd3 = SRW_CSS(n, 1, 3);
    QVector<float> gfd4 = SRW_CSS(n, 2, 4);
    QVector<float> gfd5 = SRW_CSS(n, 2, 5);

    QVector<float> gfd;
    gfd.append(gfd3);
    gfd.append(gfd4);
    gfd.append(gfd5);
    assert(gfd.size() == ALL_GRAPHLET);

    return gfd;
}

EMkStat EMkChain::update(int k, int d)
{
    MkStat xl =  getCurrentState();
    MkStat xl_1 = xl;

    switch (d) {
    case 1:
    {
        Node* nodeu = m_graph.GetNode(xl[0]);
        int du = nodeu->degree();

        if (du > 0)
        {
            // Uniformly choose one of u's neighbor
            int w = nodeu->childs[rand()%du];
            xl_1 = MkStat({w});
        }

        break;
    }
    case 2:
    {
        int u, v;
        u = xl[0];
        v = xl[1];

        Node *nodeU, *nodeV;
        nodeU = m_graph.GetNode(u);
        nodeV = m_graph.GetNode(v);

        int du, dv;
        du = nodeU->degree();
        dv = nodeV->degree();

        if (du != 1 || dv != 1)
        {
            float pu = (float)du/(du + dv); // winning probability of u
            while (true)
            {
                int winner, loser;

                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    //            qDebug() << "p(u) = " << pu << " rand r = " << r;
                if (r <= pu)
                {
                    winner = u;
                    loser = v;
                }
                else
                {
                    winner = v;
                    loser = u;
                }

                Node* nodeWinner = m_graph.GetNode(winner);
                int w = nodeWinner->childs[rand() % nodeWinner->degree()];
                if (w != loser)
                {
                    // to keep the order of nodes in Markov state
                    if (winner == u)
                    {
                        xl_1 = MkStat({winner, w});
                    }
                    else
                    {
                        xl_1 = MkStat({w, winner});
                    }

                    break;
                }
            }
        }

        break;
    }
    default:
        assert(false);
        break;
    }

    if (m_estate.size() == k-d+1)
    {
        m_estate.pop_front();
    }

    m_estate.push_back(xl_1);

    return m_estate;
}

MkStat EMkChain::getCurrentState()
{
    return m_estate[m_estate.size()-1];
}

EMkStat EMkChain::getCurrentExpandState()
{
    return m_estate;
}

int EMkChain::dgree(const MkStat &xi)
{
    int d = xi.size();
    if (d == 1)
    {
        return m_graph.GetNode(xi[0])->degree();
    }
    else if (d == 2)
    {
        return m_graph.GetNode(xi[0])->degree()-1 + m_graph.GetNode(xi[1])->degree() -1;
    }
    else
    {
        assert(false);
    }
}

int EMkChain::NB_dgree(const MkStat &xi)
{
    return std::max(dgree(xi)-1, 1);
}
