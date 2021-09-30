#include "graph.h"
#include "expandedmarkovchain.h"
#include <QDebug>

Graph::Graph()
{

}

Graph::Graph(int n)
{
    m_nodes = QVector<Node>(n);
}



GuiseGraphlet Graph::GuiseInit(int gletSize, int sid)
{
    QVector<bool> visited(nodeNum(), false);
    visited[sid] = true;

    QQueue<int> q;
    q.push_back(sid);

    GuiseGraphlet glet0;
    // bfs
    while(!q.empty())
    {
        int tmp = q.front();
        q.pop_front();

        glet0.push_back(tmp);
        if (glet0.size() == gletSize)
        {
            break;
        }

        Node* tmpNode = GetNode(tmp);
        for (int i = 0; i < tmpNode->degree(); i++)
        {
            if (!visited[tmpNode->childs[i]])
            {
                visited[tmpNode->childs[i]] = true;
                q.push_back(tmpNode->childs[i]);
            }
        }
    }

    if (glet0.size() < gletSize)
    {
        qDebug() << "Initialize failure.";
        qDebug() << glet0.size();
    }

    return glet0;
}

QVector<GuiseGraphlet> Graph::addVertex(const GuiseGraphlet &gt)
{
    GuiseGraphlet gt_ = gt;     // make a variable copy

    // Union set of adjacent nodes of remaining nodes
    // Use hash table is fast
    QSet<int> U;
    for (int i = 0; i < gt_.size(); i++)
    {
        Node* tmp = GetNode(gt_[i]);
        for (int j = 0; j < tmp->degree(); j++)
        {
            U.insert(tmp->childs[j]);
        }
    }

    QVector<GuiseGraphlet> neighbors;
    // Replace with a new node
    foreach(int u, U)
    {
        if (!gt_.contains(u))
        {
            gt_.append(u);
            neighbors.append(gt_);
            gt_.removeLast();
        }
    }

    return neighbors;
}

QVector<GuiseGraphlet> Graph::popNeighborGuise(const GuiseGraphlet &gt)
{
    QVector<GuiseGraphlet> neighbors;
    GuiseGraphlet gt_ = gt;     // make a variable copy

    int originalSize = gt_.size();
    for (int i = 0; i < originalSize; i++)
    {
        // remove the ith node
        int t = gt_[i];
        gt_.remove(i);

        if (connected(gt_))
        {
            //// A valid 3-node or 4-node Graphlet
            if (originalSize > 3)
            {
                neighbors.append(gt_);
            }

            // replacing ith node with a new node
            neighbors.append(addVertex(gt_));
        }

        // store
        gt_.insert(i, t);
    }

    if (originalSize < 5)
    {
        neighbors.append(addVertex(gt_));
    }

    return neighbors;
}

void Graph::GUISE(int sCount, int sid)
{
    qDebug() << "Apply GUISE at node: " << sid << ", sCount= " << sCount;

    GuiseGraphlet gx = GuiseInit(4, sid);
    if (gx.size() < 3)
    {
        return;
    }
    QVector<GuiseGraphlet> dgx = popNeighborGuise(gx);

    int sampled = 0;
    while (sampled < sCount)
    {
        // uniformly choose a neighbor gy of fx
        int randN = rand()%dgx.size();
        GuiseGraphlet gy = dgx[randN];

        QVector<GuiseGraphlet> dgy = popNeighborGuise(gy);

        float acProb = std::min((float)dgx.size()/dgy.size(), 1.f);
        float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (r < acProb)     // accept it
        {
            gx = gy;
            dgx = dgy;
        }
        sampled += 1;

        // insert graphlet found on current state
        GraphLet glet = inducedGraphlet(gx);
        int type = gletType(glet);

        // Add this graphlet to each node
        for(int i = 0; i < glet.size(); i++)
        {
            int nodeId = glet[i].first;
            m_gfds[nodeId][type] += 1;
        }
    }
}

void Graph::preGUISE()
{
    // initialize
    m_gfds = QVector<QVector<int>>(nodeNum(), QVector<int>(ALL_GRAPHLET, 0));

    // Get all connected components, each connected component has a vector of node
     QVector<QVector<int>> ccs = CCs();
//    QVector<int> ccs = CCs();
     qDebug() << "Connected components: " << ccs.size();

     // For each connect component, apply random walk
     for (int i = 0; i < ccs.size(); i++)
     {
         // don't need so much 1000 500
         int sCount = ccs[i].size()*1000; // 1000 * 500 *
         GUISE(sCount, ccs[i][0]);

//         int sCount = 100000;
//         GUISE(sCount, ccs[i]); // just use a node
     }
}

QVector<int> Graph::BfsLayersID(int sid, int neighborSize)
{
    // Reset layer id of all nodes
    for (int i = 0; i < nodeNum(); i++)
    {
        GetNode(i)->layerId = -1;
    }

    // first is level, second is node
    typedef QPair<int, int> tuple;

    // ensure each node is visited once
    QVector<bool> visited(nodeNum(), false);
    visited[sid] = true;

    QQueue<tuple> q;
    q.enqueue(tuple(0, sid));


    QVector<int> neighbors;
    // BFS starting from sid
    while(!q.isEmpty())
    {
        // top
        tuple t = q.head();
        int clevel = t.first;   // record this level
        int cnode = t.second;

        if (clevel >= neighborSize)
        {
            break;
        }

        neighbors.append(cnode);
        Node* n = GetNode(cnode);
        n->layerId = clevel;

        // pop
        q.dequeue();

        // push
        for (int i = 0; i< n->childs.size(); i++)
        {
            if (!visited[n->childs[i]])
            {
                visited[n->childs[i]] = true;
                q.enqueue(tuple(clevel+1, n->childs[i]));
            }
        }
    }

    return neighbors;
}

//QVector<float> Graph::GuiseGFD()
//{
//    preGUISE();
//    QVector<float> gfd(ALL_GRAPHLET, 0.f);

//    foreach (GraphLet glet, m_graphletSetGuise) {
//        gfd[gletType(glet)] += 1;
//    }

//    math_utils::normalize(gfd);
//    return gfd;
//}

bool Graph::isInduced(const GraphLet &glet)
{
    for (int i = 0; i < glet.size(); i++)
    {
        int u = glet[i].first;
        for (int j = i+1; j < glet.size(); j++)// loop for each pair of nodes in the graphlet
        {
            int v = glet[j].first;
            if (hasEdge(u, v) &&                // graph has this edge
            !(glet[i].second.contains(v)))      // but the graphlet don't
            {
//                printGraphlet(glet);
                return false;
            }
        }
    }

    return true;
}

//QVector<float> Graph::countAllGraphlets()
//{
//    QVector<QSet<GraphLet>> allGraphlets(ALL_GRAPHLET);

//    QVector<QVector<int>> ccs = CCs();
//    qDebug() << "Connected components: " << ccs.size();

//    for (int i = 0; i < ccs.size(); i++)    // for each connected components
//    {
//        QVector<int> cc = ccs[i];

//        for (int j = 0; j < cc.size(); j++) // for each node in one cc
//        {
//            QVector<QSet<GraphLet>> graphlets = GetGraphLets(cc[j]);

//            for (int k = 0; k < allGraphlets.size(); k++)   // for 29 types
//            {
//                allGraphlets[k] += graphlets[k];            // union set of each type
//            }
//        }
//    }



//    QVector<float> gfd(ALL_GRAPHLET);
//    for (int i = 0; i < gfd.size(); i++)
//    {
//        gfd[i] = allGraphlets[i].size();
//    }

//    utils::normalize(gfd);
//    return gfd;
//}

void Graph::findAllComb(QVector<QVector<bool>>& flags, QVector<bool> &flag, int i)
{
    if (i < flag.size())
    {
        flag[i] = false;
        if (i == flag.size()-1)
        {
            flags.push_back(flag);
        }
        findAllComb(flags, flag, i+1);

        flag[i] = true;
        if (i == flag.size()-1)
        {
            flags.push_back(flag);
        }
        findAllComb(flags, flag, i+1);
    }
}

Graph Graph::inducedSubgraph(const QVector<int> &nodeIds)
{
    Graph sub;
    for (int i = 0; i < nodeIds.size(); i++)
    {
        Node* node = GetNode(nodeIds[i]);
        sub.addNode(Node(node->x, node->y));// loss original node id
    }

    for (int i = 0; i < nodeIds.size(); i++)
    {
        for (int j = i + 1; j < nodeIds.size(); j++)
        {
            if (hasEdge(nodeIds[i], nodeIds[j]))
            {
//                sub.addEdge(nodeIds[i], nodeIds[j]);
                sub.addEdge(i, j);
            }
        }
    }
    return sub;
}

// without losing node id information
GraphLet Graph::inducedGraphlet(const QVector<int> &nodeIds)
{
    GraphLet glet;

    // add nodes to graphlet
    for (int i = 0; i < nodeIds.size(); i++)
    {
        GraphLetNode node;
        node.first = nodeIds[i];
        glet.append(node);
    }

    // add edges to graphlet
    for (int i = 0; i < nodeIds.size(); i++)
    {
        for (int j = i + 1; j < nodeIds.size(); j++)
        {
            if (hasEdge(nodeIds[i], nodeIds[j]))
            {
                glet[i].second.append(nodeIds[j]);
                glet[j].second.append(nodeIds[i]);
            }
        }
    }


    // sort here
    SortGraphLet(glet);
    return glet;
}

QVector<Graph> Graph::connectedSubgraphs(QVector<int> nodeIds)
{
    Graph induced = inducedSubgraph(nodeIds);
    QVector<QPair<int, int>> edges = induced.GetEdges();

    int n = induced.nodeNum();
    int m = edges.size();
    int all = int(pow(2, m));

    QVector<QVector<bool>> flags;
    QVector<bool> flag(m, false);
    // Loop for all possibilties of subgraphs
    findAllComb(flags, flag, 0);
    assert(flags.size() == all);

    QVector<Graph> subgraphs;
    for (int i = 0; i < flags.size(); i++)
    {
        Graph tmp;
        for (int j = 0; j < n; j++)
        {
            Node* n = induced.GetNode(j);
            tmp.addNode(Node(n->x, n->y));
        }

        // loop for each edge to see whether add this edge or not
        for (int j = 0; j < flags[i].size(); j++)
        {
            // edge j should be added
            if (flags[i][j])
            {
                tmp.addEdge(edges[j].first, edges[j].second);
            }
        }
        // graphlets must be connected
        if (tmp.connected())
        {
            subgraphs.append(tmp);
        }
    }

    qDebug() << "all connected spanning subgraphs of s(nodeIds)" << subgraphs.size();
    return subgraphs;
}

// This function can only be used in the original graph
QVector<QVector<int>> Graph::connectedInducedGraph(int d)
{
    // all node ids
    QVector<int> ids;
    for (int i = 0; i < m_nodes.size(); i++)
    {
        ids.append(i);
    }

    QVector<QVector<int>> dgraphs = combine(ids, d);
    QVector<QVector<int>> cdgraphs;
    for (int i = 0; i < dgraphs.size(); i++)
    {
        Graph tmp = inducedSubgraph(dgraphs[i]);
        if (tmp.connected())
        {
            cdgraphs.push_back(dgraphs[i]);
        }
    }
    return cdgraphs;
}


void Graph::dfs(int nodeId, int label, QVector<int> &labels, QVector<int> &cc)
{
    // set flag
    labels[nodeId] = label;

    // store in the connected component
    cc.push_back(nodeId);


//    qDebug() << nodeId;
    Node* node = GetNode(nodeId);
    // for each child
    for (int i = 0; i < node->degree(); i++)
    {
        if (labels[node->childs[i]] < 0)    // not visited yet
        {
            dfs(node->childs[i], label, labels, cc);
        }
    }
}

//void Graph::dfs(int nodeId, int label, QVector<int> &labels)
//{
//    // set flag
//    labels[nodeId] = label;


//    Node* node = GetNode(nodeId);
//    // for each child
//    for (int i = 0; i < node->degree(); i++)
//    {
//        if (labels[node->childs[i]] < 0)    // not visited yet
//        {
//            dfs(node->childs[i], label, labels);
//        }
//    }
//}


QVector<QVector<int> > Graph::CCs()
{
    // store all connect components
    QVector<QVector<int>> ccs;

    // not visited yet
    QVector<int> labels(nodeNum(), -1);

    int label = 0;
    // Find which one has most number of nodes
    for (int i = 0; i < nodeNum(); i++)
    {
       if (labels[i] < 0) {
           QVector<int> cc;
           dfs(i, label++, labels, cc);
           ccs.push_back(cc);
       }
    }

    return ccs;
}

//QVector<int> Graph::CCs()
//{
//    // store all connect components
//    QVector<int> ccs;

//    // not visited yet
//    QVector<int> labels(nodeNum(), -1);

//    int label = 0;
//    // Find which one has most number of nodes
//    for (int i = 0; i < nodeNum(); i++)
//    {
//       if (labels[i] < 0) {

//           // store current node
//           ccs.push_back(i);

//           dfs(i, label++, labels);
//       }
//    }

//    return ccs;
//}

//QVector<int> Graph::LCC()
//{
//    QVector<QVector<int>> ccs = CCs();

//    int maxSize = 0;
//    int maxIdx = -1;
//    for (int i = 0; i < ccs.size(); i++)
//    {
//        if (ccs[i].size() > maxSize)
//        {
//            maxSize = ccs[i].size();
//            maxIdx = i;
//        }
//    }

//    return ccs[maxIdx];
//}


bool Graph::connected()
{
//    QPair<int, QVector<int>> pair = CCLabels();
//    // only one connected components
//    return pair.first == 1;


    QVector<int> labels(nodeNum(), -1);

    QVector<int> cc;
    dfs(0, 0, labels, cc);
//    dfs(0, 0, labels);

    for (int i = 0; i < nodeNum(); i++)
    {
       if (labels[i] < 0) {
          return false;
       }
    }

    return true;
}

bool Graph::connected(const QVector<int> &nodeIds)
{
    Graph subgraph = inducedSubgraph(nodeIds);
    return subgraph.connected();
}

void Graph::setNodeNum(const int n)
{
    m_nodes = QVector<Node>(n);
}


int Graph::nodeNum() const
{
    return m_nodes.size();
}

QVector<QPair<int, int> > Graph::GetEdges()
{
    QVector<QPair<int, int>> edges;
    // for each edge
    for (int i = 0; i < nodeNum(); i++)
    {
        Node* node = GetNode(i);
        for (int j = 0; j < node->degree(); j++)
        {
            // only append one direction of edge
            QPair<int, int> edge1(i, node->childs[j]);
            QPair<int, int> edge2(node->childs[j], i);

            if (!edges.contains(edge1) && !edges.contains(edge2))
            {
                edges.append(edge1);
            }
        }
    }

    return edges;
}

int Graph::edgeNum(int d)
{
    int num = 0;

    if (d == 1)
    {
        // loop over the adjacency list
        for (int i = 0; i < nodeNum(); i++)
        {
            Node* tmp = GetNode(i);
            num += tmp->degree();
        }
    }
    else if (d == 2)
    {
       // loop over each edge
       for (int i = 0; i < nodeNum(); i++)
       {
           Node* u = GetNode(i);
           for (int j = 0; j < u->degree(); j++)
           {
               Node* v = GetNode(u->childs[j]);
               num += u->degree() + v->degree() -2;
           }
       }
    }
    else {
        qDebug() << "improper d!";
    }

#ifdef UNDIRECT_GRAPH
    num /= 2;
#endif

    return num;
}

Node* Graph::GetNode(int id)
{
    assert(id >= 0 && id < m_nodes.size());

    return &this->m_nodes[id];
}

void Graph::addNode(const Node &n)
{
    this->m_nodes.push_back(n);
}

void Graph::addEdge(int a, int b)
{
    assert(a >= 0 && a<m_nodes.size());
    assert(b >= 0 && b<m_nodes.size());

    m_nodes[a].addEdge(b);

#ifdef UNDIRECT_GRAPH
    m_nodes[b].addEdge(a);
#endif
}


bool Graph::hasEdge(int a, int b)
{
    if (a < 0 || a >= m_nodes.size() || b < 0 || b>= m_nodes.size())
    {
        return false;
    }
#ifdef UNDIRECT_GRAPH
    return m_nodes[a].hasEdge(b) && m_nodes[b].hasEdge(a);
#else
    return m_nodes[a].hasEdge(b);
#endif
}

void Graph:: renderGraph(QPainter *painter)
{
    painter->setPen(QColor(60,60,60));
    painter->setBrush(QBrush(QColor(90,60,220)));

    for(int i=0; i<m_nodes.size(); ++i)
    {
        painter->drawEllipse(QPoint(m_nodes[i].x,m_nodes[i].y),6,6);
        painter->drawText(m_nodes[i].x - 5, m_nodes[i].y-8,QString::number(i));

        for(int t=0; t<m_nodes[i].childs.size(); ++t)
        {
            int cid = m_nodes[i].childs[t];
            painter->drawLine(m_nodes[i].x,m_nodes[i].y, m_nodes[cid].x, m_nodes[cid].y);
        }
    }
}

//void Graph::renderCurGraphLets(QPainter *painter)
//{
//    if(m_cur_graphlets.isEmpty())
//        return;
//    painter->setPen(QPen(QColor(255,0,0),3));

//    painter->setBrush(QBrush(QColor(255,0,0)));
//    //    qDebug()<<"更新";
//    for(int i=0; i<m_cur_graphlets[cur_ID].size(); ++i)
//    {
//        int n0 = m_cur_graphlets[cur_ID][i].first;
//        //        qDebug()<<m_cur_graphlets[cur_ID][i].first;

//        painter->drawEllipse(QPoint(m_nodes[n0].x,m_nodes[n0].y),6,6);

//        for(int t=0; t<m_cur_graphlets[cur_ID][i].second.size(); ++t)
//        {
//            int n2 = m_cur_graphlets[cur_ID][i].second[t];
//            //            qDebug() << " " << n2;
//            painter->drawLine(m_nodes[n0].x,m_nodes[n0].y, m_nodes[n2].x, m_nodes[n2].y);
//        }
//    }

//}

void Graph::renderCurGraphLets(QPainter *painter)
{
    painter->setPen(QPen(QColor(255,0,0),3));
    painter->setBrush(QBrush(QColor(255,0,0)));

    if (m_graphlet != NULL)
    {
        for(int i=0; i<m_graphlet->nodeNum(); ++i)
        {
            Node* n1 = m_graphlet->GetNode(i);
            painter->drawEllipse(QPoint(n1->x,n1->y),6,6);

            for(int t=0; t< n1->degree(); ++t)
            {
                int j = n1->childs[t];
                Node* n2 = m_graphlet->GetNode(j);

                painter->drawLine(n1->x,n1->y, n2->x, n2->y);
            }
        }
    }
}

int Graph::HasNodeAt(QPoint p, float radius)
{
    int id = -1;
    float dist = INT_MAX;
    for(int i=0; i<m_nodes.size(); ++i)
    {
        float _dist = pow(m_nodes[i].x - p.x(),2) + pow(m_nodes[i].y-p.y(),2);

        if(_dist < dist && _dist< radius)
        {
            id = i;
            dist = _dist;
        }
    }
    return id;
}

void Graph::ShowNextGraphlet()
{
    if (!m_cur_graphlets.empty())
        cur_ID = (cur_ID+1) % m_cur_graphlets.size();
}

void Graph::clear()
{
    m_nodes.clear();
    m_cur_graphlets.clear();
}


QSet<GraphLet> Graph::SearchGraphLet(int gid, int sid)
{
    QTime time;
    time.start();

    QSet<GraphLet> glets;
    switch (gid)
    {
    case 0:
        glets = SearchGraphLet1(sid);
        break;
    case 1:
        glets = SearchGraphLet2(sid);
        break;
    case 2:
        glets = SearchGraphLet3(sid);
        break;
    case 3:
        glets = SearchGraphLet4(sid);
        break;
    case 4:
        glets = SearchGraphLet5(sid);
        break;
    case 5:
        glets = SearchGraphLet6(sid);
        break;
    case 6:
        glets = SearchGraphLet7(sid);
        break;
    case 7:
        glets = SearchGraphLet8(sid);
        break;
    case 8:
        glets = SearchGraphLet9(sid);
        break;
    case 9:
        glets = SearchGraphLet10(sid);
        break;
    case 10:
        glets = SearchGraphLet11(sid);
        break;
    case 11:
        glets = SearchGraphLet12(sid);
        break;
    case 12:
        glets = SearchGraphLet13(sid);
        break;
    case 13:
        glets = SearchGraphLet14(sid);
        break;
    case 14:
        glets = SearchGraphLet15(sid);
        break;
    case 15:
        glets = SearchGraphLet16(sid);
        break;
    case 16:
        glets = SearchGraphLet17(sid);
        break;
    case 17:
        glets = SearchGraphLet18(sid);
        break;
    case 18:
        glets = SearchGraphLet19(sid);
        break;
    case 19:
        glets = SearchGraphLet20(sid);
        break;
    case 20:
        glets = SearchGraphLet21(sid);
        break;
    case 21:
        glets = SearchGraphLet22(sid);
        break;
    case 22:
        glets = SearchGraphLet23(sid);
        break;
    case 23:
        glets = SearchGraphLet24(sid);
        break;
    case 24:
        glets = SearchGraphLet25(sid);
        break;
    case 25:
        glets = SearchGraphLet26(sid);
        break;
    case 26:
        glets = SearchGraphLet27(sid);
        break;
    case 27:
        glets = SearchGraphLet28(sid);
        break;
    case 28:
        glets = SearchGraphLet29(sid);
        break;
    default:
        assert(0);
        glets = SearchGraphLet1(sid);
        break;
    }
    return glets;
}

QSet<GraphLet> Graph::GetGraphLets(int sid)
{
    QSet<GraphLet> allGlets;
    for (int i = 0; i < ALL_GRAPHLET; i++)  // for each graphlet type
    {
        allGlets += SearchGraphLet(i, sid);
    }

    return allGlets;
}


//// TODO
//void Graph::preCount()
//{
//    qDebug() << "Start enumeration.";

//    m_graphletSet.clear();
//    // find graphlets on each node
//    for (int i = 0; i < nodeNum(); i++)
//    {
//        QSet<GraphLet> glets = GetGraphLets(i);
//        m_graphletSet += glets;
//    }
//}



void Graph::preCount()
{
    qDebug() << "Start enumeration.";

    m_gfds_enumer = QVector<QVector<int>>(nodeNum(), QVector<int>(ALL_GRAPHLET, 0));
    // find graphlets on each node
    for (int i = 0; i < nodeNum(); i++)
    {
        QSet<GraphLet> glets = GetGraphLets(i);
        foreach(const GraphLet& glet, glets)
        {
            int type = gletType(glet);
            m_gfds_enumer[i][type] ++;
        }
    }
}

//QVector<float > Graph::GetFeatureVectorAtLayer(int sid, int layerId)
//{
//    typedef QPair<int, int> tuple;  // first is level, second is node

//    QVector<bool> visited(m_nodes.size(), false);// ensure each node is visited once
//    visited[sid] = true;

//    QQueue<tuple> q;
//    q.enqueue(tuple(0, sid));

//    QVector<QSet<GraphLet>> layerGlets(ALL_GRAPHLET);
//    // BFS, limited in neighbor size of MAX_SEARCH_RANGE
//    while(!q.isEmpty())
//    {
//        // top
//        tuple t = q.head();
//        int clevel = t.first;
//        int cnode = t.second;

//        if (clevel == layerId)
//        {
//            // Accumulate graphlets found on this node
//            QVector<QSet<GraphLet>> glets = m_graphletSet[cnode];
//            for (int i = 0; i < glets.size(); i++)
//            {
//                layerGlets[i] += glets[i];
//            }
//        }
//        else if (clevel > layerId)
//        {
//            break;
//        }

//        // pop
//        q.dequeue();

//        // push
//        Node* n = GetNode(cnode);
//        for (int i = 0; i< n->childs.size(); i++)
//        {
//            if (!visited[n->childs[i]])
//            {
//                visited[n->childs[i]] = true;
//                q.enqueue(tuple(clevel+1, n->childs[i]));
//            }
//        }
//    }



//    QVector<float> f;
//    for (int i = 0; i < layerGlets.size(); i++)
//    {
//        f.append(layerGlets[i].size());
//    }
//    utils::normalize(f);

//    return f;
//}


//QVector<QSet<GraphLet> > Graph::GetNeighborGraphlets(int sid, int neighborSize)
//{
//    typedef QPair<int, int> tuple;  // first is level, second is node

//    QVector<bool> visited(m_nodes.size(), false);// ensure each node is visited once
//    visited[sid] = true;

//    QQueue<tuple> q;
//    q.enqueue(tuple(0, sid));

//    int llevel = 0;
//    QVector<QSet<GraphLet>> neighborGlets(ALL_GRAPHLET * neighborSize);
//    QVector<QSet<GraphLet>> cLevelGlets(ALL_GRAPHLET);
//    // BFS, limited in neighbor size of MAX_SEARCH_RANGE
//    while(!q.isEmpty())
//    {
//        // top
//        tuple t = q.head();
//        int clevel = t.first;
//        int cnode = t.second;

//        if (clevel != llevel)
//        {
//            // store graphlets found on last level
//            std::copy(cLevelGlets.begin(), cLevelGlets.end(), neighborGlets.begin() + llevel*ALL_GRAPHLET);
//            for (int i = 0; i < cLevelGlets.size(); i++)
//            {
//                cLevelGlets[i].clear();
//            }

//            llevel = clevel;
//            if (clevel >= neighborSize)
//            {
//                break;
//            }
//        }


//        // Accumulate graphlets found on this node
//        QVector<QSet<GraphLet>> glets = m_graphletSet[cnode];
//        for (int i = 0; i < glets.size(); i++)
//        {
//            cLevelGlets[i] += glets[i];
//        }

//        // pop
//        q.dequeue();

//        // push
//        Node* n = GetNode(cnode);
//        for (int i = 0; i< n->childs.size(); i++)
//        {
//            if (!visited[n->childs[i]])
//            {
//                visited[n->childs[i]] = true;
//                q.enqueue(tuple(clevel+1, n->childs[i]));
//            }
//        }
//    }

//    return neighborGlets;
//}

//QVector<QSet<GraphLet> > Graph::GetNeighborGraphletsAll(int sid, int neighborSize)
//{
//    typedef QPair<int, int> tuple;  // first is level, second is node


//    QVector<bool> visited(m_nodes.size(), false);// ensure each node is visited once
//    visited[sid] = true;

//    QQueue<tuple> q;
//    q.enqueue(tuple(0, sid));

//    int llevel = 0;
//    QVector<QSet<GraphLet>> allGlets(ALL_GRAPHLET);
//    // BFS, limited in neighbor size of MAX_SEARCH_RANGE
//    while(!q.isEmpty())
//    {
//        // top
//        tuple t = q.head();
//        int clevel = t.first;
//        int cnode = t.second;

//        if (clevel != llevel)
//        {
//            llevel = clevel;
//            if (clevel >= neighborSize)
//            {
//                break;
//            }
//        }

//        // accumulate graphlets found in this node
//        QVector<QSet<GraphLet>> glets = m_graphletSet[cnode];
//        for (int i = 0; i < glets.size(); i++)
//        {
//            allGlets[i] += glets[i];
//        }

//        // pop
//        q.dequeue();

//        // push
//        Node* n = GetNode(cnode);
//        for (int i = 0; i< n->childs.size(); i++)
//        {
//            if (!visited[n->childs[i]])
//            {
//                visited[n->childs[i]] = true;
//                q.enqueue(tuple(clevel+1, n->childs[i]));
//            }
//        }
//    }

//    return allGlets;
//}

QVector<float> Graph::GetfeatureVector(int sid, int neighborSize)
{
    // First get nodes that lies within the neighborhoods
    QVector<int> neighbors = BfsLayersID(sid, neighborSize);

    // Add gfd of each node to corresponding layers
    QVector<QVector<float>> layers(neighborSize, QVector<float>(ALL_GRAPHLET, 0));
    for (int i = 0; i < neighbors.size(); i++)
    {
        int nodeId = neighbors[i];
        Node* node = GetNode(nodeId);

        for (int j = 0; j < ALL_GRAPHLET; j++)
        {
            layers[node->layerId][j] += m_gfds_enumer[nodeId][j];
        }
    }

    // concat all layers
    QVector<float> featureVector;
    for (int i = 0; i < layers.size(); i++)
    {
        QVector<float> tmp = layers[i];
        math_utils::normalize(tmp);
        featureVector.append(tmp);
    }
//    utils::normalize(featureVector);

    assert(featureVector.size() == ALL_GRAPHLET*neighborSize);
    return featureVector;
}



//QVector<float> Graph::GetfeatureVectorAll(int sid, int neighborSize)
//{
//    QVector<QSet<GraphLet> > neighborGlets = GetNeighborGraphletsAll(sid, neighborSize);
//    assert(neighborGlets.size() == ALL_GRAPHLET);

//    GFD gfd(ALL_GRAPHLET);
//    int count = 0;
//    for (int i = 0; i < gfd.size(); i++)
//    {
//        int counti = neighborGlets[i].size();
//        gfd[i] = counti;
//        count += counti;
//    }


//    if (count != 0)
//    {
//        for (int i = 0; i < gfd.size(); i++)
//        {
//            gfd[i] /= (float)count;
//        }
//    }
//    return gfd;
//}


//void Graph::preRW()
//{
//    m_graphletSetGuise.clear();

//    QPair<int, QVector<int>> pair = CCLabels();
//    int labelNum = pair.first;
//    qDebug() << "number of connected components: " << labelNum;


//    EMkChain emc(*this);
//    int walkSteps = nodeNum();                   // node number
//    m_graphletSetGuise = emc.NB_SRW_GLETS(walkSteps);

//    qDebug() << "graphlets totally found: " << m_graphletSetGuise.size();
//}



QVector<float> Graph::GetfeatureVectorRW(int sid, int neighborSize)
{
    // / start BFS from s，which labels layerID for nodes which lie within layerNum.
    QVector<int> neighbors = BfsLayersID(sid, neighborSize);

    // Add gfd of each node to corresponding layers
    QVector<QVector<float>> layers(neighborSize, QVector<float>(ALL_GRAPHLET, 0));
    for (int i = 0; i < neighbors.size(); i++)
    {
        int nodeId = neighbors[i];
        Node* node = GetNode(nodeId);

        for (int j = 0; j < ALL_GRAPHLET; j++)
        {
            layers[node->layerId][j] += m_gfds[nodeId][j];
        }
    }

    // Concat all layers
    QVector<float> featureVector;
    for (int i = 0; i < layers.size(); i++)
    {
        QVector<float> tmp = layers[i];
        math_utils::normalize(tmp);
        featureVector.append(tmp);
    }
//    utils::normalize(featureVector);

    // The size of feature vector is L * 26
    assert(featureVector.size() == ALL_GRAPHLET*neighborSize);
    return featureVector;
}

//QVector<float> Graph::GetfeatureVectorRWAtLayer(int sid, int layerId)
//{
//    // reset layer id of all nodes
//    for (int i = 0; i < nodeNum(); i++)
//    {
//        GetNode(i)->layerId = -1;
//    }

//    // first is level, second is node
//    typedef QPair<int, int> tuple;

//    // ensure each node is visited once
//    QVector<bool> visited(nodeNum(), false);
//    visited[sid] = true;

//    QQueue<tuple> q;
//    q.enqueue(tuple(0, sid));

//    int neighborSize = layerId + 1;
//    // BFS starting from sid, only visit levels up to neighborSize
//    while(!q.isEmpty())
//    {
//        // top
//        tuple t = q.head();
//        int clevel = t.first;   // record this level
//        int cnode = t.second;

//        if (clevel >= neighborSize)
//        {
//            break;
//        }

//        Node* n = GetNode(cnode);
//        n->layerId = clevel;

//        // pop
//        q.dequeue();

//        // push
//        for (int i = 0; i< n->childs.size(); i++)
//        {
//            if (!visited[n->childs[i]])
//            {
//                visited[n->childs[i]] = true;
//                q.enqueue(tuple(clevel+1, n->childs[i]));
//            }
//        }
//    }


//    QVector<float> layer(ALL_GRAPHLET, 0);
//    foreach (const GraphLet& glet, m_graphletSetGuise) // m_graphletSet中有很多graphlet的节点都不正确, GetNode都只用到了那么几个节点（从0-connected components size）
//    {
//        int gType = gletType(glet);

//        for (int i = 0; i < glet.size(); i++)
//        {
//            if(GetNode(glet[i].first)->layerId == layerId) //当前层就没有加到应该加的graphlet
//            {
//                layer[gType] += 1;         // might take some time here
//            }
//        }
//    }

//    math_utils::normalize(layer);

//    assert(layer.size() == ALL_GRAPHLET);
//    return layer;
//}


//QVector<float> Graph::GetfeatureVectorAllRW(int sid, int neighborSize)
//{
//    typedef QPair<int, int> tuple;  // first is level, second is node


//    QVector<bool> visited(nodeNum(), false);// ensure each node is visited once
//    visited[sid] = true;

//    QQueue<tuple> q;
//    q.enqueue(tuple(0, sid));

//    int llevel = 0;
//    QVector<int> neighborNodes;
//    // BFS, limited in neighbor size of MAX_SEARCH_RANGE
//    while(!q.isEmpty())
//    {
//        // top
//        tuple t = q.head();
//        int clevel = t.first;
//        int cnode = t.second;


//        if (clevel != llevel)
//        {
//            llevel = clevel;
//            if (clevel >= neighborSize)
//            {
//                break;
//            }
//        }


//        // accumulate graphlets found in this node
//        neighborNodes.append(cnode);

//        // pop
//        q.dequeue();

//        // push
//        Node* n = GetNode(cnode);
//        for (int i = 0; i< n->childs.size(); i++)
//        {
//            if (!visited[n->childs[i]])
//            {
//                visited[n->childs[i]] = true;
//                q.enqueue(tuple(clevel+1, n->childs[i]));
//            }
//        }
//    }

//    // get subgraph induced by neighborhood
//    Graph g = inducedSubgraph(neighborNodes);


//    int walkSteps = (1-pow(3, neighborSize))/-2;                 // 1*(1-3^n/-2)
////    qDebug() << "walk steps:" << walkSteps;

//    EMkChain chain(g);
//    QVector<float> featureVector = chain.NB_SRW(walkSteps);
//    return featureVector;
//}




QSet<GraphLet> Graph::SearchGraphLet1(int sid)   // 输入的当前以那个点为起点找graphlet
{
    // g1: n0->n1, n1->n2
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);
    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];

        Node* t1 = GetNode(n1);
        for(int k=0; k<t1->childs.size(); ++k)
        {
            int n2 = t1->childs[k];
            if(n2 == n0)
                continue;

            GraphLet tmp;
            GraphLetNode m0,m1,m2;
            m0.first = n0;
            m1.first = n1;
            m2.first = n2;

            m0.second.append(n1);
            m1.second.append(n2);

            // for undirected graph
#ifdef UNDIRECT_GRAPH
            m1.second.append(n0);
            m2.second.append(n1);
#endif

            tmp.append(m0);
            tmp.append(m1);
            tmp.append(m2);

            if (isInduced(tmp)) // A real graphlet
            {
                SortGraphLet(tmp);
                glets.insert(tmp);
            }
        }
    }
    return glets;
}


QSet<GraphLet> Graph::SearchGraphLet2(int sid)
{
    // g1: n0->n1, n1->n2, n2->n0
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int k=0; k<t1->childs.size(); ++k)
        {
            int n2 = t1->childs[k];
            if (n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for (int m=0; m < t2->childs.size(); ++m) {
                int n3 = t2->childs[m];

                if(n3 == n0) {
                    GraphLet tmp;
                    GraphLetNode m0,m1,m2;
                    m0.first = n0;
                    m1.first = n1;
                    m2.first = n2;
                    m0.second.append(n1);
                    m1.second.append(n2);
                    m2.second.append(n0);

#ifdef UNDIRECT_GRAPH
                    m1.second.append(n0);
                    m2.second.append(n1);
                    m0.second.append(n2);
#endif
                    tmp.append(m0);
                    tmp.append(m1);
                    tmp.append(m2);

                    if (isInduced(tmp))
                    {
                        SortGraphLet(tmp);
                        glets.insert(tmp);
                    }
                }
            }
        }
    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet3(int sid)
{
    // g1: n0->n1, n1->n2, n2->n3
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int k=0; k<t1->childs.size(); ++k)
        {
            int n2 = t1->childs[k];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int m=0; m<t2->childs.size(); ++m)
            {
                int n3 = t2->childs[m];
                if (n3 == n1 || n3 == n0)
                    continue;

                GraphLet tmp;
                GraphLetNode m0, m1, m2, m3;
                m0.first = n0;
                m1.first = n1;
                m2.first = n2;
                m3.first = n3;

                m0.second.append(n1);
                m1.second.append(n2);
                m2.second.append(n3);

#ifdef UNDIRECT_GRAPH
                m1.second.append(n0);
                m2.second.append(n1);
                m3.second.append(n2);
#endif
                tmp.append(m0);
                tmp.append(m1);
                tmp.append(m2);
                tmp.append(m3);

                if (isInduced(tmp))
                {
                    SortGraphLet(tmp);
                    glets.insert(tmp);
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet4(int sid)
{
    // g4: n0->n1, n1->n2, n1->n3
    QSet<GraphLet> glets;

    int n0 = sid;   //
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int k=0; k<t1->childs.size(); ++k)
        {
            int n2 = t1->childs[k];
            if(n2 == n0)
                continue;
//            for(int m=0; m<t1->childs.size(); ++m)
            for(int m=k+1; m<t1->childs.size(); ++m)
            {
                int n3 = t1->childs[m];
                if(n3 == n0 || n3 == n2)
                    continue;
                GraphLet tmp;
                GraphLetNode m0,m1,m2,m3;
                m0.first = n0;
                m1.first = n1;
                m2.first = n2;
                m3.first = n3;
                m0.second.append(n1);
                m1.second.append(n2);
                m1.second.append(n3);

#ifdef UNDIRECT_GRAPH
                m1.second.append(n0);
                m2.second.append(n1);
                m3.second.append(n1);
#endif
                tmp.append(m0);
                tmp.append(m1);
                tmp.append(m2);
                tmp.append(m3);

                if (isInduced(tmp))
                {
                    SortGraphLet(tmp);
                    glets.insert(tmp);
                }
            }
        }
    }

    return glets;
}


QSet<GraphLet> Graph::SearchGraphLet5(int sid)
{
    // g1: n0->n1, n1->n2, n2->n3, n3->n0
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int k=0; k<t1->childs.size(); ++k)
        {
            int n2 = t1->childs[k];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int m=0; m<t2->childs.size(); ++m)
            {
                int n3 = t2->childs[m];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for(int n=0; n<t3->childs.size(); ++n)
                {
                    int n4 = t3->childs[n];
                    if (n4 == n0)
                    {
                        GraphLet tmp;
                        GraphLetNode m0, m1, m2, m3;
                        m0.first = n0;
                        m1.first = n1;
                        m2.first = n2;
                        m3.first = n3;

                        m0.second.append(n1);
                        m1.second.append(n2);
                        m2.second.append(n3);
                        m3.second.append(n0);

#ifdef UNDIRECT_GRAPH
                        m1.second.append(n0);
                        m2.second.append(n1);
                        m3.second.append(n2);
                        m0.second.append(n3);
#endif
                        tmp.append(m0);
                        tmp.append(m1);
                        tmp.append(m2);
                        tmp.append(m3);

                        if (isInduced(tmp))
                        {
                            SortGraphLet(tmp);
                            glets.insert(tmp);
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet6(int sid)
{
    // g4: n0->n1, n1->n2, n1->n3, n2->n3
    QSet<GraphLet> glets;

    int n0 = sid;   //
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int k=0; k<t1->childs.size(); ++k)
        {
            int n2 = t1->childs[k];
            if(n2 == n0)
                continue;
            Node* t2 = GetNode(n2);
            for(int m=0; m<t2->childs.size(); ++m)
            {
                int n3 = t2->childs[m];
                if(n3 == n0 || n3 == n1)
                    continue;

                Node* t3 = GetNode(n3);
                for (int n = 0; n < t3->childs.size(); ++n) {
                    int n4 = t3->childs[n];
                    if (n4 == n1) {
                        GraphLet tmp;
                        GraphLetNode m0,m1,m2,m3;
                        m0.first = n0;
                        m1.first = n1;
                        m2.first = n2;
                        m3.first = n3;
                        m0.second.append(n1);
                        m1.second.append(n2);
                        m2.second.append(n3);
                        m3.second.append(n1);

#ifdef UNDIRECT_GRAPH
                        m1.second.append(n0);
                        m2.second.append(n1);
                        m3.second.append(n2);
                        m1.second.append(n3);
#endif
                        tmp.append(m0);
                        tmp.append(m1);
                        tmp.append(m2);
                        tmp.append(m3);

                        if (isInduced(tmp))
                        {
                            SortGraphLet(tmp);
                            glets.insert(tmp);
                        }
                    }
                }
            }
        }
    }

    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet7(int sid)
{
    // g1: n0->n1, n1->n2, n2->n3, n3->n0, n3->n1
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int k=0; k<t1->childs.size(); ++k)
        {
            int n2 = t1->childs[k];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int m=0; m<t2->childs.size(); ++m)
            {
                int n3 = t2->childs[m];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for(int n=0; n<t3->childs.size(); ++n)
                {
                    int n4 = t3->childs[n];
                    if (n4 == n0)
                    {
                        for (int j = 0; j<t3->childs.size(); ++j) {
                            int n5 = t3->childs[j];
                            if (n5 == n1)
                            {

                                GraphLet tmp;
                                GraphLetNode m0, m1, m2, m3;
                                m0.first = n0;
                                m1.first = n1;
                                m2.first = n2;
                                m3.first = n3;

                                m0.second.append(n1);
                                m1.second.append(n2);
                                m2.second.append(n3);
                                m3.second.append(n0);
                                m3.second.append(n1);
#ifdef UNDIRECT_GRAPH
                                m1.second.append(n0);
                                m2.second.append(n1);
                                m3.second.append(n2);
                                m0.second.append(n3);
                                m1.second.append(n3);
#endif
                                tmp.append(m0);
                                tmp.append(m1);
                                tmp.append(m2);
                                tmp.append(m3);

                                if (isInduced(tmp))
                                {
                                    SortGraphLet(tmp);
                                    glets.insert(tmp);
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet8(int sid)
{
    // g8: n0->n1, n1->n2, n2->n0,
    // n2->n3,
    // n3->n0, n3->n1
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if (n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for (int k=0; k < t2->childs.size(); ++k) {
                int n3 = t2->childs[k];
                if(n3 == n0) {  // A loop
                    for (int m = 0; m < t2->childs.size(); ++m) {
                        n3 = t2->childs[m];
                        if (n3 == n0 || n3 == n1)
                            continue;

                        Node* t3 = GetNode(n3);
                        for (int n = 0; n < t3->childs.size(); ++n) {
                            int n4 = t3->childs[n];
                            if (n4 == n0) {
                                for (int o = 0; o < t3->childs.size(); ++o) {
                                    int n5 = t3->childs[o];
                                    if (n5 == n1) {
                                        GraphLet tmp;
                                        GraphLetNode m0,m1,m2, m3;
                                        m0.first = n0;
                                        m1.first = n1;
                                        m2.first = n2;
                                        m3.first = n3;
                                        m0.second.append(n1);
                                        m1.second.append(n2);
                                        m2.second.append(n0);
                                        m2.second.append(n3);
                                        m3.second.append(n0);
                                        m3.second.append(n1);
#ifdef UNDIRECT_GRAPH
                                        m1.second.append(n0);
                                        m2.second.append(n1);
                                        m0.second.append(n2);
                                        m3.second.append(n2);
                                        m0.second.append(n3);
                                        m1.second.append(n3);
#endif
                                        tmp.append(m0);
                                        tmp.append(m1);
                                        tmp.append(m2);
                                        tmp.append(m3);

                                        if (isInduced(tmp))
                                        {
                                            SortGraphLet(tmp);
                                            glets.insert(tmp);
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }
    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet9(int sid)
{
    // g1: n0->n1, n1->n2, n2->n3, n3->n4
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int k=0; k<t1->childs.size(); ++k)
        {
            int n2 = t1->childs[k];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int m=0; m<t2->childs.size(); ++m)
            {
                int n3 = t2->childs[m];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for (int n = 0; n<t3->childs.size(); ++n)
                {
                    int n4 = t3->childs[n];
                    if (n4 == n2 || n4 == n1 || n4 == n0)
                        continue;

                    GraphLet tmp;
                    GraphLetNode m0, m1, m2, m3, m4;
                    m0.first = n0;
                    m1.first = n1;
                    m2.first = n2;
                    m3.first = n3;
                    m4.first = n4;

                    m0.second.append(n1);
                    m1.second.append(n2);
                    m2.second.append(n3);
                    m3.second.append(n4);
#ifdef UNDIRECT_GRAPH
                    m1.second.append(n0);
                    m2.second.append(n1);
                    m3.second.append(n2);
                    m4.second.append(n3);
#endif
                    tmp.append(m0);
                    tmp.append(m1);
                    tmp.append(m2);
                    tmp.append(m3);
                    tmp.append(m4);

                    if (isInduced(tmp))
                    {
                        SortGraphLet(tmp);
                        glets.insert(tmp);
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet10(int sid)
{
    // g10: n0->n1, n1->n2,
    // n2->n3 n2->n4
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);
    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for (int j = 0; j<t1->childs.size(); j++) {
            int n2 = t1->childs[j];
            if (n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if(n3 == n0 || n3 == n1)
                    continue;

//                for(int m=0; m<t2->childs.size(); ++m)
                for(int m=k+1; m<t2->childs.size(); ++m)
                {
                    int n4 = t2->childs[m];
                    if(n4 == n0 || n4 == n1 || n4 == n3)
                        continue;

                    GraphLet tmp;
                    GraphLetNode m0,m1,m2,m3,m4;
                    m0.first = n0;
                    m1.first = n1;
                    m2.first = n2;
                    m3.first = n3;
                    m4.first = n4;

                    m0.second.append(n1);
                    m1.second.append(n2);
                    m2.second.append(n3);
                    m2.second.append(n4);
#ifdef UNDIRECT_GRAPH
                    m1.second.append(n0);
                    m2.second.append(n1);
                    m3.second.append(n2);
                    m4.second.append(n2);
#endif
                    tmp.append(m0);
                    tmp.append(m1);
                    tmp.append(m2);
                    tmp.append(m3);
                    tmp.append(m4);

                    if (isInduced(tmp))
                    {
                        SortGraphLet(tmp);
                        glets.insert(tmp);
                    }
                }
            }
        }
    }

    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet11(int sid)
{
    // g11
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];

//        for (int j = 0; j<t0->childs.size(); ++j)
        for (int j = i+1; j<t0->childs.size(); ++j)
        {
            int n2 = t0->childs[j];
            if (n2 == n0 || n2 == n1)
                continue;

//            for (int k = 0; k < t0->childs.size(); ++k)
            for (int k = j+1; k < t0->childs.size(); ++k)
            {
                int n3 = t0->childs[k];
                if (n3 == n0 || n3 == n1 || n3 == n2)
                    continue;

//                for (int m = 0; m < t0->childs.size(); ++m)
                for (int m = k+1; m < t0->childs.size(); ++m)
                {
                    int n4 = t0->childs[m];
                    if (n4 == n0 || n4 == n1 || n4 == n2 || n4 == n3)
                        continue;

                    GraphLet tmp;
                    GraphLetNode m0,m1,m2,m3,m4;
                    m0.first = n0;
                    m1.first = n1;
                    m2.first = n2;
                    m3.first = n3;
                    m4.first = n4;

                    m0.second.append(n1);
                    m0.second.append(n2);
                    m0.second.append(n3);
                    m0.second.append(n4);
#ifdef UNDIRECT_GRAPH
                    m1.second.append(n0);
                    m2.second.append(n0);
                    m3.second.append(n0);
                    m4.second.append(n0);
#endif
                    tmp.append(m0);
                    tmp.append(m1);
                    tmp.append(m2);
                    tmp.append(m3);
                    tmp.append(m4);

                    if (isInduced(tmp))
                    {
                        SortGraphLet(tmp);
                        glets.insert(tmp);
                    }
                }
            }
        }
    }

    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet12(int sid)
{
    // g12: n0->n1, n1->n2, n2->n3, n3->n4
    // n1->n3
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for (int l = 0; l<t3->childs.size(); ++l)
                {
                    int n4 = t3->childs[l];
                    if (n4 == n2 || n4 == n1 || n4 == n0)
                        continue;

                    for (int m = 0; m < t1->childs.size(); m++)
                    {
                        int n5 = t1->childs[m];
                        if (n5 == n3)
                        {
                            GraphLet tmp;
                            GraphLetNode m0, m1, m2, m3, m4;
                            m0.first = n0;
                            m1.first = n1;
                            m2.first = n2;
                            m3.first = n3;
                            m4.first = n4;

                            m0.second.append(n1);
                            m1.second.append(n2);
                            m1.second.append(n3);
                            m2.second.append(n3);
                            m3.second.append(n4);
#ifdef UNDIRECT_GRAPH
                            m1.second.append(n0);
                            m2.second.append(n1);
                            m3.second.append(n1);
                            m3.second.append(n2);
                            m4.second.append(n3);
#endif
                            tmp.append(m0);
                            tmp.append(m1);
                            tmp.append(m2);
                            tmp.append(m3);
                            tmp.append(m4);

                            if (isInduced(tmp))
                            {
                                SortGraphLet(tmp);
                                glets.insert(tmp);
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet13(int sid)
{
    // g10: n0->n1, n1->n2, n2->n3, n3->n4, n4->n2
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for (int j = 0; j<t1->childs.size(); j++)
        {
            int n2 = t1->childs[j];
            if (n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if(n3 == n0 || n3 == n1)
                    continue;

                Node* t3 = GetNode(n3);
                for(int m=0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if(n4 == n0 || n4 == n1 || n4 == n2)
                        continue;

                    Node* t4 = GetNode(n4);
                    for (int n=0; n < t4->childs.size(); ++n)
                    {
                        int n5 = t4->childs[n];
                        if(n5 == n2)
                        {
                            GraphLet tmp;
                            GraphLetNode m0,m1,m2,m3,m4;
                            m0.first = n0;
                            m1.first = n1;
                            m2.first = n2;
                            m3.first = n3;
                            m4.first = n4;

                            m0.second.append(n1);
                            m1.second.append(n2);
                            m2.second.append(n3);
                            m3.second.append(n4);
                            m4.second.append(n2);
#ifdef UNDIRECT_GRAPH
                            m1.second.append(n0);
                            m2.second.append(n1);
                            m3.second.append(n2);
                            m4.second.append(n3);
                            m2.second.append(n4);
#endif
                            tmp.append(m0);
                            tmp.append(m1);
                            tmp.append(m2);
                            tmp.append(m3);
                            tmp.append(m4);

                            if (isInduced(tmp))
                            {
                                SortGraphLet(tmp);
                                glets.insert(tmp);
                            }
                        }
                    }
                }
            }
        }
    }

    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet14(int sid)
{
    // g4: n0->n1, n1->n2, n2->n3, n3->n1
    // n1 -> n4
    QSet<GraphLet> glets;

    int n0 = sid;   //
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if(n3 == n0 || n3 == n1)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m < t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n1)
                    {
                        for (int n = 0; n<t1->childs.size(); ++n)
                        {
                            // n1's next child
                            n4 = t1->childs[n];
                            if (n4 == n0 || n4 == n2 || n4 == n3)
                                continue;

                            GraphLet tmp;
                            GraphLetNode m0,m1,m2,m3,m4;
                            m0.first = n0;
                            m1.first = n1;
                            m2.first = n2;
                            m3.first = n3;
                            m4.first = n4;

                            m0.second.append(n1);
                            m1.second.append(n2);
                            m1.second.append(n4);
                            m2.second.append(n3);
                            m3.second.append(n1);
#ifdef UNDIRECT_GRAPH
                            m1.second.append(n0);
                            m2.second.append(n1);
                            m4.second.append(n1);
                            m3.second.append(n2);
                            m1.second.append(n3);
#endif
                            tmp.append(m0);
                            tmp.append(m1);
                            tmp.append(m2);
                            tmp.append(m3);
                            tmp.append(m4);

                            if (isInduced(tmp))
                            {
                                SortGraphLet(tmp);
                                glets.insert(tmp);
                            }
                        }
                    }
                }
            }
        }
    }

    return glets;
}


QSet<GraphLet> Graph::SearchGraphLet15(int sid)
{
    // g15: n0->n1, n1->n2, n2->n3, n3->n4, n4->n0
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n2 || n4 == n1 || n4 == n0)
                        continue;

                    Node* t4 = GetNode(n4);
                    for (int n = 0; n<t4->childs.size(); ++n) {
                        int n5 = t4->childs[n];
                        if (n5 == n0) {
                            GraphLet tmp;
                            GraphLetNode m0, m1, m2, m3, m4;
                            m0.first = n0;
                            m1.first = n1;
                            m2.first = n2;
                            m3.first = n3;
                            m4.first = n4;

                            m0.second.append(n1);
                            m1.second.append(n2);
                            m2.second.append(n3);
                            m3.second.append(n4);
                            m4.second.append(n0);
#ifdef UNDIRECT_GRAPH
                            m1.second.append(n0);
                            m2.second.append(n1);
                            m3.second.append(n2);
                            m4.second.append(n3);
                            m0.second.append(n4);
#endif
                            tmp.append(m0);
                            tmp.append(m1);
                            tmp.append(m2);
                            tmp.append(m3);
                            tmp.append(m4);

                            if (isInduced(tmp))
                            {
                                SortGraphLet(tmp);
                                glets.insert(tmp);
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet16(int sid)
{
    // g16: n0->n1, n1->n2, n2->n3, n3->n4, n4->n1
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n2 || n4 == n1 || n4 == n0)
                        continue;

                    Node* t4 = GetNode(n4);
                    for (int n = 0; n<t4->childs.size(); ++n)
                    {
                        int n5 = t4->childs[n];
                        if (n5 == n1)
                        {
                            GraphLet tmp;
                            GraphLetNode m0, m1, m2, m3, m4;
                            m0.first = n0;
                            m1.first = n1;
                            m2.first = n2;
                            m3.first = n3;
                            m4.first = n4;

                            m0.second.append(n1);
                            m1.second.append(n2);
                            m2.second.append(n3);
                            m3.second.append(n4);
                            m4.second.append(n1);
#ifdef UNDIRECT_GRAPH
                            m1.second.append(n0);
                            m2.second.append(n1);
                            m3.second.append(n2);
                            m4.second.append(n3);
                            m1.second.append(n4);
#endif
                            tmp.append(m0);
                            tmp.append(m1);
                            tmp.append(m2);
                            tmp.append(m3);
                            tmp.append(m4);

                            if (isInduced(tmp))
                            {
                                SortGraphLet(tmp);
                                glets.insert(tmp);
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet17(int sid)
{
    // g17: n0->n1, n1->n2, n2->n3, n3->n4, n4->n1, n1->n3
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n2 || n4 == n1 || n4 == n0)
                        continue;

                    Node* t4 = GetNode(n4);
                    for (int n = 0; n<t4->childs.size(); ++n)
                    {
                        int n5 = t4->childs[n];
                        if (n5 == n1)
                        {
                            // Find n3
                            for (int o = 0; o < t1->childs.size(); ++o)
                            {
                                n5 = t1->childs[o];
                                if (n5 == n3)
                                {
                                    GraphLet tmp;
                                    GraphLetNode m0, m1, m2, m3, m4;
                                    m0.first = n0;
                                    m1.first = n1;
                                    m2.first = n2;
                                    m3.first = n3;
                                    m4.first = n4;

                                    m0.second.append(n1);
                                    m1.second.append(n2);
                                    m1.second.append(n3);
                                    m2.second.append(n3);
                                    m3.second.append(n4);
                                    m4.second.append(n1);
#ifdef UNDIRECT_GRAPH
                                    m1.second.append(n0);
                                    m2.second.append(n1);
                                    m3.second.append(n1);
                                    m3.second.append(n2);
                                    m4.second.append(n3);
                                    m1.second.append(n4);
#endif
                                    tmp.append(m0);
                                    tmp.append(m1);
                                    tmp.append(m2);
                                    tmp.append(m3);
                                    tmp.append(m4);

                                    if (isInduced(tmp))
                                    {
                                        SortGraphLet(tmp);
                                        glets.insert(tmp);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet18(int sid)
{
    // g18: n0->
    QSet<GraphLet> glets;

    int n0 = sid;   //
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if(n3 == n0 || n3 == n1)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m < t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n1)
                    {
                        for (int n = 0; n<t1->childs.size(); ++n)
                        {
                            // n1's next child
                            n4 = t1->childs[n];
                            if (n4 == n0 || n4 == n3 || n4 == n2)
                                continue;

                            Node* t4 = GetNode(n4);
                            for (int o = 0; o < t4->childs.size(); ++o)
                            {
                                int n5 = t4->childs[o];
                                if (n5 == n0)
                                {
                                    GraphLet tmp;
                                    GraphLetNode m0,m1,m2,m3,m4;
                                    m0.first = n0;
                                    m1.first = n1;
                                    m2.first = n2;
                                    m3.first = n3;
                                    m4.first = n4;

                                    m0.second.append(n1);
                                    m1.second.append(n2);
                                    m1.second.append(n4);
                                    m2.second.append(n3);
                                    m3.second.append(n1);
                                    m4.second.append(n0);
#ifdef UNDIRECT_GRAPH
                                    m1.second.append(n0);
                                    m2.second.append(n1);
                                    m4.second.append(n1);
                                    m3.second.append(n2);
                                    m1.second.append(n3);
                                    m0.second.append(n4);
#endif
                                    tmp.append(m0);
                                    tmp.append(m1);
                                    tmp.append(m2);
                                    tmp.append(m3);
                                    tmp.append(m4);

                                    if (isInduced(tmp))
                                    {
                                        SortGraphLet(tmp);
                                        glets.insert(tmp);
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }
    }

    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet19(int sid)
{
    // g16: n0->n1, n1->n2, n2->n3, n3->n4, n4->n1
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n2 || n4 == n1 || n4 == n0)
                        continue;

                    Node* t4 = GetNode(n4);
                    for (int n = 0; n<t4->childs.size(); ++n)
                    {
                        int n5 = t4->childs[n];
                        if (n5 == n1)
                        {
                            for (int o = 0; o < t2->childs.size(); ++o)
                            {
                                n5 = t2->childs[o];
                                if (n5 == n4)
                                {
                                    GraphLet tmp;
                                    GraphLetNode m0, m1, m2, m3, m4;
                                    m0.first = n0;
                                    m1.first = n1;
                                    m2.first = n2;
                                    m3.first = n3;
                                    m4.first = n4;

                                    m0.second.append(n1);
                                    m1.second.append(n2);
                                    m2.second.append(n3);
                                    m2.second.append(n4);
                                    m3.second.append(n4);
                                    m4.second.append(n1);
#ifdef UNDIRECT_GRAPH
                                    m1.second.append(n0);
                                    m2.second.append(n1);
                                    m3.second.append(n2);
                                    m4.second.append(n2);
                                    m4.second.append(n3);
                                    m1.second.append(n4);
#endif
                                    tmp.append(m0);
                                    tmp.append(m1);
                                    tmp.append(m2);
                                    tmp.append(m3);
                                    tmp.append(m4);

                                    if (isInduced(tmp))
                                    {
                                        SortGraphLet(tmp);
                                        glets.insert(tmp);
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet20(int sid)
{
    // g20: n0->n1, n1->n2, n2->n3, n3->n0, n0->n4, n4->n2
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for(int m=0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n0)
                    {
                        for (int n = 0; n < t0->childs.size(); ++n)
                        {
                            n4 = t0->childs[n];
                            if (n4 == n1 || n4 == n2 || n4 == n3)
                                continue;

                            Node* t4 = GetNode(n4);
                            for (int o = 0; o < t4->childs.size(); ++o)
                            {
                                int n5 = t4->childs[o];
                                if (n5 == n2)
                                {
                                    GraphLet tmp;
                                    GraphLetNode m0, m1, m2, m3, m4;
                                    m0.first = n0;
                                    m1.first = n1;
                                    m2.first = n2;
                                    m3.first = n3;
                                    m4.first = n4;

                                    m0.second.append(n1);
                                    m0.second.append(n4);
                                    m1.second.append(n2);
                                    m2.second.append(n3);
                                    m3.second.append(n0);
                                    m4.second.append(n2);
#ifdef UNDIRECT_GRAPH
                                    m1.second.append(n0);
                                    m4.second.append(n0);
                                    m2.second.append(n1);
                                    m3.second.append(n2);
                                    m0.second.append(n3);
                                    m2.second.append(n4);
#endif
                                    tmp.append(m0);
                                    tmp.append(m1);
                                    tmp.append(m2);
                                    tmp.append(m3);
                                    tmp.append(m4);

                                    if (isInduced(tmp))
                                    {
                                        SortGraphLet(tmp);
                                        glets.insert(tmp);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet21(int sid)
{
    // g21: n0->n1, n1->n2, n2->n3, n3->n4, n4->n0,
    // n0->n2
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n2 || n4 == n1 || n4 == n0)
                        continue;

                    Node* t4 = GetNode(n4);
                    for (int n = 0; n<t4->childs.size(); ++n)
                    {
                        int n5 = t4->childs[n];
                        if (n5 == n0)
                        {
                            //                            qDebug() << "loop 5";
                            for (int o = 0; o<t0->childs.size(); ++o)
                            {
                                n5 = t0->childs[o];
                                if (n5 == n2)
                                {
                                    GraphLet tmp;
                                    GraphLetNode m0, m1, m2, m3, m4;
                                    m0.first = n0;
                                    m1.first = n1;
                                    m2.first = n2;
                                    m3.first = n3;
                                    m4.first = n4;

                                    m0.second.append(n1);
                                    m0.second.append(n2);
                                    m1.second.append(n2);
                                    m2.second.append(n3);
                                    m3.second.append(n4);
                                    m4.second.append(n0);
#ifdef UNDIRECT_GRAPH
                                    m1.second.append(n0);
                                    m2.second.append(n0);
                                    m2.second.append(n1);
                                    m3.second.append(n2);
                                    m4.second.append(n3);
                                    m0.second.append(n4);
#endif
                                    tmp.append(m0);
                                    tmp.append(m1);
                                    tmp.append(m2);
                                    tmp.append(m3);
                                    tmp.append(m4);

                                    if (isInduced(tmp))
                                    {
                                        SortGraphLet(tmp);
                                        glets.insert(tmp);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}


QSet<GraphLet> Graph::SearchGraphLet22(int sid)
{
    // g22: n0->n1, n1->n2, n2->n3, n3->n0,
    // n0->n4, n4->n2,
    // n0->n2
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for(int m=0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n0)
                    {
                        for (int n = 0; n < t0->childs.size(); ++n)
                        {
                            n4 = t0->childs[n];
                            if (n4 == n2 || n4 == n3 || n4 == n1)
                                continue;
                            Node* t4 = GetNode(n4);
                            for (int o = 0; o < t4->childs.size(); ++o)
                            {
                                int n5 = t4->childs[o];
                                if (n5 == n2)
                                {
                                    for (int p = 0; p < t0->childs.size(); ++p)
                                    {
                                        n5 = t0->childs[p];
                                        if (n5 == n2)
                                        {
                                            GraphLet tmp;
                                            GraphLetNode m0, m1, m2, m3, m4;
                                            m0.first = n0;
                                            m1.first = n1;
                                            m2.first = n2;
                                            m3.first = n3;
                                            m4.first = n4;

                                            m0.second.append(n1);
                                            m1.second.append(n2);
                                            m2.second.append(n3);
                                            m3.second.append(n0);
                                            m0.second.append(n4);
                                            m4.second.append(n2);
                                            m0.second.append(n2);
#ifdef UNDIRECT_GRAPH
                                            m1.second.append(n0);
                                            m2.second.append(n1);
                                            m3.second.append(n2);
                                            m0.second.append(n3);
                                            m4.second.append(n0);
                                            m2.second.append(n4);
                                            m2.second.append(n0);
#endif
                                            tmp.append(m0);
                                            tmp.append(m1);
                                            tmp.append(m2);
                                            tmp.append(m3);
                                            tmp.append(m4);

                                            if (isInduced(tmp))
                                            {
                                                SortGraphLet(tmp);
                                                glets.insert(tmp);
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet23(int sid)
{
    // g23: n0->n1, n1->n2, n2->n3, n3->n1
    // n1->n4
    // n4->n2, n4->n3
    QSet<GraphLet> glets;

    int n0 = sid;   //
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;
            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if(n3 == n0 || n3 == n1)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m < t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n1)
                    {
                        for (int n = 0; n < t1->childs.size(); ++n)
                        {
                            n4 = t1->childs[n];
                            if (n4 == n2 || n4 == n0 || n4 == n3)
                                continue;

                            Node* t4 = GetNode(n4);
                            for (int o = 0; o < t4->childs.size(); ++o)
                            {
                                int n5 = t4->childs[o];
                                if (n5 == n2)
                                {
                                    for (int p = 0; p < t4->childs.size(); ++p)
                                    {
                                        int n6 = t4->childs[p];
                                        if (n6 == n3)
                                        {
                                            GraphLet tmp;
                                            GraphLetNode m0,m1,m2,m3, m4;
                                            m0.first = n0;
                                            m1.first = n1;
                                            m2.first = n2;
                                            m3.first = n3;
                                            m4.first = n4;

                                            m0.second.append(n1);
                                            m1.second.append(n2);
                                            m1.second.append(n4);
                                            m2.second.append(n3);
                                            m3.second.append(n1);
                                            m4.second.append(n2);
                                            m4.second.append(n3);
#ifdef UNDIRECT_GRAPH
                                            m1.second.append(n0);
                                            m2.second.append(n1);
                                            m4.second.append(n1);
                                            m3.second.append(n2);
                                            m1.second.append(n3);
                                            m2.second.append(n4);
                                            m3.second.append(n4);
#endif
                                            tmp.append(m0);
                                            tmp.append(m1);
                                            tmp.append(m2);
                                            tmp.append(m3);
                                            tmp.append(m4);

                                            if (isInduced(tmp))
                                            {
                                                SortGraphLet(tmp);
                                                glets.insert(tmp);
                                            }
                                        }
                                    }
                                }
                            }

                        }

                    }
                }
            }
        }
    }

    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet24(int sid)
{
    // g24: n0->n1, n1->n2, n2->n3, n3->n4, n4->n0
    // n0->n2, n0->n3
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n2 || n4 == n1 || n4 == n0)
                        continue;

                    Node* t4 = GetNode(n4);
                    for (int n = 0; n<t4->childs.size(); ++n) {
                        int n5 = t4->childs[n];
                        if (n5 == n0) {
                            for (int o = 0; o < t0->childs.size(); ++o)
                            {
                                int n6 = t0->childs[o];
                                if (n6 == n2)
                                {
                                    for (int p = 0; p < t0->childs.size(); ++p)
                                    {
                                        int n7= t0->childs[p];
                                        if (n7==n3)
                                        {
                                            GraphLet tmp;
                                            GraphLetNode m0, m1, m2, m3, m4;
                                            m0.first = n0;
                                            m1.first = n1;
                                            m2.first = n2;
                                            m3.first = n3;
                                            m4.first = n4;

                                            m0.second.append(n1);
                                            m0.second.append(n2);
                                            m0.second.append(n3);
                                            m1.second.append(n2);
                                            m2.second.append(n3);
                                            m3.second.append(n4);
                                            m4.second.append(n0);
#ifdef UNDIRECT_GRAPH
                                            m1.second.append(n0);
                                            m2.second.append(n0);
                                            m3.second.append(n0);
                                            m2.second.append(n1);
                                            m3.second.append(n2);
                                            m4.second.append(n3);
                                            m0.second.append(n4);
#endif
                                            tmp.append(m0);
                                            tmp.append(m1);
                                            tmp.append(m2);
                                            tmp.append(m3);
                                            tmp.append(m4);

                                            if (isInduced(tmp))
                                            {
                                                SortGraphLet(tmp);
                                                glets.insert(tmp);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet25(int sid)
{
    // g27: n0->n1, n1->n2, n2->n3, n3->n0, n0->n4,
    // n4->n1
    // n4->n3
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for(int m=0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n0)
                    {// A loop
                        for (int n = 0; n < t0->childs.size(); ++n)
                        {
                            n4 = t0->childs[n];
                            if (n4 == n1 || n4 == n2 || n4 == n3)
                                continue;

                            Node* t4 = GetNode(n4);
                            for (int o = 0; o < t4->childs.size(); ++o)
                            {
                                int n5 = t4->childs[o];
                                if (n5==n1)
                                {
                                    for (int p = 0; p < t4->childs.size(); ++p)
                                    {
                                        int n6 = t4->childs[p];
                                        if (n6==n3)
                                        {
                                            GraphLet tmp;
                                            GraphLetNode m0, m1, m2, m3, m4;
                                            m0.first = n0;
                                            m1.first = n1;
                                            m2.first = n2;
                                            m3.first = n3;
                                            m4.first = n4;

                                            m0.second.append(n1);
                                            m0.second.append(n4);
                                            m1.second.append(n2);
                                            m2.second.append(n3);
                                            m3.second.append(n0);
                                            m4.second.append(n1);
                                            m4.second.append(n3);
    #ifdef UNDIRECT_GRAPH
                                            m1.second.append(n0);
                                            m4.second.append(n0);
                                            m2.second.append(n1);
                                            m3.second.append(n2);
                                            m0.second.append(n3);
                                            m1.second.append(n4);
                                            m3.second.append(n4);
    #endif
                                            tmp.append(m0);
                                            tmp.append(m1);
                                            tmp.append(m2);
                                            tmp.append(m3);
                                            tmp.append(m4);

                                            if (isInduced(tmp))
                                            {
                                                SortGraphLet(tmp);
                                                glets.insert(tmp);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet26(int sid)
{
    // g26: n0->n1, n1->n2, n2->n3, n3->n0,
    // n0->n4,n4->n2
    // n0->n2,
    // n4->n3
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for(int m=0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n0)
                    {
                        for (int n = 0; n < t0->childs.size(); ++n)
                        {
                            n4 = t0->childs[n];
                            if (n4 == n2 || n4 == n3 || n4 == n1)
                                continue;
                            Node* t4 = GetNode(n4);
                            for (int o = 0; o < t4->childs.size(); ++o)
                            {
                                int n5 = t4->childs[o];
                                if (n5 == n2)
                                {
                                    for (int p = 0; p < t0->childs.size(); ++p)
                                    {
                                        n5 = t0->childs[p];
                                        if (n5 == n2)
                                        {
                                            for (int q = 0; q < t4->childs.size(); ++q)
                                            {
                                                n5 = t4->childs[q];
                                                if (n5 == n3)
                                                {
                                                    GraphLet tmp;
                                                    GraphLetNode m0, m1, m2, m3, m4;
                                                    m0.first = n0;
                                                    m1.first = n1;
                                                    m2.first = n2;
                                                    m3.first = n3;
                                                    m4.first = n4;

                                                    m0.second.append(n1);
                                                    m1.second.append(n2);
                                                    m2.second.append(n3);
                                                    m3.second.append(n0);
                                                    m0.second.append(n4);
                                                    m4.second.append(n2);
                                                    m4.second.append(n3);
                                                    m0.second.append(n2);
#ifdef UNDIRECT_GRAPH
                                                    m1.second.append(n0);
                                                    m2.second.append(n1);
                                                    m3.second.append(n2);
                                                    m0.second.append(n3);
                                                    m4.second.append(n0);
                                                    m2.second.append(n4);
                                                    m3.second.append(n4);
                                                    m2.second.append(n0);
#endif
                                                    tmp.append(m0);
                                                    tmp.append(m1);
                                                    tmp.append(m2);
                                                    tmp.append(m3);
                                                    tmp.append(m4);

                                                    if (isInduced(tmp))
                                                    {
                                                        SortGraphLet(tmp);
                                                        glets.insert(tmp);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet27(int sid)
{
    // g27: n0->n1, n1->n2, n2->n3, n3->n0, n0->n4,
    // n4->n1
    // n4->n2
    // n4->n3
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for(int m=0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n0)
                    {// A loop
                        for (int n = 0; n < t0->childs.size(); ++n)
                        {
                            n4 = t0->childs[n];
                            if (n4 == n1 || n4 == n2 || n4 == n3)
                                continue;

                            Node* t4 = GetNode(n4);
                            for (int o = 0; o < t4->childs.size(); ++o)
                            {
                                int n5 = t4->childs[o];
                                if (n5==n1)
                                {
                                    for (int p = 0; p < t4->childs.size(); ++p)
                                    {
                                        int n6 = t4->childs[p];
                                        if (n6==n2)
                                        {
                                            for (int q = 0; q < t4->childs.size(); ++q)
                                            {
                                                int n7 = t4->childs[q];
                                                if (n7==n3)
                                                {
                                                    GraphLet tmp;
                                                    GraphLetNode m0, m1, m2, m3, m4;
                                                    m0.first = n0;
                                                    m1.first = n1;
                                                    m2.first = n2;
                                                    m3.first = n3;
                                                    m4.first = n4;

                                                    m0.second.append(n1);
                                                    m0.second.append(n4);
                                                    m1.second.append(n2);
                                                    m2.second.append(n3);
                                                    m3.second.append(n0);
                                                    m4.second.append(n1);
                                                    m4.second.append(n2);
                                                    m4.second.append(n3);
        #ifdef UNDIRECT_GRAPH
                                                    m1.second.append(n0);
                                                    m4.second.append(n0);
                                                    m2.second.append(n1);
                                                    m3.second.append(n2);
                                                    m0.second.append(n3);
                                                    m1.second.append(n4);
                                                    m2.second.append(n4);
                                                    m3.second.append(n4);
        #endif
                                                    tmp.append(m0);
                                                    tmp.append(m1);
                                                    tmp.append(m2);
                                                    tmp.append(m3);
                                                    tmp.append(m4);

                                                    if (isInduced(tmp))
                                                    {
                                                        SortGraphLet(tmp);
                                                        glets.insert(tmp);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet28(int sid)
{
    // g28: n0->n1, n1->n2, n2->n3, n3->n0,
    // n0->n4,n4->n2
    // n0->n2,
    // n4->n3, n4->n1
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for(int m=0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n0)   // A quadra loop
                    {
                        for (int n = 0; n < t0->childs.size(); ++n)
                        {
                            n4 = t0->childs[n];
                            if (n4 == n2 || n4 == n3 || n4 == n1)
                                continue;
                            Node* t4 = GetNode(n4);
                            for (int o = 0; o < t4->childs.size(); ++o)
                            {
                                int n5 = t4->childs[o];
                                if (n5 == n2)   // A bow
                                {
                                    for (int p = 0; p < t0->childs.size(); ++p)
                                    {
                                        n5 = t0->childs[p];
                                        if (n5 == n2)   // A line
                                        {

                                            // 4-3, 4-1
                                            for (int q = 0; q < t4->childs.size(); ++q)
                                            {
                                                n5 = t4->childs[q];
                                                if (n5 == n3)
                                                {
                                                    for (int r = 0; r < t4->childs.size(); ++r)
                                                    {
                                                        n5 = t4->childs[r];
                                                        if (n5 == n1)
                                                        {
                                                            GraphLet tmp;
                                                            GraphLetNode m0, m1, m2, m3, m4;
                                                            m0.first = n0;
                                                            m1.first = n1;
                                                            m2.first = n2;
                                                            m3.first = n3;
                                                            m4.first = n4;

                                                            m0.second.append(n1);
                                                            m0.second.append(n4);
                                                            m0.second.append(n2);
                                                            m1.second.append(n2);
                                                            m2.second.append(n3);
                                                            m3.second.append(n0);
                                                            m4.second.append(n1);
                                                            m4.second.append(n2);
                                                            m4.second.append(n3);
#ifdef UNDIRECT_GRAPH
                                                            m1.second.append(n0);
                                                            m4.second.append(n0);
                                                            m2.second.append(n0);
                                                            m2.second.append(n1);
                                                            m3.second.append(n2);
                                                            m0.second.append(n3);
                                                            m1.second.append(n4);
                                                            m2.second.append(n4);
                                                            m3.second.append(n4);
#endif

                                                            tmp.append(m0);
                                                            tmp.append(m1);
                                                            tmp.append(m2);
                                                            tmp.append(m3);
                                                            tmp.append(m4);

                                                            if (isInduced(tmp))
                                                            {
                                                                SortGraphLet(tmp);
                                                                glets.insert(tmp);
                                                            }
                                                        }
                                                    }

                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }

    }
    return glets;
}

QSet<GraphLet> Graph::SearchGraphLet29(int sid)
{
    // g21: n0->n1, n1->n2, n2->n3, n3->n4, n4->n0
    // n0->n2, n0->n3
    // n1->n3, n1->n4
    // n2->n4
    QSet<GraphLet> glets;

    int n0 = sid;
    Node* t0 = GetNode(n0);

    for(int i = 0; i<t0->childs.size(); ++i)
    {
        int n1 = t0->childs[i];
        Node* t1 = GetNode(n1);

        for(int j=0; j<t1->childs.size(); ++j)
        {
            int n2 = t1->childs[j];
            if(n2 == n0)
                continue;

            Node* t2 = GetNode(n2);
            for(int k=0; k<t2->childs.size(); ++k)
            {
                int n3 = t2->childs[k];
                if (n3 == n1 || n3 == n0)
                    continue;

                Node* t3 = GetNode(n3);
                for (int m = 0; m<t3->childs.size(); ++m)
                {
                    int n4 = t3->childs[m];
                    if (n4 == n2 || n4 == n1 || n4 == n0)
                        continue;

                    Node* t4 = GetNode(n4);
                    for (int n = 0; n<t4->childs.size(); ++n) {
                        int n5 = t4->childs[n];
                        if (n5 == n0) { // A penta loop


                            for (int o = 0; o < t0->childs.size(); ++o)
                            {
                                if (t0->childs[o] == n2)
                                {
                                    for (int p = 0; p < t0->childs.size(); ++p)
                                    {
                                        if (t0->childs[p] == n3)
                                        {
                                            for (int q = 0; q < t1->childs.size(); ++q)
                                            {
                                                if (t1->childs[q] == n3)
                                                {
                                                    for (int r = 0; r < t1->childs.size(); ++r)
                                                    {
                                                        if (t1->childs[r] == n4)
                                                        {
                                                            for (int s = 0; s < t2->childs.size(); ++s)
                                                            {
                                                                if (t2->childs[s] == n4)
                                                                {
                                                                    GraphLet tmp;
                                                                    GraphLetNode m0, m1, m2, m3, m4;
                                                                    m0.first = n0;
                                                                    m1.first = n1;
                                                                    m2.first = n2;
                                                                    m3.first = n3;
                                                                    m4.first = n4;

                                                                    m0.second.append(n1);
                                                                    m0.second.append(n2);
                                                                    m0.second.append(n3);
                                                                    m0.second.append(n4);
                                                                    m1.second.append(n2);
                                                                    m1.second.append(n3);
                                                                    m1.second.append(n4);
                                                                    m2.second.append(n3);
                                                                    m2.second.append(n4);
                                                                    m3.second.append(n4);
                #ifdef UNDIRECT_GRAPH
                                                                    m1.second.append(n0);
                                                                    m2.second.append(n0);
                                                                    m3.second.append(n0);
                                                                    m4.second.append(n0);
                                                                    m2.second.append(n1);
                                                                    m3.second.append(n1);
                                                                    m4.second.append(n1);
                                                                    m3.second.append(n2);
                                                                    m4.second.append(n2);
                                                                    m4.second.append(n3);
                #endif
                                                                    tmp.append(m0);
                                                                    tmp.append(m1);
                                                                    tmp.append(m2);
                                                                    tmp.append(m3);
                                                                    tmp.append(m4);

                                                                    if (isInduced(tmp))
                                                                    {
                                                                        SortGraphLet(tmp);
                                                                        glets.insert(tmp);
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    return glets;
}

int Graph::DedupGraphLets(QVector<GraphLet> &sortedGraphlets)
{
    QTime time;
    time.start();

    // Fist sort between graphlets
    std::sort(sortedGraphlets.begin(), sortedGraphlets.end(),
              [](GraphLet& glet1, GraphLet& glet2)
    {
        int nodeNum1 = glet1.size();
        int nodeNum2 = glet2.size();
//        if (nodeNum1 != nodeNum2)
//        {
//            return nodeNum1 < nodeNum2;
//        }
        // compare node id
        for (int i = 0; i < nodeNum1; i++)
        {
            if (glet1[i].first != glet2[i].first)
            {
                return glet1[i].first < glet2[i].first;
            }
        }
        // compare edge
        for (int i = 0; i < nodeNum1; i++)
        {
            int edgeNum1 = glet1[i].second.size();
            int edgeNum2 = glet2[i].second.size();
            if (edgeNum1 != edgeNum2)
            {
                return edgeNum1 < edgeNum2;
            }

            for (int j = 0; j < edgeNum1; j++)
            {
                if (glet1[i].second[j] != glet2[i].second[j])
                {
                    return glet1[i].second[j] < glet2[i].second[j];
                }
            }
        }

        // two graphlets are the same, no matter what
        return false;
    });

//#ifdef _DEBUG
//    qDebug() << "sort graphlets cost " << time.elapsed() << "ms";
//#endif

    time.start();
    // Then find consecutive duplicated elements
//    int n = std::unique(allGraphlets.begin(), allGraphlets.end()) - allGraphlets.begin();
//    QVector<GraphLet>::iterator iter = std::unique(allGraphlets.begin(), allGraphlets.end());
    QVector<GraphLet>::iterator iter = std::unique(sortedGraphlets.begin(), sortedGraphlets.end(),
                                                   [](GraphLet glet1, GraphLet glet2)
    {
        int nodeNum = glet1.size();
        // compare node id
        for (int i = 0; i < nodeNum; i++)
        {
            if (glet1[i].first != glet2[i].first)
            {
                return false;
            }
        }
        // compare edge
        for (int i = 0; i < nodeNum; i++)
        {
            int edgeNum1 = glet1[i].second.size();
            int edgeNum2 = glet2[i].second.size();
            if (edgeNum1 != edgeNum2)
            {
                return false;
            }

            for (int j = 0; j < edgeNum1; j++)
            {
                if (glet1[i].second[j] != glet2[i].second[j])
                {
                    return false;
                }
            }
        }

        return true;
    });


//#ifdef _DEBUG
//    qDebug() << "unique operation costs " << time.elapsed() << "ms";
//#endif
    int dupN = sortedGraphlets.end() - iter;
    sortedGraphlets.erase(iter, sortedGraphlets.end());

    return dupN;
}


void Graph::SortGraphLet(GraphLet& glet) {

#ifndef UNDIRECT_GRAPH
    // if the graph is directed, we synmetric the graphlet for de-duplication
    for (int i = 0; i < glet.size(); i++)
    {
        int ni = glet[i].first;
        for (int j = 0; j < glet[i].second.size(); j++)
        {
            int  nj = glet[i].second[j];
            // if (ni, nj) exists
            for (int k = 0; k < glet.size(); k++)
            {
                // insert (nj, ni)
                if (glet[k].first == nj)
                {
                    if(!glet[nj].second.contains(ni))
                    {
                        glet[nj].second.append(ni);
                    }
                }
            }
        }
    }
#endif

    // sort nodes within the graphlet
    std::sort(glet.begin(), glet.end(),
              [](GraphLetNode& gnode1, GraphLetNode&gnode2)
    {
        // sort
        return gnode1.first < gnode2.first;
    }
    );

    // sort the out edges of each node
    for (int i = 0; i < glet.size(); i++)
    {
        std::sort(glet[i].second.begin(), glet[i].second.end(), [](int a, int b)
        {
            return a < b;
        });
    }

}

bool Graph::GraphletEqual(const GraphLet &glet1, const GraphLet &glet2)
{

    if (glet1.size() != glet2.size())
    {
        return false;
    }
    for (int i = 0; i < glet1.size(); i++)
    {
        if (glet1[i].first != glet2[i].first)
        {
            return false;
        }
        for (int j = 0; j < glet1[i].second.size(); j++)
        {
            if (glet1[i].second[j] != glet2[i].second[j])
            {
                return false;
            }
        }
    }

    return true;
}


int graphType(Graph glet)
{
    // get degrees of each node
    QVector<int> signature;
    for (int i = 0; i < glet.nodeNum(); i++)
    {
        signature.append(glet.GetNode(i)->degree());
    }

    // sort the degrees and return sorted index
    QVector<int> indices = math_utils::sortIdx(signature);


    // 3 nodes
    if (signature == QVector<int>({1, 1, 2}))
    {
        return 0;
    }
    else if (signature == QVector<int>({2, 2, 2}))
    {
        return 1;
    }


    // 4 nodes
    else if (signature == QVector<int>({1, 1, 2, 2}))
    {
        return 0;
    }
    else if (signature == QVector<int>({1, 1, 1, 3}))
    {
        return 1;
    }
    else if (signature == QVector<int>({2, 2, 2, 2}))
    {
        return 2;
    }
    else if (signature == QVector<int>({1, 2, 2, 3}))
    {
        return 3;
    }
    else if (signature == QVector<int>({2, 2, 3, 3}))
    {
        return 4;
    }
    else if (signature == QVector<int>({3, 3, 3, 3}))
    {
        return 5;
    }


    // 5 nodes
    else if (signature == QVector<int>({1, 1, 2, 2, 2}))// g51
    {
        return 0;
    }
    else if (signature == QVector<int>({1, 1, 1, 2, 3}))// g52
    {
        return 1;
    }
    else if (signature == QVector<int>({1, 1, 1, 1, 4}))// g53
    {
        return 2;
    }
    else if (signature == QVector<int>({1, 1, 2, 3, 3}))// g54
    {
        return 3;
    }

    else if (signature == QVector<int>({1, 2, 2, 2, 3}))
    {
        int index1 = indices[0];
        Node* node1 = glet.GetNode(index1);

        if (glet.GetNode(node1->childs[0])->degree() == 2)
        {
            return 4;// g55
        }
        else
        {
            return 7;// g58
        }
    }
    else if (signature == QVector<int>({1, 1, 2, 2, 4}))// g56
    {
        return 5;
    }
    else if (signature == QVector<int>({2, 2, 2, 2, 2}))// g57
    {
        return 6;
    }
    else if (signature == QVector<int>({1, 2, 2, 3, 4}))// g59
    {
        return 8;
    }
    else if (signature == QVector<int>({2, 2, 2, 2, 4}))// g510
    {
        return 9;
    }
    else if (signature == QVector<int>({1, 2, 3, 3, 3}))//g511
    {
        return 10;
    }

    else if (signature == QVector<int>({2, 2, 2, 3, 3}))
    {
        int index5 = indices[4];
        Node* node5 = glet.GetNode(index5);

        for (int i = 0; i < node5->degree(); i++)
        {
            int index = node5->childs[i];
            if (glet.GetNode(index)->degree() == 3)
            {
                return 12;//g513
            }
        }
        return 11;//g512
    }

    else if (signature == QVector<int>({2, 2, 2, 4, 4}))//g514
    {
        return 13;
    }

    else if (signature == QVector<int>({1, 3, 3, 3, 4}))//g515
    {
        return 14;
    }
    else if (signature == QVector<int>({2, 2, 3, 3, 4}))//g516
    {
        return 15;
    }
    else if (signature == QVector<int>({2, 3, 3, 3, 3}))//g517
    {
        return 16;
    }
    else if (signature == QVector<int>({2, 3, 3, 4, 4}))//g518
    {
        return 17;
    }
    else if (signature == QVector<int>({3, 3, 3, 3, 4}))//g519
    {
        return 18;
    }
    else if (signature == QVector<int>({3, 3, 4, 4, 4}))// g520
    {
        return 19;
    }
    else if (signature == QVector<int>({4, 4, 4, 4, 4}))// g521
    {
        return 20;
    }
    else
    {
        qDebug() << "undesired node number: " << glet.nodeNum();
        return -1;
    }
}


int gletTypeNum(int k)
{
    if (k == 3)
    {
        return 2;
    }
    else if (k == 4)
    {
        return 6;
    }
    else if (k == 5)
    {
        return 21;
    }
}


//QVector<int> gletDegs(Graph glet)
//{
//    QVector<int> degs;
//    for (int i = 0; i < glet.nodeNum(); i++)
//    {
//        degs.append(glet.GetNode(i)->degree());
//    }
//    return degs;
//}

void printGraphlet(const GraphLet &graphlet)
{
    // print adjacency list
    for (int i = 0; i < graphlet.size(); i++)
    {
        QString str;
        str += QString::number(graphlet[i].first) + ": ";
        for (int j = 0; j < graphlet[i].second.size(); j++)
        {
            str += QString::number(graphlet[i].second[j]) + "=>";
        }
        qDebug() << str;
    }
    qDebug() << "---------------------------";
}

int gletType(GraphLet glet)
{
    // get degrees of each node
    QVector<int> signature;
    for (int i = 0; i < glet.size(); i++)
    {
        signature.append(glet[i].second.size());
    }

    // sort the degrees and return sorted index
    QVector<int> indices = math_utils::sortIdx(signature);


    // 3 nodes
    if (signature == QVector<int>({1, 1, 2}))
    {
        return 0;
    }
    else if (signature == QVector<int>({2, 2, 2}))
    {
        return 1;
    }


    // 4 nodes
    else if (signature == QVector<int>({1, 1, 2, 2}))
    {
        return 2;
    }
    else if (signature == QVector<int>({1, 1, 1, 3}))
    {
        return 3;
    }
    else if (signature == QVector<int>({2, 2, 2, 2}))
    {
        return 4;
    }
    else if (signature == QVector<int>({1, 2, 2, 3}))
    {
        return 5;
    }
    else if (signature == QVector<int>({2, 2, 3, 3}))
    {
        return 6;
    }
    else if (signature == QVector<int>({3, 3, 3, 3}))
    {
        return 7;
    }


    // 5 nodes
    else if (signature == QVector<int>({1, 1, 2, 2, 2}))// g51
    {
        return 8;
    }
    else if (signature == QVector<int>({1, 1, 1, 2, 3}))// g52
    {
        return 9;
    }
    else if (signature == QVector<int>({1, 1, 1, 1, 4}))// g53
    {
        return 10;
    }
    else if (signature == QVector<int>({1, 1, 2, 3, 3}))// g54
    {
        return 11;
    }

    else if (signature == QVector<int>({1, 2, 2, 2, 3}))
    {
        int index1 = indices[0];
        GraphLetNode node1 = glet[index1];
        int index2 = findGletNode(glet, node1.second[0]);

        // degree == 2
        if (glet[index2].second.size() == 2)
        {
            return 12;// g55
        }
        else
        {
            return 15;// g58
        }
    }
    else if (signature == QVector<int>({1, 1, 2, 2, 4}))// g56
    {
        return 13;
    }
    else if (signature == QVector<int>({2, 2, 2, 2, 2}))// g57
    {
        return 14;
    }
    else if (signature == QVector<int>({1, 2, 2, 3, 4}))// g59
    {
        return 16;
    }
    else if (signature == QVector<int>({2, 2, 2, 2, 4}))// g510
    {
        return 17;
    }
    else if (signature == QVector<int>({1, 2, 3, 3, 3}))//g511
    {
        return 18;
    }

    else if (signature == QVector<int>({2, 2, 2, 3, 3}))
    {
        int index5 = indices[4];
        GraphLetNode node5 = glet[index5];

        for (int i = 0; i < node5.second.size(); i++) // 看度数为3的节点的邻接节点
        {
            int index = findGletNode(glet, node5.second[i]);// see wether node5 has an adjacient node with degree 3
            if (glet[index].second.size() == 3)
            {
                return 20;//g513
            }
        }
        return 19;//g512
    }

    else if (signature == QVector<int>({2, 2, 2, 4, 4}))//g514
    {
        return 21;
    }

    else if (signature == QVector<int>({1, 3, 3, 3, 4}))//g515
    {
        return 22;
    }
    else if (signature == QVector<int>({2, 2, 3, 3, 4}))//g516
    {
        return 23;
    }
    else if (signature == QVector<int>({2, 3, 3, 3, 3}))//g517
    {
        return 24;
    }
    else if (signature == QVector<int>({2, 3, 3, 4, 4}))//g518
    {
        return 25;
    }
    else if (signature == QVector<int>({3, 3, 3, 3, 4}))//g519
    {
        return 26;
    }
    else if (signature == QVector<int>({3, 3, 4, 4, 4}))// g520
    {
        return 27;
    }
    else if (signature == QVector<int>({4, 4, 4, 4, 4}))// g521
    {
        return 28;
    }
    else
    {
        qDebug() << "undesired node number: " << glet.size();
        return -1;
    }
}

int findGletNode(const GraphLet& glet, int nodeId)
{
    for (int i = 0; i < glet.size(); i++)
    {
        if (glet[i].first == nodeId)
        {
            return i;
        }
    }
    return -1;
}
