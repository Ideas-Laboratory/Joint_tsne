#include "graph.h"
#include <QDebug>

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


void Graph::SortGraphLet(GraphLet& glet) {
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
    m_gfds = QVector<QVector<int>>(nodeNum(), QVector<int>(ALL_GRAPHLET, 0));

    // Get all connected components, each connected component has a vector of nodes
     QVector<QVector<int>> ccs = CCs();
     qDebug() << "Connected components: " << ccs.size();

     // For each connect component, apply random walk
     for (int i = 0; i < ccs.size(); i++)
     {
         int sCount = ccs[i].size()*1000; // * 500
         GUISE(sCount, ccs[i][0]);
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


bool Graph::connected()
{
    QVector<int> labels(nodeNum(), -1);

    QVector<int> cc;
    dfs(0, 0, labels, cc);

    for (int i = 0; i < nodeNum(); i++)
    {
       if (labels[i] < 0) {
          return false;
       }
    }

    return true;
}


Graph Graph::inducedSubgraph(const QVector<int> &nodeIds)
{
    Graph sub;
    for (int i = 0; i < nodeIds.size(); i++)
    {
        Node* node = GetNode(nodeIds[i]);
        sub.addNode(Node(node->x, node->y));
    }

    for (int i = 0; i < nodeIds.size(); i++)
    {
        for (int j = i + 1; j < nodeIds.size(); j++)
        {
            if (hasEdge(nodeIds[i], nodeIds[j]))
            {
                sub.addEdge(i, j);
            }
        }
    }
    return sub;
}

bool Graph::connected(const QVector<int> &nodeIds)
{
    Graph subgraph = inducedSubgraph(nodeIds);
    return subgraph.connected();
}


void Graph::setNodeNum(const int n)
{
    m_nodes.clear();
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

    num /= 2;
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
    m_nodes[b].addEdge(a);
}


bool Graph::hasEdge(int a, int b)
{
    if (a < 0 || a >= m_nodes.size() || b < 0 || b>= m_nodes.size())
    {
        return false;
    }
    return m_nodes[a].hasEdge(b) && m_nodes[b].hasEdge(a);
}


void Graph::clear()
{
    m_nodes.clear();
}


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

    // The size of feature vector is L * 26
    assert(featureVector.size() == ALL_GRAPHLET*neighborSize);
    return featureVector;
}


int Graph::graphType(Graph glet)
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


int Graph::gletTypeNum(int k)
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


void Graph::printGraphlet(const GraphLet &graphlet)
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

int Graph::gletType(GraphLet glet)
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

        for (int i = 0; i < node5.second.size(); i++) // the neighbor of node which has degree of 3
        {
            int index = findGletNode(glet, node5.second[i]);// see wether node5 has a neighbor with degree 3
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

int Graph::findGletNode(const GraphLet& glet, int nodeId)
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
