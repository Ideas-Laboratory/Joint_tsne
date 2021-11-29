#ifndef GRAPH_H
#define GRAPH_H
#define ALL_GRAPHLET    29
#include <QVector>
#include <cassert>
#include <QQueue>
#include <QTime>
#include <QPoint>

#include "math_utils.h"

struct Node    // 图
{
    float x,y;
    QVector<int> childs;

    // used to store which layer the node is at from certain center
    int layerId;

    Node():x(0),y(0){}
    Node(float x, float y): x(x),y(y) {}

    int degree() const
    {
        return childs.size();
    }

    void addEdge(int id)
    {
        if(!childs.contains(id))
            childs.push_back(id);
    }

    bool hasEdge(int id) const
    {
        return childs.contains(id);
    }
};

typedef QPair<int,QVector<int>> GraphLetNode;   // node and its adjacient edges, which is an undirected graph
typedef QVector<GraphLetNode> GraphLet;         // graphlet representation
typedef QVector<int> GuiseGraphlet;             // the nodes of graphlet

class Graph
{
public:
    Graph() {

    }
    Graph(int n):m_nodes(QVector<Node>(n)) {

    }

    void setNodeNum(const int);
    int nodeNum() const;

    int edgeNum(int d);

    void addNode(const Node& n);
    void addEdge(int a, int b);

    bool hasEdge(int a, int b);

    Node *GetNode(int id);
    QVector<QPair<int, int>> GetEdges();

    void clear();

    QVector<QVector<int>> CCs();

    void preGUISE();
    QVector<float> GetfeatureVectorRW(int sid, int neighborSize);

private:
    void dfs(int nodeId, int label, QVector<int>& labels, QVector<int>& cc);
    bool connected();

    Graph inducedSubgraph(const QVector<int>& nodeIds);
    bool connected(const QVector<int>& nodeIds);

    void SortGraphLet(GraphLet& glet);
    GraphLet inducedGraphlet(const QVector<int>& nodeIds);

    // apply guise algorithm on sid with sCount as # of samples
    void GUISE(int sCount, int sid = 0);
    GuiseGraphlet GuiseInit(int gletSize, int sid = 0);
    QVector<GuiseGraphlet> addVertex(const GuiseGraphlet& gt);
    QVector<GuiseGraphlet> popNeighborGuise(const GuiseGraphlet& gt);
    QVector<int> BfsLayersID(int sid, int neighborSize);

    // utility functions
    int graphType(Graph glet);
    int gletType(GraphLet glet);
    int findGletNode(const GraphLet& glet, int nodeId);
    void printGraphlet(const GraphLet& graphlet);


    // members
    QVector<Node> m_nodes;
    QVector<QVector<int>> m_gfds;
};

#endif // GRAPH_H
