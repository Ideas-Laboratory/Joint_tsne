#ifndef GRAPH_H
#define GRAPH_H
#define MAX_SEARCH_RANGE 30
#define ALL_GRAPHLET    29

#include <qdebug.h>
#include <QVector>
#include <QPainter>
#include <cassert>
#include <QQueue>
#include <QTime>
#include "math_utils.h"
#define _DEBUG

struct Node    // 图
{
    float x,y;
    QVector<int> childs;

    // Updated by Yinqiao 3.24
    int layerId;

    Node():x(0),y(0){}
    Node(float x, float y): x(x),y(y) {}

    int degree() const
    {
        return childs.size();
    }

    void addEdge(int id)  //
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
typedef QVector<float> GFD;                     // graphlet frequency distribution


void printGraphlet(const GraphLet& graphlet);

// Undirected graph or directed graph
#define UNDIRECT_GRAPH

typedef QVector<int> GuiseGraphlet;

// Forward declaration
class Graph
{
public:
    Graph();
    Graph(int n);


    int cur_ID = 0;                                      // 当前显示的graphlet编号
    QVector<GraphLet> m_cur_graphlets;                      // 当前找到的


    Graph* m_graphlet = NULL;                            // current graphlet that RW is visiting
    QVector<Node> m_nodes;                              // 结点，没有id信息


    bool isInduced(const GraphLet& glet);          // 判断graphlet是不是Graph的诱导子图


    void findAllComb(QVector<QVector<bool>>& flags, QVector<bool>& flag, int i);
    Graph inducedSubgraph(const QVector<int>& nodeIds);            // subgraph induced by nodeIds. This function will lose original id information
    GraphLet inducedGraphlet(const QVector<int>& nodeIds);

    QVector<Graph> connectedSubgraphs(QVector<int> nodeIds);// spanning subgraphs of induced subgraphs
    QVector<QVector<int>> connectedInducedGraph(int d); // all d-node connected subgraphs

    void dfs(int nodeId, int label, QVector<int>& labels, QVector<int>& cc);
//    void dfs(int nodeId, int label, QVector<int>& labels);

    QVector<QVector<int>> CCs();
//    QVector<int> CCs();     // a node for each connected component
//    QVector<int> LCC();                                 // Largest connected components


//    void dfs(int nodeId, int label, QVector<int> &labels);
    bool connected();
    bool connected(const QVector<int>& nodeIds);

    void setNodeNum(const int);
    int nodeNum() const;


    int edgeNum(int d);

    void addNode(const Node& n);
    void addEdge(int a, int b);

    bool hasEdge(int a, int b);

    void renderGraph(QPainter* painter);
    void renderCurGraphLets(QPainter* painter);

    Node *GetNode(int id);
    int HasNodeAt(QPoint p, float radius = 100);

    QVector<QPair<int, int>> GetEdges();

    // get all graphlets from current graph. 是为了演示效果
    void ShowNextGraphlet();
    void clear();


    // Search graphlet type gid starting from node sid
    QSet<GraphLet> SearchGraphLet(int gid, int sid);

    //  return all graphlets found on given node
    QSet<GraphLet> GetGraphLets(int sid);

    ///1. My numeration based methods
//    QSet<GraphLet> m_graphletSet;
    QVector<QVector<int> > m_gfds_enumer;
    void preCount();

//    QVector<float> GetFeatureVectorAtLayer(int sid, int layerId);

     // concat gfd of different levels
//    QVector<QSet<GraphLet>> GetNeighborGraphlets(int sid, int neighborSize);
    QVector<float> GetfeatureVector(int sid, int neighborSize);

    // simply count all graphlets in the neighborhood, which is not used
//    QVector<QSet<GraphLet>> GetNeighborGraphletsAll(int sid, int neighborSize);
//    QVector<float> GetfeatureVectorAll(int sid, int neighborSize);




    ///2. GUISE
    GuiseGraphlet GuiseInit(int gletSize, int sid = 0);
    QVector<GuiseGraphlet> addVertex(const GuiseGraphlet& gt);
    QVector<GuiseGraphlet> popNeighborGuise(const GuiseGraphlet& gt);
    void GUISE(int sCount, int sid = 0);    // apply guise algorithm on sid with sCount as # of samples
    void preGUISE();
//    QVector<float> GuiseGFD();      // graphlet frequency distribution


//    QSet<GraphLet> m_graphletSetGuise;
    QVector<QVector<int>> m_gfds;
//    QVector<GraphLet> m_graphletSet;
    ///3. Expanded Markov chain, which is not suitable for our algorithm
//    void preRW();


    QVector<int> BfsLayersID(int sid, int neighborSize);
    QVector<float> GetfeatureVectorRW(int sid, int neighborSize);
//    QVector<float> GetfeatureVectorRWAtLayer(int sid, int layerId);

private:
    // A family of functions for searching graphlets
    // Search graphlet from sid
    QSet<GraphLet> SearchGraphLet1(int sid);
    QSet<GraphLet> SearchGraphLet2(int sid);
    QSet<GraphLet> SearchGraphLet3(int sid);
    QSet<GraphLet> SearchGraphLet4(int sid);
    QSet<GraphLet> SearchGraphLet5(int sid);
    QSet<GraphLet> SearchGraphLet6(int sid);
    QSet<GraphLet> SearchGraphLet7(int sid);
    QSet<GraphLet> SearchGraphLet8(int sid);
    QSet<GraphLet> SearchGraphLet9(int sid);
    QSet<GraphLet> SearchGraphLet10(int sid);
    QSet<GraphLet> SearchGraphLet11(int sid);
    QSet<GraphLet> SearchGraphLet12(int sid);
    QSet<GraphLet> SearchGraphLet13(int sid);
    QSet<GraphLet> SearchGraphLet14(int sid);
    QSet<GraphLet> SearchGraphLet15(int sid);
    QSet<GraphLet> SearchGraphLet16(int sid);
    QSet<GraphLet> SearchGraphLet17(int sid);
    QSet<GraphLet> SearchGraphLet18(int sid);
    QSet<GraphLet> SearchGraphLet19(int sid);
    QSet<GraphLet> SearchGraphLet20(int sid);
    QSet<GraphLet> SearchGraphLet21(int sid);
    QSet<GraphLet> SearchGraphLet22(int sid);
    QSet<GraphLet> SearchGraphLet23(int sid);
    QSet<GraphLet> SearchGraphLet24(int sid);
    QSet<GraphLet> SearchGraphLet25(int sid);
    QSet<GraphLet> SearchGraphLet26(int sid);
    QSet<GraphLet> SearchGraphLet27(int sid);
    QSet<GraphLet> SearchGraphLet28(int sid);
    QSet<GraphLet> SearchGraphLet29(int sid);

    // de-duplicate graphlets of certain type. require that each graphlet is already sorted
    int DedupGraphLets(QVector<GraphLet>& sortedGraphlets);
    // sort graphlet node in a graphlet and symmetrization
    void SortGraphLet(GraphLet& glet);
    // determine whether two graphlets are the same.
    // Before calling this function, glet1 and glet2 must be sorted.
    bool GraphletEqual(const GraphLet& glet1, const GraphLet& glet2);

};

QVector<int> sortIdx(QVector<int>& vec);

int graphType(Graph glet);
int gletType(GraphLet glet);
int gletTypeNum(int k);

int findGletNode(const GraphLet& glet, int nodeId);


template <typename T>
void combine_inner(T &data, int start, int n, int m, int depth, T temp, QVector<T> &result)
{
    if (depth == m - 1)
    {
        //最内层循环 将temp加入result
        for (int i = start; i < n - (m - depth - 1); ++i)
        {
            temp[depth] = data[i];
            result.push_back(temp);
        }
    }
    else
        for (int i = start; i < n - (m - depth - 1);++i)
    {
        temp[depth] = data[i];//每层输出一个元素
        combine_inner(data,i + 1, n, m, depth+1,temp,result);
    }
}

//T可以调入vector<int>, string等，需要支持下标[]操作及size()函数
template <typename T>
QVector<T> combine(T &data,int m)
{
    if (m <= 0)
        return{};
    int depth = 0;
    QVector<T> result;
    T temp(m);//,0
    combine_inner(data,0, data.size(), m, depth,temp,result);
    return result;
}

#endif // GRAPH_H
