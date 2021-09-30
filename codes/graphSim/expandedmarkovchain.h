#ifndef MARKOVCHAIN_H
#define MARKOVCHAIN_H
#include "graph.h"
#include <QSet>

typedef QVector<int> MkStat;            // A Markov state stores all nodeIds of corresponding node on G(d)
typedef QPair<MkStat, MkStat> NBMkStat; // expand Markov state to Non backtracing
typedef struct ExpandedMarkovState{
    ExpandedMarkovState()
    {

    }

    ExpandedMarkovState(const QVector<MkStat>& states):states(states)
    {
    }

    void clear()
    {
        states.clear();
    }

    void printStatus() const
    {
        QString str;
        for (int i = 0; i < states.size(); i++)
        {
            str += "(";
            for (int j= 0; j < states[i].size(); j++)
            {
                str += QString::number(states[i][j]);
                str += ", ";
            }
            str += ") => ";
        }

        qDebug() << str;
        qDebug() << "--------------------------------";
    }

    int size() const
    {
        return states.size();
    }
    MkStat& operator[](const int i)
    {
        return states[i];
    }
    MkStat operator[](const int i) const
    {
        return states[i];
    }

    void pop_front()
    {
        states.pop_front();
    }

    void push_back(const MkStat& state)
    {
        states.push_back(state);
    }


    QVector<int> nodeSet() const
    {
        QVector<int> set;
        for (int i = 0; i < states.size(); i++)
        {
            for (int j = 0; j < states[i].size(); j++)
            {
                int n = states[i][j];
                if(!set.contains(n))
                {
                    set.append(n);
                }
            }
        }
        return set;
    }

//    QSet<int> nodeSet() const
//    {
//        QSet<int> set;
//        for (int i = 0; i < states.size(); i++)
//        {
//            for (int j = 0; j < states[i].size(); j++)
//            {
//                int n = states[i][j];
//                set.insert(n);
//            }
//        }
//        return set;
//    }

    static QVector<int> commonNodes(const ExpandedMarkovState& xa, const ExpandedMarkovState& xb)
    {
        QVector<int> nodeSetA = xa.nodeSet();
        QVector<int> nodeSetB = xb.nodeSet();

        std::sort(nodeSetA.begin(), nodeSetA.end());
        std::sort(nodeSetB.begin(), nodeSetB.end());

        QVector<int> nodeSetC;
        std::set_intersection(nodeSetA.begin(), nodeSetA.end(), nodeSetB.begin(), nodeSetB.end(), std::back_inserter(nodeSetC));

        return nodeSetC;
    }

    QVector<MkStat> states;
} EMkStat;

class EMkChain
{
public:
    EMkChain();
    EMkChain(const Graph& originalGraph);

    EMkStat update(int k, int d);               // update state in Expanded Markov Chain
    EMkStat init(int d, int sid = 0);                        // start from initial node
    int coff(int k, int i, int d);
    float stableDistr(const EMkStat& xl, int d);
    QVector<float> SRW(int n, int d, int k, bool norm);
    QVector<float> SRW(int n);


    EMkStat NBInit(int k ,int d, int sid = 0);
    EMkStat NBUpdate(int k, int d);
    float NB_stableDistr(const EMkStat& xl, int d);
    QVector<float> NB_SRW(int n, int d, int k, bool norm);
    QVector<float> NB_SRW(int n);

    QSet<GraphLet> NB_SRW_GLETS(int n, int d, int k);
    QSet<GraphLet> NB_SRW_GLETS(int n);


    // CSS seems not work well...
    int isoIndicator(const EMkStat& xl, int k, int i);
    int sizeIndicator(const EMkStat& xl, int k);
    float sampleProb(const EMkStat& xl, int k, int d);
    QVector<float> SRW_CSS(int n, int d, int k);
    QVector<float> SRW_CSS(int n);
    //

    MkStat getCurrentState();
    EMkStat getCurrentExpandState();


    Graph m_graph;          // the underlying graph to exploit
private:
    int dgree(const MkStat& xi);
    int NB_dgree(const MkStat& xi);

//    bool statEqual(const MkStat& xi, const MkStat& xj);

    QVector<MkStat> popNeighbors(MkStat xi, int d);   // populate neighborhood from


    int R(int d);                               // edge number in d-node subgraph relationship graph
    int R1;
    int R2;

    EMkStat m_estate;        // current expanded status

};


#endif // MARKOVCHAIN_H
