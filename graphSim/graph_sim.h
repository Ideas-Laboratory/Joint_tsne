// This class is responsible for computing the topology feature of nodes and the similarity scores
#ifndef GRAPHSIMILARITY_H
#define GRAPHSIMILARITY_H
#include "math_utils.h"
#include "graph.h"
#include <QFile>

class GraphSimilarity
{
public:
    GraphSimilarity(QString, int);

    // compute the similarities between two graphs
    QVector<float> calcPointSims(Graph& g1, Graph& g2);


    // read in graph
    static Graph readGraph(const QString &fileName);
    // save similarities to file
    static void savePointSims(const QVector<float>& matchScores, const QString& fileName);

private:
    math_utils::KernelFunc m_kernel = math_utils::KernelFunc::COS;
    math_utils::GFDCalc m_gfd_cal = math_utils::GFDCalc::CONCAT;
    int m_bfs_level = 6;

    int m_rw_flag = 0;
};

#endif // GRAPHSIMILARITY_H
