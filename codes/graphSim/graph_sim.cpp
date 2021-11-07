#include "graph_sim.h"

Graph GraphSimilarity::readGraph(const QString &fileName)
{
    qDebug() << "read graph: " << fileName;

    Graph g;
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qDebug() << "Fails on reading file.";
        return g;
    }
    QTextStream ts(&file);

    int nodeNum = -1;
    // add edges
    while(!ts.atEnd())
    {
        QStringList line= ts.readLine().split(QRegExp("\\s+"));
        line.removeAll("");

        if( line.size() == 0)
            continue;

        if (nodeNum < 0)
        {
            nodeNum = line[0].toInt();
            g.setNodeNum(nodeNum);
        }
        else {
            int id = line[0].toInt();
            for(int i=1; i<line.size(); i+=2)
            {
                g.addEdge(id, line[i].toInt());
            }
        }

    }

    return g;
}

QVector<float> GraphSimilarity::calcPointSims(Graph &g1, Graph &g2)
{
    qDebug() << "(" << g1.nodeNum() << ", " << g1.edgeNum(1) << ")";
    qDebug() << "(" << g2.nodeNum() << ", " << g2.edgeNum(1) << ")";

    QTime time;
    time.start();

    g1.preGUISE();
    g2.preGUISE();

    qDebug() << "counting graphlets costs " << time.elapsed() << "ms";

    QVector<float> matchScores;
    for (int i = 0; i < g1.nodeNum(); i++)
    {
        QVector<float> vi_1;
        QVector<float> vi_2;

        vi_1 = g1.GetfeatureVectorRW(i, m_bfs_level);
        vi_2 = g2.GetfeatureVectorRW(i, m_bfs_level);
        float s = applyKernel(vi_1, vi_2, m_kernel);

        matchScores.push_back(s);
    }


    return matchScores;
}


void GraphSimilarity::savePointSims(const QVector<float> &matchScores, const QString &fileName)
{
    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        qDebug() << "Fails on reading file.";
        return;
    }
    else
    {
        QTextStream textStream(&file);
        for (int i = 0; i < matchScores.size(); i++)
        {
            textStream << QString::number(i);
            textStream << "\t";
            textStream << matchScores[i];
            textStream << "\n";
        }
        file.close();
    }

    qDebug() << "save file: " << fileName << " successful.";
}
