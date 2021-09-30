#include "graph_sim.h"

GraphSimilarity::GraphSimilarity(QString gfd_cal, int bfs_level):m_bfs_level(bfs_level)
{
    if (gfd_cal == "Concatenate") {
        gfd_cal = math_utils::GFDCalc::CONCAT;
    }
    else if (gfd_cal == "Accumulate") {
        gfd_cal = math_utils::GFDCalc::ACCUM;
    }
}

Graph GraphSimilarity::readGraph(const QString &fileName)
{
    qDebug() << "read graph: " << fileName;

    Graph g;
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qDebug() << "打开文件失败";
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

    if (m_rw_flag == 0)     // use random walk
    {
        g1.preGUISE();
        g2.preGUISE();
    }
    else
    {
        g1.preCount();
        g2.preCount();
    }


    #ifdef _DEBUG
        qDebug() << "counting graphlets costs " << time.elapsed() << "ms";
    #endif

    QVector<float> matchScores;
    for (int i = 0; i < g1.nodeNum(); i++)
    {
        QVector<float> vi_1;
        QVector<float> vi_2;

        if (m_rw_flag == 0)     // use random walk
        {
            vi_1 = g1.GetfeatureVectorRW(i, m_bfs_level);
            vi_2 = g2.GetfeatureVectorRW(i, m_bfs_level);
        }
        else if (m_gfd_cal == math_utils::GFDCalc::CONCAT)   // concatenate way
        {
            vi_1 = g1.GetfeatureVector(i, m_bfs_level);
            vi_2 = g2.GetfeatureVector(i, m_bfs_level);
        }

        float s = applyKernel(vi_1, vi_2, m_kernel);
//        qDebug() <<  i << "," << s;

        matchScores.push_back(s);
    }


    return matchScores;
}


void GraphSimilarity::savePointSims(const QVector<float> &matchScores, const QString &fileName)
{
    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        qDebug() << "打开文件失败";
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


//float GraphSimilarity::dist(Graph &g)
//{
//    g.preGUISE();
//    g.preCount();

//    for (int i = 0; i < g.nodeNum(); i++)
//    {
//        QVector<float> vi_1;
//        QVector<float> vi_2;

//        vi_1 = g.GetfeatureVectorRW(i, m_bfs_level);
//        vi_2 = g.GetfeatureVector(i, m_bfs_level);


//    }
//}
