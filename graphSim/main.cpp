#include "mainwindow.h"
#include <QApplication>
#include <QJsonValue>
#include <QJsonArray>
#include <QJsonObject>
#include <QJsonDocument>
#include <QJsonParseError>
#include <QIODevice>
#include <QDir>
//#include <assert.h>
#include "graph_sim.h"


int main(int argc, char *argv[])
{
//    QApplication app(argc, argv);
    //    MainWindow w;
    //    w.show();

    //    return app.exec();

    //    qDebug() << "current application path:" << QCoreApplication::applicationDirPath();
    //    qDebug() << "current path:" << QDir::currentPath();

    // parse from argv

//    assert(argc == 2);
    QString config_path = argv[1];
//    QString config_path = "/Users/joe/Codes/PythonProjects/joint_tsne_experiments/config/5_cluster_trans_split_overlap_contract.json";

//    QFileInfo info(__FILE__);
//    QString sourceDir = info.absolutePath();
//    qDebug() << sourceDir;
//    qDebug() << QDir::currentPath();
//    qDebug() << QCoreApplication::applicationDirPath();

//    QString rootDir = "../../../../../../";
    QString rootDir = "/Users/joe/Codes/PythonProjects/joint_tsne_experiments/";

    QFile configFile(config_path);
    configFile.open(QIODevice::ReadOnly | QIODevice::Text);
    QString value = configFile.readAll();
    configFile.close();

    QJsonParseError parseJsonErr;
    QJsonDocument document = QJsonDocument::fromJson(value.toUtf8(), &parseJsonErr);
    if (! (parseJsonErr.error == QJsonParseError::NoError))
    {
        qDebug() << "配置文件错误！";
        return -1;
    }

    QJsonObject jsonObject = document.object();
    qDebug() << jsonObject;

    QJsonObject algoObject = jsonObject["algo"].toObject();
    QJsonObject dataObject = jsonObject["thesne"].toObject();

    // Only access algorithm parameters
    //    qDebug() << jsonObject["k_closest_count"].toInt();
    int bfs_level = algoObject["bfs_level"].toInt();

    qDebug() << "bfs_level: " << bfs_level;
    QString gfd_calc = algoObject["gfd_calc"].toString();
    QString data_name = dataObject["data_name"].toString();

    QJsonArray data_ids = dataObject["data_ids"].toArray();



    GraphSimilarity gs(gfd_calc, bfs_level);


    QString graphDirStr = rootDir + "knn graph/" + data_name;
    QDir graphDir(graphDirStr);


    QString simDirStr = rootDir + "graphSim/" + data_name;
    QDir simDir(simDirStr);
    if (!simDir.exists())
    {
        simDir.mkpath(simDirStr);
    }
    else {
        // clear the directory
        QStringList items = simDir.entryList();
        for (int i = 0; i < items.size(); i++) {
            QString item = items[i];
            simDir.remove(item);
        }
    }
    qDebug() << simDirStr;

    // sort by name
    graphDir.setSorting(QDir::Name);
    // filter by suffix
//    QStringList nameFilter("*.txt");
//    QFileInfoList list = graphDir.entryInfoList(nameFilter, QDir::Files);

    QStringList list;
    for (int i = 0; i < data_ids.size(); i++) {
        QString data_id = "g_" + QString::number(data_ids[i].toInt());
        list.append(data_id);
    }



    QTime time;
    time.start();
    for (int i = 0; i < list.size()-1; ++i) {
//        QFileInfo fileInfo0 = list.at(i);
        QString fileName0 = graphDir.absoluteFilePath(list[i]+".txt");
        QString fileName1 = graphDir.absoluteFilePath(list[i+1]+".txt");

        // read filenames
        Graph g0 = GraphSimilarity::readGraph(fileName0);
        Graph g1 = GraphSimilarity::readGraph(fileName1);

        QVector<float> sims = gs.calcPointSims(g0, g1);

        QString simName = list[i] + "_" + list[i+1] + ".txt";
        QString simPath = simDir.absoluteFilePath(simName);

        GraphSimilarity::savePointSims(sims, simPath);
    }
    qDebug() << "Computing similairy costs " << float(time.elapsed())/1000.f << "s";

    return 0;
}
