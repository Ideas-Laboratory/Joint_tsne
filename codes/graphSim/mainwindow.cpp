#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    for(int i=0; i<ALL_GRAPHLET; ++i)
    {
        ui->comboBox->addItem(QString("Graphlet:%1").arg(i+1));
    }
    for (int i = 0; i < MAX_SEARCH_RANGE; i++)
    {
        ui->cbxRange->addItem(QString("neighborhood:%1").arg(i+1));
    }

    ui->cbxRange->setCurrentIndex(m_search_range-1);
    ui->cbxGFDCalc->setCurrentIndex(m_gfd_cal);

    //    ui->comboBox->setCurrentIndex(19);

    connect(ui->btnOpen, SIGNAL(clicked(bool)), this, SLOT(slot_openGraph()));
    connect(ui->btnSave,SIGNAL(clicked(bool)),this,SLOT(slot_saveGraph()));

    //    connect(ui->btnCalculate,SIGNAL(clicked(bool)),this,SLOT(slot_calculateGraphLet()));
    connect(ui->btnNext,SIGNAL(clicked(bool)),this,SLOT(slot_btnNext()));

    m_filenames.clear();


    srand(time(0));
}

void MainWindow::slot_openGraph()
{
    QString filename = QFileDialog::getOpenFileName(this,"Open Graph","/Users/joe/Codes/QtProjects/t-sne for comparison/GraphSimilarity/data/", "Graph (*.graph)");
    if(filename == "")
    {
        return;
    }

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QMessageBox::warning(this,tr("错误"),tr("打开文件失败"));
        return;
    }

    QTextStream in(&file);

    Graph g;

    bool flag = false;
    QString line = in.readLine();
    while (!line.isNull())
    {
        qDebug() << line;

        if (line == "#")
        {
            qDebug() << g.nodeNum();
            flag = true;
            line = in.readLine();
            continue;
        }
        QStringList list = line.split(QRegExp("\\s+"));
        if (!flag)
        {
            // We assume that nodes are stored in file sequentially
            g.addNode(Node(list[1].toFloat(), list[2].toFloat()));
        }
        else
        {
            for (int i = 1; i < list.size()-1; i++)
            {
                g.addEdge(list[0].toInt(), list[i].toInt());
                qDebug() << QString("Add edge(%1,%2)").arg(list[0], list[i]);
            }
        }
        line = in.readLine();
    }

    ui->m_widget_0->m_graph = g;
    ui->m_widget_0->update();
}

void MainWindow::slot_saveGraph()
{
    QFileDialog fileDialog;
    QString fileName = fileDialog.getSaveFileName(this,tr("Save Graph"),".",tr("Graph (*.graph)"));
    if(fileName == "")
    {
        return;
    }

    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::warning(this,tr("错误"),tr("打开文件失败"));
        return;
    }
    else
    {
        QTextStream textStream(&file);

        QString str;
        // Store positions of nodes
        for (int i = 0; i < ui->m_widget_0->m_graph.nodeNum(); i++)
        {
            str += QString::number(i);
            str += " ";

            Node* s = ui->m_widget_0->m_graph.GetNode(i);
            str += QString("%1").arg(s->x);
            str += " ";

            str += QString("%1").arg(s->y);
            str += "\n";
        }

        // Delimeter
        str += "#\n";

        // edges of the graph
        for (int i = 0; i < ui->m_widget_0->m_graph.nodeNum(); i++)
        {
            str += QString::number(i);
            str += " ";

            Node* s = ui->m_widget_0->m_graph.GetNode(i);
            for (int j = 0; j < s->childs.size(); j++)
            {
                str += QString::number(s->childs[j]);
                str += " ";
            }
            str += "\n";
        }

        textStream<<str;
        QMessageBox::warning(this,tr("提示"),tr("保存文件成功"));
        file.close();
    }
}

//void MainWindow::slot_calculateGraphLet()  // 点击计算按钮
//{
//    ui->m_widget_0->m_graph.GetGraphletFromGraph(m_graphlet_id);
//}

void MainWindow::slot_btnNext()
{
    ui->m_widget_0->m_graph.ShowNextGraphlet();
    ui->m_widget_0->update();

    ui->graphletId->setText("graphlet:" + QString::number(ui->m_widget_0->m_graph.cur_ID));
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_btnGen_clicked()
{
    int nodeNum = ui->nodeNumEdit->text().toInt();
    int edgeNum = ui->edgeNumEdit->text().toInt();

    ui->m_widget_0->m_graph = genGraph(nodeNum, edgeNum);//(200, 1000) cost 0-100ms
    ui->m_widget_0->update();
}

Graph MainWindow::genGraph(int nodeNum, int edgeNum)
{
    Graph g;
    int maxWidth = ui->m_widget_0->width();
    int maxHeight = ui->m_widget_0->height();

    //    srand(time(0));
    // Random nodes
    for (int i = 0; i < nodeNum; i++)
    {
        int nx = rand() % (maxWidth+1);
        int ny = rand() % (maxHeight+1);
        g.addNode(Node(nx, ny));
    }

    // Random edges
    for (int i = 0; i < edgeNum; i++)
    {
        int sid = rand() % nodeNum;
        int tid = rand() % nodeNum;
        g.addEdge(sid, tid);
    }

    qDebug() << QString("(%1, %2)").arg(g.nodeNum()).arg(g.edgeNum(1));

    return g;
}

void MainWindow::on_comboBox_activated(int index)
{
    m_graphlet_id = index;
}

void MainWindow::on_btnClear_clicked()
{
    ui->m_widget_0->m_graph.clear();
    ui->m_widget_0->update();
}

void MainWindow::on_btnGFD_clicked()
{
    // Pick the last file
    Graph g = read_graph(m_filenames[m_filenames.size()-1]);
    qDebug() << "read file name: " << m_filenames[m_filenames.size()-1];
    //    Graph g = ui->m_widget_0->m_graph;
    //    Graph g = readGraph("/Users/joe/Codes/QtProjects/t-sne for comparison/data/Epinion/soc-Epinions1.txt");

    // compute gfd on the whole graph
    //    EMkChain emc(lcc);
    //    // Sample size
    //    int walkSteps = 20 * g.nodeNum();
    //    // compare the gfd between SRW and out methods
    //    QVector<float> v1 = emc.NB_SRW(walkSteps);//NB_SRW    SRW
    //    qDebug() << "SRW";

    //    //// 验证Guise计算gfd的正确性
    //    QVector<float> v1 = g.GuiseGFD();
    //    math_utils::printVector(v1);
    //    qDebug() << "------------------";

    //    QVector<float> v2 = g.countAllGraphlets();
    //    qDebug() << "numeration";
    //    math_utils::printVector(v2);
    //    qDebug() << "------------------";


    //    //// 验证Guise每一层
    //    int layerId = m_search_range-1;

    //    g.preGUISE();
    //    g.preCount();
    //    for (int i = 0; i < g.nodeNum(); i++)
    //    {
    //        QVector<float> v1 = g.GetfeatureVectorRWAtLayer(i, layerId);
    //        QVector<float> v2 = g.GetFeatureVectorAtLayer(i, layerId);

    ////        math_utils::printVector(v1);
    ////        qDebug() << "------------------";
    ////        math_utils::printVector(v2);
    //        qDebug() << "node " << i << " cos similarity"<< math_utils::cosine(v1, v2);
    //        qDebug() << "------------------";
    //    }
}

//Graph MainWindow::read_fm_data(const QString& fileName)
//{
//    qDebug() << "read fm data: " << fileName;

//    QFile file(fileName);
//    //    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
//    //    {
//    //        QMessageBox::warning(this,tr("错误"),tr("打开文件失败"));
//    //        return;
//    //    }
//    file.open(QIODevice::ReadOnly);

//    Graph g;
//    QTextStream ts(&file);
//    // process TRUNCK
//    while(!ts.atEnd())
//    {
//        QStringList line= ts.readLine().split(QRegExp("\\s+"));
//        line.removeAll("");

//        if( line.size() == 0)
//            continue;

//        if(line[0][0] == "G")  // 接下来要读取Graph Adjacency List;
//            break;


//        g.addNode(Node());
//    }

//    qDebug()<<"加载完毕, 顶点数"<<g.nodeNum();

//    // add edges
//    while(!ts.atEnd())
//    {
//        QStringList line= ts.readLine().split(QRegExp("\\s+"));
//        line.removeAll("");

//        if( line.size() == 0)
//            continue;

//        int id = line[0].toInt();
//        for(int i=1; i<line.size(); i+=2)
//        {
//            g.addEdge(id, line[i].toInt());
//        }
//    }

//    return g;
//}

void MainWindow::savePointSims(const QVector<float> &matchScores, const QString &fileName)
{
    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::warning(this,tr("错误"),tr("打开文件失败"));
        return;
    }
    else
    {
        QTextStream textStream(&file);
        for (int i = 0; i < matchScores.size(); i++)
        {
            textStream << QString::number(i);
            textStream << "\t";
            textStream << QString::number(i);
            textStream << "\t";
            textStream << matchScores[i];
            textStream << "\n";
        }
        file.close();
    }

    qDebug() << "save file: " << fileName << " successful.";
}

void MainWindow::on_btnSim_clicked()
{

    QString fileName0, fileName1;
    Graph g0, g1;
    for (int i = 0; i < m_filenames.size()-1; i++)
    {
        fileName0 = m_filenames[i];
        fileName1 = m_filenames[i+1];

        // Pick the last two files
        g0 = read_graph(fileName0);
        g1 = read_graph(fileName1);


        QVector<float> pointSims = calcPointSims(g0, g1);
        // we leave it to python to compute edge similarity
        //    MatchEdgeList matchEdges = calcEdgeSims(g0, g1,pointSims);

        QFileInfo finfo0(fileName0);
        QFileInfo finfo1(fileName1);

        QStringList dirList0 = finfo0.dir().path().split("/");
        QStringList dirList1 = finfo1.dir().path().split("/");

        // /Users/joe/Codes/QtProjects/t-sne for comparison/data/data1/highdims/dim3/size100/2nn
        QString knn = dirList0[dirList0.size()-1];
        QString sizem = dirList0[dirList0.size()-2];
        QString dimn_0 = dirList0[dirList0.size()-3];
        QString data_type = dirList0[dirList0.size()-5];

        QString dimn_1 = dirList1[dirList0.size()-3];

        int d0 = (finfo0.baseName().split("_")[1]).toInt();
        int d1 = (finfo1.baseName().split("_")[1]).toInt();

        // we also consider searching range into similarities
        QString dir = QString("/Users/joe/Codes/QtProjects/t-sne for comparison/data/%1/qt_sim_1000/%2_%3/%4/%5/level%6/%7").arg(data_type).arg(dimn_0).arg(dimn_1).arg(sizem).arg(knn).arg(m_search_range).arg(ui->cbxGFDCalc->itemText(m_gfd_cal));
        QDir fileDir(dir);
        if (!fileDir.exists())
        {
            fileDir.mkpath(dir);
        }

        qDebug() << fileDir.path();
        savePointSims(pointSims,
                      (fileDir.path() + "/similar_points_%1_%2.txt").arg(d0).arg(d1));
    }

}

Graph MainWindow::read_graph(const QString &fileName)
{
    qDebug() << "read graph: " << fileName;

    Graph g;
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QMessageBox::warning(this,tr("错误"),tr("打开文件失败"));
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

//QVector<float> MainWindow::calcPointSims(Graph &g1, Graph &g2)
//{
//    QTime time;
//    time.start();

//    g1.preCalcGlets();
//    g2.preCalcGlets();

//#ifdef _DEBUG
//    qDebug() << "Computing graphlets on each node costs " << time.elapsed() << "ms";
//#endif



//    QVector<float> matchScores;
//    for (int i = 0; i < g1.nodeNum(); i++)
//    {
//        time.start();

//        QVector<float> vi_1;
//        QVector<float> vi_2;

//        if (m_gfd_cal == math_utils::GFDCalc::CONCAT)
//        {
//            vi_1 = g1.GetfeatureVector(i);
//            vi_2 = g2.GetfeatureVector(i);
//        }
//        else
//        {
//            vi_1 = g1.GetfeatureVectorAll(i);
//            vi_2 = g2.GetfeatureVectorAll(i);
//        }

//        float s = applyKernel(vi_1, vi_2, m_kernel);
//        qDebug() <<  i << "," << s;

//        matchScores.push_back(s);


//#ifdef _DEBUG
//    qDebug() << "Computing similarity on the node costs " << time.elapsed() << "ms";
//#endif
//    }


//    return matchScores;
//}

QVector<float> MainWindow::calcPointSims(Graph &g1, Graph &g2)
{
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
    qDebug() << "Count graphlets on each node costs " << time.elapsed() << "ms";
#endif



    QVector<float> matchScores;
    for (int i = 0; i < g1.nodeNum(); i++)
    {
        time.start();

        QVector<float> vi_1;
        QVector<float> vi_2;

        if (m_rw_flag == 0)     // use random walk
        {
            vi_1 = g1.GetfeatureVectorRW(i, m_search_range);
            vi_2 = g2.GetfeatureVectorRW(i, m_search_range);
        }
        else if (m_gfd_cal == math_utils::GFDCalc::CONCAT)   // concatenate way
        {
            vi_1 = g1.GetfeatureVector(i, m_search_range);
            vi_2 = g2.GetfeatureVector(i, m_search_range);
        }
        //        else                                            // accumulate way
        //        {
        //            vi_1 = g1.GetfeatureVectorAll(i, m_search_range);
        //            vi_2 = g2.GetfeatureVectorAll(i, m_search_range);
        //        }

        float s = applyKernel(vi_1, vi_2, m_kernel);
        qDebug() <<  i << "," << s;

        matchScores.push_back(s);

        //#ifdef _DEBUG
        //    qDebug() << "Computing similarity on the node costs " << time.elapsed() << "ms";
        //#endif
    }


    return matchScores;
}



void MainWindow::on_btnLoadFmData_clicked()
{
    QStringList filenames = QFileDialog::getOpenFileNames(this,"Open input sequences","/Users/joe/Codes/QtProjects/t-sne for comparison/data/highdims/", "text file(*.txt)");

    for (int i = 0; i < filenames.size(); i++)
    {
        m_filenames.push_back(filenames[i]);
        ui->fileListBox->addItem(filenames[i]);
    }
}

void MainWindow::on_cbxRange_activated(int index)
{
    m_search_range = index + 1;
    qDebug() << "set searching range: " << m_search_range;
}

void MainWindow::on_cbxKernel_activated(int index)
{
    m_kernel = math_utils::KernelFunc(index);
    qDebug() << "set kernel function: " << m_kernel << " " << ui->cbxKernel->itemText(index);
}

void MainWindow::on_cbxGFDCalc_activated(int index)
{
    m_gfd_cal = math_utils::GFDCalc(index);
}

//void MainWindow::on_btnNRW_clicked()
//{
//    // test on graph g0 to see whether graphlets gained from random walk is correct
//    EMkStat xl = m_emc.getCurrentExpandState();
//#ifdef _DEBUG
//    xl.printStatus();
//#endif

//    // subgraph induced by nodes in Xl
//    QVector<int> V = xl.nodeSet();
//    if (V.size() == m_k)//collect k distinct nodes
//    {
//        Graph s = m_emc.m_graph.inducedSubgraph(V);
//        ui->m_widget_0->m_graph.m_graphlet = new Graph(s);// might bug here

//        int i = graphType(s);
//        qDebug() << "graphlet type: " << i;
//    }

//    // draw graphlets on widget
//    ui->m_widget_0->update();

//    // uniformly choose neighbor and update
//    xl = m_emc.NBUpdate(m_k, m_d);
//}

//void MainWindow::on_btnIRW_clicked()
//{
//    m_emc = EMkChain(ui->m_widget_0->m_graph);
//    EMkStat xl = m_emc.NBInit(m_k, m_d);

//#ifdef _DEBUG
//    qDebug() << QString("k:%1, d:%2").arg(m_k).arg(m_d);
//#endif
//}


Graph MainWindow::readGraph(const QString &filePath)
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QMessageBox::warning(this,tr("错误"),tr("打开文件失败"));
        return Graph();
    }

    QTextStream in(&file);



    QString line = in.readLine();
    // first line tells us how many node and edges
    QStringList list = line.split(QRegExp("\\s+"));
    int node_num = list[0].toInt();
    int edge_num = list[1].toInt();

    qDebug() << "nodeNum: " << node_num << ", edgeNum: " << edge_num;

    QVector<int> graph_nodes;
    QVector<QPair<int,int>> graph_edges;

    line = in.readLine();
    while (!line.isNull())
    {
        QStringList list = line.split(QRegExp("\\s+"));
        int u = list[0].toInt();
        int v = list[1].toInt();

        if (!graph_nodes.contains(u))
        {
            graph_nodes.push_back(u);
        }
        if (!graph_nodes.contains(v))
        {
            graph_nodes.push_back(v);
        }
        graph_edges.push_back(qMakePair(u, v));

        line = in.readLine();
    }

    // 将原来的序号映射到(0-nodeNum)
    Graph g(node_num);
    for (int i = 0; i < graph_edges.size(); i++)
    {
        QPair<int, int> edge = graph_edges[i];

        int u = graph_nodes.indexOf(edge.first);
        int v = graph_nodes.indexOf(edge.second);
        g.addEdge(u, v);
    }

    return g;
}

void MainWindow::on_cbxRwFlag_activated(int index)
{
    m_rw_flag = index;
    if (m_rw_flag == 0)
    {
        qDebug() << "Use Random Walk.";
    }
    else
    {
        qDebug() << "Use Enumeration.";
    }
}


float MainWindow::meanVectorDist(Graph g, math_utils::DistMeasure m)
{
    float dist = 0;

    int K = 1;//10
    for (int k = 0; k < K; k++)
    {
        g.preCount();
        g.preGUISE();

        float sum = 0;
        for (int i = 0; i < g.nodeNum(); i++)
        {
            QVector<float> v1 = g.GetfeatureVectorRW(i, m_search_range);
            QVector<float> v2 = g.GetfeatureVector(i, m_search_range);

            if (m == math_utils::L1)
            {
                sum += math_utils::L1Distance(v1, v2);
            }
            else if (m == math_utils::L2)
            {
                sum += math_utils::L2Distance(v1, v2);
            }
            else
            {
                sum += math_utils::cosine(v1, v2);
            }
        }

        sum /= g.nodeNum();
        dist += sum;
    }


    return dist/K;
}

QPair<float, float> MainWindow::meanPreTime(Graph g)
{
    int K = 1;//10
    QTime timer;

    qDebug() << "Test time on enumeration...";

    timer.start();
    for (int k = 0; k < K; k++)
    {
        g.preCount();
    }
    float t1 = timer.elapsed()/K;


    qDebug() << "Test time on Guise...";

    timer.start();
    for (int k = 0; k < K; k++)
    {
        g.preGUISE();
    }
    float t2 = timer.elapsed()/K;


    return qMakePair(t1, t2);
}

QPair<float, float> MainWindow::meanTotalTime(Graph g, math_utils::DistMeasure m)
{
    int K = 1;//10
    QTime timer;


    qDebug() << "Test time on enumeration...";
    timer.start();
    for (int k = 0; k < K; k++)
    {
        g.preCount();
        for (int i = 0; i < g.nodeNum(); i++)
        {
            QVector<float> vi = g.GetfeatureVector(i, m_search_range);
        }
    }
    float t1 = timer.elapsed()/K;



    qDebug() << "Test time on Guise...";
    timer.start();
    for (int k = 0; k < K; k++)
    {
        g.preGUISE();
        for (int i = 0; i < g.nodeNum(); i++)
        {
            QVector<float> vi = g.GetfeatureVectorRW(i, m_search_range);
        }
    }
    float t2 = timer.elapsed()/K;

    return qMakePair(t1, t2);
}


void MainWindow::on_btnCalcMeanDist_clicked()
{
    // 当前选择的measure
    math_utils::DistMeasure m = math_utils::DistMeasure(ui->cbxDistMeasure->currentIndex());

    //    // 当前选择的graph
    //    QString filename = m_filenames[ui->fileListBox->currentIndex()];
    //    Graph g = read_fm_data(filename);
    Graph g = ui->m_widget_0->m_graph;

    qDebug() << QString("(%1, %2)").arg(g.nodeNum()).arg(g.edgeNum(1));

    qDebug() << "Test on L1-distance between enumeration and Guise...";
    float d = meanVectorDist(g, m);
    qDebug() << "mean distance: " << d;
}


void MainWindow::on_btnCalcMeanTime_clicked()
{
    // 当前选择的graph
    //    QString filename = m_filenames[ui->fileListBox->currentIndex()];
    //    Graph g = read_fm_data(filename);
    Graph g = ui->m_widget_0->m_graph;

    qDebug() << QString("(%1, %2)").arg(g.nodeNum()).arg(g.edgeNum(1));

    qDebug() << "Test time spent on enumeration and Guise...";

    // 当前选择的measure
    math_utils::DistMeasure m = math_utils::DistMeasure(ui->cbxDistMeasure->currentIndex());
    QPair<float, float> pair = meanTotalTime(g, m);
    qDebug() << "Enumeration vs. GUISE (ms): " << pair.first << ", " << pair.second;
}

void MainWindow::on_btnClearFileList_clicked()
{
    ui->fileListBox->clear();
    m_filenames.clear();
}
