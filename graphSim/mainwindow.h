#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMessageBox>
#include <QTextStream>
#include <qdebug.h>
#include "widget.h"
#include "math_utils.h"
#include "expandedmarkovchain.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void slot_openGraph();
    void slot_saveGraph();
//    void slot_calculateGraphLet();
    void slot_btnNext();

    void on_btnGen_clicked();

    Graph genGraph(int nodeNum = 500, int edgeNum = 499);

    void on_comboBox_activated(int index);

    void on_btnClear_clicked();

    void on_btnGFD_clicked();

    void on_btnSim_clicked();

    // read graph structure from fm_data generated from pytho
//    Graph read_fm_data(const QString& fileName);
    Graph read_graph(const QString& fileName);

    // save similarities to file
    void savePointSims(const QVector<float>& matchScores, const QString& fileName);

    // compute the similarities between two graphs
    QVector<float> calcPointSims(Graph& g1, Graph& g2);

    void on_btnLoadFmData_clicked();

    void on_cbxRange_activated(int index);

    void on_cbxKernel_activated(int index);

    void on_cbxGFDCalc_activated(int index);

//    void on_btnNRW_clicked();

//    void on_btnIRW_clicked();

    // read this kind of format...
    Graph readGraph(const QString& filePath);


    void on_cbxRwFlag_activated(int index);

    float meanVectorDist(Graph g, math_utils::DistMeasure m);
    QPair<float, float> meanPreTime(Graph g);
    QPair<float, float> meanTotalTime(Graph g, math_utils::DistMeasure m);

    void on_btnCalcMeanDist_clicked();


    void on_btnCalcMeanTime_clicked();

    void on_btnClearFileList_clicked();

private:
    int m_k = 5;
    int m_d = 2;
    EMkChain m_emc;

    Ui::MainWindow *ui;
    int m_graphlet_id = 0;

    QVector<Graph> m_k_glets;
    int m_k_glets_i = 0;

    math_utils::KernelFunc m_kernel = math_utils::KernelFunc::COS;
    math_utils::GFDCalc m_gfd_cal = math_utils::GFDCalc::CONCAT;
    int m_search_range = 6;
    int m_rw_flag = 0;

    QStringList m_filenames;
};
#endif // MAINWINDOW_H
