#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include <QMouseEvent>
#include "graph.h"

namespace Ui {
class Widget;
}

class Widget : public QWidget
{
    Q_OBJECT

public:
    explicit Widget(QWidget *parent = nullptr);
    ~Widget();

    Graph m_graph;

protected:
    void paintEvent(QPaintEvent *event);

    QPoint m_lastPos;
    int m_startID = -1;
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent *event);
private:
    Ui::Widget *ui;
};

#endif // WIDGET_H
