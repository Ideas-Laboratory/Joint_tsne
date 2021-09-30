#include "widget.h"
#include "ui_widget.h"
#include <QDebug>

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);
}

void Widget::paintEvent(QPaintEvent *event)
{
    QPainter* painter = new QPainter(this);
    // render background.
    painter->setPen(Qt::white);
    painter->setBrush(QColor(Qt::white));
    painter->drawRect(0,0,width(),height());
    painter->setFont(QFont("Calibri",10));
    painter->setRenderHint(QPainter::Antialiasing);

    // render graph
    m_graph.renderGraph(painter);

    // render graphlet
    m_graph.renderCurGraphLets(painter);

    // 如果正处于绘制边的模式
    if(m_startID != -1)
    {
        painter->drawLine(m_graph.GetNode(m_startID)->x,m_graph.GetNode(m_startID)->y,
                          m_lastPos.x(), m_lastPos.y());
    }

    delete painter;
}

void Widget::mousePressEvent(QMouseEvent *event)
{
    m_lastPos = event->pos();

    // 添加一个点
    if(event->modifiers() & Qt::ControlModifier && event->button() == Qt::LeftButton)
    {
        m_graph.addNode(Node(m_lastPos.x(),m_lastPos.y()));
        update();
        return;
    }
    // 如果选择一个现有的点，那么作为边的起始id
    if(event->button() == Qt::LeftButton)
    {
        m_startID = m_graph.HasNodeAt(m_lastPos);
    }
}

void Widget::mouseMoveEvent(QMouseEvent *event)
{
    if(m_startID != -1)
    {
        m_lastPos = event->pos();
        update();
    }
}

void Widget::mouseReleaseEvent(QMouseEvent *event)
{
    if(m_startID != -1 && event->button() == Qt::LeftButton)
    {
        int m_endID = m_graph.HasNodeAt(event->pos());
        if(m_endID != -1)
            m_graph.addEdge(m_startID,m_endID);
    }
    update();
    m_startID = -1;
}

Widget::~Widget()
{
    delete ui;
}
