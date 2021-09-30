#ifndef UTILS_H
#define UTILS_H
#include <QVector>
#include <cmath>
#include "qdebug.h"

namespace math_utils {
enum KernelFunc{
    COS,
    RBF,
    LAPLACIAN
} ;

enum GFDCalc{
    ACCUM,
    CONCAT
};

enum DistMeasure{
    L1,
    L2,
    CS
};

float applyKernel(const QVector<float>& vec1, const QVector<float>& vec2, KernelFunc func);

float cosine(const QVector<float>& vec1, const QVector<float>& vec2);
float rbf(const QVector<float>& vec1, const QVector<float>& vec2, float delta = 1.0);
float laplacian(const QVector<float>& vec1, const QVector<float>& vec2, float delta = 1.0);

float L1Distance(const QVector<float>& vec1, const QVector<float>& vec2);
float L2Distance(const QVector<float>& vec1, const QVector<float>& vec2);

void normalize(QVector<float>& vec);

//template <class T>
//void printVector(const QVector<T> &vec);

template<class T>
void printVector(const QVector<T> &vec)
{
    QString str;
    str += "[";
    for (int i = 0; i < vec.size(); i++)
    {
        str += QString::number(vec[i]);
        str += ", ";
    }
    str += "]";
    qDebug() << str;
}

// sort vec and return the sorted indices
QVector<int> sortIdx(QVector<int> &vec);
}





#endif // UTILS_H
