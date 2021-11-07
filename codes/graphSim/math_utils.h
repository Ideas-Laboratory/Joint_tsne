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

    float applyKernel(const QVector<float>& vec1, const QVector<float>& vec2, KernelFunc func);

    float cosine(const QVector<float>& vec1, const QVector<float>& vec2);
    float rbf(const QVector<float>& vec1, const QVector<float>& vec2, float delta = 1.0);
    float laplacian(const QVector<float>& vec1, const QVector<float>& vec2, float delta = 1.0);

    float L1Distance(const QVector<float>& vec1, const QVector<float>& vec2);
    float L2Distance(const QVector<float>& vec1, const QVector<float>& vec2);

    void normalize(QVector<float>& vec);

    // sort vec and return the sorted indices
    QVector<int> sortIdx(QVector<int> &vec);
}





#endif // UTILS_H
