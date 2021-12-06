#include "math_utils.h"

float applyKernel(const QVector<float> &vec1, const QVector<float> &vec2, KernelFunc func)
{
    assert(vec1.size() == vec2.size());
    switch (func) {
    case COS:
        return cosine(vec1, vec2);
    case RBF:
        return rbf(vec1, vec2);
    case LAPLACIAN:
        return laplacian(vec1, vec2);
    }
}

float cosine(const QVector<float> &vec1, const QVector<float> &vec2)
{
    assert(vec1.size() == vec2.size());

    float sum = 0.f;
    float sum1 = 0.f, sum2 = 0.f;
    for(int k = 0; k < vec1.size(); k++)
    {
        // dot product bewteen v1 and v2
        sum += vec1[k]*vec2[k];
        // dot product bewteen v1 and v1
        sum1 += vec1[k]*vec1[k];
        // dot product bewteen v2 and v2
        sum2 += vec2[k]*vec2[k];
    }
    float len1 = sqrtf(sum1);
    float len2 = sqrtf(sum2);

    float s = 0;
    if (len1 > 0 && len2 > 0)
    {
        s = sum / (len1 * len2);
    }

    return s;
}

float rbf(const QVector<float> &vec1, const QVector<float> &vec2, float delta)
{
    //
    return exp(-L2Distance(vec1, vec2)/(2*delta*delta));
}

float laplacian(const QVector<float> &vec1, const QVector<float> &vec2, float delta)
{
    return exp(-L1Distance(vec1, vec2)/delta);
}

float L1Distance(const QVector<float> &vec1, const QVector<float> &vec2)
{
    int size = vec1.size();

    float sum = 0.f;
    for (int i = 0; i < size; i++)
    {
        sum += fabs(vec1[i]-vec2[i]);
    }
    return sum;
}

float L2Distance(const QVector<float> &vec1, const QVector<float> &vec2)
{
    int size = vec1.size();

    float sum = 0.f;
    for (int i = 0; i < size; i++)
    {
        sum += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
    }
    return sum;
}

void normalize(QVector<float> &vec)
{
    float sum = 0.f;
    for (int i = 0; i < vec.size(); i++)
    {
        sum += vec[i];
    }

    if (sum != 0.f)
    {
        for (int i = 0; i < vec.size(); i++)
        {
            vec[i]/=sum;
        }
    }
}

QVector<int> sortIdx(QVector<int> &vec)
{
    // all indices
    QVector<int> indices(vec.size());
    for (int i = 0; i < indices.size(); i++)
    {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [vec](int a, int b)
    {
        return vec[a] < vec[b];
    });

    QVector<int> tmp(vec.size());
    for (int i = 0; i < tmp.size(); i++)
    {
        tmp[i] = vec[indices[i]];
    }
    // Reassign the vector
    vec = tmp;

    return indices;
}
