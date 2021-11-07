#include "math_utils.h"

float math_utils::applyKernel(const QVector<float> &vec1, const QVector<float> &vec2, KernelFunc func)
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

float math_utils::cosine(const QVector<float> &vec1, const QVector<float> &vec2)
{
    float sum = 0.f;
    float len1 = 0.f, len2 = 0.f;
    for(int k = 0; k < vec1.size(); k++)
    {
        // dot product bewteen v1 and v2
        sum += vec1[k]*vec2[k];
        // dot product bewteen v1 and v1
        len1 += vec1[k]*vec1[k];
        // dot product bewteen v2 and v2
        len2 += vec2[k]*vec2[k];
    }
    len1 = sqrtf(len1);
    len2 = sqrtf(len2);

    float s = 0;
    if (len1 != 0 && len2 != 0)
    {
        s = sum / (len1 * len2);
    }

    return s;
}

float math_utils::rbf(const QVector<float> &vec1, const QVector<float> &vec2, float delta)
{
    //
    return exp(-L2Distance(vec1, vec2)/(2*delta*delta));
}

float math_utils::laplacian(const QVector<float> &vec1, const QVector<float> &vec2, float delta)
{
    return exp(-L1Distance(vec1, vec2)/delta);
}

float math_utils::L1Distance(const QVector<float> &vec1, const QVector<float> &vec2)
{
    int size = vec1.size();

    float sum = 0.f;
    for (int i = 0; i < size; i++)
    {
        sum += fabs(vec1[i]-vec2[i]);
    }
    return sum;
}

float math_utils::L2Distance(const QVector<float> &vec1, const QVector<float> &vec2)
{
    int size = vec1.size();

    float sum = 0.f;
    for (int i = 0; i < size; i++)
    {
        sum += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
    }
    return sum;
}

void math_utils::normalize(QVector<float> &vec)
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

QVector<int> math_utils::sortIdx(QVector<int> &vec)
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
