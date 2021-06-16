#include <client.h>
#include <boundary.h>
#include <modal.h>

#include <random>

namespace EPP
{

    ModalClustering::ModalClustering()
    {
        generate = new std::mt19937(random());
    };

    ModalClustering::~ModalClustering()
    {
        delete generate;
    };

    int clusters = 0;
    int ModalClustering::findClusters(float *density)
    {
        // contiguous set starts empty
        std::fill(_contiguous, _contiguous + (N + 3) * (N + 3), false);
        // most points start undefined
        std::fill(_cluster, _cluster + (N + 3) * (N + 3), -1);
        // out of bounds points are boundary points
        for (int i = 0; i <= N; i++)
        {
            cluster(-1, i) = 0;
            cluster(N + 1, i) = 0;
            cluster(i, -1) = 0;
            cluster(i, N + 1) = 0;
        }

        // collect all the grid points
        grid_vertex *pv = vertex;
        for (int i = 0; i <= N; i++)
            for (int j = 0; j <= N; j++)
            {
                pv->f = density[i + (N + 1) * j];
                pv->i = i;
                pv->j = j;
                pv++;
            }

        // get all comparisons out of the way early and efficiently
        std::sort(vertex, vertex + (N + 1) * (N + 1), decreasing_density);

        // P(outliers) is ~3 SD from zero, i.e., 99% confidence it's not zero.
        int i = (N + 1) * (N + 1);
        double outliers = 0;
        while (outliers < 3)
            outliers += vertex[--i].f / 4 / N / N;

        // this should really be the significance threshold but this won't deadlock for now
        // float threshold = vertex[(50 * (N + 1) * (N + 1)) / 100].f;

        // for the points that are above threshold, i.e., cluster points
        clusters = 0;
        for (pv = vertex; pv < vertex + i; pv++)
        {
            outliers += pv->f / 4 / N / N;
            // if (pv->f < threshold)
            //     break;
            // visit the neighbors to see what clusters they belong to
            int result = -1;
            visit(result, pv->i - 1, pv->j);
            visit(result, pv->i + 1, pv->j);
            visit(result, pv->i, pv->j - 1);
            visit(result, pv->i, pv->j + 1);
            // if we didn't find one this is a new mode
            if (result < 0)
                result = ++clusters;
            cluster(pv->i, pv->j) = result;
            // if this point belongs to a cluster mark the neighbors as being contiguous
            if (result > 0)
            {
                contiguous(pv->i - 1, pv->j) = true;
                contiguous(pv->i + 1, pv->j) = true;
                contiguous(pv->i, pv->j - 1) = true;
                contiguous(pv->i, pv->j + 1) = true;
            }
        }
        // we don't trust these small densities so we take the rest
        // randomly so the border will grow approximately uniformly
        std::shuffle(pv, vertex + (N + 1) * (N + 1), *generate);
        for (; pv < vertex + N * N; pv++)
        {
            // find the next unassigned point that is contiguous with those already classified
            grid_vertex *pw;
            for (pw = pv; pw < vertex + (N + 1) * (N + 1); pw++)
                if (contiguous(pw->i, pw->j))
                    break;
            // if necessary swap it into position
            if (pw != pv)
                std::swap(*pv, *pw);
            // visit the neighbors and then allocate it as a background point
            int result = -1;
            visit(result, pv->i - 1, pv->j);
            visit(result, pv->i + 1, pv->j);
            visit(result, pv->i, pv->j - 1);
            visit(result, pv->i, pv->j + 1);
            cluster(pv->i, pv->j) = result;
            // if this point belongs to a cluster mark the neighbors as being contiguous
            if (result > 0)
            {
                contiguous(pv->i - 1, pv->j) = true;
                contiguous(pv->i + 1, pv->j) = true;
                contiguous(pv->i, pv->j - 1) = true;
                contiguous(pv->i, pv->j + 1) = true;
            }
        }

        return clusters;
    };

    void ModalClustering::getBoundary(float *density, ClusterBoundary &bounds)
    {
        short neighbor[8];

        bounds.clear();
        for (pv = vertex; pv < vertex + (N + 1) * (N + 1); pv++)
            // if this is a boundary point
            if (cluster(pv->i, pv->j) == 0)
            {
                // traverse the neighborhood clockwise
                neighbor[0] = cluster(pv->i, pv->j + 1);
                neighbor[1] = cluster(pv->i + 1, pv->j + 1);
                neighbor[2] = cluster(pv->i + 1, pv->j);
                neighbor[3] = cluster(pv->i + 1, pv->j - 1);
                neighbor[4] = cluster(pv->i, pv->j - 1);
                neighbor[5] = cluster(pv->i - 1, pv->j - 1);
                neighbor[6] = cluster(pv->i - 1, pv->j);
                neighbor[7] = cluster(pv->i - 1, pv->j + 1);
                int rank = 0;
                for (int i = 0; i < 8; i++)
                {
                    short left = neighbor[(i - 1) & 7];
                    short right = neighbor[(i + 1) & 7];
                    float weight = density[pv->i + (N + 1) * pv->j];
                    // if we've found a good edge create the appropriate segment
                    if (left > 0 && right > 0 && neighbor[i] == 0)
                    {
                        const double sqrt2 = sqrt(2);
                        switch (i)
                        {
                        case 0:
                            weight += density[pv->i + (N + 1) * pv->j + (N + 1)];
                            bounds.addSegment(ColoredVertical, pv->i, pv->j, left, right, weight);
                            break;
                        case 1:
                            weight += density[pv->i + 1 + (N + 1) * pv->j + (N + 1)];
                            bounds.addSegment(ColoredRight, pv->i, pv->j, left, right, weight * sqrt2);
                            break;
                        case 2:
                            weight += density[pv->i + (N + 1) + 1 * pv->j];
                            bounds.addSegment(ColoredHorizontal, pv->i, pv->j, left, right, weight);
                            break;
                        case 3:
                            weight += density[pv->i + 1 + (N + 1) * pv->j - (N + 1)];
                            bounds.addSegment(ColoredLeft, pv->i, pv->j - 1, right, left, weight * sqrt2);
                            break;
                        default:
                            // we are only responsible for the half plane head > tail
                            break;
                        };
                        // but we need to look at all of them to see if we have a vertex
                        ++rank;
                    }
                    if (rank != 2)
                    {
                        bounds.addVertex(ColoredPoint<short>(pv->i, pv->j));
                    }
                }
            };

        bounds.setColorful(clusters + 1);
    }
}