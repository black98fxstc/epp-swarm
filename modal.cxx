#include <client.h>
#include <modal.h>
#include <random>

namespace EPP
{
    ModalClustering::ModalClustering(){};

    int ModalClustering::cluster(float *density)
    {
        // contiguous set starts empty
        std::fill(_contiguous, _contiguous + (N + 2) * (N + 2), false);
        // most points start undefined
        std::fill(_cluster, _cluster + (N + 2) * (N + 2), -1);
        // out of bounds points are boundary points
        for (int i = 0; i < N; i++)
        {
            cluster(-1, i) = 0;
            cluster(N, i) = 0;
            cluster(i, -1) = 0;
            cluster(i, N) = 0;
        }

        // collect all the grid points
        grid_vertex *pv = vertex;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                pv->f = density[i + N * j];
                pv->i = i;
                pv->j = j;
                pv++;
            }

        // get all comparisons out of the way early and efficiently
        std::sort(vertex, vertex + N * N, decreasing_density);

        // this should really be the significance threshold but this won't deadlock for now
        float threshold = vertex[N * N / 2].f;

        // for the points that are above threshold, i.e., cluster points
        int clusters = 0;
        for (pv = vertex; pv < vertex + N * N; pv++)
        {
            if (pv->f < threshold)
                break;
            // visit the neighbors to see what clusters they belong to
            int result = -1;
            visit(result, pv->i - 1, pv->j);
            visit(result, pv->i + 1, pv->j);
            visit(result, pv->i, pv->j - 1);
            visit(result, pv->i, pv->j + 1);
            // if we didn't find one this is a new mode
            if (result < 0)
                cluster(pv->i, pv->j) = ++clusters;
            else
                cluster(pv->i, pv->j) = result;
        }
        // we don't trust these small densities so we take the rest
        // randomly so the border will grow approximately uniformly
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(pv, vertex + N * N, g);
        for (; pv < vertex + N * N; pv++)
        {
            // find the next unassigned point that is contiguous with those already classified
            grid_vertex *pw;
            for (pw = pv; pw < vertex + N * N; pw++)
                if (contiguous(pw->i, pw->j))
                    break;
            // if necessary swap it into position
            if (pw != pv)
            {
                grid_vertex t;
                t = *pv;
                *pv = *pw;
                *pw = t;
            }
            // visit the neighbors and then allocate it as a background point
            int result = -1;
            visit(result, pv->i - 1, pv->j);
            visit(result, pv->i + 1, pv->j);
            visit(result, pv->i, pv->j - 1);
            visit(result, pv->i, pv->j + 1);
            if (result < 0)
                cluster(pv->i, pv->j) = 0;
            else
                cluster(pv->i, pv->j) = result;
        }

        return clusters;
    };
}