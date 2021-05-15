#include <random>

namespace EPP
{
    typedef ColoredBoundary<short, short> ClusterBoundary;
    typedef ColoredMap<short, short> ClusterMap;

    class ModalClustering
    {
        // everthing is inline because we want the compiler
        // to pare down the inner loop as much as possible

        // accessors with +/-1 slop to avoid bounds checks
        short _cluster[(N + 2) * (N + 2)];
        inline short &cluster(const short &i, const short &j)
        {
            return _cluster[(i + 1) * (N + 2) * (j + 1)];
        };

        bool _contiguous[(N + 2) * (N + 2)];
        inline bool &contiguous(const short &i, const short &j)
        {
            return _contiguous[(i + 1) * (N + 2) * (j + 1)];
        };

        inline void visit(
            int &result,
            short i,
            short j)
        {
            // mark this point as contiguous with a classified point
            contiguous(i, j) = true;

            // if this point has been assigned to a cluster
            if (cluster(i, j) > 0)
                // and our starting point is unassigned or assigned to boundary neighbors
                if (result < 1)
                    // assign it to our cluster
                    result = cluster(i, j);
                else if (result != cluster(i, j))
                    // if we found something different it's a boundary point
                    result = 0;
        };

        struct grid_vertex
        {
            float f;
            short i, j;
        } vertex[N * N], *pv = vertex;

        struct
        {
            bool operator()(grid_vertex a, grid_vertex b) const { return a.f > b.f; }
        } decreasing_density;

        std::random_device rd;

    public:
        ModalClustering();
        int cluster(float *density);
        ClusterBoundary boundary();
    };
}
