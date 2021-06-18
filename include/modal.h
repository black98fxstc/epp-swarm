#include <random>

namespace EPP
{
    typedef unsigned int booleans; // at most 32 clusters or edges
    typedef ColoredBoundary<short, short> ClusterBoundary;
    typedef ColoredMap<short, short> ClusterMap;
    typedef ColoredGraph<booleans> DualGraph;

    class ModalClustering
    {
        // everything is inline because we want the compiler
        // to pare down the inner loop as much as possible

        // accessors with +/-1 slop to avoid bounds checks
        short _cluster[(N + 3) * (N + 3)];
        inline short &cluster(const int &i, const int &j)
        {
            return _cluster[(i + 1) * (N + 3) + (j + 1)];
        };

        bool _contiguous[(N + 3) * (N + 3)];
        inline bool &contiguous(const int &i, const int &j)
        {
            return _contiguous[(i + 1) * (N + 3) + (j + 1)];
        };

        inline void visit(
            int &result,
            const int &i,
            const int &j)
        {
            // if this point has been assigned to a cluster
            if (cluster(i, j) > 0)
                // and our starting point is unassigned or assigned to boundary neighbors
                if (result < 0)
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
        } vertex[(N + 1) * (N + 1)], *pv = vertex;

        struct
        {
            bool operator()(grid_vertex a, grid_vertex b) const { return a.f > b.f; }
        } decreasing_density;

        std::random_device random;
        std::mt19937 *generate;

    public:
        ModalClustering();
        ~ModalClustering();
        int findClusters(const float *density);
        void getBoundary(const float *density, ClusterBoundary &boundary);
    };
}
