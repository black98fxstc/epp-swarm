#ifndef _EPP_BOUNDARY_H
#define EPP_BOUNDARY_H 1

#include <vector>
#include <algorithm>
#include <exception>

namespace EPP
{
    /*
     * Utilities for Colored Maps
     * 
     * The map is composed of directed edges that are color labeled on each side.
     * Primary design goal is speed of the lookup function. Secondary but still
     * important, speed of pulling out a point list of the graph edges. Includes
     * support for weighing the various graph edges for EPP
     */

    // four orientations of the edge within a grid square
    enum ColoredSlope
    {
        ColoredHorizontal,
        ColoredVertical,
        ColoredRight,
        ColoredLeft
    };

    // points are given a raster order so we can search them quickly
    template <typename coordinate>
    class ColoredPoint
    {
    public:
        coordinate i;
        coordinate j;

        inline bool operator<(const ColoredPoint<coordinate> &cp) const
        {
            if (i < cp.i)
                return true;
            if (i > cp.i)
                return false;
            return j < cp.j;
        };

        inline bool operator>(const ColoredPoint<coordinate> &cp) const
        {
            if (i > cp.i)
                return true;
            if (i < cp.i)
                return false;
            return j > cp.j;
        };

        inline bool operator==(const ColoredPoint<coordinate> &cp) const
        {
            return i == cp.i && j == cp.j;
        };

        inline bool adjacent(const ColoredPoint<coordinate> cp) const
        {
            return abs(i - cp.i) <= 1 && abs(j - cp.j) <= 1;
        }

        inline ColoredPoint<coordinate>(
            coordinate i,
            coordinate j)
            : i(i), j(j){};

        inline ColoredPoint<coordinate>() = default;
        ;
    };

    /**
     * a segment is a directed edge with the two sides labeled by color.
     * each segment represents one grid square. for efficiency of lookup 
     * always stored by the lower left coordinate of the square, which 
     * means in the case of the left diagonal we are traversing it backwards.
     */
    template <typename coordinate, typename color>
    class ColoredSegment
    {
    public:
        float weight{};
        coordinate i;
        coordinate j;
        color clockwise;
        color widdershins;
        ColoredSlope slope;

        // to save space, compute the head and tail
        inline ColoredPoint<coordinate> tail() const
        {
            switch (slope)
            {
            case ColoredLeft:
                return ColoredPoint<coordinate>(i + 1, j);
            case ColoredRight:
                return ColoredPoint<coordinate>(i, j);
            case ColoredHorizontal:
                return ColoredPoint<coordinate>(i, j);
            case ColoredVertical:
                return ColoredPoint<coordinate>(i, j);
            }
            throw std::runtime_error("shouldn't happen");
        }

        inline ColoredPoint<coordinate> head() const
        {
            switch (slope)
            {
            case ColoredLeft:
                return ColoredPoint<coordinate>(i, j + 1);
            case ColoredRight:
                return ColoredPoint<coordinate>(i + 1, j + 1);
            case ColoredHorizontal:
                return ColoredPoint<coordinate>(i + 1, j);
            case ColoredVertical:
                return ColoredPoint<coordinate>(i, j + 1);
            }
            throw std::runtime_error("shouldn't happen");
        }

        // the head of one edge connects to the tail of the other
        bool adjacent(ColoredSegment<coordinate, color> ce) const
        {
            return head() == ce.head() || head() == ce.tail() || tail() == ce.tail() || tail() == ce.head();
        };

        // the point is the head or tail of the edge
        bool adjacent(ColoredPoint<coordinate> cp) const
        {
            return head() == cp || tail() == cp;
        };

        // segments are sorted for fast access
        inline bool operator<(const ColoredSegment &cs) const
        {
            if (i < cs.i)
                return true;
            if (i > cs.i)
                return false;
            return j < cs.j;
        };

        inline bool operator>(const ColoredSegment &cs) const
        {
            if (i > cs.i)
                return true;
            if (i < cs.i)
                return false;
            return j > cs.j;
        };
        //		inline bool operator<(const ColoredSegment &ce) const
        //		{
        //			return tail() < ce.tail();
        //		};
        //
        //		inline bool operator>(const ColoredSegment &ce) const
        //		{
        //			return tail() > ce.tail();
        //		};

        ColoredSegment<coordinate, color>(
            ColoredSlope slope,
            coordinate i,
            coordinate j,
            color clockwise,
            color widdershins,
            float weight)
            : slope(slope), i(i), j(j), clockwise(clockwise), widdershins(widdershins), weight(weight > 0 ? weight : std::numeric_limits<float>::min()){};

        ColoredSegment<coordinate, color>(
            ColoredSlope slope,
            coordinate i,
            coordinate j,
            color clockwise,
            color widdershins)
            : slope(slope), i(i), j(j), clockwise(clockwise), widdershins(widdershins), weight(std::numeric_limits<float>::min()){};

        ColoredSegment<coordinate, color>() = default;
        ;
    };

    // an ordered list of pointers to adjacent segments
    template <typename coordinate, typename color>
    class ColoredChain : public std::vector<ColoredSegment<coordinate, color> *>
    {
    public:
        // figure out the head and tail of the chain
        ColoredPoint<coordinate> tail()
        {
            if (this->size() == 1)
                return this->front()->tail();
            ColoredSegment<coordinate, color> *first = this->at(0);
            ColoredSegment<coordinate, color> *second = this->at(1);
            if (first->head() == second->tail() || first->head() == second->head())
                return first->tail();
            else
                return first->head();
        };

        ColoredPoint<coordinate> head()
        {
            if (this->size() == 1)
                return this->back()->head();
            ColoredSegment<coordinate, color> *ultimate = this->at(this->size() - 1);
            ColoredSegment<coordinate, color> *penultimate = this->at(this->size() - 2);
            if (penultimate->tail() == ultimate->tail() || penultimate->head() == ultimate->tail())
                return ultimate->head();
            else
                return ultimate->tail();
        };
    };

    /**
     * A directed and labeled edge longer than one grid square
     * always equivalent to a segment chain though
     */
    template <typename coordinate, typename color>
    class ColoredEdge
    {
    public:
        std::vector<ColoredPoint<coordinate>> points;
        float weight;
        color clockwise;
        color widdershins;

        ColoredEdge(
            std::vector<ColoredPoint<coordinate>> points,
            color clockwise,
            color widdershins,
            double weight)
            : points(points), clockwise(clockwise), widdershins(widdershins), weight((float)weight){};

        ColoredEdge(
            std::vector<ColoredPoint<coordinate>> points,
            color clockwise,
            color widdershins)
            : points(points), clockwise(clockwise), widdershins(widdershins), weight(0){};

        ColoredEdge(const ColoredEdge &that)
        {
            this->points = that.points;
            this->clockwise = that.clockwise;
            this->widdershins = that.widdershins;
            this->weight = that.weight;
        }

        ColoredEdge() = default;
        ;

        ~ColoredEdge() = default;
        ;

        ColoredEdge &operator=(const ColoredEdge &that)
        {
            this->points = that.points;
            this->clockwise = that.clockwise;
            this->widdershins = that.widdershins;
            this->weight = that.weight;
            return *this;
        }

        ColoredEdge &operator=(ColoredEdge &&that) noexcept
        {
            this->points = that.points;
            this->clockwise = that.clockwise;
            this->widdershins = that.widdershins;
            this->weight = that.weight;
        }
    };

    // utility class for rapid lookup of color by map position
    template <typename coordinate, typename color>
    class ColoredMap
    {
        int segments;
        ColoredSegment<coordinate, color> *boundary;
        ColoredSegment<coordinate, color> *index[N];
        color edge_color[N];

    public:
        // this is the money shot, the innermost loop
        // everything is designed to make this fast
        inline color colorAt(
            const double x,
            const double y) const
        {
            int i = (int)(x * N);
            int j = (int)(y * N);
            double dx = x * N - i;
            double dy = y * N - j;
            // jump to the first element for this i
            ColoredSegment<coordinate, color> *segment = index[i];
            for (; segment < boundary + segments; segment++)
            {
                if (segment->j == j)
                    switch (segment->slope)
                    { // we've found it so dispatch
                    case ColoredLeft:
                        if (dx >= 1 - dy)
                            return segment->clockwise;
                        else
                            return segment->widdershins;
                    case ColoredRight:
                        if (dy <= dx)
                            return segment->clockwise;
                        else
                            return segment->widdershins;
                    case ColoredHorizontal:
                        return segment->widdershins;
                    case ColoredVertical:
                        return segment->clockwise;
                    }
                if (segment->j > j || segment->i > i)
                    // definitely not here so give up and use the bottom edge
                    return edge_color[i];
                // might be another one so go around again
            }
            return edge_color[i];
        }

        explicit ColoredMap(std::vector<ColoredSegment<coordinate, color>> bounds)
        {
            segments = bounds.size();
            boundary = new ColoredSegment<coordinate, color>[segments];
            std::copy(bounds.begin(), bounds.end(), boundary);

            // sort the segments and initialize the jump table for quick lookup
            std::sort(boundary, boundary + segments);
            ColoredSegment<coordinate, color> *segment = boundary;
            color outside;
            if (segment->slope == ColoredLeft || segment->slope == ColoredVertical)
                outside = segment->widdershins;
            else
                outside = segment->clockwise;
            for (int i = 0; i < N; ++i)
            {
                if (segment->i == i)
                {
                    if (segment->slope == ColoredLeft)
                        outside = segment->widdershins;
                    else
                        outside = segment->clockwise;
                    edge_color[i] = outside;
                    index[i] = segment++;
                }
                else
                {
                    edge_color[i] = outside;
                    index[i] = boundary + segments;
                }
                for (; segment < boundary + segments; segment++)
                    if (segment->i != i)
                        break;
            }
        };

        ~ColoredMap()
        {
            delete[] boundary;
        }
    };

    /*
   The dual graph exchanges vertices and faces while inverting the meaning of edges. The initial dual points
   are the original clusters. Not clear the dual graph is planar or what the dual faces mean. Each original point 
   is connected to some others by an edge. We can simplify the graph by removing one edge and merging two clusters.
   Lather rinse repeat. Eventually we get to a simple case of two populations and one edge. There's some gotcha's 
   if things get multiply connected but basically all of these operations can be efficiently implemented as 
   boolean vectors of appropriate size.
    */
    template <typename booleans = unsigned int>
    class ColoredGraph
    {

    public:
        struct DualEdge
        {
        public:             // bits are
            booleans left;  // clusters in the left set
            booleans right; // clusters in the right set
            booleans edge;  // edges in the boundary between

            DualEdge(
                booleans left,
                booleans right,
                booleans edge)
            {
                // order is well defined although meaningless
                // except that it makes comparisons faster
                // since the edges are not directed
                if (left < right)
                {
                    this->left = left;
                    this->right = right;
                }
                else
                {
                    this->left = right;
                    this->right = left;
                }
                this->edge = edge;
            };

            inline bool same_as(const DualEdge &de) const
            {
                return left == de.left && right == de.right;
            }

            DualEdge() = default;
            ;
        };

        std::vector<booleans> nodes;
        std::vector<DualEdge> duals;

        ColoredGraph() = default;

        // implement move semantics
        ColoredGraph(std::vector<booleans> &nodes,
                     std::vector<DualEdge> &duals) : nodes(nodes), duals(duals){};

        ColoredGraph(ColoredGraph &&other) : nodes(other.nodes), duals(other.duals){};

        ColoredGraph &operator=(ColoredGraph &&other) noexcept
        {
            if (this != other)
            {
                this->nodes = other.nodes;
                this->duals = other.duals;
            }
            return *this;
        }

        // copy constructor
        ColoredGraph(const ColoredGraph &other) : nodes(other.nodes), duals(other.duals){};

        inline bool isSimple() const
        {
            return duals.size() == 1;
        }

        inline booleans left() const
        {
            return duals[0].left;
        }

        inline booleans right() const
        {
            return duals[0].right;
        }

        inline booleans edge() const
        {
            return duals[0].edge;
        }

        std::vector<ColoredGraph> simplify() const
        {
            std::vector<booleans> nodes;
            nodes.reserve(this->nodes.size() - 1);
            std::vector<DualEdge> duals;
            duals.reserve(this->duals.size() - 1);
            std::vector<ColoredGraph> graphs;
            graphs.reserve(this->duals.size());

            for (int i = 0; i < this->duals.size(); i++)
            {
                nodes.clear();
                duals.clear();
                // construct a simpler graph by removing the indicated edge
                // since left and right are disjoint this is pretty easy for the nodes
                DualEdge remove = this->duals[i];
                booleans new_node = remove.left | remove.right; // the merged result
                for (auto np : this->nodes)
                    if (!(np & new_node)) // skip the two we're merging
                        nodes.push_back(np);
                nodes.push_back(new_node); // add the merged node

                // for the edges we have to see if two or more edges collapsed into one
                for (int j = 0; j < this->duals.size(); j++)
                {
                    // skip the one we're removing
                    if (i == j)
                        continue;
                    DualEdge de = this->duals[j];
                    // this is a rapid test for the interesting cases. because of the disjunction of the
                    // nodes this is equivalent to (de.left == remove.left || de.left == remove.right)
                    int k;
                    if (de.left & new_node)
                    {
                        DualEdge nde{de.right, new_node, de.edge};
                        for (k = 0; k < duals.size(); ++k)
                            if (nde.same_as(duals[k]))
                            {
                                duals[k].edge |= nde.edge; // found it OR it in
                                break;
                            }
                        if (k == duals.size())
                            duals.push_back(nde); // new edge
                    }
                    else if (de.right & new_node)
                    {
                        DualEdge nde{de.left, new_node, de.edge};
                        for (k = 0; k < duals.size(); ++k)
                        {
                            duals[k].edge |= nde.edge; // found it OR it in
                            break;
                        }
                        if (k == duals.size())
                            duals.push_back(nde);
                    }
                    else
                    { // nothing to see here copy it forward
                        duals.push_back(de);
                    }
                }
                graphs.push_back(ColoredGraph(nodes, duals));
            }
            return graphs;
        };
    };

    template <typename coordinate, typename color, typename booleans = unsigned int>
    class ColoredBoundary
    {
        std::vector<ColoredSegment<coordinate, color>> boundary;
        std::vector<ColoredEdge<coordinate, color>> edges;
        std::vector<ColoredPoint<coordinate>> vertices;
        friend class ColoredMap<coordinate, color>;
        color colorful;

    public:
        void setColorful(const int colors)
        {
            this->colorful = colors;
            std::sort(vertices.begin(), vertices.end());
        }

        color getColorful() const { return colorful; };

        inline void addSegment(ColoredSegment<coordinate, color> segment)
        {
            boundary.push_back(segment);
        };

        void addSegment(
            ColoredSlope slope,
            int i,
            int j,
            color clockwise,
            color widdershins,
            double weight)
        {
            addSegment(ColoredSegment<coordinate, color>(slope, (coordinate)i, (coordinate)j, clockwise, widdershins, (float)weight));
        }

        void addSegment(
            ColoredSlope slope,
            int i,
            int j,
            color clockwise,
            color widdershins)
        {
            addSegment(ColoredSegment<coordinate, color>(slope, (coordinate)i, (coordinate)j, clockwise, widdershins, 0));
        }

        void addSegment(
            ColoredPoint<coordinate> tail,
            ColoredPoint<coordinate> head,
            color clockwise,
            color widdershins,
            float weight)
        {
            if (!head.adjacent(tail))
                throw std::runtime_error("points are not adjacent");

            if (head < tail)
            {
                std::swap(tail, head);
                std::swap(clockwise, widdershins);
            }

            if (tail.i == head.i)
            {
                addSegment(ColoredSegment<coordinate, color>(ColoredVertical, tail.i, tail.j, clockwise, widdershins, weight));
                return;
            }
            switch (head.j - tail.j)
            {
            case 1:
                addSegment(ColoredSegment<coordinate, color>(ColoredRight, tail.i, tail.j, clockwise, widdershins, weight));
                return;

            case 0:
                addSegment(ColoredSegment<coordinate, color>(ColoredHorizontal, tail.i, tail.j, clockwise, widdershins, weight));
                return;

            case -1:
                addSegment(ColoredSegment<coordinate, color>(ColoredLeft, tail.i, head.j, widdershins, clockwise, weight));
                return;
            }
        };

        void addSegment(
            ColoredPoint<coordinate> tail,
            ColoredPoint<coordinate> head,
            color clockwise,
            color widdershins)
        {
            addSegment(tail, head, clockwise, widdershins, 0);
        };

        void addVertex(ColoredPoint<coordinate> vertex) { vertices.push_back(vertex); };

        std::vector<ColoredPoint<coordinate>> &getVertices() { return vertices; }

        bool isVertex(ColoredPoint<coordinate> vertex)
        {
            return std::binary_search(vertices.begin(), vertices.end(), vertex);
        };

        void addEdge(ColoredEdge<coordinate, color> &edge)
        {
            auto point = edge.points.begin();
            ColoredPoint<coordinate> head, tail = *point++;
            double weight = edge.weight / (edge.points.size() - 1);
            while (point < edge.points.end())
            {
                head = *point++;
                addSegment(tail, head, edge.clockwise, edge.widdershins, weight);
                tail = head;
            }
            edges.push_back(edge);
        };

        void addEdge(std::vector<ColoredPoint<coordinate>> points, color clockwise, color widdershins, double weight)
        {
            ColoredEdge<coordinate, color> nce(points, clockwise, widdershins, weight);
            addEdge(nce);
        };

        void addEdge(std::vector<ColoredPoint<coordinate>> points, color clockwise, color widdershins)
        {
            addEdge(points, clockwise, widdershins, 0.0);
        };

        void addEdge(ColoredChain<coordinate, color> &chain)
        {
            std::vector<ColoredPoint<coordinate>> points;
            points.reserve(chain.size() + 1);
            ColoredSegment<coordinate, color> *segment = chain.front();

            color clockwise = segment->clockwise;
            color widdershins = segment->widdershins;
            double weight = 0;

            ColoredPoint<coordinate> point = chain.tail();
            if (segment->head() == point)
                std::swap(clockwise, widdershins);
            points.push_back(point);
            for (auto csp = chain.begin(); csp < chain.end(); ++csp)
            {
                segment = *csp;
                if (point == segment->tail())
                {
                    if (clockwise != segment->clockwise || widdershins != segment->widdershins)
                        throw std::runtime_error("segment chain is not consistent");
                    point = segment->head();
                }
                else
                {
                    if (clockwise != segment->widdershins || widdershins != segment->clockwise)
                        throw std::runtime_error("segment chain is not consistent");
                    point = segment->tail();
                }
                weight += segment->weight;
                points.push_back(point);
            }
            ColoredEdge<coordinate, color> edge(points, clockwise, widdershins, weight);

            edges.push_back(edge);
        };

        std::vector<bool> *done;
        // find next segment adjacent to a point
        ColoredSegment<coordinate, color> *find_next_segment(
            ColoredPoint<coordinate> point)
        {
            // we know the next adjacent segment can't be far away
            ColoredSegment<coordinate, color> low, high;
            low.i = point.i - 1;
            low.j = point.j - 1;
            high.i = point.i + 1;
            high.j = point.j + 2; // strict upper bound
            // so it's contained in a small interval
            // that we can find quickly since they are sorted
            auto lower = std::lower_bound(boundary.begin(), boundary.end(), low);
            auto upper = std::upper_bound(lower, boundary.end(), high);
            // and then we use brute force
            for (auto cp = lower; cp != upper; ++cp)
                if (!(*done)[cp - boundary.begin()])
                {
                    ColoredSegment<coordinate, color> *candidate = &(*cp);
                    if (!candidate->adjacent(point))
                        continue;
                    (*done)[cp - boundary.begin()] = true;
                    return candidate;
                }
            return NULL;
        }

        // find any segment we haven't considered yet
        ColoredSegment<coordinate, color> *find_next_segment()
        {
            for (auto csp = boundary.begin(); csp != boundary.end(); ++csp)
                if (!(*done)[csp - boundary.begin()])
                {
                    ColoredSegment<coordinate, color> *candidate = &(*csp);
                    (*done)[csp - boundary.begin()] = true;
                    return candidate;
                }
            return NULL;
        }

        // this is the other hard problem but uses
        // much less total time than the lookup
        std::vector<ColoredEdge<coordinate, color>> &getEdges()
        {
            done = new std::vector<bool>(boundary.size());
            edges.clear();
            std::sort(boundary.begin(), boundary.end());
            std::sort(vertices.begin(), vertices.end());

            ColoredChain<coordinate, color> chain;
            ColoredSegment<coordinate, color> *segment;

            // look for open edges starting and ending at a vertex
            for (auto vp = vertices.begin(); vp < vertices.end();)
            {
                ColoredPoint<coordinate> vertex = *vp;
                ColoredSegment<coordinate, color> *segment = find_next_segment(vertex);
                if (!segment)
                {
                    vp++;
                    continue;
                }
                chain.clear();
                chain.push_back(segment);
                ColoredPoint<coordinate> point;
                if (segment->tail() == vertex)
                    point = segment->head();
                else
                    point = segment->tail();
                while (!isVertex(point))
                {
                    segment = find_next_segment(point);
                    if (segment == NULL)
                        break;
                    chain.push_back(segment);
                    if (segment->tail() == point)
                        point = segment->head();
                    else
                        point = segment->tail();
                }
                addEdge(chain);
            }

            // now look for closed edges
            while ((segment = find_next_segment()))
            {
                chain.clear();
                chain.push_back(segment);
                ColoredPoint<coordinate> tail = segment->tail();
                ColoredPoint<coordinate> head = segment->head();
                while ((segment = find_next_segment(head)))
                {
                    chain.push_back(segment);
                    head = chain.head();
                    if (head == tail)
                        break;
                }
                addEdge(chain);
            }

            delete done;
            return edges;
        }

        std::unique_ptr<ColoredMap<coordinate, color>> getMap()
        {
            return std::unique_ptr<ColoredMap<coordinate, color>>(new ColoredMap<coordinate, color>(boundary));
        }

        std::unique_ptr<ColoredGraph<booleans>> getDualGraph()
        {
            std::vector<booleans> nodes(colorful);
            std::vector<typename ColoredGraph<booleans>::DualEdge> duals(edges.size());
            for (int i = 0; i < colorful; i++)
            {
                nodes[i] = 1 << i;
            }
            for (int i = 0; i < edges.size(); i++)
            {
                typename ColoredGraph<booleans>::DualEdge dual(1 << edges[i].widdershins, 1 << edges[i].clockwise, 1 << i);
                duals[i] = dual;
            }

            ColoredGraph<booleans> *graph = new ColoredGraph<booleans>(nodes, duals);
            return std::unique_ptr<ColoredGraph<booleans>>(graph);
        }

        void clear()
        {
            boundary.clear();
            edges.clear();
            vertices.clear();
            colorful = (color)0;
        };

        explicit ColoredBoundary(std::vector<ColoredEdge<coordinate, color>> edges){};

        explicit ColoredBoundary(std::vector<ColoredChain<coordinate, color>> edges){};

        explicit ColoredBoundary(std::vector<ColoredSegment<coordinate, color>> segments){};

        ColoredBoundary() = default;
        ~ColoredBoundary() = default;
    };
}
#endif /* boundary.h */