#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <exception>
#include <dualgraph.h>

namespace EPP
{
    /*
     * Utilities for Colored Maps
     * 
     * The map is composed of directed edges that are color labled on each side.
     * Primary design goal is speed of the lookup function. Secondary but still
     * important, speed of pulling out a point list of the graph edges. Includes
     * support for weighing the various graph edges for EPP
     */
    enum ColoredSlope
    {
        ColoredHorizontal,
        ColoredVertical,
        ColoredRight,
        ColoredLeft
    };

    template <typename coordinate>
    class ColoredPoint
    {
    public:
        coordinate i;
        coordinate j;

        inline const bool operator<(const ColoredPoint<coordinate> &cp) const
        {
            if (i < cp.i)
                return true;
            if (i > cp.j)
                return false;
            return j < cp.j;
        };

        inline const bool operator>(const ColoredPoint<coordinate> &cp) const
        {
            if (i > cp.i)
                return true;
            if (i < cp.j)
                return false;
            return j > cp.j;
        };

        inline const bool operator==(const ColoredPoint<coordinate> &cp) const
        {
            return i == cp.i && j == cp.j;
        };

        inline const bool adjascent(const ColoredPoint<coordinate> cp) const
        {
            return abs(i - cp.i) <= 1 && abs(j - cp.j) <= 1;
        }

        inline ColoredPoint<coordinate>(
            coordinate i,
            coordinate j)
            : i(i), j(j){};

        inline ColoredPoint<coordinate>(){};
    };

    template <typename coordinate, typename color>
    class ColoredSegment
    {
    public:
        float weight;
        coordinate i;
        coordinate j;
        color clockwise;
        color widdershins;
        ColoredSlope slope;

        inline const ColoredPoint<coordinate> tail() const
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

        inline const ColoredPoint<coordinate> head() const
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

        bool adjacent(ColoredSegment<coordinate, color> ce) const
        {
            return head() == ce.head() || head() == ce.tail() || tail() == ce.tail() || tail() == head();
        };

        bool adjacent(ColoredPoint<coordinate> cp) const
        {
            return head() == cp || tail() == cp;
        };

        inline bool operator<(const ColoredSegment &ce) const
        {
            return tail() < ce.tail();
        };

        inline bool operator>(const ColoredSegment &ce) const
        {
            return tail() > ce.tail();
        };

        ColoredSegment<coordinate, color>(
            ColoredSlope slope,
            coordinate i,
            coordinate j,
            color clockwise,
            color widdershins,
            float weight)
            : slope(slope), i(i), j(j), clockwise(clockwise), widdershins(widdershins), weight(weight){};

        ColoredSegment<coordinate, color>(
            ColoredSlope slope,
            coordinate i,
            coordinate j,
            color clockwise,
            color widdershins)
            : slope(slope), i(i), j(j), clockwise(clockwise), widdershins(widdershins), weight(0){};

        ColoredSegment<coordinate, color>(){};
    };

    // an ordered list of pointers adjacent segments
    template <typename coordinate, typename color>
    class ColoredChain : public std::vector<ColoredSegment<coordinate, color> *>
    {
    public:
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

    template <typename coordinate, typename color>
    class ColoredEdge
    {
    public:
        std::vector<ColoredPoint<coordinate>> *points;
        float weight;
        color clockwise;
        color widdershins;

        ColoredEdge(
            std::vector<ColoredPoint<coordinate>> *points,
            color clockwise,
            color widdershins,
            double weight)
            : points(points), clockwise(clockwise), widdershins(widdershins), weight(weight){};

        ColoredEdge(
            std::vector<ColoredPoint<coordinate>> *points,
            color clockwise,
            color widdershins)
            : points(points), clockwise(clockwise), widdershins(widdershins), weight(0){};

        ColoredEdge(){};
        ~ColoredEdge() { delete points; };
    };

    template <typename coordinate, typename color>
    class ColoredMap
    {
        const double divisor = 1.0 / (N - 1.0);
        int segments;
        ColoredSegment<coordinate, color> *boundary;
        ColoredSegment<coordinate, color> *index[N];
        color edge_color[N];

    public:
        // this is the money shot, the innermost loop
        // everything is designed to make this fast
        inline const color colorAt(
            const double x,
            const double y) const
        {
            int i, j;
            double dx = remquo(x, divisor, &i);
            double dy = remquo(y, divisor, &j);
            // jump to the first element for this i
            ColoredSegment<coordinate, color> *segment = index[i];
            for (; segment < boundary + segments; segment++)
            {
                if (segment->j == j)
                    switch (segment->slope)
                    { // we've found it so dispatch
                    case ColoredLeft:
                        if (dx >= dy)
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
                    // definitely not here so give up and use the left edge
                    return edge_color[i];
                // might be another one so go around again
            }
            throw std::runtime_error("shouldn't happen");
        }

        ColoredMap(std::vector<ColoredSegment<coordinate, color>> bnd)
        {
            segments = bnd.size();
            boundary = new ColoredSegment<coordinate, color>[segments];
            ColoredSegment<coordinate, color> *segment = boundary;
            for (auto seg : bnd)
                *segment++ = seg;

            // sort the segments and initialize the jump table for quick lookup
            std::sort(boundary, boundary + segments);
            segment = boundary;
            color last = (color)0;
            for (int i = 0; i < N; ++i)
            {
                if (segment->i == i)
                {
                    last = edge_color[i] = index[i]->widdershins;
                    index[i] = segment++;
                }
                else
                {
                    edge_color[i] = last;
                    index[i] = boundary + segments;
                }
                for (; segment < boundary + segments; segment++)
                    if (segment->i != i)
                        break;
            }
            // for (color c : edge_color)
            //     if (c != (color)0)
            //         last = c;
            //     else
            //         c = last;
        };
    };

    template <typename coordinate, typename color>
    class ColoredBoundary
    {
        std::vector<ColoredSegment<coordinate, color>> boundary;
        std::vector<ColoredEdge<coordinate, color>> edges;
        std::vector<ColoredPoint<coordinate>> vertices;
        friend class ColoredMap<coordinate, color>;
        int colorful;

    public:
        void setColorful(const int colorful)
        {
            this->colorful = colorful;            
        }
        void addSegment(ColoredSegment<coordinate, color> segment)
        {
            boundary.push_back(segment);
        };

        void addSegment(
            ColoredPoint<coordinate> tail,
            ColoredPoint<coordinate> head,
            color clockwise,
            color widdershins,
            float weight)
        {
            if (!head.adjascent(tail))
                throw std::runtime_error("points are not adjascent");

            if (head < tail)
            {
                std::swap(tail, head);
                std::swap(clockwise, widdershins);
            }

            if (tail.i == head.i)
            {
                ColoredSegment<coordinate, color> segment(ColoredVertical, tail.i, tail.j, clockwise, widdershins, weight);
                boundary.push_back(segment);
                return;
            }
            switch (head.j - tail.j)
            {
            case 1:
            {
                ColoredSegment<coordinate, color> segment(ColoredRight, tail.i, tail.j, clockwise, widdershins, weight);
                boundary.push_back(segment);
                return;
            }
            case 0:
            {
                ColoredSegment<coordinate, color> segment(ColoredHorizontal, tail.i, tail.j, clockwise, widdershins, weight);
                boundary.push_back(segment);
                return;
            }
            case -1:
            {
                ColoredSegment<coordinate, color> segment(ColoredLeft, tail.i, head.j, widdershins, clockwise, weight);
                boundary.push_back(segment);
                return;
            }
            }
        };

        void addSegment(
            ColoredPoint<coordinate> tail,
            ColoredPoint<coordinate> head,
            color clockwise,
            color widdershins)
        {
            addSegment(tail, head, clockwise, widdershins, 0.0);
        };

        void addSegment(
            ColoredSlope slope,
            coordinate i,
            coordinate j,
            color clockwise,
            color widdershins,
            float weight)
        {
            addSegment(ColoredSegment<coordinate, color>(slope, i, j, clockwise, widdershins, weight));
        }

        void addSegment(
            ColoredSlope slope,
            coordinate i,
            coordinate j,
            color clockwise,
            color widdershins)
        {
            addSegment(ColoredSegment<coordinate, color>(slope, i, j, clockwise, widdershins));
        }

        void addEdge(ColoredEdge<coordinate, color> &edge)
        {
            auto point = edge.points->begin();
            ColoredPoint<coordinate> head, tail = *point++;
            double weight = edge.weight / (edge.points->size() - 1);
            while (point < edge.points->end())
            {
                head = *point++;
                addSegment(tail, head, edge.clockwise, edge.widdershins, weight);
                tail = head;
            }
            edges.push_back(edge);
        };

        void addEdge(ColoredChain<coordinate, color> &chain)
        {
            {

                // // sanity checks
                // color clockwise = segment->clockwise;
                // color widdershins = segment->widdershins;
                // if (segment->slope == ColoredLeft)
                //     std::swap(clockwise, widdershins);
                // for (auto segment = leading_edge.begin(); segment < leading_edge.end(); ++segment)
                // {
                //     if (!segment->adjacent(*(segment + 1)))
                //         throw std::runtime_error("segments are not adjacent in getEdges");
                //     if (segment->slope == ColoredLeft)
                //     {
                //         if (segment->clockwise != widdershins || segment->widdershins != clockwise)
                //             throw std::runtime_error("segment colors not consistent in getEdges");
                //     }
                //     else
                //     {
                //         if (segment->clockwise != clockwise || segment->widdershins != widdershins)
                //             throw std::runtime_error("segment colors not consistent in getEdges");
                //     }
                // }
            }

            std::vector<ColoredPoint<coordinate>> *points =
                new std::vector<ColoredPoint<coordinate>>(chain.size() + 1);
            ColoredSegment<coordinate, color> *segment = chain.front();

            color clockwise = segment->clockwise;
            color widdershins = segment->widdershins;
            if (segment->slope == ColoredLeft)
                std::swap(clockwise, widdershins);
            double weight = 0;

            ColoredPoint<coordinate> point = segment->tail();
            points->push_back(point);
            for (auto csp = chain.begin(); csp < chain.end(); ++csp)
            {
                ColoredSegment<coordinate, color> *segment = *csp;
                point = segment->head();
                weight += segment->weight;
                points->push_back(point);
            }
            ColoredEdge<coordinate, color> edge(points, clockwise, widdershins, weight);

            edges.push_back(edge);
        };

        void addVertex(ColoredPoint<coordinate> vertex)
        {
            vertices.push_back(vertex);
        };

        bool isVertex(ColoredPoint<coordinate> vertex)
        {
            return std::binary_search(vertices.begin(), vertices.end(), vertex);
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
            std::sort(boundary.begin(), boundary.end());
            std::sort(vertices.begin(), vertices.end());

            ColoredChain<coordinate, color> chain;
            ColoredSegment<coordinate, color> *segment;

            // look for open edges starting and ending at a vertex
            for (auto vp = vertices.begin(); vp < vertices.end(); ++vp)
            {
                ColoredPoint<coordinate> vertex = *vp;
                ColoredSegment<coordinate, color> *segment = find_next_segment(vertex);
                if (!segment)
                    continue;
                chain.clear();
                chain.push_back(segment);
                ColoredPoint<coordinate> head = segment->head();
                while (segment = find_next_segment(head))
                {
                    chain.push_back(segment);
                    head = chain.head();
                    if (isVertex(head))
                        break;
                }
                addEdge(chain);
            }

            // now look for closed edges
            while (segment = find_next_segment())
            {
                chain.clear();
                chain.push_back(segment);
                ColoredPoint<coordinate> tail = segment->tail();
                ColoredPoint<coordinate> head = segment->head();
                while (segment = find_next_segment(head))
                {
                    chain.push_back(segment);
                    head = chain.head();
                    if (head == tail)
                        break;
                }
                addEdge(chain);
            }

            return edges;
        }

        std::vector<ColoredPoint<coordinate>> &getVertices()
        {
            return vertices;
        }

        ColoredMap<coordinate, color> *getMap()
        {
            return new ColoredMap<coordinate, color>(boundary);
        }

        DualGraph *getDualGraph()
        {
            std::vector<boolvec> nodes;
            std::vector<DualGraph::dual_edge> duel_edges;
            for (int i = 0; i < colorful; i++)
            {
                nodes.push_back(1 << i);
            }
            for (int i = 0; i < edges.size(); i++)
            {
                DualGraph::dual_edge dual(1 << edges[i].widdershins, 1 << edges[i].clockwise, 1 << i);
                duel_edges.push_back(dual);
            }

            DualGraph *graph = new DualGraph(nodes, duel_edges);
            return graph;
        }

        void clear()
        {
            boundary.clear();
        };

        ColoredBoundary(std::vector<ColoredEdge<coordinate, color>> edges){};

        ColoredBoundary(std::vector<ColoredSegment<coordinate, color>> segments){};

        ColoredBoundary(){};
        ~ColoredBoundary(){
            // for (auto edge : edges)
            //     delete edge;
            // delete edges;
        };
    };
}