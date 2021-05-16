#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <exception>

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
        std::vector<ColoredPoint<coordinate>> vertices;
        friend class ColoredMap<coordinate, color>;

    public:
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

        void addEdge(ColoredEdge<coordinate, color> edge)
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
        };

        void addVertex(ColoredPoint<coordinate> vertex)
        {
            vertices.push_back(vertex);
        };

        bool isVertex(ColoredPoint<coordinate> vertex)
        {
            return std::binary_search(vertices.begin(), vertices.end(), vertex);
        };

        // this is the other hard problem but uses
        // much less total time than the lookup
        std::vector<ColoredEdge<coordinate, color> *> *getEdges()
        {
            std::vector<bool> done(boundary.size());
            std::sort(boundary.begin(), boundary.end());
            std::sort(vertices.begin(), vertices.end());
            
            std::vector<ColoredSegment<coordinate, color>> leading_edge;

            ColoredSegment<coordinate, color> *data = leading_edge.data();
            ColoredSegment<coordinate, color> *segment;

            ColoredPoint<coordinate> head;
            ColoredSegment<coordinate, color> low, high;

            // we know the next point can't be far away
            low.i = head.i - 1;
            low.j = head.j - 1;
            high.i = head.i + 1;
            high.j = head.j + 2; // strict upper bound
            // so it's contained in a small interval
            // that we can find quickly since they sre sorted
            auto lower = std::lower_bound(boundary.begin(), boundary.end(), low);
            auto upper = std::upper_bound(lower, boundary.end(), high);
            for (auto candidate = lower; candidate != upper; ++candidate)
            {
                if (done[candidate - boundary.begin()])
                    continue;
                // relatively cheap prequalifier
                if (!segment->adjacent(*candidate))
                    continue;
                if (head == candidate->tail())
                {
                    head = candidate->head();
                }
                else if (head == candidate->head())
                {
                    head = candidate->tail();
                }
                else
                    continue;
                done[candidate - boundary.begin()] = true;
                leading_edge.push_back(*candidate);
                break;
            }
            if (isVertex(head))
            {
            }

            // sanity checks
            color clockwise = segment->clockwise;
            color widdershins = segment->widdershins;
            if (segment->slope == ColoredLeft)
                std::swap(clockwise, widdershins);
            for (segment = data; segment < data + leading_edge.size() - 1; segment++)
            {
                if (!segment->adjacent(*(segment + 1)))
                    throw std::runtime_error("segments are not adjacent in getEdges");
                if (segment->slope == ColoredLeft)
                {
                    if (segment->clockwise != widdershins || segment->widdershins != clockwise)
                        throw std::runtime_error("segment colors not consistent in getEdges");
                }
                else
                {
                    if (segment->clockwise != clockwise || segment->widdershins != widdershins)
                        throw std::runtime_error("segment colors not consistent in getEdges");
                }
            }

            std::vector<ColoredEdge<coordinate, color> *> *edges =
                new std::vector<ColoredEdge<coordinate, color> *>();

            std::vector<ColoredPoint<coordinate>> *points =
                new std::vector<ColoredPoint<coordinate>>(leading_edge.size() + 1);
            ColoredPoint<coordinate> point;
            point = segment->tail();
            points->push_back(point);
            double weight = 0;
            for (segment = data; segment < data + leading_edge.size(); segment++)
            {
                point = segment->head();
                weight += segment->weight;
                points->push_back(point);
            }
            ColoredEdge<coordinate, color> *ce =
                new ColoredEdge<coordinate, color>(points, clockwise, widdershins, weight);
            edges->push_back(ce);

            return edges;
        }

        std::vector<ColoredPoint<coordinate>> getVertices()
        {
            return NULL;
        }

        ColoredMap<coordinate, color> *getMap()
        {
            return new ColoredMap<coordinate, color>(boundary);
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