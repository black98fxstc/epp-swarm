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
     * importand, speed of pulling out a point list of the graph edges. Includes
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

    template <typename coordinate>
    class ColoredNeighborhood
    {
    public:
        const ColoredPoint<coordinate> lower;
        const ColoredPoint<coordinate> upper;

        ColoredNeighborhood<coordinate>(
            ColoredPoint<coordinate> lower,
            ColoredPoint<coordinate> upper) : lower(lower), upper(upper){};
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

        inline bool adjacent(ColoredSegment<coordinate, color> ce)
        {
            return head() == ce.head() || head() == ce.tail() || tail() == ce.tail() || tail() == head();
        };

        inline bool adjacent(ColoredPoint<coordinate> cp)
        {
            return head() == cp || tail() == cp;
        };

        inline bool operator<(const ColoredSegment &ce)
        {
            return tail() < ce.tail();
        };

        inline bool operator>(const ColoredSegment &ce)
        {
            return tail() > ce.tail();
        };

        ColoredNeighborhood<coordinate> neighborhood(
            coordinate i1,
            coordinate j1,
            coordinate i2,
            coordinate j2)
        {
            return ColoredNeighborhood<coordinate>(
                ColoredPoint<coordinate>(i1, j1),
                ColoredPoint<coordinate>(i2, j2));
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
    class ColoredEdge : public std::vector<ColoredPoint<coordinate>>
    {
    public:
        double weight;
        color clockwise;
        color widdershins;
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
                // definitely not here so give up and use the left edge
                if (segment->j > j || segment->i > i)
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
            auto point = edge.begin();
            ColoredPoint<coordinate> head, tail = *point++;
            while (point < edge.end())
            {
                head = *point++;
                addSegment(tail, head, edge.clockwise, edge.widdershins, edge.weight / (edge.size() - 1));
                tail = head;
            }
        };

        void addVertex(ColoredPoint<coordinate> vertex){

        };

        std::vector<std::vector<ColoredEdge<coordinate, color>>> getEdges()
        {
            return (std::vector<std::vector<ColoredEdge<coordinate, color>>>)0;
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
    };
}