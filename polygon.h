/*
 * Developer: Wayne Moore <wmoore@stanford.edu> 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
#ifndef _EPP_POLYGON_H
#define _EPP_POLYGON_H 1

#include "client.h"

namespace EPP
{
    void Candidate::close_clockwise(
        Polygon &polygon) const noexcept
    {
        Point tail = polygon.front();
        Point head = polygon.back();
        int edge; // which edge to start with
        if (head.j == 0)
            edge = 0;
        if (head.i == 0)
            edge = 1;
        if (head.j == Parameters::N)
            edge = 2;
        if (head.i == Parameters::N)
            edge = 3;
        while (!(head == tail))
        {
            switch (edge++ & 3) // rotate through edges clockwise
            {
            case 0: // bottom
                if (tail.j == 0 && tail.i < head.i) // is the tail on this edge
                    head = tail;                    // and clockwise from head
                else                                // no take the whole thing
                    head = Point(0, 0);             // to the corner
                break;
            case 1: // left
                if (tail.i == 0 && tail.j > head.j)
                    head = tail;
                else
                    head = Point(0, Parameters::N);
                break;
            case 2: // top
                if (tail.j == Parameters::N && tail.i > head.i)
                    head = tail;
                else
                    head = Point(Parameters::N, Parameters::N);
                break;
            case 3: // right
                if (tail.i == Parameters::N && tail.j < head.j)
                    head = tail;
                else
                    head = Point(Parameters::N, 0);
                break;
            }
            polygon.push_back(head);
        }
    }

    // Ramer–Douglas–Peucker algorithm
    void Candidate::simplify(
        const double tolerance,
        Polygon &simplified,
        const size_t lo,
        const size_t hi) const noexcept
    {
        if (lo + 1 == hi)   // empty
            return;

        double x = separatrix[hi].i - separatrix[lo].i;
        double y = separatrix[hi].j - separatrix[lo].j;
        double theta = atan2(y, x); // angle of the line from lo to hi
        double c = cos(theta);
        double s = sin(theta);
        double max = 0;
        size_t keep;
        for (size_t mid = lo + 1; mid < hi; mid++)
        {   // rotate each vector lo to mid around lo then the Y coordinate is the
            // perpendicular distance of mid from the line from lo to hi
            double d = std::abs(c * (separatrix[mid].j - separatrix[lo].j) - s * (separatrix[mid].i - separatrix[lo].i));
            if (d > max)
            {
                keep = mid;
                max = d;
            }
        }
        if (max > tolerance) // significant, so some point we must keep in here
        {
            simplify(tolerance, simplified, lo, keep);
            simplified.push_back(separatrix[keep]);
            simplify(tolerance, simplified, keep, hi);
        }
        // but if not, we don't need any of the points between lo and hi
    }

    // convenience routines
    Polygon Candidate::simplify(
        const double tolerance) const noexcept
    {
        Polygon polygon;
        polygon.reserve(separatrix.size());

        polygon.push_back(separatrix[0]);
        simplify(tolerance * Parameters::N, polygon, 0, separatrix.size() - 1);
        polygon.push_back(separatrix[separatrix.size() - 1]);

        return polygon;
    }

    Polygon Candidate::in_polygon() const noexcept
    {
        Polygon polygon;
        polygon.reserve(separatrix.size() + 4);

        for (auto &point : separatrix)
            polygon.push_back(point);

        close_clockwise(polygon);
        return polygon;
    }

    Polygon Candidate::in_polygon(
        double tolerance) const noexcept
    {
        Polygon polygon;
        polygon.reserve(separatrix.size() + 4);

        polygon.push_back(separatrix[0]);
        simplify(tolerance * Parameters::N, polygon, 0, separatrix.size() - 1);
        polygon.push_back(separatrix[separatrix.size() - 1]);

        close_clockwise(polygon);
        return polygon;
    }

    Polygon Candidate::out_polygon() const noexcept
    {
        Polygon polygon;
        polygon.reserve(separatrix.size() + 4);

        for (auto point = separatrix.rbegin(); point != separatrix.rend(); point++)
            polygon.push_back(*point);

        close_clockwise(polygon);
        return polygon;
    }

    Polygon Candidate::out_polygon(
        double tolerance) const noexcept
    {
        Polygon polygon;
        polygon.reserve(separatrix.size() + 4);

        polygon.push_back(separatrix[0]);
        simplify(tolerance * Parameters::N, polygon, 0, separatrix.size() - 1);
        polygon.push_back(separatrix[separatrix.size() - 1]);
        std::reverse(polygon.begin(), polygon.end());

        close_clockwise(polygon);
        return polygon;
    }
}

#endif /* _EPP_POLYGON_H */
