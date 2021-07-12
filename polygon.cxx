#include "client.h"

namespace EPP
{
    void Candidate::close_clockwise(
        std::vector<Point> &polygon)
    {
        Point tail = polygon.front();
        Point head = polygon.back();
        int edge;
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
            switch (edge++ & 3)
            {
            case 0:
                if (tail.j == 0 && tail.i < head.i)
                    head = tail;
                else
                    head = Point(Parameters::N, 0);
                break;
            case 1:
                if (tail.i == 0 && tail.j > head.j)
                    head = tail;
                else
                    head = Point(0, 0);
                break;
            case 2:
                if (tail.j == Parameters::N && tail.i > head.i)
                    head = tail;
                else
                    head = Point(0, Parameters::N);
                break;
            case 3:
                if (tail.i == Parameters::N && tail.j < head.j)
                    head = tail;
                else
                    head = Point(Parameters::N, Parameters::N);
                break;
            }
            polygon.push_back(head);
        }
    }

    // Ramer–Douglas–Peucker algorithm
    void Candidate::simplify(
        const double tolerance,
        std::vector<Point> &simplified,
        const unsigned short int lo,
        const unsigned short int hi)
    {
        if (lo + 1 == hi)
            return;

        double x = separatrix[hi].i - separatrix[lo].i;
        double y = separatrix[hi].j - separatrix[lo].j;
        double theta = atan2(y, x);
        double c = cos(theta);
        double s = sin(theta);
        double max = 0;
        unsigned short int keep;
        for (int mid = lo + 1; mid < hi; mid++)
        { // distance of mid from the line from lo to hi
            double d = std::abs(c * (separatrix[mid].j - separatrix[lo].j) - s * (separatrix[mid].i - separatrix[lo].i));
            if (d > max)
            {
                keep = mid;
                max = d;
            }
        }
        if (max > tolerance) // significant, so something we must keep in here
        {
            simplify(tolerance, simplified, lo, keep);
            simplified.push_back(separatrix[keep]);
            simplify(tolerance, simplified, keep, hi);
        }
        // but if not, we don't need any of the points between lo and hi
    }

    std::vector<Point> Candidate::simplify(
        const double tolerance)
    {
        std::vector<Point> polygon;
        polygon.reserve(separatrix.size());

        polygon.push_back(separatrix[0]);
        simplify(tolerance * Parameters::N, polygon, 0, separatrix.size() - 1);
        polygon.push_back(separatrix[separatrix.size() - 1]);

        return polygon;
    }

    std::vector<Point> Candidate::in_polygon()
    {
        std::vector<Point> polygon;
        polygon.reserve(separatrix.size() + 4);

        for (auto & point : separatrix)
            polygon.push_back(point);

        close_clockwise(polygon);
        return polygon;
    }

    std::vector<Point> Candidate::in_polygon(
        double tolerance)
    {
        std::vector<Point> polygon;
        polygon.reserve(separatrix.size() + 4);

        polygon.push_back(separatrix[0]);
        simplify(tolerance * Parameters::N, polygon, 0, separatrix.size() - 1);
        polygon.push_back(separatrix[separatrix.size() - 1]);

        close_clockwise(polygon);
        return polygon;
    }

    std::vector<Point> Candidate::out_polygon()
    {
        std::vector<Point> polygon;
        polygon.reserve(separatrix.size() + 4);

        for (auto point = separatrix.rbegin(); point != separatrix.rend(); point++)
            polygon.push_back(*point);

        close_clockwise(polygon);
        return polygon;
    }

    std::vector<Point> Candidate::out_polygon(
        double tolerance)
    {
        std::vector<Point> polygon;
        polygon.reserve(separatrix.size() + 4);

        polygon.push_back(separatrix[0]);
        simplify(tolerance * Parameters::N, polygon, 0, separatrix.size() - 1);
        polygon.push_back(separatrix[separatrix.size() - 1]);
        std::reverse(polygon.begin(), polygon.end());

        close_clockwise(polygon);
        return polygon;
    }
}
