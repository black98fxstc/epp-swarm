#ifndef _EPP_MODAL_H
#define _EPP_MODAL_H 1

#include <random>
#include <algorithm>

#include "boundary.h"

namespace EPP
{
	typedef unsigned int booleans;
	const int max_booleans = sizeof(booleans) * 8; // max clusters or edges in dual graph
	typedef ColoredGraph<booleans> DualGraph;

	typedef ColoredBoundary<short, short, booleans> ClusterBoundary;
	typedef ColoredMap<short, short> ClusterMap;
	typedef ColoredEdge<short, short> ClusterEdge;
	typedef std::vector<ColoredEdge<short, bool>> ClusterSeparatrix;
	typedef ColoredPoint<short> ClusterPoint;

	class ModalClustering
	{
		int clusters;

		// everything is inline because we want the compiler
		// to pare down the inner loop as much as possible

		// accessors with +/-1 slop to avoid bounds checks
		short _cluster[(N + 3) * (N + 3)];
		inline short &cluster(const int &i, const int &j) noexcept
		{
			return _cluster[(i + 1) * (N + 3) + (j + 1)];
		};

		bool _contiguous[(N + 3) * (N + 3)];
		inline bool &contiguous(const int &i, const int &j) noexcept
		{
			return _contiguous[(i + 1) * (N + 3) + (j + 1)];
		};

		inline void visit(
			int &result,
			const int &i,
			const int &j) noexcept
		{
			// if this point has been assigned to a cluster
			if (cluster(i, j) > 0)
				// and our starting point is unassigned
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
		} vertex[(N + 1) * (N + 1)], *pv;

		struct
		{
			bool operator()(grid_vertex a, grid_vertex b) const noexcept { return a.f > b.f; }
		} decreasing_density;

		std::random_device random;
		std::mt19937 *generate;

	public:
		ModalClustering() noexcept;
		~ModalClustering();
		int findClusters(const float *density) noexcept;
		void getBoundary(const float *density, ClusterBoundary &boundary) noexcept;
	};

	ModalClustering::ModalClustering() noexcept
	{
		generate = new std::mt19937(random());
	}

	ModalClustering::~ModalClustering()
	{
		delete generate;
	}

	int ModalClustering::findClusters(const float *density) noexcept
	{
		clusters = 0;
		// contiguous set starts empty
		std::fill(_contiguous, _contiguous + (N + 3) * (N + 3), false);
		// most points start undefined
		std::fill(_cluster, _cluster + (N + 3) * (N + 3), -1);

		// collect all the grid points
		grid_vertex *pv = vertex;
		for (short i = 0; i <= N; i++)
			for (short j = 0; j <= N; j++)
			{
				pv->f = density[i + (N + 1) * j];
				pv->i = i;
				pv->j = j;
				pv++;
			}

		// get all comparisons out of the way early and efficiently
		std::sort(vertex, vertex + (N + 1) * (N + 1), decreasing_density);

		int A = 3.1415926 * W * W * N * N;	// area of kernel
		if (A < 8)
			A = 8;
		double count = 0;
		int i = (N + 1) * (N + 1);
		for (int a = 0; a < A; a++)
			count += vertex[--i].f / 4 / N / N; // approximate with filter unnormalized
		int j = (N + 1) * (N + 1);
		while (count < 25 && i > 0)	// count is two standard deviations away from zero
		{
			count += vertex[--i].f / 4 / N / N;
			count -= vertex[--j].f / 4 / N / N;
		}
		assert(i > 0);
		for (pv = vertex; pv < vertex + i; pv++)
		{
			// visit the neighbors to see what clusters they belong to
			int result = -1;
			visit(result, pv->i + 1, pv->j);
			visit(result, pv->i - 1, pv->j);
			visit(result, pv->i, pv->j + 1);
			visit(result, pv->i, pv->j - 1);
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
		// std::shuffle(pv, vertex + (N + 1) * (N + 1), *generate);
		for (; pv < vertex + (N + 1) * (N + 1); pv++)
		{
			// find the next unassigned point that is contiguous with those already classified
			grid_vertex *pw;
			for (pw = pv; pw < vertex + (N + 1) * (N + 1); pw++)
				if (contiguous(pw->i, pw->j))
					break;
			// if necessary swap it into position
			if (pw != pv && pw < vertex + (N + 1) * (N + 1))
				std::swap(*pv, *pw);
			// visit the neighbors and then allocate it as a background point
			int result = -1;
			visit(result, pv->i - 1, pv->j);
			visit(result, pv->i + 1, pv->j);
			visit(result, pv->i, pv->j - 1);
			visit(result, pv->i, pv->j + 1);

			if (result < 0)
				result = 0;
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
				// for (int i = 0; i <= N; i++)
				// {
				// 	for (int j = 0; j <= N; j++)
				// 	{
				// 		char ctr;
				// 		int c = cluster(i, j);
				// 		if (c == 0)
				// 			ctr = '+';
				// 		else if (c > 9)
				// 			ctr = 'A' + c - 10;
				// 		else
				// 			ctr = '0' + c;
				// 		std::cout << ctr;
				// 	}
				// 	std::cout << std::endl;
				// }

		return clusters;
	}

	void ModalClustering::getBoundary(const float *density, ClusterBoundary &bounds) noexcept
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
				bool on_edge = false;
				for (int i = 0; i < 8; i++)
				{
					short left = neighbor[(i - 1) & 7];
					short right = neighbor[(i + 1) & 7];
					if (left > 0 && neighbor[i] == 0)
					{
						if (right == 0)
						{
							right = neighbor[(i + 2) & 7];
							if (right == 0)
							{
								right = neighbor[(i + 3) & 7];
								if (right > 0)
								{ // when there are multiple choices for head,
									// take the shortest one, i.e., the one that
									// aligns with the axes
									if (i & 1)
										i++;
								}
								else
									continue;
							}
							else if (right < 0)
								continue;
							else
							{
								if (i & 1)
									i++;
							}
						}
						else if (right < 0)
							continue;
					}
					else
					{
						if (left < 0)
							on_edge = true;
						continue;
					}

					float weight = density[pv->i + (N + 1) * pv->j];
					const double sqrt2 = sqrt(2);
					switch (i & 7)
					{
					case 0:
						weight += density[pv->i + (N + 1) * pv->j + (N + 1)];
						bounds.addSegment(ColoredVertical, pv->i, pv->j, right, left, weight);
						break;
					case 1:
						weight += density[pv->i + 1 + (N + 1) * pv->j + (N + 1)];
						bounds.addSegment(ColoredRight, pv->i, pv->j, right, left, weight * sqrt2);
						break;
					case 2:
						weight += density[pv->i + 1 + (N + 1) * pv->j];
						bounds.addSegment(ColoredHorizontal, pv->i, pv->j, right, left, weight);
						break;
					case 3:
						weight += density[pv->i + 1 + (N + 1) * pv->j - (N + 1)];
						bounds.addSegment(ColoredLeft, pv->i, pv->j - 1, left, right, weight * sqrt2);
						break;
					default:
						// we are only responsible for the half plane head > tail
						break;
					}
					// but we need to look at all of them to see if we have a vertex
					++rank;
				}
				if (rank != 2 || on_edge)
				{
					bounds.addVertex(ColoredPoint<short>(pv->i, pv->j));
				}
				if (neighbor[0] == 0 && neighbor[1] == 0 && neighbor[2] == 0)
				{
					bounds.addVertex(ColoredPoint<short>(pv->i, pv->j));
					bounds.addVertex(ColoredPoint<short>(pv->i + 1, pv->j + 1));
					bounds.addVertex(ColoredPoint<short>(pv->i + 1, pv->j));
					bounds.addVertex(ColoredPoint<short>(pv->i, pv->j + 1));
				}
			}
		bounds.setColorful(clusters + 1);

		// std::cout << std::endl;
		// for (int i = 0; i <= N; i++)
		// {
		// 	for (int j = 0; j <= N; j++)
		// 	{
		// 		char ctr;
		// 		int c = cluster(i, j);
		// 		if (c == 0)
		// 			if (bounds.isVertex(ColoredPoint<short>(i, j)))
		// 				ctr = '*';
		// 			else
		// 				ctr = '+';
		// 		else if (c > 9)
		// 			ctr = 'A' + c - 10;
		// 		else
		// 			ctr = '0' + c;
		// 		std::cout << ctr;
		// 	}
		// 	std::cout << std::endl;
		// }
	}
}
#endif /* _EPP_MODAL_H */