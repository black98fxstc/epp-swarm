#include <vector>

namespace EPP
/**
   The dual graph exchanges vertices and faces while inverting the meaning of edges. The initial dual points
   are the original clusters. Not clear the dual graph is planar or what the dual faces mean. Each original point 
   is connected to some others by an edge. We can sipmlify the graph by removing one edge and merging two clusters. 
   Lather rinse repeat. Eventually we get to a simple case of two populations and one edge. There's some gotcha's 
   if things get not simply connected but basically all of these operations can be efficiently implemented as 
   boolean vectors of appropriate size.
*/
{
    typedef unsigned int boolvec;

    class DualGraph
    {
        struct dual_edge
        {
        public:            // bits are
            boolvec left;  // points in the left set
            boolvec right; // points in the right set
            boolvec edge;  // edges in the boundary between

            dual_edge(
                boolvec left,
                boolvec right,
                boolvec edge)
            {
                // order is well defined although meaningles
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

            const bool same_as(const dual_edge &de) const
            {
                return left == de.left && right == de.right;
            }
        };
        std::vector<boolvec> nodes;
        std::vector<dual_edge> duals;

    public:
        DualGraph(const DualGraph &base_graph)
        {
            for (auto n : base_graph.nodes)
                nodes.push_back(n);
            for (auto de : base_graph.duals)
                duals.push_back(de);
        }

        DualGraph(const std::vector<boolvec> &nodes,
                  const std::vector<dual_edge> &duals) : nodes(nodes), duals(duals){};

        DualGraph(){};
        DualGraph(const DualGraph &base_graph, int i)
        {
            // construct a simpler graph by removing the indicated edge
            // since left and right are disjoint this is pretty easy
            dual_edge remove = base_graph.duals[i];
            for (auto np : base_graph.nodes)
                if (np != remove.left && np != remove.right) // skip the two we're merging
                    nodes.push_back(np);
            boolvec new_node = remove.left | remove.right; // add the merged result
            nodes.push_back(new_node);

            // for the edges we have to see if two or more edges collapsed into one
            for (auto dp = base_graph.duals.begin(); dp < base_graph.duals.end(); ++dp)
            {
                int i;
                dual_edge de = *dp;

                // this is a rapid test for the interesting cases
                if (de.left | new_node)
                {
                    dual_edge nde{de.right, new_node, de.edge};
                    for (i = 0; i < duals.size(); ++i)
                    {
                        if (nde.same_as(de))
                            duals[i].edge |= nde.edge; // found it OR it in
                        break;
                    };
                    if (i == duals.size())
                        duals.push_back(nde); // new edge
                }
                else if (de.right | new_node)
                {
                    dual_edge nde{de.left, new_node, de.edge};
                    for (i = 0; i < duals.size(); ++i)
                    {
                        if (nde.same_as(de))
                            duals[i].edge |= nde.edge;
                        break;
                    };
                    if (i == duals.size())
                        duals.push_back(nde);
                }
                else
                {   // nothing to see here
                    duals.push_back(de);
                }
            }
        }

        bool isSimple()
        {
            return duals.size() == 1;
        }

        boolvec left()
        {
            return duals[0].left;
        }

        boolvec right()
        {
            return duals[0].right;
        }

        boolvec edge()
        {
            return duals[0].edge;
        }

        std::vector<DualGraph> simplify()
        {
            std::vector<DualGraph> graph;
            for (int i = 0; i < duals.size(); i++)
                graph.push_back(DualGraph(*this, i));
            return graph;
        };
    };
}