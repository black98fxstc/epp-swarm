#include <vector>

namespace EPP
/**
   The dual graph exchanges vertices and faces while inverting the meaning of edges. The initial dual points
   are the original clusters. Not clear the dual graph is planar or what the dual faces mean. Each original point 
   is connected to some others by an edge. We can sipmlify the graph by removing one edge and merging two clusters. 
   Lather rinse repeat. Eventually we get to a simple case of two populations and one edge. There's some gotcha's 
   if things get multiply connected but basically all of these operations can be efficiently implemented as 
   boolean vectors of appropriate size.
*/
{
    typedef unsigned int boolvec;

    class DualGraph
    {
        struct DualEdge
        {
        public:            // bits are
            boolvec left;  // points in the left set
            boolvec right; // points in the right set
            boolvec edge;  // edges in the boundary between

            DualEdge(
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

            inline const bool same_as(const DualEdge &de) const
            {
                return left == de.left && right == de.right;
            }

            DualEdge(){};
        };

        std::vector<boolvec> nodes;
        std::vector<DualEdge> duals;

    public:
        DualGraph(){};

        DualGraph(std::vector<boolvec> &nodes,
                  std::vector<DualEdge> &duals) : nodes(nodes), duals(duals){};

        inline const bool isSimple() const
        {
            return duals.size() == 1;
        }

        inline const boolvec left() const
        {
            return duals[0].left;
        }

        inline const boolvec right() const
        {
            return duals[0].right;
        }

        inline const boolvec edge() const
        {
            return duals[0].edge;
        }

        std::vector<DualGraph> simplify() const
        {
            std::vector<boolvec> nodes(this->nodes.size() - 1);
            std::vector<DualEdge> duals(this->duals.size() - 1);
            std::vector<DualGraph> graphs;

            for (int i = 0; i < this->duals.size(); i++)
            {
                nodes.clear();
                duals.clear();
                // construct a simpler graph by removing the indicated edge
                // since left and right are disjoint this is pretty easy for the nodes
                DualEdge remove = this->duals[i];
                for (auto np : this->nodes)
                    if (np != remove.left && np != remove.right) // skip the two we're merging
                        nodes.push_back(np);
                boolvec new_node = remove.left | remove.right; // add the merged result
                nodes.push_back(new_node);

                // for the edges we have to see if two or more edges collapsed into one
                for (auto dp = this->duals.begin(); dp < this->duals.end(); ++dp)
                {
                    int i;
                    DualEdge de = *dp;

                    // this is a rapid test for the interesting cases
                    if (de.left & new_node)
                    {
                        DualEdge nde{de.right, new_node, de.edge};
                        for (i = 0; i < duals.size(); ++i)
                        {
                            if (nde.same_as(de))
                                duals[i].edge |= nde.edge; // found it OR it in
                            break;
                        };
                        if (i == duals.size())
                            duals.push_back(nde); // new edge
                    }
                    else if (de.right & new_node)
                    {
                        DualEdge nde{de.left, new_node, de.edge};
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
                    { // nothing to see here copy it forward
                        duals.push_back(de);
                    }
                }
                graphs.push_back(DualGraph(nodes, duals));
            }
            return graphs;
        };
    };
}