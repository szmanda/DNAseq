using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DNAseq
{
    internal class AlgorythmSTSP
    {

        /*private List<Arc> CalculateMinimumSpanningTree(Graph graph)
        {
            List<Arc> minimumSpanningTree = new List<Arc>();
            // Sort the Arcs in ascending order of their costs
            List<Arc> sortedArches = graph.Arches.OrderBy(e => e.Cost).ToList();
            // Create a disjoint set for the vertices
            DisjointSet<Oligonucleotide> disjointSet = new DisjointSet<Oligonucleotide>(graph.Vertices);
            foreach (Arc arc in sortedArches)
            {
                if (!disjointSet.AreInSameSet(arc.StartVertex, arc.EndVertex))
                {
                    disjointSet.Union(arc.StartVertex, arc.EndVertex);
                    minimumSpanningTree.Add(arc);
                }
            }
            return minimumSpanningTree;
        }*/
    }
}
