using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Channels;
using System.Threading.Tasks;

namespace DNAseq
{
    internal class CycleBreakerAlgorythm {
        public Instance instance;
        public List<Vertex> vertices;
        public List<int> cycles;
        public int[,] offsets;

        public CycleBreakerAlgorythm(Instance instance) {
            int s = instance.oligonucleotides.Count;
            this.instance = instance;
            this.cycles = new List<int>();
            this.offsets = new int[s,s];
        }
        public void run() {
            generateOffsets(instance.oligonucleotides);
            int count = instance.oligonucleotides.Count;
            List<Oligonucleotide> oligs = instance.oligonucleotides;
            vertices = new List<Vertex>();

            foreach (Oligonucleotide olig in instance.oligonucleotides) {
                List<Arc> arcs = generateSuccessors(olig, 100);
                List<Arc> ordered = arcs.OrderBy(a => a.cost).ToList();
                vertices.Add(new Vertex(olig.id, olig, ordered));
                /*Console.WriteLine($"\tSuccessors of {olig.text}");
                foreach (Arc arc in ordered)
                {
                    Console.WriteLine($"{oligs[arc.from].text} --{arc.cost}--> {oligs[arc.to].text}");
                }*/
            }
            foreach (Vertex v in vertices)
            {
                Console.WriteLine($"== {v.olig.text} ==");
                (v.pathLength, v.cycleId) = findPath(v, new List<int>());
                if (v.cycleId != null && !cycles.Contains(v.cycleId ?? 0))
                {
                    cycles.Add(v.cycleId ?? 0);
                }
                Console.WriteLine($"path length: {v.pathLength}, cycle: {v.cycleId}");
                // mark non heads, thi important for merges later
                vertices[v.arcs[0].to].incoming += 1;
            }
            Console.WriteLine(vertices[0].pathLength);
            int withLongestPath = vertices.OrderByDescending(a => a.pathLength).First().id;
            //int cycleId = vertices[withLongestPath].cycleId ?? 0;

            // simple search for one of the maximum costs.
            // for now ignoring other potential vertices for breaking the cycle
            /*int current = vertices[cycleId].arcs[0].to;
            int withLargestCost = cycleId;
            while (true)
            {
                Console.WriteLine($"{current} \t {vertices[current].olig.text} \t {vertices[current].arcs[0].cost}");
                current = vertices[current].arcs[0].to;
                if (current == cycleId) { break; }
                if (vertices[withLargestCost].arcs[0].cost < vertices[current].arcs[0].cost)
                {
                    withLargestCost = vertices[current].arcs[0].to;
                }
            }
            Console.WriteLine($"with largest cost {withLargestCost}");
            Console.WriteLine($"\tSuccessors: ");
            foreach (Arc arc in vertices[withLargestCost].arcs)
            {
                Console.WriteLine($"{oligs[arc.from].text} --{arc.cost}--> {oligs[arc.to].text}");
            }*/

            foreach (int c in cycles)
            {
                // set vertices' cycleOffset property
                int counter = 0;
                foreach (Vertex v in findVerticesInCycle(c))
                {
                    v.cycleOffset = counter++;
                }

                // potential cycle breakpoints
                /*Console.WriteLine($"\tcycle {c}");
                List<Vertex> potentialBreakpoints = findVerticesInCycle(c)
                    .OrderByDescending(v => v.arcs[0].cost).Take(5).ToList();
                foreach (Vertex v in potentialBreakpoints)
                {
                    Console.WriteLine($"{v.olig.text}, current cost {v.arcs[0].cost}");
                    foreach (Arc a in v.arcs.Where(a => vertices[a.to].incoming == 0 && vertices[a.to].cycleId != c).Take(5))
                    {
                        Console.WriteLine($"{v.olig.text} --{a.cost}--> {vertices[a.to].olig.text} c: {vertices[a.to].cycleId} len: {vertices[a.to].pathLength}");
                    }
                }*/
            }

            /// Breaking the cycle and joining in to anoter vertex
            /// Simplest version joining only to vertices:
            /// - without predacessors,
            /// - leading to a different cycle
            /// - not checking successors for a better match.
            while (cycles.Count > 0)
            {
                int cycleId = cycles[0];
                Console.WriteLine($"\t Breaking the cycle {cycleId}");
                List<Vertex> potentialBreakpoints = findVerticesInCycle(cycleId)
                    .OrderByDescending(v => v.arcs[0].cost).Take(5).ToList(); 

                List<Vertex> potentialJoin = vertices.Where(a => a.incoming == 0 && a.cycleId != cycleId).ToList();
                foreach (Vertex vertex in potentialJoin)
                {
                    Console.WriteLine($"join? : {vertex.olig.text} cycle: {vertex.cycleId} length: {vertex.pathLength}");
                    // evaluating the connection againts potential preakpoints
                    /*foreach (Vertex br in potentialBreakpoints)
                    {
                        int? offset = calculateOffset(br.olig, vertex.olig);
                        if (offset != null && offset > 0)
                        {
                            Console.WriteLine($"eval: [ connection_cost: {offset}, length: +{vertex.pathLength} ]");
                            
                        }
                    }*/
                    List<Vertex> brpts = potentialBreakpoints.OrderBy(a => offsets[a.id, vertex.id]).Take(5).ToList();
                    foreach (Vertex br in brpts)
                    {
                        Console.WriteLine($"{br.id}\t {br.olig.text} offset: {offsets[br.id, vertex.id]}");
                    }
                }
                cycles.Remove(cycleId);
            }
            Console.WriteLine("Finished.");

        }

        public List<Arc> generateSuccessors(Oligonucleotide olig, int maxCost = 100) {
            List<Arc> result = new List<Arc>();
            foreach (Oligonucleotide o in instance.oligonucleotides) {
                int? offset = calculateOffset(olig, o);
                if (offset == null) continue;
                else if (offset > 0 && offset < maxCost) {
                    result.Add(new Arc(olig.id, o.id, (int)offset, (int)offset));
                }
            }
            return result;
        }

        public int findPathLength(Vertex v, List<int> visited)
        {
            visited.Add(v.id);
            if (v.arcs.Count == 0) {
                Console.WriteLine("Path ended without vertices");
                return 1;
            }
            int next = v.arcs[0].to;
            if (visited.Contains(next))
            {
                /*Console.WriteLine("Path ended with a cycle");*/
                return 0;
            }
            return findPathLength(vertices[next], visited) + 1;
        }

        // returns: length, cycleId
        public (int, int?) findPath(Vertex v, List<int> visited)
        {
            visited.Add(v.id);
            if (v.arcs.Count == 0)
            {
                Console.WriteLine("Path ended without vertices");
                return (1, null);
            }
            int next = v.arcs[0].to;
            if (visited.Contains(next))
            {
                //Console.WriteLine("Path ended with a cycle");
                return (0, findCycleId(v));
            }
            (int length, int? cycleId) = findPath(vertices[next], visited);
            return (length + 1, cycleId);
        }

        // smallest id in the cycle, null if not a cycle
        public int? findCycleId(Vertex v, List<int>? visited = null)
        {
            if (visited == null) { visited = new List<int>(); }
            if (visited.Contains(v.id)) { return v.id; }
            if (v.arcs.Count == 0)
            {
                Console.WriteLine("err: This is not a cycle!");
                return null;
            }
            visited.Add(v.id);
            int next = v.arcs[0].to;
            int? qqq = findCycleId(vertices[next], visited);
            if (qqq == null) { return null; }
            return Math.Min(qqq ?? 0, v.id);
        }

        public List<Vertex> findVerticesInCycle(int id)
        {
            List<Vertex> visited = new List<Vertex>();
            int i = id;
            while (true)
            {
                if (visited.Contains(vertices[i]))
                {
                    return visited;
                }
                else
                {
                    visited.Add(vertices[i]);
                }
                i = vertices[i].arcs[0].to;
            }
        }

        /// Calculates offset between the oligonucleotides:
        /// 
        /// `(AAAA, AATT) -> 2, (AATT, AAAA) -> -2, (AAAA, TTTT) -> null`
        public int? calculateOffset(Oligonucleotide olig1, Oligonucleotide olig2) {
            int n = olig1.text.Length;
            List<char> o1 = new List<char>(olig1.text);
            List<char> o2 = new List<char>(olig2.text);

            for (int i = 0; i < n; i++) {
                bool matchPositive = true, matchNegative = true;
                for (int j = 0; j < n - i; j++) {
                    if (o1[i+j] != o2[j]) {
                        matchPositive = false;
                    }
                    if (o2[i+j] != o1[j]) {
                        matchNegative = false;
                    }
                }
                if (matchPositive || matchNegative) {
                    return matchPositive ? i : -i;
                }
            }
            return null;
        }

        // calculates a matrix of all offsets, this will hopefully be a more efficient method of accessing the offset
        public void generateOffsets(List<Oligonucleotide> oligs)
        {
            for (int i = 0; i < oligs.Count; i++)
            {
                for (int j = 0; j < oligs.Count; j++)
                {
                    (int pos, int neg) = calculatePositiveNegativeOffset(oligs[i], oligs[j]);
                    {
                        offsets[j, i] = pos;
                        offsets[i, j] = neg;
                    }
                }
            }
        }

        public (int, int) calculatePositiveNegativeOffset(Oligonucleotide olig1, Oligonucleotide olig2)
        {
            int n = olig1.text.Length;
            List<char> o1 = new List<char>(olig1.text);
            List<char> o2 = new List<char>(olig2.text);
            int pos = n;
            int neg = n;

            for (int i = 0; i < n; i++)
            {
                bool matchPositive = true, matchNegative = true;
                for (int j = 0; j < n - i; j++)
                {
                    if (o1[i + j] != o2[j])
                    {
                        matchPositive = false;
                    }
                    if (o2[i + j] != o1[j])
                    {
                        matchNegative = false;
                    }
                }
                if (matchPositive) pos = i;
                if (matchNegative) neg = i;
                if (pos != n && neg != n) return (pos, neg);
            }
            return (pos, neg);
        }
    }
}
