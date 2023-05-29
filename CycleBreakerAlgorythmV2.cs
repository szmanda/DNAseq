using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using System.Text;
using System.Threading.Channels;
using System.Threading.Tasks;

namespace DNAseq
{
    internal class CycleBreakerAlgorythmV2
    {
        public Instance instance;
        public List<int> cycleList;
        public int[,] offsets;
        public int[] successors;
        public int[] lengths;
        public int[] cycles;
        public int s;

        public CycleBreakerAlgorythmV2(Instance instance)
        {
            s = instance.oligonucleotides.Count;
            this.instance = instance;
            cycleList = new List<int>();
            offsets = new int[s, s];
            successors = new int[s];
            lengths = new int[s];
            cycles = new int[s];
        }
        public void run()
        {
            int count = instance.oligonucleotides.Count;
            List<Oligonucleotide> oligs = instance.oligonucleotides;

            generateOffsets(oligs);
            generateSuccessors();
            for (int i = 0; i < count; i++) {
                generateLengths(i);
                // Console.WriteLine($"{instance.oligonucleotides[i].text} --- len: {lengths[i]}, cycle: {cycles[i]} ({calculateCycleLength(cycles[i]).Item1}), suc: {instance.oligonucleotides[successors[i]].text}");
            }

            //printSuccessors(324);
            foreach (int c in cycleList)
            {
                Console.WriteLine($"CYCLE: {c}");
                printSuccessors(c);
                Console.WriteLine("Potential breakpoints");
                for (int i = 0; i < lengths[i]; i++)
                {
                    //offsets[i,successors[i]]
                    i = successors[i];
                }
            }

            


            Console.WriteLine("Finished.");
        }

        public List<int> getCycle(int id)
        {
            List<int> visited = new List<int>();
            while (!visited.Contains(id))
            {
                visited.Add(id);
                id = successors[id];
            }
            return visited;
        }

        public void printSuccessors(int id)
        {
            Console.WriteLine($"Printing successors of {id}: {instance.oligonucleotides[id].text}");
            List<int> visited = new List<int>();
            int i = 1;
            while (!visited.Contains(id))
            {
                Console.WriteLine($"{i++}.\t{instance.oligonucleotides[id].text} {id}");
                visited.Add(id);
                id = successors[id];
            }
        }


        // calculate length and cycle for a vertex with given id (also calculates all successors)
        // (length, cycleId, isCycle)
        public (int, int, bool) generateLengths(int id)
        {
            if (lengths[id] > 0) return (lengths[id], cycles[id], false);
            if (lengths[id] == -1)
            {
                lengths[id] = -2; // mark as last in cycle
                (int cycleLength, int cycleId) = calculateCycleLength(id);
                if (!cycleList.Contains(cycleId)) cycleList.Add(cycleId);  
                return (cycleLength, cycleId, true);
            }
            lengths[id] = -1; // mark as visited
            
            (int l, int c, bool isCycle) = generateLengths(successors[id]);
            bool exitingCycle = lengths[id] == -2;
            lengths[id] = isCycle ? l : l + 1;
            if (exitingCycle) isCycle = false;
            cycles[id] = c;

            return (lengths[id], cycles[id], isCycle);
        }

        // (length, minId)
        public (int, int) calculateCycleLength(int id)
        {
            List<int> visited = new List<int>();
            int length = 0;
            int minId = id;
            while (!visited.Contains(id))
            {
                visited.Add(id);
                length++;
                if (id < minId) minId = id;
                id = successors[id];
            }
            return (length, minId);
        }
        public void generateSuccessors()
        {
            for (int i = 0; i < s; i++)
            {
                int minId = 0;
                int min = 100;
                for (int j = 0; j < s; j++)
                {
                    if (i == j) continue;
                    if (offsets[i, j] < min)
                    {
                        min = offsets[i, j];
                        minId = j;
                    }
                        
                }
                successors[i] = minId;
            }
        }

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


        /// Calculates offset between the oligonucleotides:
        /// 
        /// `(AAAA, AATT) -> 2, (AATT, AAAA) -> -2, (AAAA, TTTT) -> null`
        public int? calculateOffset(Oligonucleotide olig1, Oligonucleotide olig2)
        {
            int n = olig1.text.Length;
            List<char> o1 = new List<char>(olig1.text);
            List<char> o2 = new List<char>(olig2.text);

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
                if (matchPositive || matchNegative)
                {
                    return matchPositive ? i : -i;
                }
            }
            return null;
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
