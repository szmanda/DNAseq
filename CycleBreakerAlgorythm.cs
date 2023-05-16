using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DNAseq
{
    internal class CycleBreakerAlgorythm {
        public Instance instance;

        public CycleBreakerAlgorythm(Instance instance) {
            this.instance = instance;
        }
        public void run() {
            int count = instance.oligonucleotides.Count;
            List<Oligonucleotide> oligs = instance.oligonucleotides;
            foreach (Oligonucleotide olig in instance.oligonucleotides) {
                List<Arc> arcs = generateSuccessors(olig, 5);
                List<Arc> ordered = arcs.OrderBy(a => a.cost).ToList();
                Console.WriteLine($"\tSuccessors of {olig.text}");
                foreach (Arc arc in ordered)
                {
                    Console.WriteLine($"{oligs[arc.from].text} --{arc.cost}--> {oligs[arc.to].text}");
                }
            }
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
    }
}
