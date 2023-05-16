using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DNAseq
{
    public class Oligonucleotide
    {
        public int id { get; set; }
        public string text { get; set; }

        public Oligonucleotide(int id = -1, string text = "")
        {
            this.id = id;
            this.text = text;
        }
    }

    /*public class Arc
    {
        public Oligonucleotide StartVertex { get; set; }
        public Oligonucleotide EndVertex { get; set; }
        public int Cost { get; set; }
    }*/

    public class Graph
    {
        public List<Oligonucleotide> Vertices { get; set; }
        public List<Arc> Arches { get; set; }
    }

    /// <summary>
    /// For CycleBreaker algorythm
    /// </summary>
    public class Vertex
    {
        public int id { get; set; }
        public Oligonucleotide olig { get; set; }
        public List<Arc> arcs { get; set; }
    }

    public class Arc
    {
        public Arc(int from = 0, int to = 0, int cost = 0, int offset = 0) {
            this.from = from;
            this.to = to;
            this.cost = cost;
            this.offset = offset;
        }
        public int from { get; set; }
        public int to { get; set; }
        public int cost { get; set; }
        public int offset { get; set; }
    }

}
