using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DNAseq
{
    public class Instance
    {
        public enum ErrorType
        {
            NONE,
            NEGATIVE_RANDOM,
            NEGATIVE_REPEAT,
            POSITIVE_RANDOM,
            POSITIVE_WRONG_ENDING
        }

        public string filepath;
        public ErrorType errorType = ErrorType.NONE;
        public int numErrors = 0;
        public int n = 0; // dna sequence length
        public int s = 0; // number of oligonucleotides
        public int l = 0; // oligonucleotide length
        public string name = string.Empty;
        public List<Oligonucleotide> oligonucleotides = new List<Oligonucleotide>();

        public Instance(string filepath)
        {
            this.filepath = filepath;
            ExtractInstanceInfo();
        }

        public Instance()
        {
        }

        private void ExtractInstanceInfo()
        {
            name = Path.GetFileName(filepath);
            l = 10;

            int dotPos = name.IndexOf('.');
            int plusPos = name.IndexOf('+');
            int minusPos = name.IndexOf('-');
            int numErrorsOffset = 0;
            if (plusPos != -1)
                numErrorsOffset = plusPos;
            else if (minusPos != -1)
                numErrorsOffset = minusPos;
            else
                throw new InvalidOperationException("Unexpected filename");

            s = int.Parse(name.Substring(dotPos + 1, numErrorsOffset - dotPos - 1));
            n = s + l - 1;
            numErrors = int.Parse(name.Substring(numErrorsOffset + 1, name.Length - numErrorsOffset - 1));

            if (minusPos != -1 && numErrors >= 40)
                errorType = ErrorType.NEGATIVE_RANDOM;
            else if (minusPos != -1)
                errorType = ErrorType.NEGATIVE_REPEAT;
            else if (plusPos != -1 && numErrors >= 80)
                errorType = ErrorType.POSITIVE_RANDOM;
            else if (plusPos != -1)
                errorType = ErrorType.POSITIVE_WRONG_ENDING;
        }

        public string toString()
        {
            return $"Instance {name}: {n} {s} {l} : {oligonucleotides.Count} olinonucleotides, first: {oligonucleotides.First().text}";
        }
    }
}
