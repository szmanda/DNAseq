using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DNAseq
{
    public class Utils
    {
        private static string projectDirectory = Directory.GetParent(Environment.CurrentDirectory).Parent.Parent.FullName;

        public static List<Instance> LoadInstances()
        {
            List<Instance> instances = new List<Instance>();
            string path = Path.Combine(projectDirectory, "instances");
            Console.WriteLine($"Loading instances from {path}");
            foreach (string entry in Directory.EnumerateFiles(path))
            {
                Instance instance = new Instance(entry);
                int counter = 0;
                using (StreamReader file = new StreamReader(entry))
                {
                    string line;
                    while ((line = file.ReadLine()) != null)
                    {
                        Oligonucleotide olig = new Oligonucleotide(counter++, line.Trim());
                        if (olig.text.Length == instance.l)
                            instance.oligonucleotides.Add(olig);
                    }
                }

                instances.Add(instance);
            }
            return instances;
        }

        public static List<Instance> LoadMockInstances()
        {
            List<Instance> instances = new List<Instance>();
            Instance instance = new Instance();
            instance.n = 5;
            instance.s = 4;
            instance.l = 3;
            // Mock instance: ACTCTGG, S = {ACT, CTC, TCT, TGG}
            // string[] oligsString = { "ACTACTACT", "TTACTACTA", "TCTTCTTCT", "TGGTGGTGG" };
            string[] oligsString = { "ACT", "CTC", "TCT", "TGG" };
            for (int i = 0; i < oligsString.Length; i++)
            {
                instance.oligonucleotides.Add(new Oligonucleotide(i, oligsString[i]));
            }
            instances.Add(instance);
            return instances;
        }
    }
}
