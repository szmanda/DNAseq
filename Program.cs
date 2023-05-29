// See https://aka.ms/new-console-template for more information
using DNAseq;

Console.WriteLine("Hello, World! Testing Testing");

CycleBreakerAlgorythm alg = new CycleBreakerAlgorythm(new Instance());
Oligonucleotide olig1 = new Oligonucleotide();
olig1.text = "TTTT";
Oligonucleotide olig2 = new Oligonucleotide();
olig2.text = "AAAT";
Console.WriteLine($"Offset is {alg.calculateOffset(olig1, olig2) ?? 99}");

/*Console.WriteLine("Loading mock instance");
List<Instance> mockInstances = Utils.LoadMockInstances();
foreach (Instance instance in mockInstances)
{
    Console.WriteLine(instance.toString());
    alg = new CycleBreakerAlgorythm(instance);
    alg.run();
}
*/

Console.WriteLine("\n\nLoading instances");
List<Instance> instances = Utils.LoadInstances();
//for (int i = 0; i < instances.Count; i++)
//int i = 30; // instance 55.300-120 (many negative errors)
//int i = 25; // instance 53.500+50 (few positive errors)
int i = 24;
{
    Console.WriteLine(instances[i].toString());
    alg = new CycleBreakerAlgorythm(instances[i]);
    alg.run();
    //if (Console.ReadKey().KeyChar == 'q') break;
}