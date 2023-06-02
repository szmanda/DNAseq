// See https://aka.ms/new-console-template for more information
using DNAseq;
using System.Diagnostics;

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

float[,] times = new float[53, 20];

Console.WriteLine("\n\nLoading instances");
List<Instance> instances = Utils.LoadInstances();
for (int t = 0; t < 1; t++)
{
    for (int i = 0; i < instances.Count; i++)
    //int i = 30; // instance 55.300-120 (many negative errors)
    //int i = 25; // instance 53.500+50 (few positive errors)
    {
        //Console.WriteLine(instances[i].toString());
        alg = new CycleBreakerAlgorythm(instances[i]);

        Console.WriteLine($"Running instance {instances[i].name}");
        System.Console.SetOut(new System.IO.StreamWriter(System.IO.Stream.Null));
        Stopwatch stopWatch = new Stopwatch();
        stopWatch.Start();

        string result = alg.run();

        stopWatch.Stop();
        var standardOutput = new StreamWriter(Console.OpenStandardOutput());
        standardOutput.AutoFlush = true;
        Console.SetOut(standardOutput);
        Console.WriteLine($"result length: {result.Length}, includes {alg.solutionOligCount} oligs from input");

        TimeSpan ts = stopWatch.Elapsed;
        /*string elapsedTime = String.Format("{0:00}.{1:000}", ts.Seconds, ts.Milliseconds);
        Console.WriteLine(elapsedTime, "elapsed time");*/
        times[i, t] = ts.Seconds + ts.Milliseconds / (float)1000;

        Console.WriteLine(result);

        //Utils.SaveToFile("result.txt", result);
        // if (Console.ReadKey().KeyChar == 'q') break;
    }
    /*var standardOutput = new StreamWriter(Console.OpenStandardOutput());
    standardOutput.AutoFlush = true;
    Console.SetOut(standardOutput);
    Console.WriteLine(t);*/

}
for (int i = 0; i < instances.Count; i++)
{
    for (int t = 0; t < 20; t++)
    {
        Console.Write($"{times[i, t]}; ");
    }
    Console.WriteLine();
}