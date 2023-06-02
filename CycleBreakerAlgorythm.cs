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
        public int solutionOligCount = 0;

        public CycleBreakerAlgorythm(Instance instance) {
            int s = instance.oligonucleotides.Count;
            this.instance = instance;
            this.cycles = new List<int>();
            this.offsets = new int[s,s];
        }
        public String run() {
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
                generatePathCosts(v);
            }
            foreach (Vertex v in vertices)
            {
                // Console.WriteLine($"== {v.olig.text} ==");
                (v.pathLength, v.cycleId) = findPath(v, new List<int>());
                if (v.cycleId != null && !cycles.Contains(v.cycleId ?? 0))
                {
                    cycles.Add(v.cycleId ?? 0);
                }
                // Console.WriteLine($"path length: {v.pathLength}, cycle: {v.cycleId}");
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

            /// testing input data
            /// 
            Vertex longest = vertices.Where(a => a.incoming == 0).OrderByDescending(a => a.pathLength).First();
            Console.WriteLine($"oligCount: {longest.pathLength}, cost: {10 + longest.pathCost}");
            // printSuccessors(longest.id);

            Console.WriteLine($"============starting-with-{cycles.Count}-cycles============");

            /// Breaking the cycle and joining in to anoter vertex
            /// Simplest version joining only to vertices:
            /// - without predacessors,
            /// - leading to a different cycle
            /// - not checking successors for a better match.
            int iterationsLimiter = 3;
            while (cycles.Count >= 1 && iterationsLimiter > 0)
            {
                int cycleId = cycles[0];
                int cycleLength = findPathLength(vertices[cycleId], new List<int>());
                if (cycleLength == 0) {
                    Console.WriteLine($"ERROR: Cycle {cycleId} has no length!!!");
                    printSuccessors(cycleId);
                    Console.WriteLine("Removing this cycle!");
                    cycles.RemoveAt(0);
                    Console.ReadKey();
                    continue;
                }

                Console.WriteLine($"\t Breaking the cycle {cycleId} of length {cycleLength}");
                printSuccessors(cycleId);
                List<Vertex> potentialBreakpoints = findVerticesInCycle(cycleId)
                    .OrderByDescending(v => v.arcs[0].cost).Take(5).ToList();

                List<Vertex> startVerticesInThisCycle = vertices.Where(a => a.incoming == 0 && a.cycleId == cycleId).OrderByDescending(a => a.pathLength).Take(1).ToList();
                List<Vertex> potentialJoin = vertices.Where(a => a.incoming == 0 && a.cycleId != cycleId).ToList();
                List<Vertex> potentialSelfJoin = vertices.Where(a => a.incoming == 0 && a.cycleId == cycleId).ToList();
                List<Arc> newArcs = new List<Arc>();
                List<Arc> newSelfArcs = new List<Arc>();
                foreach (Vertex br in potentialBreakpoints)
                {
                    // Console.WriteLine($"? Breaking {br.id}: {br.olig.text} --{br.arcs[0].cost}--> {vertices[br.arcs[0].to].olig.text} ?");
                    // Console.WriteLine($"? cycle: {br.cycleId} cycleOffset: {br.cycleOffset}");
                    int bestBr = 0;
                    int bestNew = 0;
                    int costNew = 100;

                    // Evaluating the connection againts potential breakpoints to connect to other subgraph
                    // Finding the best arch from this breakpoint
                    foreach (Vertex vertex in potentialJoin)
                    {
                        Vertex v = vertex;
                        for (int i = 0; i < 3; i++) 
                        {
                            int? offset = calculatePositiveNegativeOffset(br.olig, v.olig).Item1;
                            if (offset != null && offset > 0)
                            {
                                // Console.WriteLine($"    ?join : {v.olig.text} cycle: {v.cycleId} length: {v.pathLength} eval: [ connection_cost: {offset}, length: +{v.pathLength} ]");
                                if (offset < costNew)
                                {
                                    bestBr = br.id;
                                    bestNew = v.id;
                                    costNew = offset ?? 100;
                                    if (costNew == 100) throw new Exception("offset is null");
                                }
                            }
                            v = vertices[v.arcs[0].to]; // next
                        }
                    }
                    newArcs.Add(new Arc(bestBr, bestNew, costNew, costNew));
                    // Evaluating self joins
                    // Finding the best arch, leading to the same cycle
                    costNew = 100;
                    foreach (Vertex vertex in potentialSelfJoin)
                    {
                        Vertex v = vertex;
                        for (int i = 0; i < 3; i++)
                        {
                            int? offset = calculateOffset(br.olig, v.olig);
                            if (offset != null && offset > 0)
                            {
                                // Console.WriteLine($"    ?join : {v.olig.text} cycle: {v.cycleId} length: {v.pathLength} eval: [ connection_cost: {offset}, length: +{v.pathLength} ]");
                                if (offset < costNew)
                                {
                                    bestBr = br.id;
                                    bestNew = v.id;
                                    costNew = offset ?? 100;
                                    if (costNew == 100) throw new Exception("offset is null");
                                }
                            }
                            v = vertices[v.arcs[0].to]; // next
                        }
                    }
                    newSelfArcs.Add(new Arc(bestBr, bestNew, costNew, costNew));
                }
                Arc bestArc = new Arc(0,1,100,100);
                float bestRatio = -1000;
                newArcs = newArcs.OrderBy(a => a.cost).Take(5).ToList();
                Console.WriteLine("5 best new arcs:");
                foreach (Arc arc in newArcs)
                {
                    Vertex from = vertices[arc.from];
                    Vertex to = vertices[arc.to];
                    int extensionLength = to.pathLength;
                    int extensionCost = to.pathCost;
                    Console.WriteLine($"{from.olig.text} --{arc.cost}--> {to.olig.text}");
                    Console.WriteLine($"    This would mean extension with a fragment of length {extensionLength}, and cost: {extensionCost}");
                    foreach (Vertex start in startVerticesInThisCycle)
                    {
                        int o = findCycleOffset(start);
                        int cutOffLength = ((-1 + o - from.cycleOffset) % cycleLength);
                        int lengthDiff = extensionLength - cutOffLength;
                        int lengthNew = start.pathLength + lengthDiff;
                        float costDiff = (float)extensionCost; // TODO: subtract the cost of the part after cut-off point
                        float costNew = start.pathCost + extensionCost;
                        float lengthToCostRatio = (costNew != 0) ? lengthNew * lengthNew / costNew : -1000; // using squared length to promote longer solutions
                        if (lengthToCostRatio > bestRatio && arc.from != arc.to)
                        {
                            bestRatio = lengthToCostRatio;
                            bestArc = arc;
                        }
                        Console.WriteLine($"\tBy joining to {start.olig.text} (with length: {start.pathLength}) with offset {o}:");
                        Console.WriteLine($"\twould increase it's length by {lengthDiff} to {lengthNew} adding cost: {costDiff}");
                        Console.WriteLine($"\tlength/cost: {lengthToCostRatio}, lengthDiff/costDiff: {lengthDiff / costDiff}, sq(length)/cost: {lengthNew * lengthNew / costNew}");
                    }
                }
                cycleId = vertices[bestArc.from].cycleId ?? cycleId;
                // TODO: Function calculating length and cost *until* a given vertex (null if loop and not found)
                // TODO: Add Self-join!!!
                // TODO: REwrite this stuff, clearly separating diffs.
                Arc bestSelfArc = new Arc(0, 1, 100, 100);
                float bestSelfRatio = -1000;
                newSelfArcs = newSelfArcs.OrderBy(a => a.cost).Take(5).ToList();
                Console.WriteLine("5 best new SELF arcs:");
                foreach (Arc arc in newSelfArcs)
                {
                    Vertex from = vertices[arc.from];
                    Vertex to = vertices[arc.to];
                    if (from.cycleId != to.cycleId) Console.WriteLine($"Self join cycles are different: {from.cycleId} != {to.cycleId}.");
                    int extensionLength = to.pathLength - cycleLength; // TODO: add an amount depending on break offset 
                        int end = to.id; // for calculating cost without redundancy (it'd probably be clearer with a function)
                        for (int i = 0; i < extensionLength - 1; i++) { end = vertices[end].arcs[0].to; }
                    int extensionCost = to.pathCost - vertices[end].pathCost; // cost without the inside of a cycle
                    Console.WriteLine($"{from.olig.text} --{arc.cost}--> {to.olig.text} (SELF)");
                    Console.WriteLine($"    This would mean increasing the length by {extensionLength} at offset {from.cycleOffset}");
                    foreach (Vertex start in startVerticesInThisCycle)
                    {
                        int o = findCycleOffset(start);
                        int cutOffLength = ((-1 + o - from.cycleOffset) % cycleLength);
                        int lengthDiff = extensionLength - cutOffLength;
                        int lengthNew = start.pathLength + lengthDiff;
                        float costDiff = (float)extensionCost; // TODO: subtract the cost of the part after cut-off point
                        float costNew = start.pathCost + extensionCost;
                        float lengthToCostRatio = (costNew != 0) ? lengthNew * lengthNew / costNew : -1000; // using squared length to promote longer solutions
                        if (lengthToCostRatio > bestSelfRatio && arc.from != arc.to)
                        {
                            bestSelfRatio = lengthToCostRatio;
                            bestSelfArc = arc;
                        }
                        Console.WriteLine($"\t{start.olig.text} (with length: {start.pathLength}) with offset {start.cycleOffset} {o}:");
                        Console.WriteLine($"\twould increase it's length by {lengthDiff} to {lengthNew} adding cost: {costDiff}");
                        Console.WriteLine($"\tlength/cost: {lengthToCostRatio}, lengthDiff/costDiff: {lengthDiff / costDiff}, sq(length)/cost: {lengthNew*lengthNew / costNew}");
                    }
                }

                bool doSelfJoin = false;
                if (cycles.Count == 1 || bestSelfRatio > bestRatio) { doSelfJoin = true; }
                int newCycleId = vertices[bestArc.to].cycleId ?? -1;
                Console.ForegroundColor = ConsoleColor.Yellow;
                if (!doSelfJoin) Console.WriteLine($"Decision made - joining {vertices[bestArc.from].olig.text} --{bestArc.cost}--> {vertices[bestArc.to].olig.text}");
                else             Console.WriteLine($"Decision made - self joining {vertices[bestSelfArc.from].olig.text} --{bestSelfArc.cost}--> {vertices[bestSelfArc.to].olig.text} (SELF)");
                Console.ForegroundColor = ConsoleColor.White;
                Console.WriteLine($"Before join, cycle {newCycleId} of length {findPathLength(vertices[newCycleId], new List<int>())} has starting points of such lengths:");
                foreach (Vertex v in vertices) {
                    if (v.cycleId == newCycleId && v.incoming == 0)
                        Console.Write($"{v.pathLength} ");
                }
                Console.WriteLine();

                ////// investigating offset:
                /*List<Vertex> cycleVertecies = findVerticesInCycle(cycleId);
                foreach(Vertex v in vertices)
                {
                    if (v.cycleId == cycleId)
                        Console.WriteLine($"{v.olig.text} offset: {v.cycleOffset}, findOffset: {findCycleOffset(v)}");
                }*/

                if (doSelfJoin)
                {
                    //////////// SELF JOIN
                    Console.WriteLine();
                    List<Vertex> cycle = findVerticesInCycle(cycleId);
                    startVerticesInThisCycle = vertices.Where(a => a.incoming == 0 && a.cycleId == cycleId).OrderByDescending(a => a.pathLength).ToList();

                    int outCycleOffset = vertices[bestSelfArc.from].cycleOffset;
                    int inCycleOffset = findCycleOffset(vertices[bestSelfArc.to]);
                    int outCycleOffsetRelative = cycleLength;
                    int inCycleOffsetRelative = ((-1 + inCycleOffset - outCycleOffset) % cycleLength);
                    Console.WriteLine($"Old offsets (this will all be part of a cycle); from: {outCycleOffset} (rel:{outCycleOffsetRelative}) to: {inCycleOffset}(rel:{inCycleOffsetRelative}) (total cycle length: {cycleLength})");

                    int newCycleLength = cycleLength - inCycleOffsetRelative + vertices[bestSelfArc.to].pathLength - cycleLength;
                    Console.WriteLine($"New cycle will have a length of {newCycleLength}");
                    if (newCycleLength <= cycleLength) {
                        Console.WriteLine("WARNING: after self join, cycle length will not increase. Unpredictable results possible! Breaking algorythm loop");
                        break;
                    }

                    List<int> visited = new List<int>();
                    // fetchnig old offsets for starting vertices
                    foreach (Vertex start in startVerticesInThisCycle)
                    {
                        start.cycleOffset = findCycleOffset(start);
                    }

                    // joining cycles:
                    Console.WriteLine($"{bestSelfArc.from} --> {bestSelfArc.to}");
                    vertices[bestSelfArc.from].arcs[0] = bestSelfArc;

                    newCycleId = findVerticesInCycle(bestSelfArc.from).OrderBy(v => v.id).First().id;
                    int newCycleCost = bestSelfArc.cost;
                    int tmpId = bestSelfArc.to;
                    for (int i = 0; i < newCycleLength; i++)
                    {
                        newCycleCost += vertices[tmpId].arcs[0].cost;
                        tmpId = vertices[tmpId].arcs[0].to; // next
                        Vertex v = vertices[tmpId];
                    }
                    int counter = 0;
                    Console.WriteLine($"Old cycle: [id: {cycleId}, length:{cycleLength}]");
                    Console.WriteLine($"New cycle: [id: {newCycleId}, length:{newCycleLength}, cost: {newCycleCost}]");
                    Vertex tmpIn = vertices[cycleId];
                    for (int i = 0; i < vertices[bestSelfArc.to].cycleOffset; i++)
                    {
                        tmpIn = vertices[tmpIn.arcs[0].to];
                    }
                    Vertex oldIn = new Vertex(tmpIn.id, tmpIn.olig, tmpIn.arcs);
                    oldIn.incoming = tmpIn.incoming;
                    oldIn.pathLength = tmpIn.pathLength;
                    oldIn.arcs = tmpIn.arcs;
                    oldIn.pathCost = tmpIn.pathCost;
                    Console.WriteLine($"Old IN {oldIn.id} {oldIn.olig.text} {oldIn.cycleOffset}    cost:{oldIn.arcs[0].cost}");

                    foreach (Vertex v in findVerticesInCycle(newCycleId))
                    {
                        v.cycleOffset = counter++;
                        v.cycleId = newCycleId;
                        v.pathLength = newCycleLength;
                        v.pathCost = newCycleCost;

                        visited.Add(v.id);
                        // Display new cycle:
                        // Console.WriteLine($"{v.id}\t{v.olig.text} {v.cycleOffset}\t cost:{v.arcs[0].cost}");
                    }

                    foreach (Vertex start in startVerticesInThisCycle)
                    {
                        int id = start.id;
                        Console.Write($"{start.id} ");
                        int oldOffset = start.cycleOffset;
                        int offsetRelative = ((-1 + oldOffset - outCycleOffset + cycleLength) % cycleLength);
                        int inOffsetRelative = ((-1 + oldIn.cycleOffset - outCycleOffset + cycleLength) % cycleLength);

                        int newOffset = findCycleOffset(start); // only works for joiners with offset >= inOffset
                        int lenDiff = start.pathLength - cycleLength - inOffsetRelative; // only works for joiners with offset > inOffset
                        Console.WriteLine($"Changing vertecies starting with {start.olig.text} (length: {start.pathLength}) in cycle {start.cycleId} with offset {oldOffset} (rel:{offsetRelative})");
                        for (int i = 0; i < start.pathLength - cycleLength; i++)
                        {
                            if (visited.Contains(id)) break;
                            visited.Add(id);
                            vertices[id].cycleId = newCycleId;
                            if (offsetRelative < inCycleOffsetRelative)
                            {
                                int pathLengthDiffUntilNewCycle = inCycleOffsetRelative - offsetRelative;
                                Console.WriteLine($"(!!! NOT IMPLEMENTED YET !!!) Path from old entrance to the cycle to the new entrance is {pathLengthDiffUntilNewCycle}");
                                // todo: recalculate new offset
                            }
                            if (offsetRelative > inCycleOffsetRelative)
                            {
                                vertices[id].cycleOffset = newOffset;
                                vertices[id].pathLength += lenDiff;
                                //vertices[id].pathCost += costDifference;
                                vertices[id].cycleId = newCycleId;
                            }
                            if (offsetRelative == inCycleOffsetRelative)
                            {
                                vertices[id].cycleOffset = newOffset;
                                vertices[id].pathLength += lenDiff - vertices[id].pathLength + cycleLength;
                                //vertices[id].pathCost += costDifference;
                                vertices[id].cycleId = newCycleId;
                            }
                            // ! Commented out ! be careful not to override things such as cycle offset, since it is used in further iterations,
                            // ! i might need to calculate it for each vertex before looping over it.
                            // vertices[id].cycleOffset = 0;// newOffset;
                            // vertices[id].pathLength += 0;// lenDiff;
                            // vertices[id].pathCost += 0;// costDifference;
                            // vertices[id].cycleId = newCycleId;
                            id = vertices[id].arcs[0].to; // next
                        }
                    }



                    cycles.Remove(cycleId);
                    cycles.Add(newCycleId);
                    Console.WriteLine("---\n\n");


                    //////////// END OF SELF JOIN
                }
                else if (!doSelfJoin)
                {
                    //// Actually break the cycle and join to another using `bestArc`
                    // the part we are joining to does not need adjustment, only the cycle currently being broken.
                    List<Vertex> cycle = findVerticesInCycle(cycleId);
                    startVerticesInThisCycle = vertices.Where(a => a.incoming == 0 && a.cycleId == cycleId).OrderByDescending(a => a.pathLength).ToList();

                    // int newCycleId = 0; //! find cycle id
                    int newOffset = findCycleOffset(vertices[bestArc.to]);
                    int costDifference = vertices[bestArc.to].pathCost; // for now - ignoring old cycle offset in this calculation
                    List<int> visited = new List<int>(); // pewnie można by to zrobić lepiej niż lista, ale taki quick fix żeby wielokroenie nie modyfikować
                    vertices[vertices[bestArc.to].arcs[0].to].incoming -= 1; // stop pointing arc to this vertex
                    foreach (Vertex v in findVerticesInCycle(bestArc.from))
                    {
                        int cutOffLength = (-1 + newOffset - v.cycleOffset) % cycleLength;
                        vertices[v.id].cycleOffset = newOffset;
                        vertices[v.id].pathLength += vertices[bestArc.to].pathLength - cutOffLength;
                        visited.Add(v.id);
                    }
                    foreach (Vertex start in startVerticesInThisCycle)
                    {
                        int id = start.id;
                        int oldOffset = findCycleOffset(start);
                        int lenDiff = -((-1 + oldOffset - vertices[bestArc.from].cycleOffset) % cycleLength) + vertices[bestArc.to].pathLength;
                        for (int i = 0; i < start.pathLength - cycleLength; i++)
                        {
                            if (visited.Contains(id)) break;
                            visited.Add(id);
                            vertices[id].cycleOffset = newOffset;
                            vertices[id].pathLength += lenDiff;
                            vertices[id].pathCost += costDifference;
                            vertices[id].cycleId = newCycleId;
                            id = vertices[id].arcs[0].to; // next
                        }
                    }

                    // actually joining graphs:
                    vertices[bestArc.from].arcs[0] = bestArc;
                    cycles.Remove(cycleId);

                    Console.WriteLine($"After join, cycle {newCycleId} of length {findPathLength(vertices[newCycleId], new List<int>())} has starting points of such lengths:");
                    foreach (Vertex v in vertices)
                    {
                        if (v.cycleId == newCycleId && v.incoming == 0)
                            Console.Write($"{v.pathLength} ");
                    }
                    Console.WriteLine();
                }
                Console.WriteLine($"========================{cycles.Count}-cycles-remaining========================");
                iterationsLimiter--;
                //break; // Breaking the loop (Testing first iteration)
            }

            Console.WriteLine($"Given instance {instance.name}");
            Console.WriteLine("Best result:");
            Vertex starting = vertices.Where(a => a.incoming == 0).OrderByDescending(a => a.pathLength).First();
            Console.WriteLine($"oligCount: {starting.pathLength}, length: {10 + starting.pathCost}");
            printSuccessors(starting.id);
            
            
            //////////// Formatting the output
            List<int> tmpVisited = new List<int>();
            int tmpIdx = starting.id;
            int oligLength = 10;
            int tmpOffset = 10;
            int totalOffset = 0;
            StringBuilder sequenceStringBuilder = new StringBuilder();
            while (!tmpVisited.Contains(tmpIdx))
            {
                tmpVisited.Add(tmpIdx);
                if (oligLength - tmpOffset < 0) { Console.WriteLine($"ERROR: Unexpected offset {oligLength - tmpOffset}! after {tmpVisited.Count} elements"); break; }
                string a = vertices[tmpIdx].olig.text.Substring(oligLength - tmpOffset);
                Console.WriteLine($"{tmpIdx}\t{a}");
                sequenceStringBuilder.Append(a);
                tmpOffset = calculateOffset(vertices[tmpIdx].olig, vertices[vertices[tmpIdx].arcs[0].to].olig) ?? 10;//vertices[tmpIdx].arcs[0].offset; // sometimes errors?
                if (tmpOffset < 0) tmpOffset = 10;
                tmpIdx = vertices[tmpIdx].arcs[0].to;
                /*for (int i = 0; i < totalOffset; i++) Console.Write(" ");
                Console.WriteLine($"{vertices[tmpIdx].olig.text}");*/
                totalOffset += tmpOffset;
                solutionOligCount++;
            }
            String result = sequenceStringBuilder.ToString();
            Console.WriteLine($"Finished instance {instance.name}.");
            return result;
        }

        public void printSuccessors(int id)
        {
            Console.WriteLine($"Printing successors of {id}: {instance.oligonucleotides[id].text}");
            List<int> visited = new List<int>();
            int i = 0;
            int totalCost = 10;
            while (!visited.Contains(id))
            {
                int offset = calculateOffset(instance.oligonucleotides[id], instance.oligonucleotides[vertices[id].arcs[0].to]) ?? 100;
                Console.WriteLine($"{i}.\t{instance.oligonucleotides[id].text} {id}\t {offset}");
                visited.Add(id);
                id = vertices[id].arcs[0].to;
                totalCost += offset;
                i++;
            }
            Console.WriteLine($"Included oligonucleotides: {i}, Total cost: {totalCost}");
        }

        // search for closest vertex with .cycleOffset set 
        public int findCycleOffset(Vertex v)
        {
            if (v.cycleOffset != -1) return v.cycleOffset;
            return findCycleOffset(vertices[v.arcs[0].to]);
        }

        // it produces correct results only for vertecies outside cycles.
        public int generatePathCosts(Vertex v)
        {
            // if (v.pathCost > 0) return v.pathCost;
            if (v.pathCost == -1) return 0;
            v.pathCost = -1; // mark as visited
            int cost = generatePathCosts(vertices[v.arcs[0].to]) + calculateOffset(vertices[v.arcs[0].from].olig, vertices[v.arcs[0].to].olig) ?? 1000;
            v.pathCost = cost;
            return cost;
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
