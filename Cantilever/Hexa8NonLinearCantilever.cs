using System;
using System.Collections.Generic;
using System.IO;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public static class Hexa8NonLinearCantilever
    {
        private const int subdomainID = 0;

        [Fact]
        private static void SolveModel()
        {
            var model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

            var load = 850;
            BuildCantileverModel(model, load);

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            // Analyzers
            int increments = 10;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-5;
            childAnalyzerBuilder.MaxIterationsPerIncrement = 100;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Setup displacements log
            var watchDof = new Dictionary<int, int[]>();
            //watchDof.Add(subdomainID, new int[1] { 46 });
            watchDof.Add(subdomainID, new int[1] { 94 });
            var log1 = new IncrementalDisplacementsLog(watchDof);
            childAnalyzer.IncrementalDisplacementsLog = log1;

            // Run the analysis 
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Write log data to file
            WriteDisplacementsToFile(log1, watchDof, load);
        }

        private static void BuildCantileverModel(Model model, double load_value)
        {
            var material1 = new ElasticMaterial3D() { PoissonRatio = 0.3, YoungModulus = 1353000 };

            double[,] nodeData = new double[,] {
            //{-0.250000,-0.250000,-1.000000}, {0.250000,-0.250000,-1.000000}, {-0.250000,0.250000,-1.000000}, {0.250000,0.250000,-1.000000},
            //{-0.250000,-0.250000,-0.500000}, {0.250000,-0.250000,-0.500000}, {-0.250000,0.250000,-0.500000}, {0.250000,0.250000,-0.500000},
            //{-0.250000,-0.250000,0.000000}, {0.250000,-0.250000,0.000000}, {-0.250000,0.250000,0.000000}, {0.250000,0.250000,0.000000},
            //{-0.250000,-0.250000,0.500000}, {0.250000,-0.250000,0.500000}, {-0.250000,0.250000,0.500000}, {0.250000,0.250000,0.500000},
            //{-0.250000,-0.250000,1.000000}, {0.250000,-0.250000,1.000000}, {-0.250000,0.250000,1.000000}, {0.250000,0.250000,1.000000},
            //};

            {-0.250000,-0.250000,-2.000000}, {0.250000,-0.250000,-2.000000}, {-0.250000,0.250000,-2.000000}, {0.250000,0.250000,-2.000000},
            {-0.250000,-0.250000,-1.500000}, {0.250000,-0.250000,-1.500000}, {-0.250000,0.250000,-1.500000}, {0.250000,0.250000,-1.500000},
            {-0.250000,-0.250000,-1.000000}, {0.250000,-0.250000,-1.000000}, {-0.250000,0.250000,-1.000000}, {0.250000,0.250000,-1.000000},
            {-0.250000,-0.250000,-0.500000}, {0.250000,-0.250000,-0.500000}, {-0.250000,0.250000,-0.500000}, {0.250000,0.250000,-0.500000},
            {-0.250000,-0.250000,0.000000}, {0.250000,-0.250000,0.000000}, {-0.250000,0.250000,0.000000}, {0.250000,0.250000,0.000000},
            {-0.250000,-0.250000,0.500000}, {0.250000,-0.250000,0.500000}, {-0.250000,0.250000,0.500000}, {0.250000,0.250000,0.500000},
            {-0.250000,-0.250000,1.000000}, {0.250000,-0.250000,1.000000}, {-0.250000,0.250000,1.000000}, {0.250000,0.250000,1.000000},
            {-0.250000,-0.250000,1.500000}, {0.250000,-0.250000,1.500000}, {-0.250000,0.250000,1.500000}, {0.250000,0.250000,1.500000},
            {-0.250000,-0.250000,2.000000}, {0.250000,-0.250000,2.000000}, {-0.250000,0.250000,2.000000}, {0.250000,0.250000,2.000000},
            };

            int[,] elementData = new int[,] {
            //{1,8,7,5,6,4,3,1,2},
            //{2,12,11,9,10,8,7,5,6},
            //{3,16,15,13,14,12,11,9,10},
            //{4,20,19,17,18,16,15,13,14},
            
            {1,8,7,5,6,4,3,1,2},
            {2,12,11,9,10,8,7,5,6},
            {3,16,15,13,14,12,11,9,10},
            {4,20,19,17,18,16,15,13,14},
            {5,24,23,21,22,20,19,17,18},
            {6,28,27,25,26,24,23,21,22},
            {7,32,31,29,30,28,27,25,26},
            {8,36,35,33,34,32,31,29,30},
            };

            // Put nodes and nodal coordinates in model
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y:  nodeData[nNode, 1], z: nodeData[nNode, 2] ));
            }

            // Define elements
            Element e1;
            int subdomainID = Hexa8NonLinearCantilever.subdomainID;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8                    
                };
                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e1.ID, e1);
            }

            // Constrain nodes at z = -1
            for (int k = 1; k <= 4; k++)
            {
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
            }

            // Top load
            Load load1;
            //for (int k = 17; k <= 20; k++)
            for (int k = 33; k <= 36; k++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[k],
                    DOF = StructuralDof.TranslationY,
                    Amount = 1 * load_value
                };
                model.Loads.Add(load1);
            }
        }

        private static void WriteDisplacementsToFile(IncrementalDisplacementsLog log1, Dictionary<int, int[]> watchDof, double load)
        {
            // Output
            string pathName = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);

            string FileName = "ForcesDisplacementsData.txt";
            string Extension = Path.GetExtension(FileName);
            string fileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(FileName));
            string File = string.Format("{0}{1}", fileNameOnly, Extension);

            var append = false;
            using (var writer = new StreamWriter(File, append))
            {
                writer.WriteLine(" IncrementalLoad     IncrementalDisplacement");
                writer.WriteLine(" ");
            }
            append = true;
            for (int i = 0; i < log1.dofDisplacementsPerIter.Count; i++)
            {
                using (var writer = new StreamWriter(File, append)) // append mode to continue from previous increment
                {
                    writer.WriteLine($"       {load * (i + 1) / log1.dofDisplacementsPerIter.Count}             {log1.dofDisplacementsPerIter[i][0][watchDof[0][0]]}");
                }
            }

        }
    }

}
