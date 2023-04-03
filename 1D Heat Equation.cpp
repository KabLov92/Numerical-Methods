// Simple One - Dimensional Steady State Heat conduction Solver
// 
//

#include <iostream>
#include <cmath>
#include <vector>
#include <matplot/matplot.h>
//include "SOR.h"

using namespace std;

int FDMHeatConduction_OneDimension(int nodes, double LeftWallBoundary, double RightWallBoundary, double tolerance,
    bool insulated)
{
    using namespace matplot;
    int iterations= 0;
    //double dx = domainLength / nodes;
    vector<double> NodesTemps = vector<double>(nodes, 0.0);
    vector<double> oldIteration = vector<double>(nodes, 0.0);
    int p = 1; 
    std::vector<std::vector<double>> lines;
    
    double mse = 1.;
    double t = 0;
    double tOld = 0;
    
    oldIteration = NodesTemps;
    while (mse > tolerance)
    {
            iterations++;
            NodesTemps[0] = LeftWallBoundary;
            if (insulated)
            {
                NodesTemps[nodes - 1] = NodesTemps[nodes - 2];
            }
            for (int i = 1; i < nodes - p; i++)
            {
                NodesTemps[i] = (NodesTemps[i - 1] + NodesTemps[i + 1]) / 2;
                t += NodesTemps[i];
                tOld+= oldIteration[i];
            }

            NodesTemps[nodes - 1] = RightWallBoundary;
            lines.push_back(NodesTemps);
            mse = abs((t / nodes) - (tOld / nodes))*100./ (t / nodes);
            oldIteration = NodesTemps;
            //cout << "MSE: " << mse << "\n";
            t = 0;
            tOld = 0;

        
    }

   
    cout <<"iteration num "  << iterations <<"\n";
    xlabel("Nodes");
    ylabel("Temperature ( C°)");
    plot(lines);
    

    show();

  return iterations;
}



int SORSolver(std::vector<std::vector<double>> a, std::vector<double> b, double e, double guess, double Omega, int iteration, bool Insulated)
{

    using namespace matplot;


    int p = b.size();
    std::vector<double> u(p, guess);
    std::vector<double> uOld(p, guess);
    double tol = e;
    double delta = tol + 1;
    int Fromlast = 1;
    //Boundary conditions at delta x 0 and delta x L
    int iterations = 0;

    auto Xx = linspace(0, p);
    std::vector<std::vector<double>> lines;
    while (delta > tol)
    {
        iterations++;
        //Loop through  iterations or exit if solution is converged
    
            u[0] = b[0];
            uOld[0] = b[0];
            for (int i = 1; i < p - 1; i++)
            {
                if (Insulated && i == p - 2)
                {
                    u[p - 1] = u[p - 2];
                    uOld[p - 1] = uOld[p - 2];
                }
                u[i] = Omega * (guess + (1 / a[i][i]) * (u[i - 1] + uOld[i + 1])) + (1 - Omega) * uOld[i];
                delta = abs(u[i] - uOld[i]) * 100. / u[i];
                uOld[i] = u[i];

            }
            u[p - 1] = b[p - 1];
            uOld[p - 1] = b[p - 1];

            lines.push_back(u);
        
    }
    cout << "SOR iterations: " << iterations << "\n";
    plot(lines);
    hold(true);
    xlabel("Nodes");
    ylabel("Temperature ( C°)");
    show();
    
    return iterations;
}



std::vector<std::vector<double>>  CreateMatrixStuctureGS(int Nodes)
{
    std::vector<std::vector<double>> Matrix = std::vector<std::vector<double>>(Nodes, std::vector<double>(Nodes, 0.0));
	Matrix[0][0] = 2;
	Matrix[Nodes - 1][Nodes - 1] = 2;
    for (int i = 1; i < Nodes - 1; i++)
    {
		Matrix[i][i - 1] = 1;
		Matrix[i][i] = 2;
		Matrix[i][i + 1] = 1;
	}
    


	return Matrix;
}


void ShowIterationNumber()
{
    using namespace matplot;
    std::vector<int> x;
    std::vector<int> xGridlines;
    int iterationTotal = 0;
    //Heat conduction using fdm
    for (int i = 1; i < 100; i++)
    {
        iterationTotal = FDMHeatConduction_OneDimension(5 * i, 20, 60, 1e-6, false);
        x.push_back(iterationTotal);
        xGridlines.push_back(5 * i);
        cout << "Total iterations for " << 5 * i << "nodes :" << iterationTotal << "\n";

    }
    plot(xGridlines, x);
    //xticklabels(xGridlines);
    show();
}


using std::cin, std::cout;

int main() {
   using namespace matplot;
   std::vector<int> x;
   ///uncomment below for FDM
   //FDMHeatConduction_OneDimension(10, 30, 30, 1e-6, true);
   
    int gridPoints = 10;
    double boundaryLeft =30;
    double boundaryRight =30;
    
    //diagonal matrix derived from Taylor expansion
    std::vector<std::vector<double>> a = CreateMatrixStuctureGS(gridPoints);
    std::vector<double> b(gridPoints, 0);
    
    //boundary conditions at start and end
    b[0] = boundaryLeft;
    b[gridPoints - 1] =  boundaryRight;
  
    std::vector<double> NodeTemps2 = {};
    //int e = SORSolver(a,b,1e-6,0,1.2, 300, true);

    
    SORSolver(a, b, 1e-6, 0, 1.2, 300, true);
       
   
   

  return 0;
}