// Simple One - Dimensional Steady State Heat conduction Solver
// 
//

#include <iostream>
#include <cmath>
#include <vector>
#include <matplot/matplot.h>
//include "SOR.h"

using namespace std;

std::vector<double> FDMHeatConduction_OneDimension(int nodes, double LeftWallBoundary, double RightWallBoundary,
    int iterations, bool insulated)
{
    
    //double dx = domainLength / nodes;
    vector<double> NodesTemps = vector<double>(nodes, 0.0);
    int p = 1;
    NodesTemps[0] = LeftWallBoundary;
    
    NodesTemps[nodes-1] = RightWallBoundary;



    for (int j = 0; j < iterations; j++)
    {
        if (insulated)
        {
            NodesTemps[nodes - 1] = NodesTemps[nodes - 2];
        }
        for (int i = 1; i < nodes - p; i++)
        {
            NodesTemps[i] = (NodesTemps[i - 1] + NodesTemps[i + 1]) / 2;
        }
        
    }
    
    return NodesTemps;


}



std::vector<double> SORSolver(std::vector<std::vector<double>> a, std::vector<double> b, double e, double guess, double Omega, int iteration, bool Insulated)
{
    int p = b.size();
    std::vector<double> u(p, guess);
    std::vector<double> uOld(p, guess);
    double tol = e;
    double delta = tol + 1;
    int Fromlast = 1;
    //Boundary conditions at delta x 0 and delta x L
    u[0]= b[0];
    uOld[0] = b[0];
    u[p - 1] = b[p - 1];
    uOld[p - 1] = b[p - 1];
    


    //Loop through  iterations or exit if solution is converged
    for (int j = 0; j < iteration; j++)
    {
        
        for (int i = 1; i < p - 1; i++)
        {            
            if (Insulated && i ==p - 2)
            {
                u[p - 1] = u[p - 2];
                uOld[p - 1] = uOld[p - 2];
            }
            u[i] = Omega * (guess +(1/a[i][i]) * (u[i-1] +  uOld[i+1])) + (1 - Omega) * uOld[i]; 
            delta = abs(u[i]- uOld[i]) * 100. / u[i];
            uOld[i] = u[i];
            if (delta <= tol)
            {
                cout << "Converged at iteration " << j << "\n";
                return u;
            }
        }
    }

    return u;
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


using std::cin, std::cout;

int main() {
   using namespace matplot;

   //Heat conduction using fdm
   vector<double> NodeTemps = FDMHeatConduction_OneDimension(10, 20, 60, 100, false);
    
    
    std::vector<std::vector<double>> C = {
        NodeTemps, NodeTemps};
    //Uncomment the below to see the FDM Solution 
    //imagesc(C);
   //colorbar();
   //grid(false);
   //xlabel("Node");
   //yticklabels({});
    
    ///
    //show();
    

    /// Below part is for SOR method
    int gridPoints = 10;
    double boundaryLeft =20;
    double boundaryRight =60;
    
    //diagonal matrix derived from Taylor expansion
    std::vector<std::vector<double>> a = CreateMatrixStuctureGS(gridPoints);
    std::vector<double> b(gridPoints, 0);
    
    //boundary conditions at start and end
    b[0] = boundaryLeft;
    b[gridPoints - 1] =  boundaryRight;
  
    std::vector<double> NodeTemps2 = SORSolver(a,b,0.000001,0,1.2, 100, false);
    std::vector<std::vector<double>> D = {
       NodeTemps2, NodeTemps2 };
    imagesc(D);
    colorbar();
    grid(false);
    xlabel("Node");
    yticklabels({});
    
   
    show();
}