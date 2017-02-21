#include <lapacke.h>
#include <functional>
#include <cstring>
#include "../Utilities/LagrangePolynomials.h"
#include "../Utilities/PolyEval.h"
#include "../Utilities/LobattoIntegration.h"

using namespace std;

#ifndef MASSMATRIX_H
#define MASSMATRIX_H

void massMatrix(float *MassMatrix,unsigned N)
{
    float Poly[N+1][N+1];
    float **poly;
    poly   =   new float*[N+1];
    lagrangePolynomials(*Poly,N);
    unsigned i,j;

    for(i=0;i<=N;i++)
    {
        poly[i] =   new float[N+1];
        memcpy(poly[i],Poly[i],(N+1)*sizeof(float));
    }

    function<float(float)> eval;
    for(i=0;i<=N;i++)
    {
        for(j=0;j<=N;j++)
        {
            eval = [&poly,&i,&j,&N](float x){return (((polyEval(poly[i],N,x))*(polyEval(poly[j],N,x))));};
            MassMatrix[i*(N+1)+j] = lobattoIntegration(-1.0,1.0,N+1,eval);
        }
    }

    for(i=0;i<=N;i++)
        delete[] poly[i];

    delete[] poly;
    return ;
}

void twoDMassMatrix(float *MassMatrix, unsigned N)
{
    float m[N+1][N+1];
    massMatrix(*m,N);
    unsigned i1,i2,j1,j2;

    for(i1=0;i1<=N;i1++)
        for(j1=0;j1<=N;j1++)
            for(i2=0;i2<=N;i2++)
                for(j2=0;j2<=N;j2++)
                    MassMatrix[(i1*(N+1)+j1)*(N+1)*(N+1)+i2*(N+1)+j2] = m[i1][i2]*m[j1][j2];
    return ;
}

#endif