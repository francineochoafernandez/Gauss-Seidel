#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
int i,j;
using namespace std;
typedef struct GaussSeidel //x= [ L+D]^(−1)  * (b− U x)
{
  float a[20][20],br[20][20], x[20],b[20];
  int n, maxit,indi;
  float err;

  void RellenaMat(void)
  {
    printf("\n\nDar la cantidad de incógnitas de la ecuación: ");
    scanf("%d",&n);
    for(i=0;i<n;i++)
    {
        printf("\n\nDa todos los coeficientes de la ecuacion %d: \n",i+1);
        for(j=0;j<n;j++)
          scanf("%f",&a[i][j]);
    }

    printf("\n\nDa los coeficientes del vector b: \n");
    for(j=0;j<n;j++)
      scanf("%f",&b[j]);

    for(j=0;j<n;j++)
      x[j]=0;

    printf("\n\nDar el error relativo y el numero de iteraciones:  \n");
    scanf("%f%d",&err,&maxit);

  }

  void ImprimeOp(void)
  {
    //Imprimiendo matrices
    for(i=0;i<n;i++)
    {
      cout << "|";
      for(j=0;j<n;j++)
      {
        printf("\t%f \t",a[i][j]);
      }
      cout << "|";
      if(i==n/2)
      {
        printf("   *   ");
      }
      else
      {
        printf("       ");
      }

      printf("| %f |\t",x[i]);

      if(i==n/2)
      {
        printf("   =   ");
      }
      else
      {
        printf("       ");
      }

      printf("| %f | \t",b[i]);
      printf("\n");
    }
    printf("\n");
  }

  void MetodoGauss(void)
  {
    int k=0;
    float suma;
    while(k<maxit)
    {
      suma=0;
      k++;
      for(i=0;i<n;i++)
      {
        suma=0;
        for(j=0;j<n;j++)
        {
          if(i!=j)
          {
            suma=suma+a[i][j]*x[j];
          }
          x[i]=(b[i]-suma)/a[i][i];
        }
      }
    }
    printf("\nDespues de %i\n",k);
    for(i=0;i<n;i++)
      printf("x[%i]: %f",i+1,x[i]);
      printf(" \n " );
  }

  void CopiaMatriz(float mat1[20][20],float mat2[20][20])
  {
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
      {
        mat2[i][j]=mat1[i][j];
      }
    }
  }

  void AcomodandoMatriz()
  {
    float renc, renc2,aux;
    int y,r=0;

    //Combirtiendola en diagonal Dominante
    for(i=0;i<n;i++)
    {
      renc=0;
      for(j=0;j<n;j++)
      {
        if (i!=j)
        {
          renc= renc + abs(br[i][j]);
        }
      }

      y=0;
      r=i;
      if(renc>abs(br[i][i]))//El renglon no cumple así que se cambia orden
      {
        for(int l=0;l<n;l++)//Este for solo es para que se repitan las lineas la n cantidad de veces.
        {//Aquí se hacen todas las posibles combinaciones para ver cual digito es el que iria en la diagonal
          renc2=0;
          aux=0;
          for(int k=0;k<n;k++)
          {
            if(k!=y)
              renc2= renc2 + abs(br[i][k]);
          }

          if(abs(br[i][y])>renc2)
          {
            for(i=0;i<n;i++)
            {
  	           aux=br[i][r];
  	           br[i][r]=br[i][y];
  	           br[i][y]=aux;
            }
          }
          y++;

        }

      }

    }

  }

  void ComparaDiago(int op)
  {
    //Comparando diagonal
    float ren;
    int c=0;

    for(i=0;i<n;i++)
    {
      ren=0;
      for(j=0;j<n;j++)
      {
        if (i!=j)
        {
          ren= ren + abs(br[i][j]);
        }
      }

      if(ren>abs(br[i][i]))
      {
        if(op==1)
        {
          printf("\nLa matriz NO es dominante, no va a convergir.\n");
          printf("\nSe tratara de mover la matriz.\n");
          AcomodandoMatriz();
        }
        else
        {
          printf("\nLa matriz NO se pudo acomodar para volverla dominante.\n");
          indi=3;
        }
        i=n;
        j=n;
        c++;
      }

    }
    if(c==0)
    {
      printf("\nLa matriz SI es dominante, sí va a convergir.\n");
      indi=1;
      if(op==2)
      {
        CopiaMatriz( br, a);
        ImprimeOp();
      }

    }
  }

}GS;

int main(void)
{
  GS sisecs;
  sisecs.RellenaMat();
  sisecs.ImprimeOp();
  sisecs.CopiaMatriz( sisecs.a, sisecs.br);
  sisecs.ComparaDiago(1);
  if(sisecs.indi!=1)
  {
    sisecs.ComparaDiago(2);
  }
  sisecs.MetodoGauss();

  return 0;
}
