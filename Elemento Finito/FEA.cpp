#include<iostream>
#include<cstdlib>
#include<fstream>
#include<chrono>
#include<unistd.h>

using namespace std;
int r=0;
int c=0;
int pot;
float factor;
float pivote;
float mult;

int main()
{
    ifstream f("stiffness_matrix.txt"); //Lectura de la matriz de rigidez
    /*El primer valor del archivo es el tamaño de la matriz que es cuadrada por definicion.
    El segundo valor es el factor de división*/
    f >> r >> pot;
    float K[r][r+1], cK[r][r+1], coefK[r][r], idenK[r][r]; //Definicion de las matrices a usar.
    float F[r]; float NodDis[r]; //Definicion de los vectores a usar.

    for (int i = 0; i < r; i++) //Formacion de la matriz de rigidez expandida con los valores del archivo.
        {
            for (int j = 0; j<=r; j++)
            {
                f >> K[i][j];
                cK[i][j]=K[i][j]; //Creacion de una copia de la matriz para realizar las operaciones sobre ella.
            }

        }
    //EXTRACCION DE LA MATRIZ DE COEFICIENTES
    for (int j=0; j<r; j++)
    {
        for (int k=0; k<r; k++)
        {
            coefK[j][k]=K[j][k];
        }
    }
    //EXTRACCION DEL VECTOR FUERZA
    for(int i=0; i<r; i++)
    {
        F[i]=K[i][r];
    }

    //METODO DE GAUSS-JORDAN
    auto start_gauss=chrono::system_clock::now(); //Inicio del cronometro
    //Proceso de creacion de matriz triangular superior
    for(int p=0; p<=r-2; p++)
		{
			pivote=cK[p][p];
			for(int g=p+1; g<=r-1; g++)
			{
				mult=cK[g][p];
				for(int a=p; a<=r; a++)
				{
					cK[g][a]=((pivote*cK[g][a])-(mult*cK[p][a]));
				}
			}
		}
    //Proceso de eliminacion de los valores superiores de la matriz
    for(int p=r-1; p>=r-2; p--)
    {
        pivote=cK[p][p];
        for (int g=0; g<=p-1; g++)
        {
            factor=((cK[g][p])/(pivote));
            for (int a=p; a<=r; a++)
            {
                cK[g][a]=((cK[g][a])-(factor*cK[p][a]));
            }

        }
    }
    //Obtencion de los valores buscados del sistema
    for (int g=0; g<r; g++)
    {
        cK[g][r]=((cK[g][r])/(cK[g][g]));
        cK[g][g]=((cK[g][g])/(cK[g][g]));
    }
    //Division del resultado entre el factor de la matriz (Multiplicacion de escalar por matriz)
    for(int i=0; i<r; i++)
    {
            cK[i][r]=(cK[i][r])/pot;
    }
    sleep(1);
    auto end_gauss=chrono::system_clock::now(); //Finalizacion de tiempo de ejecucion Gauss-Jordan
    chrono::duration<double, milli> duration_gauss = end_gauss - start_gauss; //Tiempo total de ejecucion


    //METODO POR MATRIZ INVERSA
    auto start_minv=chrono::system_clock::now(); //Inicio del cronometro
    for (int i = 0; i < r; i++) { //Se genera una matriz identidad
        for (int j = 0; j < r; j++)
            {
            if (j==i) //Se agregan 1's en las entradas de la diagonal principal
                idenK[j][i] = 1;
            else    //Se hacen 0 a las demas entradas
                idenK[j][i] = 0;
        }
    }

    for (int i=0; i<r; i++) //Operaciones para obtener la inversa
    {
        float piv = coefK[i][i]; //Pivote
        if (piv == 0) //Se verifica que el pivote sea distinto a 0
        {
            for (int j=0; j<r; j++)
            {
                for (int k=0; k<r; k++)
                {
                    if (k == i) //Nos colocamos en la posicion del pivote
                    {
                        if (j > i) //Se buscan valores a partir de la posicion del pivote
                        {
                            if (coefK != 0) //Columnas diferentes de 0
                            {
                                for (int l=0; l<r; l++)
                                {
                                    //Se opera sobre la matriz identidad generada
                                    float var_temporal = idenK[k][l];
                                    idenK[k][l] = idenK[j][l];
                                    idenK[j][l] = var_temporal;

                                    // Operacion en la matriz original

                                    var_temporal = coefK[k][l];
                                    coefK[k][l] = coefK[j][l];
                                    coefK[j][l] = var_temporal;

                                }
                            }
                        }
                    }
                }
            }
        }
        piv = coefK[i][i]; //Se actualiza el valor del pivote si se cambian filas
        for (int j=0; j<r; j++) //Se divide la fila actual entre el pivote
        {
            idenK[i][j] = (idenK[i][j]/piv);
            coefK[i][j] = (coefK[i][j]/piv);
        }

        for (int j=0; j<r; j++) //Restas entre filas para reducir la matriz
            {
                if (i != j) //Nos colocamos en una posicion diferente del pivote
                {
                    float temporal = coefK[j][i];
                    for (int l=0; l<r; l++)
                    {
                        idenK[j][l] = idenK[j][l] - temporal*idenK[i][l];
                        coefK[j][l] = coefK[j][l] - temporal*coefK[i][l];
                    }
                }
            }
    }

   //Multiplicacion matriz inversa * Vector fuerza (x=A^-1*b)

    for (int h=0; h<r; h++)
    {
        NodDis[h]=0;
        for (int l=0; l<r; l++)
        {
            NodDis[h] += idenK[h][l]*F[l];
        }
        NodDis[h]=NodDis[h]/pot; //Division del vector resultante entre el factor de la matriz
    }
    sleep(1);
    auto end_minv=chrono::system_clock::now(); //Finalizacion de tiempo de ejecucion por Matriz Inversa
    chrono::duration<double, milli> duration_minv = end_minv - start_minv; //Tiempo total de ejecucion




    //IMPRESION DE DATOS Y RESULTADOS

    cout<<"---------------------------------------------------"<<endl;
    cout<<"              Matriz original extendida "<<endl;
	cout<<"---------------------------------------------------"<<endl;
	for(int g=0; g<r; g++)
	{
		for(int a=0; a<=r; a++)
		{
			cout<<K[g][a]<<"\t";
		}
		cout<<endl;
	}
	cout<<endl;

	cout<<"---------------------------------------------------"<<endl;
	cout<<"                  Matriz Reducida "<<endl;
	cout<<"---------------------------------------------------"<<endl;
	for(int g=0; g<r; g++)
	{
		for(int a=0; a<=r; a++)
		{
		  cout<<cK[g][a]<<"\t";
		}
		cout<<endl;
	}
	cout<<endl;
	cout<<"---------------------------------------------------"<<endl;
	cout<<"                  Matriz Inversa "<<endl;
	cout<<"---------------------------------------------------"<<endl;
	for(int g=0; g<r; g++)
    {
        for(int h=0; h<r; h++)
        {
            cout<<idenK[g][h]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;

	cout<<"---------------------------------------------------"<<endl;
	cout<<"               Desplazamientos nodales "<<endl;
	cout<<"---------------------------------------------------"<<endl;
	cout<<"******************Por Gauss-Jordan*****************"<<endl;
    for (int g=0; g<r; g++)
    {
        cout<<cK[g][r]<<"\t";
    }
    cout<<endl;
    cout<<"*****************Por Matriz inversa*****************"<<endl;
    for (int g=0; g<r; g++)
    {
        cout<<NodDis[g]<<"\t";
    }
    cout<<endl;

    cout<<"---------------------------------------------------"<<endl;
	cout<<"                     Duracion "<<endl;
	cout<<"---------------------------------------------------"<<endl;
	cout<<"******************Por Gauss-Jordan*****************"<<endl;
    cout<<duration_gauss.count()<<" ms"<<endl;
    cout<<"*****************Por Matriz Inversa****************"<<endl;
    cout<<duration_minv.count()<<" ms"<<endl;


    return 0;
}

