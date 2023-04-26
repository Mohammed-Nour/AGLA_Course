#include<bits/stdc++.h>
using namespace std;

const double EPSILION = abs(1e-4); // Define a constant value for the epsilon error

// Define Matrix class
class Matrix
{
public:
    // Constructor to initialize matrix with size n x m and values from input vector
    Matrix(int n, int m, vector<vector<double>> matrix)
    {
        this->n = n;
        this->m = m;
        this->matrix= matrix;
    }
 Matrix(int n)
    {
        this->n = n;
        this->m = n;
        this->matrix.resize(n,vector<double>(n,0));
    }
    // Getter method to return number of rows in matrix
    int getN()
    {
        return this->n;
    }

    // Getter method to return number of columns in matrix
    int getM()
    {
        return this->m;
    }

    // Getter method to return matrix as a vector of vectors
    vector<vector<double>> getMatrix()
    {
        return this->matrix;
    }

    // Overload assignment operator to assign one matrix to another
    Matrix* operator=(Matrix& matrix)
    {
        this->matrix = matrix.getMatrix();
        return this;
    }

    // Overload addition operator to add two matrices together
    Matrix* operator+(Matrix &matrix)
    {
        Matrix* D;
        int m =matrix.getM() ;
        int n =matrix.getN() ;
        vector<vector<double>> newMatrix(n,vector<double>(m,0));
        for(int i = 0 ; i<n; i++)
        {
            for(int  j = 0 ; j <m; j++)
            {
                newMatrix[i][j]=this->getMatrix()[i][j] + matrix.getMatrix()[i][j];
            }
        }
        D = new Matrix(n,m,newMatrix);
        return D;
    }
    friend ostream &operator<<( ostream &output, Matrix& newOne )
    {
        output << setprecision(2) << fixed;
        output <<  newOne.matrix[newOne.n][newOne.m];
        return output;
    }

    friend istream &operator>>( istream  &input, Matrix &newOne)
    {

        input >>newOne.n >> newOne.matrix[newOne.n][newOne.m];
        return input;
    }
    // Overload subtraction operator to subtract one matrix from another
    Matrix* operator-(Matrix &matrix)
    {
        Matrix* E;
        int m =matrix.getM() ;
        int n =matrix.getN() ;
        vector<vector<double>> newMatrix(n,vector<double>(m,0));
        for(int i = 0 ; i<n; i++)
        {
            for(int  j = 0 ; j <m; j++)
            {
                newMatrix[i][j]=this->getMatrix()[i][j] - matrix.getMatrix()[i][j];
            }
        }
        E = new Matrix(n,m,newMatrix);
        return E;
    }

    // Overload multiplication operator to multiply two matrices together
    Matrix* operator*(Matrix& matrix)
    {
        Matrix* F;
        int secondM =matrix.getM() ;
        int secondN =matrix.getN() ;
        int firstM =this->getM() ;
        int firstN=this->getN() ;
        vector<vector<double>> newMatrix(firstN,vector<double>(secondM,0));
        for(int i = 0 ; i<firstN; i++)
        {
            for(int  j = 0 ; j <secondM; j++)
            {
                for(int k = 0 ; k < secondN ; k++)
                {
                    newMatrix[i][j] +=   this->getMatrix()[i][k]* matrix.getMatrix()[k][j];
                }
            }
        }
        F = new Matrix(firstN,secondM,newMatrix);
        return F;
    }

    // This function returns a new matrix that is the transpose of the original matrix.
    Matrix* transposeMatrix()
    {
        Matrix* G;
        int m = this->getM(); // get the number of columns in the original matrix
        int n = this->getN(); // get the number of rows in the original matrix
        vector<vector<double>> newMatrix(m, vector<double>(n, 0)); // create a new matrix with swapped dimensions
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < m; j++)
            {
                // swap the indices of the current element to get the transposed element
                newMatrix[j][i] = this->getMatrix()[i][j];
            }
        }
        G = new Matrix(m, n, newMatrix);
        return G;
    }
     void printMatrix()
    {
        cout << setprecision(4) << fixed;
        for(int i = 0; i < this->getN(); i++)
        {
            for(int j = 0; j < this->getM(); j++)
            {
                cout<<this->getMatrix()[i][j] << " ";
            }
            cout << endl;
        }
    }
public:
    vector<vector<double>> matrix; // The matrix itself is represented as a 2D vector
    int n, m; // The number of rows and columns in the matrix
};

class SequareMatrix: public  Matrix
{
public:
    // Constructor to initialize matrix with size n x m and values from input vector
    SequareMatrix(int n,vector<vector<double>> matrix) :  Matrix(n,n,matrix) {};
    SequareMatrix(int n) :  Matrix(n) {};
    // Overload assignment operator to assign one matrix to another
    SequareMatrix* operator=(SequareMatrix& matrix)
    {
        this->matrix = matrix.getMatrix();
        return this;
    }

    // Overload multiplication operator to multiply two matrices together
    SequareMatrix* operator*(SequareMatrix& matrix)
    {
        SequareMatrix* F;
        vector<vector<double>> newMatrix(this->getN(),vector<double>(this->getN(),0));
        for(int i = 0 ; i<this->getN(); i++)
        {
            for(int  j = 0 ; j <this->getN(); j++)
            {
                for(int k = 0 ; k < this->getN() ; k++)
                {
                    newMatrix[i][j] +=   this->getMatrix()[i][k]* matrix.getMatrix()[k][j];
                }
            }
        }
        F = new SequareMatrix(this->getN(),newMatrix);
        return F;
    }
    void diagonalNormalization(SequareMatrix* sequareMatrix)
    {
        for(int i = 0 ; i< this->getN() ;i++)
        {
            double piovt =this->getMatrix()[i][i];
            for(int j =0 ; j< this->getN() ;j++)
            {
                this->matrix[i][j] = this->matrix[i][j]/piovt;
                sequareMatrix->matrix[i][j] = sequareMatrix->matrix[i][j] /piovt;
            }
        }
    }


};

class IdentityMatrix: public  SequareMatrix
{
public:
    // Constructor to initialize matrix with size n x m and values from input vector
    IdentityMatrix(int n) :  SequareMatrix(n) {};
    IdentityMatrix(int n, vector<vector<double>> matrix) :  SequareMatrix(n)
    {
        this->matrix = matrix;
    };
       IdentityMatrix(SequareMatrix& sequareMatrix,int n) :  SequareMatrix(n)
    {
        this->matrix = sequareMatrix.getMatrix();
    };
    IdentityMatrix* identity(int n)
    {
        for(int j = 0 ; j < n ; j++)
        {
            for(int k = 0 ; k < n ; k++)
            {
                if(j==k)
                {

                    this->matrix[j][k] = 1;

                }
                else
                {
                    this->matrix[j][k] = 0;
                }
            }
        }
        return this;
    }

};
class EliminationMatrix: public  SequareMatrix
{
public:
    // Constructor to initialize matrix with size n x m and values from input vector
    EliminationMatrix(int n) :  SequareMatrix(n) {};
    EliminationMatrix(int n, vector<vector<double>> matrix) :  SequareMatrix(n)
    {
        this->matrix = matrix;
    };

    EliminationMatrix* elimination(SequareMatrix*  matrix,int index,int secondIndex,double pivot)
    {
        for(int i = 0 ; i <this->getN() ; i++)
        {
            for(int j = 0; j<this->getN() ; j++)
            {
                if(i == index &&j==secondIndex)
                {
                    this->matrix[i][j] = -(matrix->getMatrix()[i][j]/pivot);
                }
                else if( i ==j )
                {
                    this->matrix[i][j]  = 1;
                }
                else
                {
                    this->matrix[i][j]  =0;
                }
            }
        }
        return this;
    }

};
class PermutationMatrix: public  SequareMatrix
{
public:
    // Constructor to initialize matrix with size n x m and values from input vector
    PermutationMatrix(int n,vector<vector<double>> matrix) :  SequareMatrix(n)
    {
        this->matrix = matrix;
    };
    PermutationMatrix(int n) :  SequareMatrix(n) {};
    PermutationMatrix* permutation(SequareMatrix*  matrix,int index,int secondIndex)
    {
        for(int i = 0 ; i <this->getN(); i++)
        {
            for(int j = 0; j<this->getN() ; j++)
            {
                if(i==j)
                {
                    this->matrix[i][j] =1;
                }
                else
                {
                    this->matrix[i][j] =0;
                }
            }
        }
        for(int i = 0 ; i <this->getN() ; i++)
        {
            for(int j = 0; j<this->getN() ; j++)
            {
                if( i==index)
                {
                    if(this->matrix[i][j]  == 1)
                    {
                        this->matrix[i][j] =0;
                        this->matrix[i][secondIndex] =1;
                    }
                }
                else if(i ==secondIndex)
                {
                    if(this->matrix[i][j]  == 1)
                    {
                        this->matrix[i][j] =0;
                        this->matrix[i][index] =1;
                    }
                }
            }
        }
        return this;
    }
};
int m;
Matrix* A;
Matrix* ATranspose;
Matrix* B;
Matrix* ATransposeMultiA;
Matrix* ATransposeMultiB;
SequareMatrix* inverse;
SequareMatrix* AInverse;
Matrix* resultX;
vector<double>t;
int degree;

void readInput()
{
    cin >> m;
    vector<vector<double>> b(m,vector<double>(1,0));
    t.resize(m,0);
    for(int j = 0 ; j < m ; j++)
    {
        cin>>t[j]>>b[j][0];
    }
    B = new Matrix(m,1,b);

    cin>>degree;
    int index =0;
    vector<vector<double>> matrix(m,vector<double>(degree+1,0));
    for(int j = 0 ; j < m ; j++)
    {
        for(int k = 0 ; k <degree+1 ; k++)
        {
            if(k==0)
            {
                matrix[j][k] = 1;
            }
            else
            {
                matrix[j][k] = pow(t[index],k);
            }
        }
        index++;
    }

    A = new Matrix(m,degree+1,matrix);
}
void directWay()
{
    inverse = (SequareMatrix*)ATransposeMultiA;
    IdentityMatrix* identityMatrix = new IdentityMatrix(inverse->getN());
    identityMatrix = identityMatrix->identity(inverse->getN());
    AInverse =identityMatrix;
    for(int i = 0 ; i <inverse->getN() ; i++)
    {
        double maximum=abs(inverse->getMatrix()[i][i]);
        int index = i,secondIndex = i;
        for(int j = i+1 ; j< inverse->getN() ; j++)
        {
            if(abs(inverse->getMatrix()[j][i]) > maximum )
            {
                maximum = abs(inverse->getMatrix()[j][i]) ;
                index =j;
                secondIndex =i;
            }
        }
        if(maximum !=abs(inverse->getMatrix()[i][i]))
        {
            PermutationMatrix* permutationMatrix = new PermutationMatrix(inverse->getN());
            permutationMatrix =permutationMatrix->permutation(inverse,index,secondIndex);
            inverse = *permutationMatrix * *inverse;
            AInverse = *permutationMatrix * *AInverse;
        }
        double piovt = inverse->getMatrix()[i][i];
        for(int k = i+1 ; k < inverse->getN() ; k++)
        {
            if(inverse->getMatrix()[k][i]!=0)
            {
                EliminationMatrix* eliminationMatrix = new EliminationMatrix(inverse->getN());
                eliminationMatrix =eliminationMatrix->elimination(inverse,k,i,piovt);
                inverse = *eliminationMatrix * *inverse;
                AInverse = *eliminationMatrix * *AInverse;
            }
        }
    }
}
void backWay()
{
    for(int i = inverse->getN()-1 ; i >-1 ; i--)
    {
        double piovt =inverse->getMatrix()[i][i];
        for(int k = i-1 ; k > -1; k--)
        {
            if(inverse->getMatrix()[k][i]!=0)
            {
                EliminationMatrix* eliminationMatrix = new EliminationMatrix(inverse->getN());
                eliminationMatrix =eliminationMatrix->elimination(inverse,k,i,piovt);
                inverse = *eliminationMatrix * *inverse;
                AInverse =*eliminationMatrix * *AInverse;
            }
        }
    }
    inverse->diagonalNormalization(AInverse);
}

int main()
{


    readInput();
    cout<<"A:\n";
    A->printMatrix();
    ATranspose = A->transposeMatrix();
    ATransposeMultiA = *ATranspose * *A;
    cout<<"A_T*A:\n";
    ATransposeMultiA->printMatrix();
    directWay();
    backWay();
    cout<<"(A_T*A)^-1:\n";
    AInverse->printMatrix();
    cout<<"A_T*b:\n";
    ATransposeMultiB = *ATranspose * *B;
    ATransposeMultiB->printMatrix();
    cout<<"x~:\n";
    resultX= (Matrix)*AInverse * *ATransposeMultiB;
    resultX->printMatrix();

    return 0;

}
