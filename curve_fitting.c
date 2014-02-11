#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define IDEBUG 0
#define SIZE 8
int DEGREE ;
int vector_print_off (int nr, double *x)    //prints the vector with array x and row size nr
{
    int i;
    if ( nr <= 0 ) return (-1);
    for (i = 0; i < nr; i++) 
    {
        printf ("%9.4f  \n", x[i]);
    }
    printf("\n");  
    return (0);
}
void gauss(double ** a,int N,double * sol)  // does gauss elimination and stores the result in sol array
{
    int i,j,k;
    float t,det=1;  
          
    for(i=0;i<N;i++)  
    {
        for(j=0;j<N;j++)  
        {
            if(i!=j)  
            {  
                t=a[j][i]/a[i][i];  
                for(k=0;k<N+1;k++)  
                    a[j][k]-=a[i][k]*t;  
            }
        }
    }  
    for(i=0;i<N;i++)
    {
        det*=a[i][i];
    }
    if(det==0)    // if determinant is 0
    {  
        printf("\nThe matrix is singular .\n");
    }  

    if(IDEBUG==1)
    {       
        printf("\nThe Gauss-Jordan Matrix is :\n\n");  
        for(i=0;i<N;i++)  
        {  
            for(j=0;j<N+1;j++)  
                printf("%.4f ",a[i][j]);  
            printf("\n");  
       }
    }     
    for(i=0;i<N;i++)
        sol[i]=a[i][N]/a[i][i];
}

int summation_2(int r,double data_set_y[],double data_set_x[])     // function for summation[(x^i)y]
{
    int ret=0;
    int i,j;
    for(i=0;i<SIZE;i++)
    {
        ret=ret+(pow(data_set_x[i],r)*data_set_y[i]);

    }
    return ret;
}


int degree_summation(double data_set_x [] ,int R)
{
    int i,ret=0;
    for(i=0;i<SIZE;i++)
    {
        ret=ret+(pow(data_set_x[i],R));   
    }
    return ret;
}

void fitting(double data_x[],double data_y[],double * sol)              //curve fitting
{
   
    int j;
    int i=0;
  
    double** data_matrix;
    data_matrix = (double**)malloc((DEGREE+1)*sizeof(double*));
    for(i=0;i<DEGREE+1;i++)
    {
      data_matrix[i] = (double*)malloc((DEGREE+2)*sizeof(double));
    }

    for(i=0;i<DEGREE+1;i++)
    {
        for(j=0;j<DEGREE+1;j++)
        {
            data_matrix[i][j]=degree_summation(data_x,i+j);
        }
    }
    
    for(i=0;i<DEGREE+1;i++)
    {
        data_matrix[i][DEGREE+1]=summation_2(i,data_y,data_x);
    }
    if(IDEBUG==1)
    {   
        for(i=0;i<DEGREE+1;i++)
        {
            for(j=0;j<DEGREE+2;j++)
                printf("%f   ",data_matrix[i][j]);
            printf("\n");
        }
    }
    gauss(data_matrix,DEGREE+1,sol);   
 //   vector_print_off(DEGREE+1,sol);
}

double value(double *sol,double data_x)         //value of the polynomial at particular x
{
    int i=0;
    double  ret=0;
    int j=0;
    for(j=0;j<DEGREE+1;j++)
    {
        int c=pow(data_x,j);
        ret=ret+( sol[j] * c );
        
    }
    return ret;
}

void error(double data_x[],double data_y[],double sol[],double *err)   // difference between y and func(x) stored in err array 
{
    int i, ret,j=0;
    for(i=0;i<SIZE;i++)
    {
        err[i]=data_y[i]-value(sol,data_x[i]);
    }
}

double error2(double data_x[],double data_y[],double sol[])           // error = 1/2[summation(err^2) from 0 to SIZE ]
{
    int i;
    double ret=0;
    double err[SIZE];
    error(data_x,data_y,sol,err);
    for(i=0;i<SIZE;i++)
    {
        ret=ret + (err[i]*err[i]);

    }
    double c=ret/2;
    return c;
}

int main()
{
       double data_x[]={-1,0,1,2,3,5,7,9};       
       double data_y[]={-1,3,2.5,5,4,2,5,4};
       //double data_x[]={1,2,3,4,5,6,7,8};
       //double data_y[]={1,2,3,4,5,6,7,8};

      for(DEGREE=1;DEGREE<SIZE;DEGREE++)
      {
    	    int i=0;
    	    double sol[DEGREE+1];   
    	    double err[SIZE];
	    fitting(data_x,data_y,sol);         
            error(data_x,data_y,sol,err);
        
	  /*  printf("\n  error \n");
            for(i=0;i<SIZE;i++)
            {
	            printf("%f \n",err[i]); 
            }*/
            double z=error2(data_x,data_y,sol);
            printf("\n error= %f for  degree= %d  \n ",z,DEGREE);
    }
}

