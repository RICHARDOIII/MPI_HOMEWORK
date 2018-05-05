#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <mpi.h>
using std::cout;
using std::endl;

void smooth(float* x,float* y,int n,float a,float b, float c,int size)
{
	int n2=n+2;
	for( int i=1;i<n/size;i++)
	{	for( int j=1;j<n2;j++)
		{y[i*n2+j]=a * (x[(i-1)*n2+ (j-1)] + x[(i-1)*n2+ (j+1)] 
				+ x[(i+1)*n2 +(j-1)]+ x[(i+1)*n2 +(j+1)])
			+b * (x[(i-1)*n2+(j)] + x[(i+1)*n2 +(j)]
				+x[(i)*n2 +(j-1)] + x[(i)*n2+(j+1)])
			+ c*(x[(i)*n2+(j)]);
		}
	}
};
void initialize(float* x,int n)
{
	srand(time(NULL));
	int n2=n+2;
	for( int i=0;i<n2;i++)
	{	for(int j=0;j<n2;j++)
		{
			 x[i*n2 + j]= (float)rand()/(float)RAND_MAX;
		}
	}
};

	
void count(float* x,int n,float t,int& count1,int size)
{	count1=0;
	int n2=n+2;
	for (int i=1;i<n/size;i++)
	{	for(int j=1;j<n2;j++)
		{  if(x[i*n2 + j]<t)
			{ count1+=1;
			}
		}
	}
};

void main( int argc, char *argv[])
{	double begin,end,begin1,begin2,begin3,begin4, e1,e2,e3,e4,time,time1,time2,time3,time4;
	int rank,size;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int n,countx,county;
	float a,b,t,c,*x,*y;	
		
	n=pow(2,15);
	int n2=n+2;
/////////////////////////////////////////////////////////////////////////	
	//time the x allocation
	if(rank==0)
	{
		//double begin,end;
		begin=MPI_Wtime();
		x= (float *)malloc((n+2) * (n+2) * sizeof(long int));
		end = MPI_Wtime();
		 time = (end - begin);
	}
///////////////////////////////////////////////////////////////////
	// Setting constants
	a=.05;b=.1;c=.4;t=.1;


////////////////////////////////////////////////////////////////////
	if (rank==0)
	{
	// Initalize the thing with timings
		begin1 = MPI_Wtime();
		initialize(x,n);
		e1 = MPI_Wtime();
		time1= (e1 -begin1);
	}


//////////////////////////////////////////////////////////////////////
	float *X,*Y;
	int split=n/size;
	X=(float *)malloc( (split+2)*n2* sizeof(int));
	Y=(float *)malloc( (split+2)*n2 * sizeof(int));
	
	int *displ = (int *)malloc(size * sizeof(int)); 
	int *scount = (int *)malloc(size * sizeof(int));	
	for (int i=0; i<size;i++)
	{
		displ[i] = (split)*n2*i;
		scount[i] =(split+2)*n2; 
	}
	MPI_Scatterv( x ,scount,displ,MPI_FLOAT, X,scount[0], MPI_FLOAT,0, MPI_COMM_WORLD);
//////////////////////////////////////////////////////////////////////
	//time the smoothening
	if (rank==0)
	{
		begin1 =MPI_Wtime();
	}

	smooth(X,Y,n,a,b,c,size);
	if (rank==0)
	{
	e1 =MPI_Wtime();
	 time2= (e1-begin1);
///////////////////////////////////////////////////////////////////////
	}	 
	// time the count of the x array
	if (rank==0)
	{
	 begin2 = MPI_Wtime(); 
	}
	count(X,n,t,countx,size);
///////////////////////////////////////////////////////////////////////////////
	int globalx;
	MPI_Reduce(&countx, &globalx,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	if (rank == 0)
	{
	e2 = MPI_Wtime();
	time3=(e2-begin2);
//////////////////////////////////////////////////////////////		
	// time the count of the y array
	begin3 = MPI_Wtime(); 
	}
	count(Y,n,t,county,size);
//////////////////////////////////////////////////////////////
	int globaly;
	MPI_Reduce(&county,&globaly,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	if (rank==0)
	{	

	e3 = MPI_Wtime();
	 time4=(e3-begin3);
		
	//PRINT TIME
	cout<<endl;
	cout<< "Summary"<<endl;
	cout<< "-------"<<endl;
	cout << "Number of elements in a row/column         :: "<<n2<<endl;
	cout << "Number of inner elements in a row/column   :: "<<n<<endl;
	cout << "Total number of elements                   :: "<<(n2)*(n2)<<endl;
 	cout << "Total number of inner elements             :: " << n*n <<endl; 
	cout << "Memory (GB) used per array                 :: " << (n+2)*(n+2)*sizeof(float) / (float)(1024*1024*1024) << endl;
 	cout << "Threshold                                  :: " << t <<endl;
	cout << "Smoothing constants (a, b, c)              :: " << a << " " << b << " " << c <<endl;
	cout << "Number   of elements below threshold (X)   :: " << globalx << std::endl; 
	cout << "Fraction of elements below threshold       :: " << globalx / (float)(n*n) << endl;
	cout << "Number   of elements below threshold (Y)   :: " << globaly << endl;
	cout << "Fraction of elements below threshold       :: " << globaly / (float)(n*n) << endl;
	cout<<endl;
	cout << "------"<<endl;
	cout << "CPU: Alloc-X   :: "<<time<<endl;
	cout << "CPU: Init-X    :: "<<time1<<endl; 
	cout << "CPU: Smooth    :: "<<time2<<endl;
	cout << "CPU: Count-X   :: "<<time3<<endl;
	cout << "CPU: Count-Y   :: "<<time4<<endl;
	}
};	


