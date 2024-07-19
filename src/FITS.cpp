/******************************************************************************************/
/***************************    Head Files needed     **********************************/
/******************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include<math.h>
#include<vector>
#include<algorithm>
using namespace std;
/******************************************************************************************/
/********************    Data Structure and Global Variables    ************************/
/******************************************************************************************/
typedef struct Solution{
	int *p;
	double *SizeG;
	double cost;
}Solution;

typedef struct Neighborhood{
	int  type;
	int  v;
	int  g;
	int  x;
	int  y;
}Neighborhood;

typedef struct Pair{
	int i;
	int g;
}Pair;

#define Num 10					//种群大小
#define NumInner 10000			//禁忌搜索深度
#define TabuLen 10
#define UP 5

char * File_Name;
char * Output_File_Name;
char * Solution_File;

int N, K;					 // node number and group number
double f, f_best, f_best1;	//f_best1 记录IVNS_TS method IVNS阶段跑完的最好结果
Solution CS, NS, GS, OS;	//CS:current solution ,GS: global solution
Neighborhood *Neighbors;
Pair *PairSet;
int **TabuTenure;

double total_time, starting_time, Time_limit;
int * p;					// assign  every vertex to a cluster, 
int * bestp;				// patition array for each vertex 
//int * bestp ;				
double * SizeG;
double ** Delta_Matrix;		// incremental matrix 
double ** Distance;			// distance matrix between elements
double ** DistanceMulti;
double * LB;				// Lower bound of the number of elements in the group i 
double * UB;				// Upper bound of the number of elements in the group i 
double * w;					// node weights

double penalty_factor;
double ys, ys_best;
double gs;
int penalty_count;

//以下变量用于混合进化算法
Solution *Pop;
//int **Cluster1_Node;				//存储解（s1）每个cluster中顶点编号
//int **Cluster2_Node;				//存储解（s2）每个cluster中顶点编号
double *Cluster_F[2];				//存储解（s）每个cluster的object value
double *Pop_f;						//种群每个个体的object value

int **Pop_Cluster_Len;				//用于计算种群相似性的变量，存储个体每个Cluster的长度
int ***Pop_Cluster_Vec;				//用于计算种群相似性的变量，存储个体每个Cluster包含的顶点编号
//double *Cluster2_F;				//存储解（s2）每个cluster的object value
double ComputeMinF();
void RandomShake(int L, int p1[], double SizeGroup[]);
void RandomShake1(int L, int p1[], double SizeGroup[], int iter);
/******************************************************************************************/
/********************    Inputing Data and Allocate memory   ***************************/
/******************************************************************************************/
void Initializing()
{
	int i, j, k;
	int x1, x2;
	float d;
	ifstream FIC;
	FIC.open(File_Name);
	if (FIC.fail())
	{
		cout << "### Erreur open, File_Name " << File_Name << endl;
		exit(0);
	}

	FIC >> N >> K;
	char StrReading[100];
	FIC >> StrReading;
	if (FIC.eof())
	{
		cout << "### Error open, File_Name " << File_Name << endl;
		exit(0);
	}
	if (strcmp(StrReading, "ds") == 0 || strcmp(StrReading, "ss") == 0)
	{
		LB = new double[K];
		UB = new double[K];
		for (i = 0; i<K; i++)
		{
			FIC >> LB[i];
			FIC >> UB[i];
		}
		// for(i=0;i<K;i++)printf("%lf %lf\n",LB[i],UB[i]);
	}
	FIC >> StrReading;
	if (FIC.eof())
	{
		cout << "### Error open, File_Name " << File_Name << endl;
		exit(0);
	}
	if (strcmp(StrReading, "W") == 0)
	{
		w = new double[N];
		for (i = 0; i<N; i++)
			FIC >> w[i];
	}

	Distance = new double *[N];
	for (i = 0; i<N; i++)
		Distance[i] = new double[N];
	DistanceMulti = new double *[N];
	for (i = 0; i<N; i++)
		DistanceMulti[i] = new double[N];

	while (!FIC.eof())
	{
		FIC >> x1 >> x2 >> d;
		//cout << x1 <<"  "<< x2 <<"  "<<d<<" "<< endl;
		if (x1<0 || x2<0 || x1 >= N || x2 >= N)
		{
			cout << "### Error of node : x1="
				<< x1 << ", x2=" << x2 << endl;
			exit(0);
		}
		if (x1 != x2)
		{

			Distance[x2][x1] = d;
			Distance[x1][x2] = Distance[x2][x1];
			DistanceMulti[x2][x1] = 2.0*d;
			DistanceMulti[x1][x2] = DistanceMulti[x2][x1];
		}
	}
	for (i = 0; i<N; i++)
		DistanceMulti[i][i] = Distance[i][i] = 0.0;
	FIC.close();
}

//针对handover算例
void Initializing1()
{
	int i, j, k;
	int x1, x2;
	float d;
	double U;
	ifstream FIC;
	FIC.open(File_Name);
	if (FIC.fail())
	{
		cout << "### Erreur open, File_Name " << File_Name << endl;
		exit(0);
	}
	FIC >> N >> K;
	if (FIC.eof())
	{
		cout << "### Error open, File_Name " << File_Name << endl;
		exit(0);
	}

	LB = new double[K];
	UB = new double[K];

	FIC >> U;
	for (i = 0; i<K; i++){ LB[i] = 0; UB[i] = U; }
	w = new double[N];
	for (i = 0; i<N; i++) FIC >> w[i];

	Distance = new double *[N];
	for (i = 0; i<N; i++) Distance[i] = new double[N];

	DistanceMulti = new double *[N];
	for (i = 0; i<N; i++) DistanceMulti[i] = new double[N];

	int c = 0;
	for (x1 = 0; x1 < N; x1++)
		for (x2 = 0; x2 < N; x2++)
		{
		FIC >> d;
		Distance[x1][x2] = d;
		DistanceMulti[x1][x2] = 2.0*d;
		c++;
		}

	FIC.close();
	/*
	for(i=0;i<N;i++)
	{
	printf("\n");
	for(j=0;j<N;j++) printf("%d ", int (D[i][j]));
	}
	*/
}


void AllocateMemory()
{
	int i, j;
	p = new int[N];
	//bestp = new int [N];
	SizeG = new double[K];

	Delta_Matrix = new double *[N];
	for (int i = 0; i<N; i++)
		Delta_Matrix[i] = new double[K];

	CS.p = new int[N];
	NS.p = new int[N];
	GS.p = new int[N];
	OS.p = new int[N];

	CS.SizeG = new double[K];
	NS.SizeG = new double[K];
	GS.SizeG = new double[K];
	OS.SizeG = new double[K];

	Neighbors = new Neighborhood[N*(N - 1) / 2 + N*K];
	PairSet = new Pair[N*K];
	TabuTenure = new int*[N];
	for (i = 0; i < N; i++)
		TabuTenure[i] = new int[K];

	Pop = new Solution[Num];
	for (int i = 0; i < Num; i++)
	{
		Pop[i].p = new int[N];
		Pop[i].SizeG = new double[K];
	}
	//Cluster1_Node = new int*[K];
	//Cluster2_Node = new int*[K];
	for (int i = 0; i < K; i++)
	{
		//Cluster1_Node[i] = new int[N];
		//Cluster2_Node[i] = new int[N]; 
	}
	Cluster_F[0] = new double[K];
	Cluster_F[1] = new double[K];
	Pop_f = new double[Num];

	Pop_Cluster_Len = new int*[Num];
	Pop_Cluster_Vec = new int**[Num];
	for (int i = 0; i < Num; i++)
	{
		Pop_Cluster_Len[i] = new int[K];
		Pop_Cluster_Vec[i] = new int*[K];
	}


	for (i = 0; i < Num; i++)
	{
		for (j = 0; j < K; j++)
		{
			Pop_Cluster_Vec[i][j] = new int[N];
		}
	}
	//Cluster2_F = new double[K];

}

void ReleaseMemory()
{
	delete[] p; p = NULL;
	//delete [] bestp; bestp = NULL; 
	delete[] SizeG; SizeG = NULL;

	delete[] CS.p; CS.p = NULL;
	delete[] CS.SizeG; CS.SizeG = NULL;
	delete[] GS.p; GS.p = NULL;
	delete[] GS.SizeG; GS.SizeG = NULL;
	delete[] NS.p; NS.p = NULL;
	delete[] NS.SizeG; NS.SizeG = NULL;
	delete[] OS.p; OS.p = NULL;
	delete[] OS.SizeG; OS.SizeG = NULL;

	delete[] LB; LB = NULL;
	delete[] UB; UB = NULL;
	delete[] Neighbors; Neighbors = NULL;

	for (int i = 0; i<N; i++)
	{
		delete[] Delta_Matrix[i]; Delta_Matrix[i] = NULL;
		delete[] Distance[i]; Distance[i] = NULL;
		delete[] DistanceMulti[i]; DistanceMulti[i] = NULL;
	}

}
/******************************************************************************************/
/*********************************    OutPuting Results   ******************************/
/******************************************************************************************/
int Proof(Solution &S)
{
	int i, j;
	double ff;
	int flag;
	ff = 0.0;
	for (i = 0; i < N; i++)
		for (j = i + 1; j < N; j++)					//原始的计算f 的方式
		{
		if (S.p[i] == S.p[j])
		{
			ff += Distance[i][j];
		}
		}
	S.cost = ff;
	for (i = 0; i<K; i++)
		S.SizeG[i] = 0.0;
	for (i = 0; i<N; i++)
		S.SizeG[S.p[i]] += w[i];
	flag = 1;
	for (i = 0; i < K; i++)
	{
		if (S.SizeG[i] < LB[i] - 1.0e-10 || S.SizeG[i]> UB[i] + 1.0e-10)
		{
			flag = 0;
			break;
		}
	}
	return flag;
}

// 输出.sol文件
void Outputing(Solution &S, char *filename)
{
	int i, r;
	FILE *fp;
	char buff[80];
	r = rand() % 1000;
	if (Proof(S) == 0)
		return;
	sprintf(buff, "%s.fval", filename);
	fp = fopen(buff, "a+");
	fprintf(fp, "N = %d  G = %d  f = %lf\n", N, K, S.cost);
	fclose(fp);
	//for (i = 0; i<K; i++)
	//	fprintf(fp, "%lf   %lf   %lf\n", LB[i], UB[i], S.SizeG[i]);
	//printf("\n");

	sprintf(buff, "%s", filename);
	fp = fopen(buff, "a+");
	for (i = 0; i<N; i++)
		fprintf(fp, "%d ", S.p[i]);
	fprintf(fp, "\n");
	fclose(fp);
}


void Out_results(double best, double ave, double worst, double H_best, double H_ave, double H_worst, double AvgTime, char *filename, char instance[])
{
	int i;
	FILE *fp;
	char buff[80];
	sprintf(buff, "%s", filename);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s   %lf   %lf   %lf    %lf  %lf  %lf  %lf\n", instance, best, ave, worst, AvgTime, H_best, H_ave, H_worst);
	fclose(fp);
}

void Out_results1(double fi, char *output_filename, char instance[], double time, double handSum)
{
	int i;
	FILE *fp;
	char buff[80];
	sprintf(buff, "%s", output_filename);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s   %lf %lf %lf  %lf\n", instance, fi, f_best1, time, handSum);
	fclose(fp);
}




// Initial
void RandomInitiaSol(int p[], double SizeG[])
{
	int i, j;
	int p1;
	int count, c1, Nc;
	int *Flag = new int[N];
	double *SizeGroup = new double[K];
	int *G = new int[K];
	int *VN = new int[N];
	for (i = 0; i<K; i++)
		SizeGroup[i] = 0.0;
	for (i = 0; i<N; i++)
		Flag[i] = 0;

	count = 0;
	Nc = 0;
	for (i = 0; i<K; i++)
		if (SizeGroup[i] < LB[i])
			G[count++] = i;
	for (j = 0; j<N; j++)
		if (Flag[j] == 0)
			VN[Nc++] = j;

	while (count > 0)
	{
		while (1)
		{
			p1 = VN[rand() % Nc];			//random pick a node from VN;
			c1 = G[rand() % count];			//randomly pick a cluster form G
			if (SizeGroup[c1] + w[p1] <= UB[c1])
			{
				p[p1] = c1;
				Flag[p1] = 1;
				SizeGroup[c1] += w[p1];
				break;
			}
		}
		count = 0;
		Nc = 0;
		for (i = 0; i<K; i++)
			if (SizeGroup[i] < LB[i])
				G[count++] = i;
		for (j = 0; j<N; j++)
			if (Flag[j] == 0)
				VN[Nc++] = j;
	}
	count = 0;
	Nc = 0;
	for (i = 0; i<K; i++)
		if (SizeGroup[i] < UB[i])
			G[count++] = i;
	for (j = 0; j<N; j++)
		if (Flag[j] == 0)
			VN[Nc++] = j;

	while (Nc > 0)
	{
		while (1)
		{
			p1 = VN[rand() % Nc];
			c1 = G[rand() % count];
			if (SizeGroup[c1] + w[p1] <= UB[c1])
			{
				p[p1] = c1;
				Flag[p1] = 1;
				SizeGroup[c1] += w[p1];
				break;
			}
		}
		count = 0;
		Nc = 0;
		for (i = 0; i<K; i++)
			if (SizeGroup[i] < UB[i])
				G[count++] = i;
		for (j = 0; j<N; j++)
			if (Flag[j] == 0)
				VN[Nc++] = j;
	}

	for (i = 0; i<K; i++)
		SizeG[i] = SizeGroup[i];
	delete[] SizeGroup; SizeGroup = NULL;
	delete[] Flag; Flag = NULL;
	delete[] G; G = NULL;
	delete[] VN; VN = NULL;
	// printf("finish construction \n");
	// for(i=0;i<K;i++)printf("%lf   %lf   %lf\n",LB[i], UB[i], SizeG[i]);
}

//针对handover算例
void RandomInitiaSol1(int p[], double SizeG[])
{
	int i, j;
	int p1;
	int count, c1, Nc;
	int *Flag = new int[N];
	double v_max;
	int j_max;
	double *SizeGroup = new double[K];
	int *G = new int[K];
	//int *VN = new int [N]; 
	for (i = 0; i<K; i++)
		SizeGroup[i] = 0.0;
	for (i = 0; i<N; i++)
		Flag[i] = 0;

	count = 0;
	v_max = -99999.0;
	for (i = 0; i<K; i++)
		if (SizeGroup[i] < LB[i])
			G[count++] = i;
	for (j = 0; j < N; j++)			//找出weights最大的node
	{
		if (Flag[j] == 0 && w[j] > v_max)
		{
			v_max = w[j];
			j_max = j;
		}
	}
	while (count > 0)	//assigned node to cluster one by one, greedy 
	{
		while (1)
		{
			p1 = j_max;
			c1 = G[rand() % count];
			if (SizeGroup[c1] + w[p1] <= UB[c1])
			{
				p[p1] = c1;
				Flag[p1] = 1;				//p1 is assigned to c1
				SizeGroup[c1] += w[p1];
				break;
			}
		}
		count = 0;
		v_max = -99999.0;
		for (i = 0; i<K; i++)
			if (SizeGroup[i] < LB[i])
				G[count++] = i;
		for (j = 0; j < N; j++)
		{
			if (Flag[j] == 0 && w[j] > v_max)
			{
				v_max = w[j];
				j_max = j;
			}
		}
	}
	//  for(i=0;i<K;i++) printf("%lf  %lf ",LB[i], SizeGroup[i]);   
	count = 0;
	v_max = -99999.0;
	Nc = 0;
	for (i = 0; i<K; i++)
		if (SizeGroup[i] < UB[i])
			G[count++] = i;
	for (j = 0; j < N; j++)
	{
		if (Flag[j] == 0 && w[j] > v_max)
		{
			v_max = w[j];
			j_max = j;
		}
	}
	for (j = 0; j<N; j++)
		if (Flag[j] == 0)
			Nc++;						//not assigned node all the same
	while (Nc > 0)
	{
		while (1)
		{
			p1 = j_max;
			c1 = G[rand() % count];
			if (SizeGroup[c1] + w[p1] <= UB[c1])
			{
				p[p1] = c1;
				Flag[p1] = 1;
				SizeGroup[c1] += w[p1];
				break;
			}
			// for(j=0;j<K;j++) printf("%lf  ",SizeGroup[j]);  printf("%d",count);
			// for(j=0;j<N;j++) if(Flag[j]==0) printf(" %lf ", w[j]); printf("\n");
		}
		count = 0;
		v_max = -99999.0;
		Nc = 0;
		for (i = 0; i<K; i++)
			if (SizeGroup[i] < UB[i])
				G[count++] = i;
		for (j = 0; j < N; j++)
		{
			if (Flag[j] == 0 && w[j] > v_max)
			{
				v_max = w[j];
				j_max = j;
			}
		}
		for (j = 0; j<N; j++)
			if (Flag[j] == 0) Nc++;
	}

	for (i = 0; i<K; i++)
		SizeG[i] = SizeGroup[i];
	delete[] SizeGroup; SizeGroup = NULL;
	delete[] Flag; Flag = NULL;
	delete[] G; G = NULL;
	//delete [] VN; VN = NULL; 
	// printf("\n finish construction \n");
	// for(i=0;i<K;i++)printf("%lf   %lf   %lf\n",LB[i], UB[i], SizeG[i]);
}

double ComputeF(int p1[])
{
	double sum = 0.0;
	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			if (p1[i] == p1[j])
				sum += Distance[i][j];
		}
	}
	//cout << "sum=" << sum << endl;
	return sum;
}

int VerifySolution()
{
	int i, j;
	int flag = 0;
	for (i = 0; i < K; i++)
		CS.SizeG[i] = 0;
	for (i = 0; i < N; i++)
		if (CS.p[i] != -1)
			CS.SizeG[CS.p[i]] += w[i];
	for (i = 0; i < K; i++)
	{
		if (CS.SizeG[i] > UB[i] || CS.SizeG[i] < LB[i])
		{
			//cout << "i=" << i << " size=" << CS.SizeG[i] << "  not proper cluster" << endl;
			flag = 1;
		}
		//else
		//	cout << "i=" << i << " size=" << CS.SizeG[i] << "  proper cluster " << endl;
	}
	return flag;
}

void BuildNeighbors()
{
	int i, j, g;
	int count;
	int SN = N*(N - 1) / 2 + N*K;
	count = 0;
	for (i = 0; i<N; i++)
		for (g = 0; g<K; g++)
		{
		Neighbors[count].type = 1;
		Neighbors[count].v = i;
		Neighbors[count].g = g;
		count++;
		}
	for (i = 0; i<N; i++)
		for (j = i + 1; j<N; j++)
		{
		Neighbors[count].type = 2;
		Neighbors[count].x = i;
		Neighbors[count].y = j;
		count++;
		}
}

// Clear delta matrix
void Clear_Delta_Matrix()
{
	int x, g;
	f = 0.0;
	for (x = 0; x < N; x++)
		for (g = 0; g < K; g++)
			Delta_Matrix[x][g] = 0.0;
	for (x = 0; x < N; x++)
		for (g = 0; g < K; g++)
			TabuTenure[x][g] = 0;

}

// Build delta matrix
void Build_Delta_Matrix()
{
	int i, j;
	Clear_Delta_Matrix();
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			Delta_Matrix[i][p[j]] += Distance[i][j];
	f = 0.0;
	for (i = 0; i < N; i++)
		f += Delta_Matrix[i][p[i]];
	f = f / 2;
}


void Build_Delta_Matrix_Repair(int p1[], double &ff)
{
	int i, j, x, g;
	//Clear_Delta_Matrix();

	for (x = 0; x < N; x++)
		for (g = 0; g < K; g++)
			Delta_Matrix[x][g] = 0.0;

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			Delta_Matrix[i][p1[j]] += Distance[i][j];

	ff = 0.0;
	for (i = 0; i < N; i++)
		ff += Delta_Matrix[i][p1[i]];
	ff = ff / 2;
}

// Update one move delta matrix
void Update_Delta_Matrix(int i, int g0, int g1)
{
	int x, j, k;

	for (j = 0; j<N; j++)
	{
		if (j != i)
		{
			Delta_Matrix[j][g0] -= Distance[i][j];
			Delta_Matrix[j][g1] += Distance[i][j];
		}
	}

}
int num_tabu_best_1 = 0, num_best_1 = 0;												 // the number of tabu neighbors and non-tabu neighbors
int num_tabu_best_2 = 0, num_best_2 = 0;
int num_tabu_best_3 = 0, num_best_3 = 0;
int num_tabu_best_4 = 0, num_best_4 = 0;
int best_x_1[50], best_c_1[50];
int best_x_2[50], best_y_2[50];
int best_x_3[50], best_y_3[50], best_z_3[50];

int tabu_best_x_1[50], tabu_best_c_1[50];
int tabu_best_x_2[50], tabu_best_y_2[50];
int tabu_best_x_3[50], tabu_best_y_3[50], tabu_best_z_3[50];

double tabu_best_delta_1, best_delta_1;
double tabu_best_delta_2, best_delta_2;
double tabu_best_delta_3, best_delta_3;

double tabu_best_delta_111, best_delta_111;
double tabu_best_delta_222, best_delta_222;
double tabu_best_delta_333, best_delta_333;


//N1:one-Move
void Swap_Move_Tabu_Search_N1(int p1[], double SizeG[], int iter)
{
	int x, j;
	double delt;

	for (x = 0; x < N; x++)
	{
		for (j = 0; j < K; j++)
		{
			if ((p1[x] != j) && (SizeG[p1[x]] - w[x] >= LB[p1[x]]) && (SizeG[j] + w[x] <= UB[j]))
			{
				delt = Delta_Matrix[x][j] - Delta_Matrix[x][p1[x]];
				if (TabuTenure[x][j] <= iter)
				{
					if (delt > best_delta_1)
					{
						best_x_1[0] = x;
						best_c_1[0] = j;
						best_delta_1 = delt;
						num_best_1 = 1;
					}
					else if (delt == best_delta_1 && num_best_1 < 50)
					{
						best_x_1[num_best_1] = x;
						best_c_1[num_best_1] = j;
						num_best_1++;
					}
				}
				else if (TabuTenure[x][j] > iter)
				{
					if (delt > tabu_best_delta_1)
					{
						tabu_best_x_1[0] = x;
						tabu_best_c_1[0] = j;
						tabu_best_delta_1 = delt;
						num_tabu_best_1 = 1;
					}
					else if (delt == tabu_best_delta_1 && num_tabu_best_1 < 50)
					{
						tabu_best_x_1[num_tabu_best_1] = x;
						tabu_best_c_1[num_tabu_best_1] = j;
						num_tabu_best_1++;
					}
				}
				//printf("delt =%lf\n",delt);
			}
		}
	}

}

//N2:Swap-Move tabu search
void Swap_Move_Tabu_Search_N2(int p1[], double SizeG[], int iter)
{
	int x, y;
	double delt;

	for (x = 0; x < N; x++)
	{
		for (y = x + 1; y < N; y++)
		{
			if ((p1[x] != p1[y]) && (SizeG[p1[x]] + w[y] - w[x] >= LB[p1[x]])
				&& (SizeG[p1[x]] + w[y] - w[x] <= UB[p1[x]])
				&& (SizeG[p1[y]] + w[x] - w[y] >= LB[p1[y]])
				&& (SizeG[p1[y]] + w[x] - w[y] <= UB[p1[y]]))
			{
				// delt = Delta[ x ][ p1[y] ] + Delta[ y ][ p1[x] ] - DistanceMulti[x][y];   
				delt = (Delta_Matrix[x][p1[y]] - Delta_Matrix[x][p1[x]]) + (Delta_Matrix[y][p1[x]]
					- Delta_Matrix[y][p1[y]]) - DistanceMulti[x][y];

				if ((TabuTenure[x][p1[y]] <= iter) && (TabuTenure[y][p1[x]] <= iter))				 // if this is not tabued 
				{
					if (delt > best_delta_2)
					{
						best_x_2[0] = x;
						best_y_2[0] = y;

						best_delta_2 = delt;
						num_best_2 = 1;
					}
					else if (delt == best_delta_2 && num_best_2 < 50)
					{
						best_x_2[num_best_2] = x;
						best_y_2[num_best_2] = y;

						num_best_2++;
					}
				}
				else if ((TabuTenure[x][p1[y]] > iter) || (TabuTenure[y][p1[x]] > iter))			// if it is tabu 
				{
					if (delt > tabu_best_delta_2)
					{
						tabu_best_x_2[0] = x;
						tabu_best_y_2[0] = y;

						tabu_best_delta_2 = delt;
						num_tabu_best_2 = 1;
					}
					else if (delt == tabu_best_delta_2 && num_tabu_best_2 < 50)
					{
						tabu_best_x_2[num_tabu_best_2] = x;
						tabu_best_y_2[num_tabu_best_2] = y;

						num_tabu_best_2++;
					}
				}

			}
		}
	}
}

//N3:2-1 Exchange Move tabu_search
void Swap_Move_Tabu_Search_N3(int p1[], double SizeG[], int iter)
{
	double delt;
	int x, y, z;

	for (x = 0; x < N; x++)
	{
		for (y = x + 1; y < N; y++)
		{
			if (p1[y] == p1[x])								//添加的判断条件
			{
				for (z = 0; z < N; z++)
				{
					if ((p1[x] == p1[y] && p1[x] != p1[z])
						&& (SizeG[p1[x]] - w[y] - w[x] + w[z] >= LB[p1[x]])
						&& (SizeG[p1[x]] - w[y] - w[x] + w[z] <= UB[p1[x]])
						&& (SizeG[p1[z]] + w[x] + w[y] - w[z] >= LB[p1[z]])
						&& (SizeG[p1[z]] + w[x] + w[y] - w[z] <= UB[p1[z]]))
					{
						// delt = (Delta[ x ][ p[z] ] + Delta[ z ][ p[x] ] - DistanceMulti[x][z]) + Delta[ y ][ p[z] ] + DistanceMulti[y][x]- DistanceMulti[y][z];  
						delt = (Delta_Matrix[x][p1[z]] - Delta_Matrix[x][p1[x]]) +
							(Delta_Matrix[z][p1[x]] - Delta_Matrix[z][p1[z]]) - DistanceMulti[x][z];
						delt += ((Delta_Matrix[y][p1[z]] - Delta_Matrix[y][p1[y]]) + (DistanceMulti[y][x] - DistanceMulti[y][z]));

						if ((TabuTenure[x][p1[z]] <= iter) && (TabuTenure[y][p1[z]] <= iter)
							&& (TabuTenure[z][p1[x]] <= iter))											 // if this is not tabued 
						{
							if (delt > best_delta_3)
							{
								best_x_3[0] = x;
								best_y_3[0] = y;
								best_z_3[0] = z;
								best_delta_3 = delt;
								num_best_3 = 1;
							}
							else if (delt == best_delta_3 && num_best_3 < 50)
							{
								best_x_3[num_best_3] = x;
								best_y_3[num_best_3] = y;
								best_z_3[num_best_3] = z;
								num_best_3++;
							}
						}
						else if ((TabuTenure[x][p1[z]] > iter) || (TabuTenure[y][p1[z]] > iter)
							|| (TabuTenure[z][p1[x]] > iter))												// if it is tabu 
						{
							if (delt > tabu_best_delta_3)
							{
								tabu_best_x_3[0] = x;
								tabu_best_y_3[0] = y;
								tabu_best_z_3[0] = z;
								tabu_best_delta_3 = delt;
								num_tabu_best_3 = 1;
							}
							else if (delt == tabu_best_delta_3 && num_tabu_best_3 < 50)
							{
								tabu_best_x_3[num_tabu_best_3] = x;
								tabu_best_y_3[num_tabu_best_3] = y;
								tabu_best_z_3[num_tabu_best_3] = z;
								num_tabu_best_3++;
							}
						}
					}
				}
			}
		}
	}
}

//根据flag_type 进行顶点移动
void ReallyMove(int p1[], int flag_type, double &ff, double best_delta_global, double tabu_best_delta_global, int iter, double SizeG[])
{
	int select, cluster, old_group, old_group1, old_group2;
	switch (flag_type)
	{
	case 1:
		ff += best_delta_global;
		select = rand() % num_best_1;
		old_group = p1[best_x_1[select]];
		cluster = best_c_1[select];

		SizeG[old_group] -= w[best_x_1[select]];
		SizeG[cluster] += w[best_x_1[select]];
		Update_Delta_Matrix(best_x_1[select], old_group, cluster);
		p1[best_x_1[select]] = cluster;
		TabuTenure[best_x_1[select]][old_group] = TabuLen + iter;
		break;
	case 2:
		ff += tabu_best_delta_global;
		select = rand() % num_tabu_best_1;
		old_group = p1[tabu_best_x_1[select]];
		cluster = tabu_best_c_1[select];

		SizeG[old_group] -= w[tabu_best_x_1[select]];
		SizeG[cluster] += w[tabu_best_x_1[select]];
		Update_Delta_Matrix(tabu_best_x_1[select], old_group, cluster);
		p1[tabu_best_x_1[select]] = cluster;
		TabuTenure[tabu_best_x_1[select]][old_group] = TabuLen + iter;
		break;
	case 3:
		ff += best_delta_global;
		select = rand() % num_best_2;
		old_group = p1[best_x_2[select]];
		old_group1 = p1[best_y_2[select]];
		//cout << "delta3=" << Delta_Matrix[best_x_2[select]][]
		SizeG[old_group] += (w[best_y_2[select]] - w[best_x_2[select]]);
		SizeG[old_group1] += (w[best_x_2[select]] - w[best_y_2[select]]);
		Update_Delta_Matrix(best_x_2[select], old_group, old_group1);

		p1[best_x_2[select]] = p1[best_y_2[select]];
		Update_Delta_Matrix(best_y_2[select], old_group1, old_group);
		p1[best_y_2[select]] = old_group;
		TabuTenure[best_x_2[select]][old_group] = TabuLen + iter;
		TabuTenure[best_y_2[select]][old_group1] = TabuLen + iter;
		break;
	case 4:
		ff += tabu_best_delta_global;
		select = rand() % num_tabu_best_2;
		old_group = p1[tabu_best_x_2[select]];
		old_group1 = p1[tabu_best_y_2[select]];
		SizeG[old_group] += (w[tabu_best_y_2[select]] - w[tabu_best_x_2[select]]);
		SizeG[old_group1] += (w[tabu_best_x_2[select]] - w[tabu_best_y_2[select]]);
		Update_Delta_Matrix(tabu_best_x_2[select], old_group, old_group1);

		p1[tabu_best_x_2[select]] = p1[tabu_best_y_2[select]];
		Update_Delta_Matrix(tabu_best_y_2[select], old_group1, old_group);
		p1[tabu_best_y_2[select]] = old_group;
		TabuTenure[tabu_best_x_2[select]][old_group] = TabuLen + iter;
		TabuTenure[tabu_best_y_2[select]][old_group1] = TabuLen + iter;
		break;
	case 5:
		ff += best_delta_global;
		select = rand() % num_best_3;

		old_group = p1[best_x_3[select]];
		old_group1 = p1[best_y_3[select]];
		old_group2 = p1[best_z_3[select]];

		SizeG[old_group] += (w[best_z_3[select]] - w[best_y_3[select]] - w[best_x_3[select]]);
		SizeG[old_group2] += (w[best_x_3[select]] + w[best_y_3[select]] - w[best_z_3[select]]);

		Update_Delta_Matrix(best_x_3[select], old_group, old_group2);

		p1[best_x_3[select]] = p1[best_z_3[select]];
		Update_Delta_Matrix(best_z_3[select], old_group2, old_group);
		p1[best_z_3[select]] = old_group;

		Update_Delta_Matrix(best_y_3[select], old_group1, old_group2);
		p1[best_y_3[select]] = old_group2;

		TabuTenure[best_x_3[select]][old_group] = TabuLen + iter;
		TabuTenure[best_y_3[select]][old_group1] = TabuLen + iter;
		TabuTenure[best_z_3[select]][old_group2] = TabuLen + iter;
		break;
	case 6:
		ff += tabu_best_delta_global;
		select = rand() % num_tabu_best_3;

		old_group = p1[tabu_best_x_3[select]];
		old_group1 = p1[tabu_best_y_3[select]];
		old_group2 = p1[tabu_best_z_3[select]];

		SizeG[old_group] += (w[tabu_best_z_3[select]] - w[tabu_best_y_3[select]] - w[tabu_best_x_3[select]]);
		SizeG[old_group2] += (w[tabu_best_x_3[select]] + w[tabu_best_y_3[select]] - w[tabu_best_z_3[select]]);

		Update_Delta_Matrix(tabu_best_x_3[select], old_group, old_group2);

		p1[tabu_best_x_3[select]] = p1[tabu_best_z_3[select]];
		Update_Delta_Matrix(tabu_best_z_3[select], old_group2, old_group);
		p1[tabu_best_z_3[select]] = old_group;

		Update_Delta_Matrix(tabu_best_y_3[select], old_group1, old_group2);
		p1[tabu_best_y_3[select]] = old_group2;

		TabuTenure[tabu_best_x_3[select]][old_group] = TabuLen + iter;
		TabuTenure[tabu_best_y_3[select]][old_group1] = TabuLen + iter;
		TabuTenure[tabu_best_z_3[select]][old_group2] = TabuLen + iter;
		break;
	default:
		break;
	}
}

int CheckMove(int p1[], double ff, double &sum1)
{
	int i, j;
	double sum = 0;
	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			if (p1[i] == p1[j])
				sum += Distance[i][j];
		}
	}
	sum1 = sum;
	if (fabs(sum - ff)< 1.0e-10)
		return 0;
	else return 1;

}

//N1:one-Move
void Swap_Move_Tabu_Search_N1_No_Restrict(int p1[], double SizeG[], int iter)
{
	int x, j;
	double delt, delt_fs, delt1, delt2;
	double a1, a2;
	double a1_prev, a2_prev;
	for (x = 0; x < N; x++)
	{
		for (j = 0; j < K; j++)
		{
			if ((p1[x] != j))
			{
				//delt = Delta_Matrix[x][j] - Delta_Matrix[x][p1[x]];
				if (SizeG[j] + w[x]>UB[j])
					a1 = (SizeG[j] + w[x]) - UB[j];
				else if (SizeG[j] + w[x] < LB[j])
					a1 = LB[j] - (SizeG[j] + w[x]);
				else
					a1 = 0;

				if (SizeG[p1[x]] - w[x]>UB[p1[x]])
					a2 = (SizeG[p1[x]] - w[x]) - UB[p1[x]];
				else if (SizeG[p1[x]] - w[x] < LB[p1[x]])
					a2 = LB[p1[x]] - (SizeG[p1[x]] - w[x]);
				else
					a2 = 0;

				if (SizeG[j] > UB[j])
					a1_prev = SizeG[j] - UB[j];
				else if (SizeG[j]<LB[j])
					a1_prev = LB[j] - SizeG[j];
				else
					a1_prev = 0;
				if (SizeG[p1[x]]>UB[p1[x]])
					a2_prev = SizeG[p1[x]] - UB[p1[x]];
				else if (SizeG[p1[x]] < LB[p1[x]])
					a2_prev = LB[p1[x]] - SizeG[p1[x]];
				else
					a2_prev = 0;

				delt1 = Delta_Matrix[x][j] - Delta_Matrix[x][p1[x]];
				delt2 = (a1 + a2) - (a1_prev + a2_prev);
				//delt = Delta_Matrix[x][j] - Delta_Matrix[x][p1[x]] - penalty_factor*((a1+a2)-(a1_prev+a2_prev));
				delt = delt1 - penalty_factor*delt2;
				if (TabuTenure[x][j] <= iter)
				{
					if (delt > best_delta_1)
					{
						best_x_1[0] = x;
						best_c_1[0] = j;
						best_delta_1 = delt;
						best_delta_111 = delt1;
						//delt222 = delt2;
						//best_delt_fs1 = delt_fs;
						num_best_1 = 1;
					}
					else if (delt == best_delta_1 && num_best_1 < 50)
					{
						best_x_1[num_best_1] = x;
						best_c_1[num_best_1] = j;
						num_best_1++;
					}
				}
				else if (TabuTenure[x][j] > iter)
				{
					if (delt > tabu_best_delta_1)
					{
						tabu_best_x_1[0] = x;
						tabu_best_c_1[0] = j;
						tabu_best_delta_1 = delt;
						tabu_best_delta_111 = delt1;
						//tabu_best_delt_fs1 = delt;
						num_tabu_best_1 = 1;
					}
					else if (delt == tabu_best_delta_1 && num_tabu_best_1 < 50)
					{
						tabu_best_x_1[num_tabu_best_1] = x;
						tabu_best_c_1[num_tabu_best_1] = j;
						num_tabu_best_1++;
					}
				}
				//printf("delt =%lf\n",delt);
			}
		}
	}

}

//N2:Swap-Move tabu search
void Swap_Move_Tabu_Search_N2_No_Restrict(int p1[], double SizeG[], int iter)
{
	int x, y;
	double delt, delt1, delt2;
	double a1, a2;
	double a1_prev, a2_prev;

	for (x = 0; x < N; x++)
	{
		for (y = x + 1; y < N; y++)
		{
			if ((p1[x] != p1[y]))
			{

				if (SizeG[p1[y]] + w[x] - w[y] > UB[p1[y]])
					a1 = SizeG[p1[y]] + w[x] - w[y] - UB[p1[y]];
				else if (SizeG[p1[y]] + w[x] - w[y] < LB[p1[y]])
					a1 = LB[p1[y]] - (SizeG[p1[y]] + w[x] - w[y]);
				else
					a1 = 0;

				if (SizeG[p1[x]] + w[y] - w[x] > UB[p1[x]])
					a2 = SizeG[p1[x]] + w[y] - w[x] - UB[p1[x]];
				else if (SizeG[p1[x]] + w[y] - w[x]  < LB[p1[x]])
					a2 = LB[p1[x]] - (SizeG[p1[x]] + w[y] - w[x]);
				else
					a2 = 0;

				if (SizeG[p1[y]] > UB[p1[y]])
					a1_prev = SizeG[p1[y]] - UB[p1[y]];
				else if (SizeG[p1[y]] < LB[p1[y]])
					a1_prev = LB[p1[y]] - SizeG[p1[y]];
				else
					a1_prev = 0;
				if (SizeG[p1[x]] > UB[p1[x]])
					a2_prev = SizeG[p1[x]] - UB[p1[x]];
				else if (SizeG[p1[x]] < LB[p1[x]])
					a2_prev = LB[p1[x]] - SizeG[p1[x]];
				else
					a2_prev = 0;


				delt1 = (Delta_Matrix[x][p1[y]] - Delta_Matrix[x][p1[x]]) + (Delta_Matrix[y][p1[x]]
					- Delta_Matrix[y][p1[y]]) - DistanceMulti[x][y];
				delt2 = (a1 + a2) - (a1_prev + a2_prev);

				//delt = (Delta_Matrix[x][p1[y]] - Delta_Matrix[x][p1[x]]) + (Delta_Matrix[y][p1[x]]
				//	- Delta_Matrix[y][p1[y]]) - DistanceMulti[x][y] - penalty_factor*((a1+a2) - (a1_prev + a2_prev));

				delt = delt1 - penalty_factor*delt2;
				if ((TabuTenure[x][p1[y]] <= iter) && (TabuTenure[y][p1[x]] <= iter))				 // if this is not tabued 
				{
					if (delt > best_delta_2)
					{
						best_x_2[0] = x;
						best_y_2[0] = y;

						best_delta_2 = delt;
						best_delta_222 = delt1;
						num_best_2 = 1;
					}
					else if (delt == best_delta_2 && num_best_2 < 50)
					{
						best_x_2[num_best_2] = x;
						best_y_2[num_best_2] = y;
						num_best_2++;
					}
				}
				else if ((TabuTenure[x][p1[y]] > iter) || (TabuTenure[y][p1[x]] > iter))			// if it is tabu 
				{
					if (delt > tabu_best_delta_2)
					{
						tabu_best_x_2[0] = x;
						tabu_best_y_2[0] = y;

						tabu_best_delta_2 = delt;
						tabu_best_delta_222 = delt1;
						num_tabu_best_2 = 1;
					}
					else if (delt == tabu_best_delta_2 && num_tabu_best_2 < 50)
					{
						tabu_best_x_2[num_tabu_best_2] = x;
						tabu_best_y_2[num_tabu_best_2] = y;

						num_tabu_best_2++;
					}
				}

			}
		}
	}
}

//N3:2-1 Exchange Move tabu_search
void Swap_Move_Tabu_Search_N3_No_Restrict(int p1[], double SizeG[], int iter)
{
	double delt, delt1, delt2;
	int x, y, z;

	double a1, a2;
	double a1_prev, a2_prev;

	for (x = 0; x < N; x++)
	{
		for (y = x + 1; y < N; y++)
		{
			if (p1[y] == p1[x])								//添加的判断条件
			{
				for (z = 0; z < N; z++)
				{
					if ((p1[x] == p1[y] && p1[x] != p1[z]))
					{

						if (SizeG[p1[z]] + w[x] + w[y] - w[z] > UB[p1[z]])
							a1 = SizeG[p1[z]] + w[x] + w[y] - w[z] - UB[p1[z]];
						else if (SizeG[p1[z]] + w[x] + w[y] - w[z] < LB[p1[z]])
							a1 = LB[p1[z]] - (SizeG[p1[z]] + w[x] + w[y] - w[z]);
						else
							a1 = 0;

						if (SizeG[p1[x]] + w[z] - w[x] - w[y] > UB[p1[x]])
							a2 = SizeG[p1[x]] + w[z] - w[x] - w[y] - UB[p1[x]];
						else if (SizeG[p1[x]] + w[z] - w[x] - w[y]  < LB[p1[x]])
							a2 = LB[p1[x]] - (SizeG[p1[x]] + w[z] - w[x] - w[y]);
						else
							a2 = 0;

						if (SizeG[p1[z]] > UB[p1[z]])
							a1_prev = SizeG[p1[z]] - UB[p1[z]];
						else if (SizeG[p1[z]] < LB[p1[z]])
							a1_prev = LB[p1[z]] - SizeG[p1[z]];
						else
							a1_prev = 0;
						if (SizeG[p1[x]] > UB[p1[x]])
							a2_prev = SizeG[p1[x]] - UB[p1[x]];
						else if (SizeG[p1[x]] < LB[p1[x]])
							a2_prev = LB[p1[x]] - SizeG[p1[x]];
						else
							a2_prev = 0;

						// delt = (Delta[ x ][ p[z] ] + Delta[ z ][ p[x] ] - DistanceMulti[x][z]) + Delta[ y ][ p[z] ] + DistanceMulti[y][x]- DistanceMulti[y][z];  
						delt1 = (Delta_Matrix[x][p1[z]] - Delta_Matrix[x][p1[x]]) +
							(Delta_Matrix[z][p1[x]] - Delta_Matrix[z][p1[z]]) - DistanceMulti[x][z];
						delt1 += ((Delta_Matrix[y][p1[z]] - Delta_Matrix[y][p1[y]]) + (DistanceMulti[y][x] - DistanceMulti[y][z]));

						delt2 = (a1 + a2) - (a1_prev + a2_prev);

						//delt = (Delta_Matrix[x][p1[z]] - Delta_Matrix[x][p1[x]]) +
						//	(Delta_Matrix[z][p1[x]] - Delta_Matrix[z][p1[z]]) - DistanceMulti[x][z];
						//delt += ((Delta_Matrix[y][p1[z]] - Delta_Matrix[y][p1[y]]) + (DistanceMulti[y][x] - DistanceMulti[y][z]));

						delt = delt1 - penalty_factor*delt2;
						if ((TabuTenure[x][p1[z]] <= iter) && (TabuTenure[y][p1[z]] <= iter)
							&& (TabuTenure[z][p1[x]] <= iter))											 // if this is not tabued 
						{
							if (delt > best_delta_3)
							{
								best_x_3[0] = x;
								best_y_3[0] = y;
								best_z_3[0] = z;
								best_delta_3 = delt;
								best_delta_333 = delt1;
								num_best_3 = 1;
							}
							else if (delt == best_delta_3 && num_best_3 < 50)
							{
								best_x_3[num_best_3] = x;
								best_y_3[num_best_3] = y;
								best_z_3[num_best_3] = z;
								num_best_3++;
							}
						}
						else if ((TabuTenure[x][p1[z]] > iter) || (TabuTenure[y][p1[z]] > iter)
							|| (TabuTenure[z][p1[x]] > iter))												// if it is tabu 
						{
							if (delt > tabu_best_delta_3)
							{
								tabu_best_x_3[0] = x;
								tabu_best_y_3[0] = y;
								tabu_best_z_3[0] = z;
								tabu_best_delta_3 = delt;
								tabu_best_delta_333 = delt1;
								num_tabu_best_3 = 1;
							}
							else if (delt == tabu_best_delta_3 && num_tabu_best_3 < 50)
							{
								tabu_best_x_3[num_tabu_best_3] = x;
								tabu_best_y_3[num_tabu_best_3] = y;
								tabu_best_z_3[num_tabu_best_3] = z;
								num_tabu_best_3++;
							}
						}
					}
				}
			}
		}
	}
}

double ComputeExceedS(int p1[], double SizeG[])
{
	int x, g;
	double sum = 0;
	for (g = 0; g < K; g++)
	{
		if (SizeG[g] < LB[g])
			sum += fabs(SizeG[g] - LB[g]);
		else if (SizeG[g]>UB[g])
			sum += fabs(SizeG[g] - UB[g]);
	}
	gs = sum;
	return sum;
}


int IsProperSol(int p1[], double SizeG[])
{
	int i, j;
	int flag = 1;
	for (i = 0; i < K; i++)
	{
		if (SizeG[i] > UB[i] || SizeG[i] < LB[i])
		{
			//cout << "i=" << i << " size=" << CS.SizeG[i] << "  not proper cluster" << endl;
			flag = 0;
		}
	}
	return flag;
}

int non_improve;

//feasible local search procedure
void Feasible_Local_Search(int p1[], double SizeG[], int &iter, double &ff)
{
	double best_delta_global, tabu_best_delta_global;
	int flag_type;
	double sum;
	if (CheckMove(p1, ff, sum))
	{
		cout << "1111  in infeasible local search sum !=ff" << " ff=" << ff << " sum=" << sum << endl;
		getchar();
	}
	/*
	if (!IsProperSol(p1, SizeG))
	{
	cout << "is not proper" << endl;
	getchar();
	}*/
	best_delta_global = -999999;
	tabu_best_delta_global = -999999;
	best_delta_1 = -999999;
	tabu_best_delta_1 = -999999;
	best_delta_2 = -999999;
	tabu_best_delta_2 = -999999;
	best_delta_3 = -999999;
	tabu_best_delta_3 = -999999;

	num_tabu_best_1 = 0;
	num_best_1 = 0;
	num_tabu_best_2 = 0;
	num_best_2 = 0;
	num_tabu_best_3 = 0;
	num_best_3 = 0;

	Swap_Move_Tabu_Search_N1(p1, SizeG, iter);
	Swap_Move_Tabu_Search_N2(p1, SizeG, iter);
	Swap_Move_Tabu_Search_N3(p1, SizeG, iter);

	if (best_delta_1 > best_delta_global)
	{
		best_delta_global = best_delta_1;
		flag_type = 1;
	}
	if (best_delta_2 > best_delta_global)
	{
		best_delta_global = best_delta_2;
		flag_type = 3;
	}
	if (best_delta_3 > best_delta_global)
	{
		best_delta_global = best_delta_3;
		flag_type = 5;
	}

	if (num_tabu_best_1 > 0 && tabu_best_delta_1 > best_delta_global && tabu_best_delta_1 > tabu_best_delta_global &&ff + tabu_best_delta_1 > f_best)
	{
		tabu_best_delta_global = tabu_best_delta_1;
		flag_type = 2;
	}

	if (num_tabu_best_2 > 0 && tabu_best_delta_2 > best_delta_global && tabu_best_delta_2 > tabu_best_delta_global &&ff + tabu_best_delta_2 > f_best)
	{
		tabu_best_delta_global = tabu_best_delta_2;
		flag_type = 4;
	}

	if (num_tabu_best_3 > 0 && tabu_best_delta_3 > best_delta_global && tabu_best_delta_3 > tabu_best_delta_global &&ff + tabu_best_delta_3 > f_best)
	{
		tabu_best_delta_global = tabu_best_delta_3;
		flag_type = 6;
	}
	ReallyMove(p1, flag_type, ff, best_delta_global, tabu_best_delta_global, iter, SizeG);
	if (CheckMove(p1, ff, sum))
	{
		cout << "in feasible local search sum !=ff" << " ff=" << ff << " sum=" << sum << endl;
		getchar();
	}
	if (IsProperSol(p1, SizeG))
	{
		if (ff > f_best)
		{
			f_best = ff;
			total_time = (clock() - starting_time) / CLOCKS_PER_SEC;
			for (int i = 0; i < N; i++)
			{
				GS.p[i] = p1[i];
			}
			non_improve = 0;
		}
		else
			non_improve++;

	}
	iter++;
}

//infeasible local search procedure
void Infeasible_Local_Search(int p1[], double SizeG[], int &iter, double &ff)
{

	double best_delta_global, tabu_best_delta_global;
	double best_fff1 = -999999, tabu_best_fff1 = -999999;
	double f1 = -999999;
	int *p_temp = new int[N];
	double *sz = new double[K];
	int flag_type, flag_f1 = 0;
	int iter_ils = 0;
	double sum;
	penalty_count = 0;
	//penalty_factor = 2;
	memcpy(p_temp, p1, sizeof(int)*N);
	while (iter_ils < 200)
	{

		best_delta_global = -999999;
		tabu_best_delta_global = -999999;
		best_delta_1 = -999999;
		tabu_best_delta_1 = -999999;
		best_delta_2 = -999999;
		tabu_best_delta_2 = -999999;
		best_delta_3 = -999999;
		tabu_best_delta_3 = -999999;
		best_delta_111 = -999999;
		tabu_best_delta_111 = -999999;

		num_tabu_best_1 = 0;
		num_best_1 = 0;
		num_tabu_best_2 = 0;
		num_best_2 = 0;
		num_tabu_best_3 = 0;
		num_best_3 = 0;
		num_tabu_best_4 = 0;
		num_best_4 = 0;

		Swap_Move_Tabu_Search_N1_No_Restrict(p1, SizeG, iter);
		Swap_Move_Tabu_Search_N2_No_Restrict(p1, SizeG, iter);
		Swap_Move_Tabu_Search_N3_No_Restrict(p1, SizeG, iter);

		if (best_delta_1 > best_delta_global)
		{
			best_delta_global = best_delta_1;
			best_fff1 = best_delta_111;

			flag_type = 1;
		}
		if (best_delta_2 > best_delta_global)
		{
			best_delta_global = best_delta_2;
			best_fff1 = best_delta_222;
			flag_type = 3;
		}
		if (best_delta_3 > best_delta_global)
		{
			best_delta_global = best_delta_3;
			best_fff1 = best_delta_333;
			flag_type = 5;
		}

		if (num_tabu_best_1 > 0 && tabu_best_delta_1 > best_delta_global && tabu_best_delta_1 > tabu_best_delta_global &&ys + tabu_best_delta_1 > ys_best)
		{
			tabu_best_delta_global = tabu_best_delta_1;
			tabu_best_fff1 = tabu_best_delta_111;
			flag_type = 2;
		}

		if (num_tabu_best_2 > 0 && tabu_best_delta_2 > best_delta_global && tabu_best_delta_2 > tabu_best_delta_global &&ys + tabu_best_delta_2 > ys_best)
		{
			tabu_best_delta_global = tabu_best_delta_2;
			tabu_best_fff1 = tabu_best_delta_222;
			flag_type = 4;
		}

		if (num_tabu_best_3 > 0 && tabu_best_delta_3 > best_delta_global && tabu_best_delta_3 > tabu_best_delta_global &&ys + tabu_best_delta_3 > ys_best)
		{
			tabu_best_delta_global = tabu_best_delta_3;
			tabu_best_fff1 = tabu_best_delta_333;
			flag_type = 6;
		}
		//cout << "ff=" << ff << " gs="<<gs<< " ComputeF="<<ComputeF(p1)<<" flag_type="<<flag_type<<endl;
		ReallyMove(p1, flag_type, ff, best_fff1, tabu_best_fff1, iter, SizeG);
		if (CheckMove(p1, ff, sum))
		{

			cout << "here 1 in ils sum !=ff" << " ff=" << ff << " sum=" << sum << " best_fff1=" << best_fff1 << " tabu_bff1=" << tabu_best_fff1 << " CompF=" << ComputeF(p1) << endl;
			getchar();
		}

		gs = ComputeExceedS(p1, SizeG);
		ys = ff - penalty_factor*gs;
		if (ys > ys_best)
			ys_best = ys;

		if (!IsProperSol(p1, SizeG))
		{
			penalty_count++;
		}
		iter_ils++;
		//iter += iter_ils;
		iter++;
		//cout << "iter_ils=" << iter_ils  << " iter="<<iter<<endl;
		//printf("iter_ils=%d iter=%d\n", iter_ils, iter);
		if (IsProperSol(p1, SizeG))
		{
			if (ff > f1)
			{
				f1 = ff;
				flag_f1 = 1;
				memcpy(p_temp, p1, sizeof(int)*N);
				memcpy(sz, SizeG, sizeof(double)*K);
				non_improve = 0;
			}
			else
				non_improve++;
			if (ff > f_best)
			{
				f_best = ff;
				total_time = (clock() - starting_time) / CLOCKS_PER_SEC;
				for (int i = 0; i < N; i++)
				{
					GS.p[i] = p1[i];
				}

				//break;
			}
		}
		else
			non_improve++;

		if ((non_improve + 1) % 100 == 0)
		{
			//Pertubation(0.3*N, p1, SizeG);
			//RandomShake(0.1*N, p1, SizeG);
			RandomShake1(0.10*N, p1, SizeG, iter);
			Build_Delta_Matrix_Repair(p1, ff);

		}

		if ((iter_ils + 1) % UP == 0)
		{
			if (penalty_count >= UP - 1)
				penalty_factor *= 2.0;
			else if (penalty_count <= 1)
				penalty_factor /= 2.0;
			penalty_count = 0;
		}

		//cout << "iter=" << iter << " ff=" << ff << " f_best=" << f_best << " penalty_factor=" << penalty_factor << " gs=" << gs << endl;
		//	cout << "iter=" << iter << " ff=" << ff << " f_best=" << ComputeMinF() << " penalty_factor=" << penalty_factor << " gs=" << gs << endl;

	}
	if (non_improve != 0)
		non_improve++;

	if (!IsProperSol(p1, SizeG))
	{
		//if (flag_f1 == 1)
		memcpy(p1, p_temp, sizeof(int)*N);
		//else
		//	memcpy(p1, GS.p, sizeof(int)*N);
		for (int i = 0; i < K; i++)
			SizeG[i] = 0.0;
		for (int i = 0; i < N; i++)
			SizeG[p1[i]] += w[i];
		RandomShake1(0.10*N, p1, SizeG, iter);
		Build_Delta_Matrix_Repair(p1, ff);


	}

	delete p_temp; p_temp = NULL;
	delete sz; sz = NULL;
}


void Tabu_Search(int p1[], double SizeG[], double *cost)
{
	int i, j, v, g;
	int x, y, z;
	int iter, select, cluster;
	int old_group, old_group1, old_group2, swap;
	int flag_type;												//标记选择的邻域类型
	int stat_count[9];
	double best_delta_global, tabu_best_delta_global;
	double ff, sum;
	int count1 = 0, count2 = 0;
	for (i = 0; i<N; i++)
		p[i] = p1[i];
	for (i = 0; i < 9; i++)
		stat_count[i] = 0;
	Build_Delta_Matrix();
	*cost = f;
	ff = f;
	f_best = f;
	iter = 0;
	penalty_factor = 2;
	penalty_count = 0;
	non_improve = 0;
	ys_best = f;
	while (1.0*(clock() - starting_time) / CLOCKS_PER_SEC < Time_limit)
	{

		Feasible_Local_Search(p1, SizeG, iter, ff);
		if ((non_improve + 1) % 1000 == 0)							//trapped in a deep local optimum
		{
			Infeasible_Local_Search(p1, SizeG, iter, ff);
		}
		if ((non_improve + 1) % 500 == 0)
		{
			RandomShake1(0.10*N, p1, SizeG, iter);
			Build_Delta_Matrix_Repair(p1, ff);
		}
		cout << "iter=" << iter << " f_best=" << f_best << endl;

	}
	*cost = ff;

}

void RandomShake1(int L, int p1[], double SizeGroup[], int iter)
{
	int i, v, g, x, y;
	int NumberNeighbors, old_g, old_g1, swap;
	int cur_index, theta, count = 0;
	theta = L;

	NumberNeighbors = N*(N - 1) / 2 + N*K;
	do
	{
		cur_index = rand() % NumberNeighbors;
		if (Neighbors[cur_index].type == 1)
		{
			v = Neighbors[cur_index].v;
			g = Neighbors[cur_index].g;
			if ((p1[v] != g) && (SizeGroup[p1[v]] - w[v] >= LB[p1[v]]) && (SizeGroup[g] + w[v] <= UB[g]) && (TabuTenure[v][g] <= iter))
			{
				old_g = p1[v];
				SizeGroup[old_g] = SizeGroup[old_g] - w[v];
				SizeGroup[g] = SizeGroup[g] + w[v];
				p1[v] = g;
				count++;
			}
		}
		else if (Neighbors[cur_index].type == 2)
		{
			x = Neighbors[cur_index].x;
			y = Neighbors[cur_index].y;
			if ((p1[x] != p1[y]) && (SizeGroup[p1[x]] + (w[y] - w[x]) >= LB[p1[x]])
				&& (SizeGroup[p1[x]] + (w[y] - w[x]) <= UB[p1[x]])
				&& (SizeGroup[p1[y]] + (w[x] - w[y]) >= LB[p1[y]])
				&& (SizeGroup[p1[y]] + (w[x] - w[y]) <= UB[p1[y]])
				&& (TabuTenure[x][p1[y]] <= iter) && (TabuTenure[y][p1[x]] <= iter))
			{
				old_g = p1[x];
				old_g1 = p1[y];

				SizeGroup[p1[x]] += (w[y] - w[x]);
				SizeGroup[p1[y]] += (w[x] - w[y]);

				swap = p1[x];
				p1[x] = p1[y];
				p1[y] = swap;
				count++;
			}
		}
		//cout << "count=" << count << endl;
	} while (count < theta);
}
//计算handover minimation problem 目标函数值
double ComputeMinF()
{
	int *p;
	int i, j;
	double *q, sum = 0.0;
	int *flag;
	p = new int[K];
	q = new double[K];
	flag = new int[N];
	for (i = 0; i < K; i++)
	{
		p[i] = 0;
		q[i] = 0.0;

	}
	for (i = 0; i < N; i++)
		flag[i] = 0;
	for (i = 0; i < N; i++)
	{
		//printf("OS.p[%d]=%d\n", i, OS.p[i]);
		p[GS.p[i]]++;
	}

	for (i = 0; i < N; i++)
	{
		if (flag[i] != 1)
		{
			for (j = 0; j < N; j++)
			{
				if (GS.p[i] != GS.p[j])
				{
					//cout << "i=" << i << " j=" << j <<" OS.p[i]="<<OS.p[i]<< endl;
					q[GS.p[i]] += Distance[i][j];
					flag[i] = 1;
					//flag[j] = 1;
				}
			}
		}
	}
	for (i = 0; i < K; i++)
	{
		sum += q[i];
	}
	//cout << "sum=" << sum << endl;

	delete p; p = NULL;
	delete q; q = NULL;
	delete flag; flag = NULL;
	return sum;
}

void TS()
{
	double ctime = 0;
	int i, j;
	GS.cost = -9999999;
	f_best = -999999;
	starting_time = clock();
	RandomInitiaSol(CS.p, CS.SizeG);
	ctime = 1.0*(clock() - starting_time) / CLOCKS_PER_SEC;
	cout << "initial finished" << endl;
	Tabu_Search(CS.p, CS.SizeG, &CS.cost);
}
int main(int argc, char *argv[])
{
	int i, j;
	int i1, j1;
	int seed;
	double handSum;
	const int  Times = 100;
	double F[Times], H[Times];
	double Ctime[Times];
	double AvgTime = 0.0;
	double F_best = -99999999, F_worst = 999999999, F_ave = 0.0, deviation = 0;
	double H_best = -99999999, H_worst = 99999999, H_ave = 0.0;
	seed = time(NULL) % 1000000;
	srand(seed);

	File_Name = argv[1];
	Output_File_Name = argv[2];
	char * analyze_crx = argv[3];
	Initializing();
	//Initializing1();
	AllocateMemory();
	Time_limit = 0.2 * N;
	BuildNeighbors();

	OS.cost = -9999999.0;
	for (j = 0; j<Times; j++)
	{
		F[j] = 0.0;
		H[j] = 0.0;
		Ctime[j] = 0;
	}
	for (i = 0; i < Times; i++)
	{
		TS();
		if (Proof(GS))
		{
			F[i] = GS.cost;
			Ctime[i] = total_time;
			if (F[i]> OS.cost)
			{
				for (i1 = 0; i1<N; i1++)
					OS.p[i1] = GS.p[i1];
				for (j1 = 0; j1<K; j1++)
					OS.SizeG[j1] = GS.SizeG[j1];
				OS.cost = GS.cost;
			}
		}
		handSum = ComputeMinF();							//针对Handover算例
		if (Proof(GS))
			H[i] = handSum;
		printf("f[%d] = %lf h[%d] = %lf\n", i, F[i], i, H[i]);
		Out_results1(F[i], Output_File_Name, File_Name, Ctime[i], H[i]);
		Outputing(GS, analyze_crx);
	}
	for (i = 0; i<Times; i++)
	{
		if (F[i] > F_best)
			F_best = F[i];
		if (F[i] < F_worst)
			F_worst = F[i];
		F_ave += F[i];
		AvgTime += Ctime[i];
	}
	for (i = 0; i<Times; i++)
	{
		if (H[i] > H_best)
			H_best = H[i];
		if (H[i] < H_worst)
			H_worst = H[i];
		H_ave += H[i];
	}
	F_ave /= Times;
	H_ave /= Times;
	AvgTime /= Times;
	Out_results(F_best, F_ave, F_worst, H_best, H_ave, H_worst, AvgTime, Output_File_Name, File_Name);
	//Outputing(OS, File_Name);
	//ObserveNodeLocate();
	ReleaseMemory();
	getchar();
	return 0;
}
