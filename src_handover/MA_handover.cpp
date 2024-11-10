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

//the following are used parameters
int Search_depth_FLS;								//search depth of each FLS phase
int Iters_InfLS;									//maximum number of iterations of each InfLS phase
int Tabu_tenure;									//tabu tenure
int Frequency_shake;								//frequency of shake
double Shake_stren;									//shake strength
int Iters_FITS;										//number of iterations of FITS of MA
int Pop_num;										//population size


int PSIZE;							//another elment for the offspring		
#define UP 5

char * File_Name;
char * Output_stat_Name;
char * Output_sol_file;
char * Solution_File;

int N, K;											// node number and group number
double f, f_best, f_best1;							//f_best1 记录IVNS_TS method IVNS阶段跑完的最好结果
Solution CS, NS, GS, OS;							//CS:current solution ,GS: global solution
Neighborhood *Neighbors;

int **TabuTenure;
//以下变量用于tabu search N1,N2，N3,领域(N4未使用)
int num_tabu_best_1 = 0, num_best_1 = 0;
int num_tabu_best_2 = 0, num_best_2 = 0;
int num_tabu_best_3 = 0, num_best_3 = 0;
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
double total_time, starting_time, Time_limit;
int * p;								// assign  every vertex to a cluster, 
int * bestp;							// patition array for each vertex 

int non_improve;
//int * bestp ;				
double * SizeG;
double ** Delta_Matrix;					// incremental matrix 
double ** Distance;						// distance matrix between elements
double ** DistanceMulti;
double * LB;							// Lower bound of the number of elements in the group i 
double * UB;							// Upper bound of the number of elements in the group i 
double * w;								// node weights

double penalty_factor;
double ys, ys_best;
double gs;
int penalty_count;

//以下变量用于混合进化算法
Solution *Pop;
//int **Cluster1_Node;					//存储解（s1）每个cluster中顶点编号
//int **Cluster2_Node;					//存储解（s2）每个cluster中顶点编号
double *Parent_Cluster_F[2];			//存储解（s）每个cluster的object value

int *Offspring;							//父代P1，P2,杂交产生的子代O
double **Dist;							//种群中个体间的distance（n-similarity）
double *Avg_Dist;						//平均dist
double *QDF;							//quality-and-distance(QDF for short)
int *OR_FS;								//the rank of fs of the i-th individual
int *DR_Avg_Dist;						//the rank of avg_dist of the i-th individual


int *Parent_Cluster_Len[2];
int **Parent_Cluster_Vec[2];

int **Pop_Cluster_Len;					//用于计算种群相似性的变量，存储个体每个Cluster的长度
int ***Pop_Cluster_Vec;					//用于计算种群相似性的变量，存储个体每个Cluster包含的顶点编号

void RandomShake(int L, int p1[], double SizeGroup[], int iter);
void Tabu_Search(int p1[], double SizeG[], double *cost);
void Build_Delta_Matrix();
//reading handover instances
void Initializing_handover()
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
}


void AllocateMemory()
{
	int i, j;
	p = new int[N];
	Offspring = new int[N];
	//bestp = new int [N];
	SizeG = new double[K];

	Delta_Matrix = new double *[N];
	Dist = new double *[PSIZE];												//最后一个单元留给offspring
	for (int i = 0; i<N; i++)
		Delta_Matrix[i] = new double[K];
	for (int i = 0; i < PSIZE; i++)
		Dist[i] = new double[PSIZE];
	Avg_Dist = new double[PSIZE];
	QDF = new double[PSIZE];
	OR_FS = new int[PSIZE];
	DR_Avg_Dist = new int[PSIZE];

	CS.p = new int[N];
	NS.p = new int[N];
	GS.p = new int[N];
	OS.p = new int[N];

	CS.SizeG = new double[K];
	NS.SizeG = new double[K];
	GS.SizeG = new double[K];
	OS.SizeG = new double[K];

	Neighbors = new Neighborhood[N*(N - 1) / 2 + N*K];

	TabuTenure = new int*[N];
	for (i = 0; i < N; i++)
		TabuTenure[i] = new int[K];

	Pop = new Solution[PSIZE];
	for (int i = 0; i < PSIZE; i++)
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
	Parent_Cluster_F[0] = new double[K];
	Parent_Cluster_F[1] = new double[K];
	Parent_Cluster_Len[0] = new int[K];
	Parent_Cluster_Len[1] = new int[K];
	
	
	Parent_Cluster_Vec[0] = new int*[K];
	Parent_Cluster_Vec[1] = new int*[K];


	Pop_Cluster_Len = new int*[PSIZE];
	Pop_Cluster_Vec = new int**[PSIZE];

	for (int i = 0; i < PSIZE; i++)
	{
		Pop_Cluster_Len[i] = new int[K];
		Pop_Cluster_Vec[i] = new int*[K];
	}


	for (i = 0; i < PSIZE; i++)
	{
		for (j = 0; j < K; j++)
		{
			Pop_Cluster_Vec[i][j] = new int[N];
		}
	}
	for (j = 0; j < K; j++)
	{
		Parent_Cluster_Vec[0][j] = new int[N];
		Parent_Cluster_Vec[1][j] = new int[N];		
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
void Outputing(Solution &S, char *filename, char *output_file)
{
	int i;
	FILE *fp;
	char buff[80];
	if (Proof(S) == 0)
		return;
	sprintf(buff, output_file);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s N = %d  G = %d  f = %lf\n", filename, N, K, S.cost);
	for (i = 0; i<K; i++)
		fprintf(fp, "%lf   %lf   %lf\n", LB[i], UB[i], S.SizeG[i]);
	printf("\n");
	for (i = 0; i<N; i++)
		fprintf(fp, "%5.4d   %5.3d\n", i, S.p[i]);
	fclose(fp);
}


void Out_results(double best, double ave, double worst, double AvgTime, char *filename, char instance[])
{
	FILE *fp;
	char buff[80];
	sprintf(buff, "%s", filename);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s   %lf   %lf   %lf    %lf\n", instance, best, ave, worst, AvgTime);
	fclose(fp);
}

void Out_results1(double fi,char *output_filename, char instance[], double time)
{
	FILE *fp;
	char buff[80];
	sprintf(buff, "%s", output_filename);
	fp = fopen(buff, "a+");
	fprintf(fp, "%s %lf  %lf\n", instance,fi, time);
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

//计算两个个体之间的相似度
int Calculate_Sim_Between_Two(int *p1, int *p)
{
	int *flag1 = new int[N];
	int *flag2 = new int[N];
	int *one_Len = new int[K];
	int *two_Len = new int[K];
	int **one_Vec = new int*[K];
	int **two_Vec = new int*[K];
	int sum = 0, count, i, sum1, delta;
	int m1, m2, m3, m4, target1, target2;
	for (int i = 0; i < K; i++)
	{
		one_Vec[i] = new int[N];
		two_Vec[i] = new int[N];
	}

	memset(flag1, 0, sizeof(int)*K);
	memset(flag2, 0, sizeof(int)*K);
	memset(one_Len, 0, sizeof(int)*K);
	memset(two_Len, 0, sizeof(int)*K);
	for (i = 0; i < N; i++)
	{
		one_Vec[p1[i]][one_Len[p1[i]]++] = i;
		two_Vec[p[i]][two_Len[p[i]]++] = i;
	}

	for (count = 0; count < K; count++)
	{
		delta = 0;
		for (m1 = 0; m1 < K; m1++)
		{

			if (flag1[m1] != 1)
			{
				for (m2 = 0; m2 < K; m2++)
				{
					if (flag2[m2] != 1)																		//表示染色体2 当前色族未被选中
					{
						sum1 = 0;

						for (m3 = 0; m3 <one_Len[m1]; m3++)
						{
							for (m4 = 0; m4 <two_Len[m2]; m4++)
							{
								if (one_Vec[m1][m3] == two_Vec[m2][m4])
								{
									sum1++;
									//cout << "sum1=" << sum1 << endl;
									break;
								}
							}

						}
						//cout << "sum1=" << sum1 <<" m1="<<m1<<" m2="<<m2<< endl;
						if (sum1 > delta)												//最大匹配
						{
							delta = sum1;
							//cout << "delta=" << delta << endl;
							target1 = m1;
							target2 = m2;
						}
					}

				}
			}
		}
		//cout << "delta=" << delta << " target1=" << target1 << " target2=" << target2 <<" count="<<count<< endl;
		//getchar();
		flag1[target1] = 1;
		flag2[target2] = 1;
		sum += delta;
		//diversity = MaxVtx - sum;
		//sumDiversity += diversity;
	}
	delete flag1; flag1 = NULL;
	delete flag2; flag2 = NULL;
	delete one_Len; one_Len = NULL;
	delete two_Len; two_Len = NULL;
	for (i = 0; i < K; i++)
	{
		delete[]one_Vec[i]; one_Vec[i] = NULL;
		delete[]two_Vec[i]; two_Vec[i] = NULL;
	}

	return sum;

}

//初始化种群:RanReal instances
void InitialPopulation()
{
	int i;
	int flag;		
	double *sizeG = new double[K];
	double cost;
	for (i = 0; i < Pop_num;i++)
	{
		flag = 1;
		RandomInitiaSol(p, sizeG);
		Tabu_Search(p, sizeG, &cost);			
		if (flag != 0)								//可以插入种群
		{
			Build_Delta_Matrix();
			memcpy(Pop[i].p, p, sizeof(int)*N);
			memcpy(Pop[i].SizeG, sizeG, sizeof(double)*K);
			Pop[i].cost = f;			
		}		
		//cout << "i=" << i << endl;
	}
	delete sizeG; sizeG = NULL;
}
void VerifySolution()
{
	int i;
	for (i = 0; i < K; i++)
		CS.SizeG[i] = 0;
	for (i = 0; i < N; i++)
	{
		if (CS.p[i] != -1)
		{
			CS.SizeG[CS.p[i]] += w[i];
		}
		else
		{
			cout << "illegal CS.p";
			getchar();
		}
	}
	for (i = 0; i < K; i++)
	{
		if (CS.SizeG[i] > UB[i] || CS.SizeG[i] < LB[i])
		{						
			cout << "illegal CS.sizeG " << "sizeG []" << i << " " << CS.SizeG[i];
			getchar();
		}		
	}
	
}

//随机分配剩余的reserve顶点, 成功返回1，失败返回0
int Assign_ResiVec_Random(int p[])
{
	int i, j;
	int p1;
	int count, c1, Nc;
	int *Flag = new int[N];
	double *SizeGroup = new double[K];
	int *G = new int[K];
	int *VN = new int[N];
	int try_num = 0;
	memset(SizeGroup, 0, sizeof(double)*K);
	for (i = 0; i < N;i++)
	{
		if (p[i] != -1)
		{
			SizeGroup[p[i]] += w[i];
		}
	}
	for (i = 0; i<N; i++)
		Flag[i] = 1;
	for (i = 0; i < N; i++)
		if (p[i] == -1)							//顶点尚未分配
		{
		Flag[i] = 0;
		}
	count = 0;
	Nc = 0;
	for (i = 0; i<K; i++)
		if (SizeGroup[i] < LB[i])
			G[count++] = i;						//count 为size小于LB的cluster的个数
	for (j = 0; j<N; j++)
		if (Flag[j] == 0)
		{

		VN[Nc++] = j;
		}
	//cout << "count=" << count <<" Nc="<<Nc<< endl;
	while (count > 0)
	{
		while (1)
		{
			//cout << "Nc=" << Nc <<"  Count="<<count<<  endl;
			if (Nc == 0)
				return 0;
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
	//cout << "count=" << count << " Nc=" << Nc << endl;
	try_num = 0;
	while (Nc > 0)
	{
		while (1)
		{

			p1 = VN[rand() % Nc];
			c1 = G[rand() % count];
			try_num++;
			if (try_num == 100000)
				return 0;
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
		try_num = 0;
		for (i = 0; i<K; i++)
			if (SizeGroup[i] < UB[i])
				G[count++] = i;
		for (j = 0; j<N; j++)
			if (Flag[j] == 0)
				VN[Nc++] = j;
	}	
	return 1;

	delete[] SizeGroup; SizeGroup = NULL;
	delete[] Flag; Flag = NULL;
	delete[] G; G = NULL;
	delete[] VN; VN = NULL;
}

//辅助函数：用于构造Parent_Cluster_Len和Parent_Cluster_Vec 数组,	p1,p2 为两个父代
void Build_Parent_LenAndVec(int p1, int p2)
{
	int j;

	for (j = 0; j < K; j++)
	{
		Parent_Cluster_Len[0][j] = 0;
		Parent_Cluster_Len[1][j] = 0;
	}
	for (j = 0; j < N; j++)
	{		
		Parent_Cluster_Vec[0][Pop[p1].p[j]][Parent_Cluster_Len[0][Pop[p1].p[j]]++] = j;
		Parent_Cluster_Vec[1][Pop[p2].p[j]][Parent_Cluster_Len[1][Pop[p2].p[j]]++] = j;
	}
}

void Update_Parent_LenAndVec(int pa, int flag_p, int flag[])
{
	int j;
	if (flag_p == 0)
	{
		for (j = 0; j < K; j++)
			Parent_Cluster_Len[1][j] = 0;
		for (j = 0; j < N; j++)
		{
			if (flag[j] != 1)					//j仍然未分配
			{
				Parent_Cluster_Vec[1][Pop[pa].p[j]][Parent_Cluster_Len[1][Pop[pa].p[j]]++] = j;
			}
		}
	}
	else if (flag_p == 1)
	{
		for (j = 0; j < K; j++)
			Parent_Cluster_Len[0][j] = 0;
		for (j = 0; j < N; j++)
		{
			if (flag[j] != 1)
			{
				Parent_Cluster_Vec[0][Pop[pa].p[j]][Parent_Cluster_Len[0][Pop[pa].p[j]]++] = j;
			}
		}
	}
}

void Build_Parent_Cluster_F(int v1, int k1)
{
	int j, k, h;
	int m, m2;

	for (j = 0; j < K; j++)
	{
		Parent_Cluster_F[0][j] = 0;
		Parent_Cluster_F[1][j] = 0;
	}
	for (j = 0; j < K; j++)
	{
		for (k = 0; k < Parent_Cluster_Len[0][j]; k++)
		{
			m = Parent_Cluster_Vec[0][j][k];
			for (h = k + 1; h < Parent_Cluster_Len[0][j]; h++)
			{
				m2 = Parent_Cluster_Vec[0][j][h];
				if (Pop[v1].p[m] == Pop[v1].p[m2])				//肯定成立
					Parent_Cluster_F[0][j] += Distance[m][m2];

			}
		}
	}

	for (j = 0; j < K; j++)
	{
		for (k = 0; k < Parent_Cluster_Len[1][j]; k++)
		{
			m = Parent_Cluster_Vec[1][j][k];
			for (h = k + 1; h < Parent_Cluster_Len[1][j]; h++)
			{
				m2 = Parent_Cluster_Vec[1][j][h];
				if (Pop[k1].p[m] == Pop[k1].p[m2])				//肯定成立
					Parent_Cluster_F[1][j] += Distance[m][m2];

			}
		}
	}
}

void Update_Parent_Cluter_F(int pa, int flag_p)
{
	int j, k, h;
	int m, m2;


	if (flag_p == 0)
	{
		for (j = 0; j < K; j++)
			Parent_Cluster_F[1][j] = 0;
		for (j = 0; j < K; j++)
		{
			for (k = 0; k < Parent_Cluster_Len[1][j]; k++)
			{
				m = Parent_Cluster_Vec[1][j][k];
				for (h = k + 1; h < Parent_Cluster_Len[1][j]; h++)
				{
					m2 = Parent_Cluster_Vec[1][j][h];
					if (Pop[pa].p[m] == Pop[pa].p[m2])				//肯定成立
						Parent_Cluster_F[1][j] += Distance[m][m2];

				}
			}
		}
	}
	else if (flag_p == 1)
	{
		for (j = 0; j < K; j++)
			Parent_Cluster_F[0][j] = 0;
		for (j = 0; j < K; j++)
		{
			for (k = 0; k < Parent_Cluster_Len[0][j]; k++)
			{
				m = Parent_Cluster_Vec[0][j][k];
				for (h = k + 1; h < Parent_Cluster_Len[0][j]; h++)
				{
					m2 = Parent_Cluster_Vec[0][j][h];
					if (Pop[pa].p[m] == Pop[pa].p[m2])				//肯定成立
						Parent_Cluster_F[0][j] += Distance[m][m2];

				}
			}
		}
	}
}

//交叉操作：cluster 最大f值:相比Cross_Over，利用空间换时间，降低了时间复杂度，
void Cross_Over_Cluster()
{
	int i1, j1, l, sn;
	int select_clu;
	int select_vec;
	int *flag, flag_succ = 0;			
	double fff = -999, fff1;
	flag = new int[N];
	i1 = rand() % Pop_num;
	j1 = rand() % Pop_num;
	while (i1 == j1)
	{
		j1 = rand() % Pop_num;
	}
	Build_Parent_LenAndVec(i1, j1);
	Build_Parent_Cluster_F(i1, j1);		
	for (int i = 0; i < N; i++)
	{
		flag[i] = 0;
		CS.p[i] = -1;
	}
	l = 0;
	sn = 0;
	while (l < K && sn <= N )
	{
		fff = -999;
		fff1 = -999;
		if (l % 2 == 0)
		{
			for (int i = 0; i < K; i++)
			{
				if (Parent_Cluster_F[0][i] > fff)
				{
					fff = Parent_Cluster_F[0][i];
					select_clu = i;

				}
			}
			for (int i = 0; i < Parent_Cluster_Len[0][select_clu]; i++)
			{
				select_vec = Parent_Cluster_Vec[0][select_clu][i];
				CS.p[select_vec] = l;
				flag[select_vec] = 1;
				sn++;
			}

			Parent_Cluster_F[0][select_clu] = -999999;
			Parent_Cluster_Len[0][select_clu] = 0;

			Update_Parent_LenAndVec(j1, 0, flag);
			Update_Parent_Cluter_F(j1, 0);
		}
		else
		{
			for (int i = 0; i < K; i++)
			{
				if (Parent_Cluster_F[1][i] > fff1)
				{
					fff1 = Parent_Cluster_F[1][i];
					select_clu = i;
				}
			}
			for (int i = 0; i < Parent_Cluster_Len[1][select_clu]; i++)
			{
				select_vec = Parent_Cluster_Vec[1][select_clu][i];
				CS.p[select_vec] = l;
				flag[select_vec] = 1;
				sn++;
			}
			Parent_Cluster_F[1][select_clu] = -999999;
			Parent_Cluster_Len[1][select_clu] = 0;

			Update_Parent_LenAndVec(i1, 1, flag);
			Update_Parent_Cluter_F(i1, 1);
		}

		l++;
		//getchar();
	}
	//cout << "l=" << l << " sn=" << sn << endl;			
	if (sn < N )
		flag_succ = Assign_ResiVec_Random(CS.p);	
	if (flag_succ == 0)
	{
		cout << "in cross_over_cluster method, flag_succ=0" << endl;
		memcpy(CS.p, Pop[j1].p, sizeof(int)*N);
	}
	for (int i = 0; i < K; i++)
		CS.SizeG[i] = 0.0;
	for (int i = 0; i < N; i++)
		CS.SizeG[CS.p[i]] += w[i];
	VerifySolution();
	delete flag; flag = NULL;	
}



//更新种群：更新质量最差的individual（pool worst）
void UpdatePopulation()
{
	int maxm = 99999999;
	int select;
	for (int i = 0; i<Pop_num; i++)
	{
		if (Pop[i].cost<maxm)
		{

			maxm = Pop[i].cost;
			select = i;

		}
	}
	if (CS.cost > maxm)
	{
		for (int i = 0; i<N; i++)
			Pop[select].p[i] = CS.p[i];
		Pop[select].cost = CS.cost;
	}
}


//辅助函数：用于构造Pop_Cluster_Len和Pop_Cluster_Vec 数组
void BuildPopLenAndVec()
{
	int i, j;

	for (i = 0; i < Pop_num; i++)
	{
		for (j = 0; j < K; j++)
			Pop_Cluster_Len[i][j] = 0;
	}
	for (i = 0; i < Pop_num; i++)
	{
		for (j = 0; j < N; j++)
		{
			Pop_Cluster_Vec[i][Pop[i].p[j]][Pop_Cluster_Len[i][Pop[i].p[j]]++] = j;
		}
	}
}

//构造Dist矩阵
void Build_Pop_Dist()
{
	int i, j, sim;
	BuildPopLenAndVec();
	for (i = 0; i < Pop_num; i++)
	{
		for (j = i + 1; j < Pop_num; j++)
		{
			sim = Calculate_Sim_Between_Two(Pop[i].p, Pop[j].p);
			Dist[i][j] = N - sim;
			Dist[j][i] = Dist[i][j];
		}
	}

	for (i = 0; i < PSIZE; i++)
		Dist[i][i] = 0;
	for (i = 0; i < PSIZE; i++)
	{
		Dist[Pop_num][i] = 0;
		Dist[i][Pop_num] = 0;
	}
}

//辅助函数
void Build_Individual_LenAndVec(int num)
{
	int j;
	for (j = 0; j < K; j++)
		Pop_Cluster_Len[num][j] = 0;

	for (j = 0; j < N; j++)
	{
		//cout << " pop_cluster_len []"<<Pop_Cluster_Len[i][Pop[i].p[j]]<<endl;
		Pop_Cluster_Vec[num][Pop[num].p[j]][Pop_Cluster_Len[num][Pop[num].p[j]]++] = j;

	}

}

void BuildNeighbors()
{
	int i, j, g;
	int count;
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
	int j;

	for (j = 0; j<N; j++)
	{
		if (j != i)
		{
			Delta_Matrix[j][g0] -= Distance[i][j];
			Delta_Matrix[j][g1] += Distance[i][j];
		}
	}

}

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
		TabuTenure[best_x_1[select]][old_group] = Tabu_tenure + iter;
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
		TabuTenure[tabu_best_x_1[select]][old_group] = Tabu_tenure + iter;
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
		TabuTenure[best_x_2[select]][old_group] = Tabu_tenure + iter;
		TabuTenure[best_y_2[select]][old_group1] = Tabu_tenure + iter;
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
		TabuTenure[tabu_best_x_2[select]][old_group] = Tabu_tenure + iter;
		TabuTenure[tabu_best_y_2[select]][old_group1] = Tabu_tenure + iter;
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

		TabuTenure[best_x_3[select]][old_group] = Tabu_tenure + iter;
		TabuTenure[best_y_3[select]][old_group1] = Tabu_tenure + iter;
		TabuTenure[best_z_3[select]][old_group2] = Tabu_tenure + iter;
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

		TabuTenure[tabu_best_x_3[select]][old_group] = Tabu_tenure + iter;
		TabuTenure[tabu_best_y_3[select]][old_group1] = Tabu_tenure + iter;
		TabuTenure[tabu_best_z_3[select]][old_group2] = Tabu_tenure + iter;
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
	while (iter_ils < Iters_InfLS)
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

			cout << "here 1 in ils sum !=ff" << " ff=" << ff << " sum=" << sum << " best_fff1=" << best_fff1 << " tabu_bff1=" << tabu_best_fff1  << endl;
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

		if ((non_improve + 1) % Iters_InfLS == 0)
		{
			RandomShake(Shake_stren*N, p1, SizeG, iter);
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
		RandomShake(Shake_stren*N, p1, SizeG, iter);
		Build_Delta_Matrix_Repair(p1, ff);
	}
	delete p_temp; p_temp = NULL;
	delete sz; sz = NULL;
}

void Tabu_Search(int p1[], double SizeG[], double *cost)
{
	int i;
	int iter, select, cluster;
	int old_group, old_group1, old_group2, swap;
	int flag_type;												//标记选择的邻域类型
	double best_delta_global, tabu_best_delta_global;
	double ff, sum;
	int count1 = 0, count2 = 0;
	for (i = 0; i<N; i++)
		p[i] = p1[i];
	Build_Delta_Matrix();
	*cost = f;
	ff = f;
	//f_best = f;
	iter = 0;
	penalty_factor = 2;
	penalty_count = 0;
	non_improve = 0;
	ys_best = f;
	while (iter < Iters_FITS && 1.0*(clock() - starting_time) / CLOCKS_PER_SEC  <  Time_limit)
	{
		if ((non_improve + 1) % Frequency_shake == 0)
		{
			RandomShake(Shake_stren*N, p1, SizeG, iter);
			Build_Delta_Matrix_Repair(p1, ff);
		}
		Feasible_Local_Search(p1, SizeG, iter, ff);

		if ((non_improve + 1) % Search_depth_FLS == 0)							//trapped in a deep local optimum
		{
			Infeasible_Local_Search(p1, SizeG, iter, ff);
		}
		//cout << "iter=" << iter << " ff=" << ff << " f_best=" << f_best << endl;		
	}
	*cost = ff;
}

void RandomShake(int L, int p1[], double SizeGroup[], int iter)
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

//计算种群多样性
double Calculate_Pop_Diversity()
{
	double diversity = 0;
	for (int i = 0; i < Pop_num; i++)
	{
		for (int j = i + 1; j < Pop_num; j++)
		{
			diversity += Dist[i][j];
		}
	}
	return diversity / (Pop_num*(Pop_num - 1) / 2);
}

double ComputeMinF()
{
	double sum = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = i+1; j < N; j++)
		{
			sum+=Distance[i][j];
		}
	}
	sum -= f_best;
	sum *= 2;
	return sum;
}

//混合进化搜索
void Hybrid_TS()
{
	int i, j;
	int cnt = 0;
	double diversity;
	GS.cost = -9999999.0;
	f_best = -999999;
	starting_time = clock();
	InitialPopulation();
	
	cout << "here1" << endl; 


	for (i = 0; i < Pop_num; i++)
	{
		printf("pop[%d]=%lf\n", i, Pop[i].cost);
	
		if (Pop[i].cost>f_best)
		{
			f_best = Pop[i].cost;
			for (j = 0; j < N; j++)
				GS.p[j] = Pop[i].p[j];
		}
	}	
	while (1.0*(clock() - starting_time) / CLOCKS_PER_SEC  <  Time_limit)		
	{
		Cross_Over_Cluster();
		Tabu_Search(CS.p, CS.SizeG, &CS.cost);	
		UpdatePopulation();	
		Build_Pop_Dist();
		diversity = Calculate_Pop_Diversity();
		cout << "generation cnt=" << cnt++ << " fbest=" << f_best << ",diversity="<<diversity<< endl;
	}
}

int main(int argc, char *argv[])
{
	const int Times = 20;
	double F[Times], H[Times];
	double Ctime[Times];
	double AvgTime = 0.0;
	double F_best = -99999999, F_worst = 999999999, F_ave = 0.0, deviation = 0;
	double H_best = -99999999, H_worst = 99999999, H_ave = 0.0;
	
	if (argc < 4)
	{
		cout << "MA_CCP usage: input_file output_sol_file out_stat_file" << endl;
		cout << "(where input_file is the instance name, output_sol_file is a file used to store the solution information,	\
			output_stat_file stores the running information)" << endl;
		exit(-1);
	}

	File_Name = argv[1];
	Output_sol_file = argv[2];
	Output_stat_Name = argv[3];
	srand(unsigned(time(NULL)));
	
	//the following are used parameters
	Search_depth_FLS = 1000;								//search depth of each FLS phase
	Iters_InfLS = 200;									//maximum number of iterations of each InfLS phase
	Tabu_tenure = 10;									//tabu tenure
	Frequency_shake = 500;								//frequency of shake
	Shake_stren = 0.10;									//shake strength
	Iters_FITS = 10000;										//number of iterations of FITS of MA
	Pop_num = 5;										//population size

	
	PSIZE = Pop_num + 1;
	Initializing_handover();
	AllocateMemory();
	BuildNeighbors();
	

	Time_limit = 1.0 * N;
	OS.cost = -9999999.0;
	for (int j = 0; j < Times; j++)
	{
		F[j] = 0.0;
		H[j] = 0.0;
		Ctime[j] = 0;
	}
	for (int i = 0; i < Times; i++)
	{		
		Hybrid_TS();
		if (Proof(GS))
		{
			F[i] = GS.cost;
			Ctime[i] = total_time;
			if (F[i]> OS.cost)
			{
				for (int i1 = 0; i1 < N; i1++)
					OS.p[i1] = GS.p[i1];
				for (int j1 = 0; j1 < K; j1++)
					OS.SizeG[j1] = GS.SizeG[j1];
				OS.cost = GS.cost;
			}
		}	
		double handSum = ComputeMinF();
		if(Proof(GS))
			H[i] = handSum;
		printf("f[%d] = %lf \n", i, F[i]);
 	        printf(" h[%d] = %lf\n", i, H[i]);
		
		Out_results1(H[i], Output_stat_Name, File_Name, Ctime[i]);
	}
	for (int i = 0; i < Times; i++)
	{
		if (F[i] > F_best)
			F_best = F[i];
		if (F[i] < F_worst)
			F_worst = F[i];
		F_ave += F[i];
		AvgTime += Ctime[i];
	}
	for (int i = 0; i < Times; i++)
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
	Out_results(H_best, H_ave, H_worst, AvgTime, Output_stat_Name, File_Name);
	Outputing(OS, File_Name, Output_sol_file);
	//ObserveNodeLocate();
	ReleaseMemory();
	getchar();
	return 0;
}
