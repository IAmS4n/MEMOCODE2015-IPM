#include <stdio.h>
#include <cstdlib>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <set>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
//#include <xmmintrin.h>

//#pragma comment(linker, "/STACK:10000000000000")
//#pragma comment(linker,  "/HEAP:10000000000000")
using namespace std;
#define dominate_predict_mod3 1

#define rep(type,value,start,end)	for(type value=start;value<end;value++)

typedef unsigned long long ULL;
typedef unsigned int INX_TYPE;
typedef unsigned short VAL_TYPE;
typedef unsigned short TIME_TYPE;
typedef unsigned short SUM_TYPE;
//typedef unsigned int SUM_TYPE;//7*VAL_TYPE


void initValues();
void outputAndCleanup();
inline bool dominates(INX_TYPE a, INX_TYPE b);
int cmpfunc(const void *a, const void *b);

////////////////////////////////////////////////////////////////INPUT/////////////////////////////////////////////////////////////
const int dimensions=7;
const INX_TYPE datasetSize=800000;
const TIME_TYPE maxTimeSteps=40989;
const INX_TYPE maxSkylineElements=2008;
#define OUTPUT_PATH	"large.out"
#define INPUT_PATH	"large.input"
#define TIME_PATH	"large.times"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

INX_TYPE **skyline;
TIME_TYPE datasetTimes[datasetSize][2];
INX_TYPE numSkyline[datasetSize];
VAL_TYPE dataset[datasetSize][dimensions];
SUM_TYPE datasum[datasetSize];


struct mypair{
	INX_TYPE i;
	SUM_TYPE v;
	mypair(){};

	mypair(INX_TYPE init_index, SUM_TYPE init_value)
	{
		i = init_index;
		v = init_value;
	}
}
datasetEtimes[datasetSize];

struct mylongint{
	unsigned long long a,b;
}
dataset_quad[datasetSize];

#ifdef	dominate_predict_mod1
	mylongint MSB[datasetSize][2];
#elif	dominate_predict_mod2
	char MSB[datasetSize][2];
#elif	dominate_predict_mod3
	mylongint MSB[datasetSize][2];
#endif



/////////////////////////////////////////////////////////////////////////////////////////////////////////
struct compare {
	bool operator() (const mypair &firstElem, const mypair &secondElem) const
	{return ((firstElem.v!=secondElem.v)?(firstElem.v < secondElem.v):(firstElem.i < secondElem.i));}
};
bool pairCompare(const mypair &firstElem, const mypair &secondElem){
	return ((firstElem.v!=secondElem.v)?(firstElem.v < secondElem.v):(firstElem.i < secondElem.i));
}
typedef set<mypair,compare>::iterator SIT_TYPE;

inline void findsl(INX_TYPE pi,const set<mypair,compare> &sds,const set<mypair,compare> &ssl,vector<INX_TYPE>& vr){
	/*
	SIT_TYPE it=sds.lower_bound(mypair(0,datasum[pi]));
	if (it==sds.end())
		it=sds.begin();
	for (;it!=sds.end();++it)
		if(dominates(pi, it->i)==1)
			vr.push_back(it->i);
	*/

	/*
	SIT_TYPE								sit=sds.lower_bound(mypair(0,datasum[pi]));
	set<mypair,compare>::reverse_iterator	eit=sds.rbegin();
	int step=((eit->v)-(sit->v))/2;
	vector<INX_TYPE> res[2]; 
	#pragma omp parallel for
	for(int i=0;i<2;i++)
	{
		int from=sit->v+i*step;
		int to=sit->v+(i+1)*step;
		for(SIT_TYPE it=sds.lower_bound(mypair(0,from));it!=sds.end() && it->v<to;it++)
			if(dominates(pi, it->i)==1)
				res[i].push_back(it->i);
	}
	vr.reserve( res[0].size() +res[1].size());
	rep(int,i,0,2)
		vr.insert( vr.end(), res[i].begin(), res[i].end() );
	for(SIT_TYPE it=sds.lower_bound(mypair(0,sit->v+2*step));it!=sds.end();it++)
		vr.push_back(it->i);
	*/

	vr.clear();
	vector<INX_TYPE> good;
	SIT_TYPE h=sds.upper_bound(mypair(0,datasum[pi]+1));
	int emrg=0;
	for (;emrg<5 && h!=sds.end();++h){
		emrg++;
		if(dominates(pi, h->i)==1)
			good.push_back(h->i);
	}
	good.push_back(ssl.begin()->i);
	
	SIT_TYPE it=sds.lower_bound(mypair(0,datasum[pi]));
	if (it==sds.end())
		it=sds.begin();
	if(good.size()!=0){
		for (;it!=sds.end();++it)
			if(dominates(pi, it->i)==1){
				bool ok=true;
				rep(int,i,0,good.size())
					if(dominates(good[i], it->i)==1){
						ok=false;
						break;
					}
				if(ok)
					vr.push_back(it->i);
			}
	}
	else{
		for (;it!=sds.end();++it)
			if(dominates(pi, it->i)==1)	
				vr.push_back(it->i);
	}
	
}
inline void del_sl(set<mypair,compare>& ssl,INX_TYPE element_idx){
	if(1!=ssl.erase(mypair(element_idx,datasum[element_idx])))
		cout<<"\nERROR,"<<element_idx;
}
inline void del_sl(set<mypair,compare>& ssl,vector<INX_TYPE> velement_idx){
	rep(INX_TYPE,i,0,velement_idx.size())
		del_sl(ssl,velement_idx[i]);
}

inline void add_sl(set<mypair,compare> & ssl,INX_TYPE element_idx){
	SUM_TYPE sum=datasum[element_idx];
	ssl.insert(mypair(element_idx,sum));
}
inline INX_TYPE index_of_time(TIME_TYPE t){
	INX_TYPE a,b,m;
	a=0;
	b=datasetSize;
	while(b-a!=1){
		m=(a+b)/2;
		if(datasetTimes[m][0]<=t)
			a=m;
		else
			b=m;
	}
	return a;
}
inline void solve(TIME_TYPE part_start,TIME_TYPE part_end){
	
	set<mypair,compare>	sds,ssl;
	INX_TYPE i_start,i_end;
	
	i_end=i_start=0;
	if (part_start>1024){
		i_end=datasetEtimes[datasetSize-1].i;//index e bozorgtarin expire time

		//find start index for part_start time
		i_start=index_of_time( part_start-1024 );//index of first effectiv

		while(i_start<datasetSize && datasetTimes[i_start][0]<=part_start){
			if(datasetTimes[i_start][1]>part_start){
				if(datasetTimes[i_end][1]>datasetTimes[i_start][1])
					i_end=i_start;

				sds.insert(mypair(i_start,datasum[i_start]));
				bool isskyline=true;			
				for(SIT_TYPE iti=ssl.begin();iti!=ssl.end() && iti->v<=datasum[i_start];iti++)
					if(dominates(iti->i,i_start)==1){
						isskyline=false;
						break;
					}
				if (isskyline){
					vector<INX_TYPE> must_be_remove;
					for (SIT_TYPE it=ssl.lower_bound(mypair(0,datasum[i_start]))/*ssl.begin()*/;it!=ssl.end(); ++it)
						if(dominates(i_start,it->i)==1)
							must_be_remove.push_back(it->i);
					del_sl(ssl,must_be_remove);
					add_sl(ssl,i_start);
				}
			}
			i_start++;
		}
		i_end--;
	}

	rep(TIME_TYPE,time,part_start,part_end){
	//rep(TIME_TYPE,time,0,maxTimeSteps){
		while(i_start<datasetSize && datasetTimes[i_start][0]<=time){
			INX_TYPE element_idx=i_start;

			sds.insert(mypair(element_idx,datasum[element_idx]));

			bool isskyline=true;
			bool predict=false;
			
			if(!predict)
			{
				for(SIT_TYPE iti=ssl.begin();iti!=ssl.end() && iti->v<=datasum[element_idx];iti++)
				{
					if(dominates(iti->i,element_idx)==1){
						isskyline=false;
						break;
					}
				}
			}




			if (isskyline){
				vector<INX_TYPE> must_be_remove;
				for (SIT_TYPE it=ssl.lower_bound(mypair(0,datasum[element_idx]))/*ssl.begin()*/;it!=ssl.end(); ++it)
					if(dominates(element_idx,it->i)==1)
						must_be_remove.push_back(it->i);
				del_sl(ssl,must_be_remove);
				add_sl(ssl,element_idx);
			}
			i_start++;
		}
		
		vector<INX_TYPE> del_onskyline;
		while(i_end<datasetSize && datasetEtimes[i_end].v<=time){
			INX_TYPE element_idx=datasetEtimes[i_end].i;
			SIT_TYPE eraseit=sds.find(mypair(element_idx,datasum[element_idx]));
			if(eraseit!=sds.end()){
				if(ssl.find(mypair(element_idx,datasum[element_idx]))!=ssl.end()){
					del_sl(ssl,element_idx);
					del_onskyline.push_back(element_idx);
				}
				sds.erase(eraseit);
			}
			i_end++;
		}
		//cout<<endl<<i_start<<"\t"<<i_end;

		vector<INX_TYPE> tmp;
		rep(INX_TYPE,i,0,del_onskyline.size()){
			
			findsl(del_onskyline[i],sds,ssl,tmp);

			rep(INX_TYPE,j,0,tmp.size()){
				bool ok=true;
				for (SIT_TYPE k=ssl.begin();k!=ssl.end() && k->v<=datasum[tmp[j]];k++)
					if(dominates(k->i,tmp[j])==1){
						ok=false;
						break;
					}	
				if(ok){
					vector<INX_TYPE> mbd;
					for (SIT_TYPE k=ssl.lower_bound(mypair(0,datasum[tmp[j]]))/*ssl.begin()*/;k!=ssl.end();k++)
						if(dominates(tmp[j],k->i)==1)
							mbd.push_back(k->i);
					del_sl(ssl,mbd);
					add_sl(ssl,tmp[j]);
				}
			}
		}

		numSkyline[time]=ssl.size();
		INX_TYPE cur_index=0;

		for (SIT_TYPE k=ssl.begin();k!=ssl.end();k++){
			skyline[time][cur_index++]=k->i;
		}
	}
	cout<<"\nEND\t"<<part_start;
}

int main(int argc, char *argv[]) {
	initValues();
	printf("Start timing\n");
//	clock_t t;
///	t = clock();
        struct timeval t1,t2;
        printf("Start timing\n");
         gettimeofday(&t1,NULL);


	rep(INX_TYPE,i,0,datasetSize){
		unsigned int sum=0;
		datasum[i]=0;
		rep(int,d,0,dimensions)
			sum+=dataset[i][d];
		datasum[i]=((SUM_TYPE)(sum>>3));
		
		/*
		dataset_quad[i].a=0;
		rep(int,d,0,dimensions)
			dataset_quad[i].a+=dataset[i][d];
		dataset_quad[i].a>>=3;	
		rep(char,byte_inx,8,16){
			char same_bits=0;
			rep(int,d,0,dimensions)
				same_bits|=(dataset[i][d]&(1<<byte_inx))>>(byte_inx-d);
			dataset_quad[i].a|=((unsigned long long)same_bits)<<(7*(byte_inx-7));
		}

		dataset_quad[i].b=0;
		rep(int,d,0,dimensions)
			dataset_quad[i].b+=(dataset[i][d]&(1<<8-1));
		dataset_quad[i].b>>=3;	
		rep(char,byte_inx,0,8){
			char same_bits=0;
			rep(int,d,0,dimensions)
				same_bits|=(dataset[i][d]&(1<<byte_inx))>>(byte_inx-d);
			dataset_quad[i].b|=((unsigned long long)same_bits)<<(7*(byte_inx-7));
		}
		*/
		#ifdef dominate_predict_mod1
			//bug dare
			MSB[i][0].a=0;
			ULL mask	=	0x000000000000007fULL;//0111 1111
			ULL bitmask	=	0x0000000000010000ULL;//2^16
			for(int byte_inx=15;byte_inx>=8;byte_inx--){
				bitmask>>=1;
				ULL same_bits=0;
				rep(int,d,0,dimensions)
					same_bits|=(dataset[i][d]&bitmask)>>(byte_inx-d);
				same_bits&=mask;
				mask^=same_bits;
				MSB[i][0].a|=same_bits<<(7*byte_inx-48);
			}
			MSB[i][1].a=~MSB[i][0].a;
			
			MSB[i][0].b=0;
			mask	=	0x000000000000007fULL;//0111 1111
			bitmask	=	0x0000000000010000ULL;//2^16
			for(int byte_inx=7;byte_inx>=0;byte_inx--){
				bitmask>>=1;
				ULL same_bits=0;
				rep(int,d,0,dimensions)
					same_bits|=(dataset[i][d]&bitmask)>>(byte_inx-d);
				same_bits&=mask;
				mask^=same_bits;
				MSB[i][0].b|=same_bits<<(7*byte_inx-48);
			}
			MSB[i][1].b=~MSB[i][0].b;
		#elif dominate_predict_mod2
			MSB[i][0]=0;
			rep(int,d,0,dimensions)
				MSB[i][0]|=(dataset[i][d]&0x8000)>>(15-d);
			MSB[i][1]=~MSB[i][0];
		#elif dominate_predict_mod3
			MSB[i][0].b=MSB[i][0].a=0x0000000000000000ULL;
			rep(int,d,0,3){
				MSB[i][0].a<<=16;
				MSB[i][0].a|=dataset[i][d];
				MSB[i][0].a<<=1;

				MSB[i][0].b<<=16;
				MSB[i][0].b|=dataset[i][d+3];
				MSB[i][0].b<<=1;
			}
			MSB[i][0].a<<=12;
			MSB[i][0].a|=dataset[i][6]>>4;
			MSB[i][1].a=0x8000400020001000ULL|MSB[i][0].a; // 80 00 40 00 20 00 10 00

			MSB[i][0].b<<=12;
			MSB[i][0].b|=dataset[i][6]&0x000f;
			MSB[i][1].b=0x8000400020001000ULL|MSB[i][0].b; // 80 00 40 00 20 00 10 00
		#endif

		datasetEtimes[i].i=i;
		datasetEtimes[i].v=datasetTimes[i][1];
	}
	sort(datasetEtimes, datasetEtimes+datasetSize,pairCompare);

	////////////////////////////////////////////////////////////
	/*
	const TIME_TYPE mid=30000;
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			printf ("id = %d, \n", omp_get_thread_num());
			solve(mid,maxTimeSteps);
		}
		#pragma omp section
		{
			printf ("id = %d, \n", omp_get_thread_num());
			solve(0,mid);
		}
		
	}
	*/
	int core_number=1;
	if(argc>=2)
		core_number=atoi(argv[1]);
	else
		cout<<"\npart thread\n";
	if(argc>=3)
		omp_set_num_threads(atoi(argv[2]));

	int step=maxTimeSteps/core_number+1;
	#pragma omp parallel for schedule(dynamic)
		for(int i=core_number-1;i>=0;i--){
			printf ("\nid = %d, \n", omp_get_thread_num());
			solve(i*step,min((TIME_TYPE)((i+1)*step),(TIME_TYPE)maxTimeSteps));
		}
	
        gettimeofday(&t2,NULL);
        double elapsed_time = ((t2.tv_sec - t1.tv_sec)*1000.0) + ((t2.tv_usec - t1.tv_usec)/1000.0);
        printf("Elapsed Time  = %.5f ms\n",elapsed_time);

	

//	t = clock() - t;
//	printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	printf("End timing\n");
	outputAndCleanup();
}


// Allocate memory and read data from input files
void initValues() {
  int i, j, z;
  skyline = (INX_TYPE **)malloc(sizeof(*skyline) * maxTimeSteps);
 if (skyline) {
    for (i=0; i<maxTimeSteps; i++)
		skyline[i] = (INX_TYPE *)malloc(sizeof(INX_TYPE) * maxSkylineElements);
  }else
	  cout<<"\nMALLOC ERROR";
  
  FILE *f;
  if((f = fopen(INPUT_PATH, "r")) == NULL ) {
    printf("Can't open file input\n");
    assert(0);
  } 
  for (i=0; i < datasetSize; i++) {
    for (j=0; j < dimensions; j++) {
      int x;
      z = fscanf(f, "%d", &x);
      dataset[i][j] = (VAL_TYPE)x;
    }
  }
  fclose(f);

  if ((f = fopen(TIME_PATH, "r" )) == NULL ) {
    printf("Can't open file time\n");
    assert(0);
  } 
  
  for (i=0; i < datasetSize; i++) {
    int x;
    z = fscanf(f, "%d", &x);
    datasetTimes[i][0] = (TIME_TYPE)x;
    z = fscanf(f, "%d", &x);
    datasetTimes[i][1] = (TIME_TYPE)x;
  }
  fclose(f);

  for (i=0; i<maxTimeSteps; i++)
    numSkyline[i] = 0;

}
 

// Does element at dataset[a] dominate element at dataset[b]? Return 1 if yes, 
// 0 otherwise.

inline bool dominates(INX_TYPE a, INX_TYPE b) {
	#ifdef dominate_predict_mod1
		if(MSB[b][1].a & MSB[a][0].a)
			return false;
		//if(MSB[b][1].b & MSB[a][0].b)
			//return false;
	#elif dominate_predict_mod2
		if(MSB[b][1]&MSB[a][0])
			return false;
	#elif dominate_predict_mod3
		ULL def1=MSB[b][1].a-MSB[a][0].a;
		if( (def1 & 0x8000400020001000ULL)!=0x8000400020001000ULL)
			return false;
		ULL def2=MSB[b][1].b-MSB[a][0].b;
		if((def2 & 0x8000400020001000ULL)==0x8000400020001000ULL){
			return !(def1==def2 && def2==0x8000400020001000ULL);
				//return false;//mosabian daqiqan
		}
		if((def2 & 0x8000400020000000ULL)==0x8000400020000000ULL)
			return ((def1&0x0000000000000fffULL)!=0);
		return false;
	#endif	
	
	bool res=false;
	rep(char,d,0,dimensions){
		if (dataset[a][d] > dataset[b][d])
			return 0;
		else if (dataset[a][d] < dataset[b][d])
			res=true;
	}
	return res;
	
}   



void outputAndCleanup() {
  int i, l;
  FILE *f;

  if ((f = fopen(OUTPUT_PATH, "w")) == NULL) {
    printf("Can't open file output\n");
    assert(0);
  } 

  int maxSkyline = 0;
  
  for (i=0; i < (int)maxTimeSteps; i++) {
    qsort(skyline[i], numSkyline[i], sizeof(INX_TYPE), cmpfunc);
	for (l=0; l < numSkyline[i]; l++) {
		fprintf(f, "%d ", skyline[i][l]);
	}
	fprintf(f, "\n");
  }
  fclose(f);

  
}

// Just used to compare two indices
int cmpfunc(const void *a, const void *b) {
  return (*(INX_TYPE*)a - *(INX_TYPE*)b);
}
