
/*
attprior - detecting attractors in a synchronous Boolean network using prior information

Copyright: Tatsuya Akutsu, Kyoto University 

This program can be used only for acamemic purposes.
The author does not have any responsibility on problems caused by this program.
*/


#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<string.h>

/*
#define	ONEPERIODICONLY	1
*/

#define	OUTPUTALLSTATES 1

#define	MAXTRIALS	100000

#define	MAXLINESIZE	10000

#define MAXGENES	4000
#define MAXNAMELEN	10
#define	MAXDIFPROB	3

#define	MAXATTR		10000

#define	MAXPERIOD	1000
#define	MAXDEPTH	100000
#define	MAXBASINDEPTH	200

#define	MAXLIST		200000

#define	MAXCHANGEBITS	30
#define	MAXTOTALCHGBITS	MAXCHANGEBITS*MAXCHANGEBITS*MAXCHANGEBITS

char	BASINPATH[MAXDEPTH][MAXGENES];
char	ATTPOOL[MAXATTR][MAXGENES];
int	ATTPERIOD[MAXATTR];
int	AttNum,PattNum;
char	PERIODSTATES[MAXPERIOD][MAXGENES];
char	PERIODTEMP[MAXGENES];

char	GENENAMES[MAXGENES][MAXNAMELEN];
int	GeneNum;

#define	MAXDIFPERIOD	1000
int	NumDifPeriod;
int	DIFPERIODINF[MAXDIFPERIOD][2];

int	STATE0[MAXGENES];
int	STATE1[MAXGENES];
int	STATE2[MAXGENES];
int	INITSTATE[MAXGENES];
double	PROB[MAXGENES];
int	IDX[MAXGENES];
int	NUMPROB[MAXDIFPROB];
double	DIFPROB[MAXDIFPROB];

double	PROB2[MAXTOTALCHGBITS];
int	IDX2[MAXTOTALCHGBITS];
int	IDXK1[MAXTOTALCHGBITS];
int	IDXK2[MAXTOTALCHGBITS];
int	IDXK3[MAXTOTALCHGBITS];
int	CNTP;
int	N1,N2,N3;

double	ProbCertain;

int	NumTrial;


typedef	struct list_str {
	int type;
	int gene_id;
	struct list_str	*left;
	struct list_str	*right;
}	List;

#define	TYPE_GENE	0
#define	TYPE_NOT	1
#define	TYPE_AND	2
#define	TYPE_OR		3
#define	TYPE_IDENT	4
/* 2021Oct */
#define	TYPE_ZERO	5
#define	TYPE_ONE	6
/* 2021Oct */


#define	MAXTIMESTEP	2000
short	STATES[MAXGENES][MAXTIMESTEP];

List	LISTPOOL[MAXLIST];
int	ListNum;

typedef struct node_str {
	int	gene_id;
	List	*func;
}	Node;

Node	NODES[MAXGENES];
int	NodeCnt;

double	AveTrial;

int load_net(char fname[]);
int getnames(char str[]);
int add_gene(char gname[]);
int get_bool(char str[]);
List *get_bool_sub(char str[],int *idx,int len);
int get_gene_id(char str[],int *idx);
int print_func(int idx);
int print_func_sub(List *lptr,int level);
int load_init_val(char fname[]);
int print_states();
int calc_next_val(List *lptr);
int simulate(int ti,int *periodstart);
int output_time_series(int tn);
int examine_state();
int modify_init_val();
int print_states0();
int print_states2();
int search_att1(int k1,int k2,int k3);
int search_att2(int k1,int k2,int k3);
int search_att3(int k1,int k2,int k3);
int search_att_sub1(int curi,int curk,int k1,int k2,int k3);
int search_att_sub2(int curi,int curk,int k1,int k2,int k3);
int search_att_sub3(int curi,int curk,int k1,int k2,int k3);
int search_att_wp();
int sort_prob();
int add_new_att(int period);
int find_minimalstate(int period);
int fprint_states(FILE *fp);
int output_periodic_attractor(FILE *fp, int period);

int msort_real(double dat[],int n,int idx[]);
int msort_real_body(double dat[],int n,int idx[]);
int merge_real(double dat[],int m,int n,int idx[]);

/* 2021Oct */
int is_blank(char str[]);
/* 2021Oct */

int main(int argc,char *argv[])
{
	int	i,ti,gid;
	double	ex,p;


	if (argc != 3) {
		printf("Usage: AttPrior BoolNetFileName, PriorInfoFileName\n");
		exit(0);
	}
	load_net(argv[1]);
	load_init_val(argv[2]);

	for(gid = 0;gid < GeneNum;gid++) {
		INITSTATE[gid] = STATE0[gid];
	}

	printf("ListNum: %d :  %d  :  %d\n",ListNum,NodeCnt,GeneNum);

	AttNum = 0;
	PattNum = 0;

	sort_prob();
	search_att_wp();

}

int modify_init_val()
{
	int	i,hdist;
	double	r;

	hdist = 0;
	for(i = 0;i < GeneNum;i++) {
		r = drand48();
		if (r > ProbCertain) {
			hdist++;
			STATE0[i] = 1-INITSTATE[i];
		}
		else {
			STATE0[i] = INITSTATE[i];
		}
	}
}

int search_att_wp()
{
	int	i,j,h,pi,k1,k2,k3,retval;
	FILE	*fp;

	NumTrial = 0;

printf("N1: %d  N2: %d  N3: %d\n",N1,N2,N3);

	for(i = 0;i < CNTP;i++) {
		j = IDX2[i];
		k1 = IDXK1[j];
		k2 = IDXK2[j];
		k3 = IDXK3[j];
printf("k1: %d  k2: %d  k3: %d\n",k1,k2,k3);
		retval = search_att1(k1,k2,k3);
		if (retval == -1) {
			break;
		}
#ifdef	ONEPERIODICONLY
		if (retval == 1) {
			break;
		}
#endif
	}
	printf("\n------------\n");

 	if ((fp = fopen("found_attractors.txt","w"))==NULL) {
		printf("cannot open found_attractors.txt\n");
		exit(1);
	}
	fprintf(fp,"Following %d/%d singleton/periodic attractors were found\n",AttNum-PattNum,PattNum);
	for(j = 0;j < AttNum;j++) {
		for(i = 0;i < NodeCnt;i++) {
			STATE1[i] = ATTPOOL[j][i];
		}
		if (ATTPERIOD[j] > 1) {
			fprintf(fp,"%d-th (periodic) attractor with period %d ---\n",j+1,ATTPERIOD[j]);
		}
		else {
			fprintf(fp,"%d-th (singleton) attractor with period %d ---\n",j+1,ATTPERIOD[j]);
		}
#ifdef	OUTPUTALLSTATES
		if (ATTPERIOD[j] == 1) {
			fprint_states(fp);
		}
		else {
			output_periodic_attractor(fp,ATTPERIOD[j]);
		}
#else
		fprint_states(fp);
#endif
	}
	NumDifPeriod = 0;
	for(j = 0;j < AttNum;j++) {
		pi = ATTPERIOD[j];
		for(h = 0;h < NumDifPeriod;h++) {
			if (DIFPERIODINF[h][0] == pi) {
				DIFPERIODINF[h][1]++;
				break;
			}
		}
		if (h == NumDifPeriod) {
			DIFPERIODINF[h][0] = pi;
			DIFPERIODINF[h][1] = 1;
			NumDifPeriod++;
			if (NumDifPeriod == MAXDIFPERIOD) {
				printf("Error: Too many different periods\n");
				break;
			}
		}
	}
	fprintf(fp,"------ Distribution of Attractors -----\n");
		fprintf(fp,"%6s : %10s\n","PERIOD","#ATTRACTORS");
	for(h = 0;h < NumDifPeriod;h++) {
		fprintf(fp,"%6d : %10d\n",DIFPERIODINF[h][0],DIFPERIODINF[h][1]);
	} 
	fclose(fp);
	

}

int search_att1(int k1,int k2,int k3)
{
	return search_att_sub1(0,0,k1,k2,k3);
}

int search_att2(int k1,int k2,int k3)
{
	return search_att_sub2(N1,0,k1,k2,k3);
}

int search_att3(int k1,int k2,int k3)
{
	return search_att_sub3(N2,0,k1,k2,k3);
}

int search_att_sub1(int curi,int curk,int k1,int k2,int k3)
{
	int	i,gid,retval;

	if (curk == k1) {
		return search_att2(k1,k2,k3);
	}
	for(i = curi;i < N1;i++) {
		gid = IDX[i];
		STATE0[gid] = 1-STATE0[gid];
		retval = search_att_sub1(i+1,curk+1,k1,k2,k3);
		STATE0[gid] = 1-STATE0[gid];
		if (retval == -1) {
			return -1;
		}
#ifdef	ONEPERIODICONLY	
		if (retval == 1) {
			return 1;
		}
#endif
	}
	return 0;
}

int search_att_sub2(int curi,int curk,int k1,int k2,int k3)
{
	int	i,gid,retval;

	if (curk == k2) {
		return search_att3(k1,k2,k3);
	}
	for(i = curi;i < N2;i++) {
		gid = IDX[i];
		STATE0[gid] = 1-STATE0[gid];
		retval = search_att_sub2(i+1,curk+1,k1,k2,k3);
		STATE0[gid] = 1-STATE0[gid];
		if (retval == -1) {
			return -1;
		}
#ifdef	ONEPERIODICONLY	
		if (retval == 1) {
			return 1;
		}
#endif
	}
	return 0;
}
int search_att_sub3(int curi,int curk,int k1,int k2,int k3)
{
	int	i,gid,retval;

	if (curk == k3) {
		return examine_state();
	}
	for(i = curi;i < N3;i++) {
		gid = IDX[i];
		STATE0[gid] = 1-STATE0[gid];
		retval = search_att_sub3(i+1,curk+1,k1,k2,k3);
		STATE0[gid] = 1-STATE0[gid];
		if (retval == -1) {
			return -1;
		}
#ifdef	ONEPERIODICONLY	
		if (retval == 1) {
			return 1;
		}
#endif
	}
	return 0;
}
	
int examine_state()
{
	int	gid,ti,retval,pt;

	NumTrial++;

	if (NumTrial > MAXTRIALS) {
		printf("\nToo many trials\n");
		return -1;
	}

	for(gid = 0;gid < GeneNum;gid++) {
		STATE1[gid] = STATE0[gid];
	}
	for(gid = 0;gid < GeneNum;gid++) {
		STATES[gid][0] = STATE1[gid];
		BASINPATH[0][gid] = STATE1[gid];
	}
	for(ti = 1;ti < MAXBASINDEPTH;ti++) {

		retval = simulate(ti,&pt);

		if (retval == 1) {
			printf("Periodic Attractor was found at %d-th trial (period %d)\n",NumTrial,ti-pt);
			for(gid = 0;gid < GeneNum;gid++) {
				STATE1[gid] = BASINPATH[pt][gid];
			}
			print_states();
			find_minimalstate(ti-pt);
			AveTrial += NumTrial;

/*
			output_time_series(31);
*/

			return 1;
		}
		else if (retval == 2) {
/*
			printf("Singleton Attractor was found at %d-the step of %d-th trial\n",ti,NumTrial);
*/
/*
			output_time_series(1);
			print_states();
*/
			putchar('+');
			fflush(stdout);
			add_new_att(1);

			return 0;
		}
		for(gid = 0;gid < GeneNum;gid++) {
			STATES[gid][ti] = STATE1[gid];
		}
	}
	return 0;
}


int output_time_series(int tn)
{
	FILE	*fp;
	int	i,gid,ti;

	if ((fp=fopen("time_series.txt","w"))==NULL) {
		printf("Cannnot open time_series.txt\n");
		exit(0);
	}
	for(i = 0;i < NodeCnt;i++) {
		gid = NODES[i].gene_id;
		fprintf(fp,"%-10s ",GENENAMES[gid]);
		for(ti = 0;ti < tn;ti++) {
			fprintf(fp,"%d  ",STATES[gid][ti]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

int find_minimalstate(int period)
{
	int	i,j,gid,pcnt,mint,ti,flag;
	List	*lptr;


	pcnt = 0;
	for(i = 0;i < NodeCnt;i++) {
		PERIODSTATES[pcnt][i] = (char) STATE1[i];
	}
	for(pcnt = 1;pcnt < MAXPERIOD;pcnt++) {
		for(i = 0;i < NodeCnt;i++) {
			gid = NODES[i].gene_id;
			lptr = NODES[i].func;
			j = calc_next_val(lptr);
			STATE2[gid] = j;
		}
		for(i = 0;i < NodeCnt;i++) {
			PERIODSTATES[pcnt][i] = (char) STATE2[i];
			STATE1[i] = STATE2[i];
		}
		for(i = 0;i < NodeCnt;i++) {
			if (STATE2[i] != PERIODSTATES[0][i]) {
				break;
			}
		}
		if (i == NodeCnt) {
			break;
		}
	}
	if (pcnt == MAXPERIOD) {
		printf("Too long period\n");
		return 0;
	}
if (pcnt == 1 || pcnt != period) {
printf("Error in period: %d %d\n",pcnt,period);
exit(1);
}

	mint = 0;
	for(i = 0;i < NodeCnt;i++) {
		PERIODTEMP[i] = PERIODSTATES[0][i];
	}

	for(ti = 1;ti < pcnt;ti++) {
		flag = 0;

		for(i = 0;i < NodeCnt;i++) {
			if (PERIODSTATES[ti][i] < PERIODTEMP[i]) {
				flag = 1;
				break;
			}
			if (PERIODSTATES[ti][i] > PERIODTEMP[i]) {
				break;
			}
		}
		if (flag == 1) {
			for(i = 0;i < NodeCnt;i++) {
				PERIODTEMP[i] = PERIODSTATES[ti][i];
			}
			mint = ti;
		}
	}
	for(i = 0;i < NodeCnt;i++) {
		STATE1[i] = PERIODTEMP[i];
	}

	add_new_att(pcnt);
}

int add_new_att(int period)
{
	int	i,j;

	for(j = 0;j < AttNum;j++) {
		for(i = 0;i < NodeCnt;i++) {
			if (ATTPOOL[j][i] != STATE1[i]) {
				break;
			}
		}
		if (i == NodeCnt) {
/* not new attractor */
			return 0;
		}
	}


	printf("\nNew Attractor with Period %d was Found: #Att is %d\n",period,AttNum);
	print_states();

	for(i = 0;i < NodeCnt;i++) {
		ATTPOOL[AttNum][i] = STATE1[i];
	}
	ATTPERIOD[AttNum] = period;
if (AttNum >= MAXATTR) { printf("Use larger MAXATTR\n"); exit(0); }
	AttNum++;
	if (period > 1) {
if (PattNum >= MAXATTR) { printf("Use larger MAXATTR\n"); exit(0); }
		PattNum++;
	}
	return 1;
}


int print_states0()
{
	int	i,gid,iv;

	for(i = 0;i < NodeCnt && i < 500;i++) {
		gid = NODES[i].gene_id;
		iv = STATE0[gid];
		printf("%d",iv);
		if ((i % 70)==69) {
			printf("\n");
		}
	}
	printf("\n");
	printf("-----------\n");
}

int print_states()
{
	int	i,gid,iv;

	for(i = 0;i < NodeCnt && i < 500;i++) {
		gid = NODES[i].gene_id;
		iv = STATE1[gid];
		printf("%d",iv);
		if ((i % 70)==69) {
			printf("\n");
		}
	}
	printf("\n");
	printf("-----------\n");
}

int fprint_states(FILE *fp)
{
	int	i,gid,iv;

#ifdef	OUTPUTALLSTATES
	for(i = 0;i < NodeCnt;i++) {
#else
	for(i = 0;i < NodeCnt && i < 500;i++) {
#endif
		gid = NODES[i].gene_id;
		iv = STATE1[gid];
		fprintf(fp,"%d",iv);
		if ((i % 70)==69) {
			fprintf(fp,"\n");
		}
	}
	fprintf(fp,"\n");
	fprintf(fp,"-----------\n");
}

int print_states2()
{
	int	i,gid,iv;

	for(i = 0;i < NodeCnt && i < 500;i++) {
		gid = NODES[i].gene_id;
		iv = STATE2[gid];
		printf("%d",iv);
		if ((i % 70)==69) {
			printf("\n");
		}
	}
	printf("\n");
	printf("-----------\n");
}

int load_net(char fname[])
{
	FILE	*fp;
	int	i,len,maxlinesize;
static	char	buf1[MAXLINESIZE+1],buf2[MAXLINESIZE+1],buf3[MAXLINESIZE+1];

	GeneNum = 0;

	if ((fp = fopen(fname,"r"))==NULL) {
		printf("cannot open %s\n",fname);
		exit(0);
	}
	
	maxlinesize = 0;
	while(fgets(buf1,MAXLINESIZE,fp)!=NULL) {
		len = strlen(buf1);
		if (len > maxlinesize) {
			maxlinesize = len;
		}
/* 2021Oct */
		if (strncmp(buf1,"targets,",8)==0) {
			if (sscanf(&buf1[8],"%s",buf2)!=1 || strcmp(buf2,"factors")!=0) {
				printf("BoolNet file format error: %s : %s \n",buf1,buf2);
				exit(0);
			}
			continue;
		}
		if (buf1[0] == '#' || is_blank(buf1)==1) {
			continue;
		}
/* 2021Oct */
		getnames(buf1);
	}
	fclose(fp);

/*
	printf("maximum line length = %d\n",maxlinesize);
*/

	printf("#genes = %d\n",GeneNum);
	for(i = 0;i < GeneNum;i++) {
		printf("%-8s ",GENENAMES[i]);
		if ((i % 8)==7) {
			printf("\n");
		}
	}
	printf("\n");

	NodeCnt = 0;
	ListNum = 0;

	if ((fp = fopen(fname,"r"))==NULL) {
		printf("cannot open %s\n",fname);
		exit(0);
	}

	while(fgets(buf1,MAXLINESIZE,fp)!=NULL) {
/* 2021Oct */
		if (strncmp(buf1,"targets,",8)==0) {
			if (sscanf(&buf1[8],"%s",buf2)!=1 || strcmp(buf2,"factors")!=0) {
				printf("BoolNet file format error: %s : %s \n",buf1,buf2);
				exit(0);
			}
			continue;
		}
		if (buf1[0] == '#' || is_blank(buf1)==1) {
			continue;
		}
/* 2021Oct */
		get_bool(buf1);
	}

	fclose(fp);

	printf("#NODES = %d\n",NodeCnt);
	printf("#GeneNames = %d\n",GeneNum);
	printf("#LISTS = %d\n",ListNum);

/*
	for(i = 0;i < GeneNum;i++) {
		print_func(i);
	}
*/

}

int simulate(int ti,int *periodstart)
{
	int	i,j,k,gid,t;
	List	*lptr;

	for(i = 0;i < NodeCnt;i++) {
		gid = NODES[i].gene_id;
		lptr = NODES[i].func;
		j = calc_next_val(lptr);
		STATE2[gid] = j;
	}
	for(i = 0;i < NodeCnt;i++) {
		if (STATE2[i] != STATE1[i]) {
			break;
		}
	}
	if (i == NodeCnt) {
/* Singleton Attractor */
		return 2;
	}
	for(t = 0;t < ti;t++) {
		for(i = 0;i < NodeCnt;i++) {
			if (STATE2[i] != BASINPATH[t][i]) {
				break;
			}
		}
		if (i == NodeCnt) {
/* Periodic Attractor */
			*periodstart = t;
			return 1;
		}
	}
	for(i = 0;i < NodeCnt;i++) {
		STATE1[i] = STATE2[i];
		BASINPATH[ti][i] = STATE2[i];
	}
	return 0;
}

int calc_next_val(List *lptr)
{
	int	iv1,iv2,gid;

/* 2021Oct */
	if (lptr->type==TYPE_ZERO) {
		gid = lptr->gene_id;
		return 0;
	}
	else if (lptr->type==TYPE_ONE) {
		gid = lptr->gene_id;
		return 1;
	}
	else if (lptr->type==TYPE_GENE) {
/* 2021Oct */
		gid = lptr->gene_id;
		return STATE1[gid];
	}
	else if (lptr->type==TYPE_IDENT) {
		iv1 = calc_next_val(lptr->left);
		return iv1;
	}
	else if (lptr->type==TYPE_NOT) {
		iv1 = calc_next_val(lptr->left);
		return 1-iv1;
	}
	else if (lptr->type==TYPE_AND) {
		iv1 = calc_next_val(lptr->left);
		iv2 = calc_next_val(lptr->right);
		if (iv1 == 1 && iv2 == 1) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else if (lptr->type==TYPE_OR) {
		iv1 = calc_next_val(lptr->left);
		iv2 = calc_next_val(lptr->right);
		if (iv1 == 1 || iv2 == 1) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		printf("Error in CalcNextVal\n");
		print_func_sub(lptr,0);
		exit(0);
	}
}
		


int load_init_val(char fname[])
{
	FILE	*fp;
	int	i,j,id,idx;
	char	str1[5000],str2[5000];
	double	p;

	for(i = 0;i < GeneNum;i++) {
		STATE1[i] = 2;
	}
	if ((fp = fopen(fname,"r"))==NULL) {
		printf("Cannot open %s\n",fname);
		exit(0);
	}
	while(fgets(str1,1000,fp)!=NULL) {
		if (sscanf(str1,"%s %d %lf",str2,&j,&p)!=3) {
			printf("Error in %s\n",str1);
			exit(0);
		}
		idx = 0;
		id = get_gene_id(str2,&idx);
		STATE1[id] = j;
		PROB[id] = p;
	}
	fclose(fp);

	for(i = 0;i < GeneNum;i++) {
		if (STATE1[i] != 0 && STATE1[i] != 1) {
			printf("Error in state for %s : %d\n",GENENAMES[i],STATE1[i]);
			exit(0);
		}
	}

	for(i = 0;i < GeneNum && i < 20;i++) {
		id = NODES[i].gene_id;
		printf("%-10s : %d\n",GENENAMES[id],STATE1[id]);
	}

	for(i = 0;i < GeneNum;i++) {
		STATE0[i] = STATE1[i];
	}
}

		

int print_func(int idx)
{
	int	id;
	List	*lptr;

	id = NODES[idx].gene_id;
	printf("%s   ----\n",GENENAMES[id]);
	lptr = NODES[idx].func;
	print_func_sub(lptr,0);
}


int print_func_sub(List *lptr,int level)
{
	int	i,id;

/* 2021Oct */
	if (lptr->type == TYPE_ZERO) {
		for(i = 0;i < level;i++) {
			printf("   ");
		}
		printf("0\n");
	}
	else if (lptr->type == TYPE_ONE) {
		for(i = 0;i < level;i++) {
			printf("   ");
		}
		printf("1\n");
	}
	else if (lptr->type == TYPE_GENE) {
/* 2021Oct */
		id = lptr->gene_id;
		for(i = 0;i < level;i++) {
			printf("   ");
		}
		printf("%s\n",GENENAMES[id]);
	}
	else if (lptr->type == TYPE_AND) {
		for(i = 0;i < level;i++) {
			printf("   ");
		}
		printf("AND\n");
		print_func_sub(lptr->left,level+1);
		print_func_sub(lptr->right,level+1);
	}
	else if (lptr->type == TYPE_OR) {
		for(i = 0;i < level;i++) {
			printf("   ");
		}
		printf("OR\n");
		print_func_sub(lptr->left,level+1);
		print_func_sub(lptr->right,level+1);
	}
	else if (lptr->type == TYPE_NOT) {
		for(i = 0;i < level;i++) {
			printf("   ");
		}
		printf("NOT\n");
		print_func_sub(lptr->left,level+1);
	}
	else if (lptr->type == TYPE_IDENT) {
		for(i = 0;i < level;i++) {
			printf("   ");
		}
		printf("IDENT\n");
		print_func_sub(lptr->left,level+1);
	}
	else {
		printf("Error in print_func\n");
	}
}
		
		

int getnames(char str[])
{
	int	i,j,k,len;
	char	ch,gname[100];
	

	len = strlen(str);
	k = 0;
	for(i = 0;i < len;i++) {
		ch = str[i];
		if ((ch >= 'a' && ch <= 'z')||(ch >= 'A' && ch <= 'Z')||(ch >= '0' && ch <= '9')) {
			gname[k] = ch;
			k++;
		}
		else {
/* 2021Oct */
			if (k > 0 && !(k==1 && (gname[0]=='0' || gname[0]=='1'))) {
/* 2021Oct */
				gname[k] = '\0';
				j = add_gene(gname);
				k = 0;
			}
		}
	}
}


int add_gene(char gname[])
{
	int	i,glen;

	glen = strlen(gname);
	if (glen >= MAXNAMELEN-1) {
		printf("Too long gene name: %s\n",gname);
		exit(0);
	}
	for(i = 0;i < GeneNum;i++) {
		if (strcmp(GENENAMES[i],gname)==0) {
			return 0;
		}
	}
	strcpy(GENENAMES[GeneNum],gname);
	GeneNum++;

	if (GeneNum >= MAXGENES-1) {
		printf("Too many genes\n");
		exit(0);
	}

	return GeneNum;
}
	
int get_bool(char str[])
{
	int	len,idx,id,id2;
	char	ch1;
	List	*lptr,*lptr2;

	len = strlen(str);

	idx = 0;
	id = get_gene_id(str,&idx);
	if (str[idx] != ',') {
		printf("Error01 in %s : %s\n",str,str[idx]);
		exit(0);
	}
	idx++;
	idx++;
/* 2021Oct */
	if ((str[idx] == '0' || str[idx] == '1') && (str[idx+1]==' ' || str[idx+1] == '\t' || str[idx+1] == '\0' || str[idx+1] == '\n')) {
		lptr = &LISTPOOL[ListNum];
if (ListNum >= MAXLIST-1) { printf("Use larger MAXLIST\n"); exit(0); }
		ListNum++;
		if (str[idx] == '0') {
			lptr->type = TYPE_ZERO;
		}
		else {
			lptr->type = TYPE_ONE;
		}
		lptr->gene_id = id;
	}
	else if (str[idx] != '(' && str[idx] != '!') {
/* 2021Oct */
/* case of: GENEID, GENEID */
		lptr = &LISTPOOL[ListNum];
if (ListNum >= MAXLIST-1) { printf("Use larger MAXLIST\n"); exit(0); }
		ListNum++;
		lptr->type = TYPE_IDENT;
		lptr->gene_id = id;
		lptr2 = &LISTPOOL[ListNum];
if (ListNum >= MAXLIST-1) { printf("Use larger MAXLIST\n"); exit(0); }
		ListNum++;
		id2 = get_gene_id(str,&idx);
		lptr2->type = TYPE_GENE;
		lptr2->gene_id = id2;
		lptr->left = lptr2;
	}
	else if (str[idx] == '!') {
/* case of: GENEID, !(GENEID) */
		idx++;
		if (str[idx] != '(') {
			printf("Error1 in %s\n",str);
			exit(0);
		}
		lptr = &LISTPOOL[ListNum];
if (ListNum >= MAXLIST-1) { printf("Use larger MAXLIST\n"); exit(0); }
		ListNum++;
		idx++;
		lptr2 = get_bool_sub(str,&idx,len);
		lptr->left = lptr2;
		lptr->type = TYPE_NOT;
	}
	else {
		idx++;
		lptr = get_bool_sub(str,&idx,len);
	}
	
	NODES[NodeCnt].gene_id = id;
	NODES[NodeCnt].func = lptr;
	NodeCnt++;

	if (NodeCnt >= MAXGENES) {
		printf("Too many nodes\n");
		exit(0);
	}
}



List *get_bool_sub(char str[],int *idx,int len)
{
	int	id,flag;
	char	ch1,ch2,gname;
	List	*lptr1,*lptr2,*lptr3;

	flag = 0;
	while(*idx < len) {
		ch1 = str[*idx];
		if (ch1 == ' ') { 
			(*idx)++;
			continue;
		}
		else if (ch1 == ')') { 
			(*idx)++;
			return lptr1;
		}
		else if ((ch1 >= 'a' && ch1 <= 'z') || (ch1 >= 'A' && ch1 <= 'Z') || (ch1 >= '0' && ch1 <= '9')) {
			lptr1 = &LISTPOOL[ListNum];
if (ListNum >= MAXLIST-1) { printf("Use larger MAXLIST\n"); exit(0); }
			ListNum++;
			id = get_gene_id(str,idx);
			lptr1->type = TYPE_GENE;
			lptr1->gene_id = id;
		}
		else if (ch1 == '&') {
			(*idx)++;
			lptr2 = &LISTPOOL[ListNum];
if (ListNum >= MAXLIST-1) { printf("Use larger MAXLIST\n"); exit(0); }
			ListNum++;
			lptr2->type = TYPE_AND;
			lptr2->left = lptr1;
			lptr3 = get_bool_sub(str,idx,len);
			lptr2->right = lptr3;
			return lptr2;
		}
		else if (ch1 == '|') {
			(*idx)++;
			lptr2 = &LISTPOOL[ListNum];
if (ListNum >= MAXLIST-1) { printf("Use larger MAXLIST\n"); exit(0); }
			ListNum++;
			lptr2->type = TYPE_OR;
			lptr2->left = lptr1;
			lptr3 = get_bool_sub(str,idx,len);
			lptr2->right = lptr3;
			return lptr2;
		}
		else if (ch1 == '!') {
			(*idx)++;
			if (str[*idx] != '(') {
				printf("Error1 in %s\n",str);
				exit(0);
			}
			(*idx)++;
			lptr1 = &LISTPOOL[ListNum];
if (ListNum >= MAXLIST-1) { printf("Use larger MAXLIST\n"); exit(0); }
			ListNum++;
			lptr3 = get_bool_sub(str,idx,len);
			lptr1->left = lptr3;
			lptr1->type = TYPE_NOT;
		}
		else if (ch1 == '(') {
			(*idx)++;
			lptr1 = get_bool_sub(str,idx,len);
			flag = 1;
		}

	}
	if (flag == 1) {
		return lptr1;
	}
	printf("Error2 in %s\n",str);
	exit(0);
}
	


int get_gene_id(char str[],int *idx)
{
		
	int	i;
	char	ch1;
	char	gname[MAXNAMELEN];

	i = 0;
	while(1) {
		ch1 = str[*idx];
		if ((ch1 >= 'a' && ch1 <= 'z') || (ch1 >= 'A' && ch1 <= 'Z') || (ch1 >= '0' && ch1 <= '9')) {
			gname[i] = ch1;
			i++;
			(*idx)++;
			if (i >= MAXNAMELEN-1) {
				printf("Too long name: %s\n",str);
				exit(0);
			}
		}
		else {
			break;
		}
	}
	gname[i] = '\0';

	for(i = 0;i < GeneNum;i++) {
		if (strcmp(GENENAMES[i],gname)==0) {
			break;
		}
	}
	if (i == GeneNum) {
		printf("No name: %s\n",gname);
		printf("str: %s\n",str);
		exit(0);
	}
	return i;
}

int output_periodic_attractor(FILE *fp, int period)
{
	int	i,j,k,gid,t;
	List	*lptr;

	for(i = 0;i < NodeCnt;i++) {
		STATE0[i] = STATE1[i];
	}
	for(k = 0;k < period;k++) {
		fprint_states(fp);
		for(i = 0;i < NodeCnt;i++) {
			gid = NODES[i].gene_id;
			lptr = NODES[i].func;
			j = calc_next_val(lptr);
			STATE2[gid] = j;
		}
		for(i = 0;i < NodeCnt;i++) {
			STATE1[i] = STATE2[i];
		}
	}
	fprintf(fp,"------- End of periodic attractor ---------\n");
	for(i = 0;i < NodeCnt;i++) {
		if (STATE1[i] != STATE0[i]) {
			break;
		}
	}
	if (i != NodeCnt) {
		fprintf(fp,"Error was found: the above attractor is inconsistent\n");
		return 0;
	}
	return 1;
}


/************************************
Sort Probabilities
************************************/

int sort_prob()
{
	int	i,j,k,difp,n1,n2,n3,k1max,k2max,k3max,k1,k2,k3;
	double	p,p1,p2,p3;


	for(i = 0;i < GeneNum;i++) {
		IDX[i] = i;
	}
	msort_real(PROB,GeneNum,IDX);
	difp = 0;
	for(j = 0;j < MAXDIFPROB;j++) {
		NUMPROB[j] = 0;
	}
	j = 0;
	for(i=1;i < GeneNum;i++) {
		p1 = PROB[IDX[i-1]];
		p2 = PROB[IDX[i]];
		if (p2 != p1) {
			NUMPROB[difp] = i-j;
			DIFPROB[difp] = p1;
			j = i;
			difp++;
			if (difp > 3) {
				printf("Too many probability values: %d\n",difp);
				exit(1);
			}
		}
	}
	NUMPROB[difp] = GeneNum-j;
	DIFPROB[difp] = PROB[IDX[GeneNum-1]];
	difp++;
/*
if (difp < 3) {
printf("#different probabilities is less than 3: %d\n",difp);
exit(1);
}
*/

	n1 = NUMPROB[0];
	n2 = NUMPROB[1];
	n3 = NUMPROB[2];
	N1 = n1;
	N2 = n1+n2;
	N3 = n1+n2+n3;
	p1 = DIFPROB[0];
	p2 = DIFPROB[1];
	p3 = DIFPROB[2];
	if ((p3 < 0.5 && n3 > 0) || (p2 < 0.5 && n2 > 0) || p1 < 0.5) {
		printf("Too low probability was specified: %lf %lf %lf\n",p1,p2,p3);
		exit(1);
	}
	if (NUMPROB[0] < MAXCHANGEBITS) { k1max = NUMPROB[0]; }
	else { k1max = MAXCHANGEBITS; }
	if (NUMPROB[1] < MAXCHANGEBITS) { k2max = NUMPROB[1]; }
	else { k2max = MAXCHANGEBITS; }
	if (NUMPROB[2] < MAXCHANGEBITS) { k3max = NUMPROB[2]; }
	else { k3max = MAXCHANGEBITS; }
printf("k1max: %d  prob: %lf\n",k1max,p1);
printf("k2max: %d  prob: %lf\n",k2max,p2);
printf("k3max: %d  prob: %lf\n",k3max,p3);

/* assuming p1 > p2 > p3 */

	CNTP = 0;

	for(k1 = 0;k1 <= k1max;k1++) {
		if (p1 == 1.0 && k1 > 0) {
			continue;
		}
		if(k2max == 0) {
			if (p1 == 1.0) {
				p = 0.0;
			}
			else {
				p = k1*log(1.0-p1) + (n1-k1)*log(p1);
			}
			PROB2[CNTP] = p;
			IDX2[CNTP] = CNTP;
			IDXK1[CNTP] = k1;
			IDXK2[CNTP] = 0;
			IDXK3[CNTP] = 0;
			CNTP++;
			continue;
		}
		for(k2 = 0;k2 <= k2max;k2++) {
			if(k3max == 0) {
				if (p1 == 1.0) {
					p = k2*log(1.0-p2)+(n2-k2)*log(p2);
				}
				else {
					p = k1*log(1.0-p1) + (n1-k1)*log(p1) + k2*log(1.0-p2) + (n2-k2)*log(p2);
				}
				PROB2[CNTP] = p;
				IDX2[CNTP] = CNTP;
				IDXK1[CNTP] = k1;
				IDXK2[CNTP] = k2;
				IDXK3[CNTP] = 0;
				CNTP++;
				continue;
			}
			for(k3 = 0;k3 <= k3max;k3++) {
				if (p1 == 1.0) {
					p = k2*log(1.0-p2)+(n2-k2)*log(p2) + k3*log(1.0-p3)+(n3-k3)*log(p3);
				}
				else {
					p = k1*log(1.0-p1) + (n1-k1)*log(p1) + k2*log(1.0-p2) + (n2-k2)*log(p2) + k3*log(1.0-p3) + (n3-k3)*log(p3);
				}
				PROB2[CNTP] = p;
				IDX2[CNTP] = CNTP;
				IDXK1[CNTP] = k1;
				IDXK2[CNTP] = k2;
				IDXK3[CNTP] = k3;
				CNTP++;
			}
		}
	}
	msort_real(PROB2,CNTP,IDX2);
/*
printf("sort ok: %d\n",CNTP);
	for(i = 0;i < 100 && i < CNTP;i++) {
		j = IDX2[i];
		printf("%2d : %-10.5lf  :  %d %d %d\n",j,PROB2[j],IDXK1[j],IDXK2[j],IDXK3[j]);
	}
	printf("\n");
*/

}

					
/* 2021Oct */
int is_blank(char str[])
{
	int	i,flag;

	flag = 1;
	for(i = 0; i < MAXLINESIZE;i++) {
		if (str[i] == '\0' || str[i] == '\n') {
			break;
		}
		if (str[i] != ' ' && str[i] != '\t') {
			flag = 0;
		}
	}
	return flag;
}
/* 2021Oct */

	
	
/************************************************
*						*
*		Sort				*
*						*
************************************************/

int msort_real(double dat[],int n,int idx[])
{
	int	i;

	for(i = 0;i < n;i++) {
		idx[i] = i;
	}
	msort_real_body(dat,n,idx);
}

int msort_real_body(double dat[],int n,int idx[])
{
	int	m;

	if (n < 2) {
		return 0;
	}
	m = n/2;
	msort_real_body(dat,m,idx);
	msort_real_body(dat,n-m,&idx[m]);
	merge_real(dat,m,n,idx);
}

int merge_real(double dat[],int m,int n,int idx[])
{
	static	int	idxtmp[MAXGENES];
	int	i,j,k;
	double	d1,d2;

	i = 0;
	j = m;
	
	for(k = 0;k < n;k++) {
		if (i == m) {
			idxtmp[k] = idx[j];
			j++;
		}
		else if (j == n) {
			idxtmp[k] = idx[i];
			i++;
		}
		else {
			d1 = dat[idx[i]];
			d2 = dat[idx[j]];
			if (d1 >= d2) {
				idxtmp[k] = idx[i];
				i++;
			}
			else {
				idxtmp[k] = idx[j];
				j++;
			}
		}
	}
	for(i = 0;i < n;i++) {
		idx[i] = idxtmp[i];
	}
}

