#pragma once
#pragma once


#pragma once
#ifndef _pn_H_
#define _pn_H_
#include"individual.h"
#include"niches.h"
#include"CEC2013.h"
#include<iostream>
#include <iomanip> 
#include <cmath> 
using namespace std;

class evol
{
public:

	evol();
	~evol();


	int popsize;
	vector<individual> parentpop;
	vector<individual>  childpop;
	vector<individual>  mixedpop;

	vector<individual> copop;
	vector<individual>  gbest;


	vector<nich> niching;
	double  meandis;
	double varp;
	int evolflag;
	int nn;
	int num;
	int c_num;
	double pmdis;
	int  iter;            //
	double nichdis;
	double mdis;
	double  cmdis;
	int csolutions;
	int step;
	void sortfit(vector<individual>& pop, vector<int>& sort);
	void cal_distance(vector<individual>& pop, vector<int>& sort);
	void select_elite(vector<individual>& pop, vector<int>& sort, CEC2013* pfunc);
	void  report_niching(int pt);

	double gyita;

	void select_best();
	vector<individual>  localpop;
	void iteration_p(CEC2013* pfunc);
	void iteration(CEC2013* pfunc);
	void evalpop(CEC2013* pfunc, vector<individual>& pop);
	void statisticpop(vector<individual>& pop);
	int calopt(double curu, CEC2013* pFunc);
	void run(CEC2013* pfunc);

	void nevol(CEC2013* pfunc);

	double bestfit;
	double worstfit;
	void iteration_replace(CEC2013* pfunc);
	void iteration1(CEC2013* pfunc);

	void  generate(CEC2013* pfunc);

	void mergepop();
	void iteration_anich(CEC2013* pfunc);

	void clustering(int csize, vector<individual>& pop, vector<int>& sort);
	void find_n(int csize, vector<individual>& pop, vector<int>& sort);

	void sample(vector<individual>& lpop, CEC2013* pfunc);
	ofstream frank;
	ofstream fniche;
	ofstream fevol;
	void report(int i);
};



evol::evol()
{


	popsize = pops;

	step = 1;
	for (int i = 0; i < popsize; i++)
	{
		individual a;
		a.rnd_init();
		a.dis = nvar * 1000000;

		parentpop.push_back(a);

	}
	/////////////////////caldistall
	dist_all = 0;
	for (int i = 0; i < nvar; i++)
		dist_all += (ubound[i] - lbound[i]) * (ubound[i] - lbound[i]);
	dist_all = sqrt(dist_all);

	bestof = -100000000000;
	worstof = 10000000000;
	iter = 0;
	char filename1[1024];
	nn = sqrt(pops);
	meandis = 0.0;
	sprintf_s(filename1, " %d.pop.txt", ID);
	frank.open(filename1, std::ios::out);

	char filename2[1024];
	nn = sqrt(pops);

	sprintf_s(filename2, " .%dnich66.txt", ID);

	fniche.open(filename2, std::ios::out);


}
evol::~evol()
{

	parentpop.clear();
	childpop.clear();
	gbest.clear();

}


void evol::sortfit(vector<individual>& pop, vector<int>& sort)
{

	int d = 0;
	for (int i = 0; i < pop.size(); i++)
	{

		sort.push_back(i);

	}
	for (int i = 0; i < pop.size(); i++)
	{
		for (int j = i + 1; j < pop.size(); j++)
		{
			if (pop[sort[j]].obj > pop[sort[i]].obj)
			{
				int dchange = sort[i];
				sort[i] = sort[j];
				sort[j] = dchange;
			}
		}

	}

}

void evol::cal_distance(vector<individual>& pop, vector<int>& sort)
{
	pop[sort[0]].dis = dist_all;
	pop[sort[0]].parent = -1;
	for (int i = 1; i < sort.size(); i++)
	{
		double dis_r = 1000000000000.0;

		for (int j = 0; j < i; j++)
		{
			double	d = pop[sort[i]].caldist(pop[sort[j]]);
			if (d <= dis_r)
			{
				dis_r = d;
				pop[sort[i]].parent = sort[j];

			}

		}
		pop[sort[i]].dis = dis_r;

	}


}



void evol::select_elite(vector<individual>& pop, vector<int>& sort, CEC2013* pfunc)
{

	double yita;
	yita = (1 - Factor * Factor);
	yita = yita * yita * yita * yita * yita;



	if (yita < 0.001 )
		yita = 0.001;
	


	cal_distance(pop, sort);
	double av = 0;

	double maxd = 0;  //找除了初始条件之外，距离信息最优的个体；
	for (int i = 1; i < pop.size(); i++)
	{

		pop[sort[i]].l2 = pop[sort[i]].dis;
		av += pop[sort[i]].dis;
		pop[sort[i]].l1 = (pop[sort[i]].obj - worstof) / (bestof - worstof);

		//cout << avfd << "	";
		if (maxd < pop[sort[i]].dis)
			maxd = pop[sort[i]].dis;
	}
	
	pop[sort[0]].pdis = 0;
	pop[sort[0]].dis = dist_all;
	pop[sort[0]].l1 = (pop[sort[0]].obj - worstof) / (bestof - worstof);
	pop[sort[0]].l2 = pop[sort[0]].dis;
	av = (av) / (pop.size());

	if (iter == 0)
	{

		meandis = av;
	}

	cmdis = av;








	for (int i = 0; i < pop.size(); i++)
	{

	pop[sort[i]].l2 = pop[sort[i]].dis /meandis;
	//	pop[sort[i]].l2 = pop[sort[i]].dis ;
	}



	


	for (int i = 0; i < pop.size(); i++)
	{


		if (pop[i].dis == 0.0)
			pop[i].lfit = 0.0;
		else
			pop[i].lfit = (yita)*pop[i].l2 /nvar+ (1 - yita) * pop[i].l1;		
	}

	vector<int> lsort;
	lsort = sort;

	for (int i = 0; i < lsort.size(); i++)
	{

		for (int j = i + 1; j < lsort.size(); j++)
		{

			if (pop[lsort[i]].lfit < pop[lsort[j]].lfit)
			{

				int dchange = lsort[i];
				lsort[i] = lsort[j];
				lsort[j] = dchange;
			}
		}
	}
	int count = 0;




	parentpop.clear();
	for (int i = 0; i < popsize; i++)
	{

		parentpop.push_back(pop[lsort[i]]);
	}
	mdis = dist_all;

	vector<int> dsort;

	for (int i = 0; i < parentpop.size(); i++)
	{
		if (pop[lsort[i]].dis >  meandis)
			dsort.push_back(lsort[i]);    //选择gbest个体
	}
	gp = dsort.size();
	if (gp < 2*ceil(sqrt(pops))&&gp>ceil(sqrt(pops)))
	{
		
			step = 2;
		
	}
	if ( gp<=ceil(sqrt(pops)))
	{

		step =3;

	}
	else
	{
		step == 1;
	}
	
	{
		if (dsort.size() > ceil(sqrt(popsize)))

		{
			for (int i = 0; i < dsort.size(); i++)
				pop[dsort[i]].lfit = pop[dsort[i]].l2 * yita / nvar + pop[dsort[i]].l1;


			for (int i = 0; i < dsort.size(); i++)
			{
				for (int j = i + 1; j < dsort.size(); j++)
				{
					if (pop[dsort[i]].lfit < pop[dsort[j]].lfit)
					{
						int change = dsort[i];
						dsort[i] = dsort[j];
						dsort[j] = change;
					}
				}
			}
			gbest.clear();
			for (int i = 0; i < (sqrt(popsize)); i++)
				gbest.push_back(pop[dsort[i]]);
		}
		else
		{
			gbest.clear();
			for (int i = 0; i < dsort.size(); i++)
				gbest.push_back(pop[dsort[i]]);
		}

		//clustering

	
		for (int i = 0; i < pop.size(); i++)
		{
			pop[i].rank = -1;
		}

		for (int i = 0; i < gbest.size(); i++)
			pop[dsort[i]].rank = i;
		for (int i = 0; i < sort.size(); i++)
		{

			if (pop[sort[i]].rank == -1)
			{
				int p = pop[sort[i]].parent;
				pop[sort[i]].rank = pop[p].rank;
			}
		}
		//	for (int i = 0; i < sort.size(); i++)
		//	{

			//	cout << pop[sort[i]].rank << "	";
			//}
		niching.clear();
		for (int i = 0; i < gbest.size(); i++)
		{
			nich a;
			int d1 = popsize / gbest.size();
			int d2 = sqrt(pops);
			if (d1 > d2)
				a.psize = popsize / gbest.size();
			else
				a.psize = ceil(sqrt(pops));
			a.psize = ceil(sqrt(pops));
			if (a.psize < nvar)
			a.psize = nvar;
			
			a.mdis = meandis;
			for (int j = 0; j < pop.size(); j++)
			{
				if (pop[sort[j]].rank == i && pop[sort[j]].dis != 0)
				{
					a.npop.push_back(pop[sort[j]]);
				}
			}
			niching.push_back(a);
		}
		for (int i = 0; i < niching.size(); i++)
		{
			niching[i].select_fitness(niching[i].psize, niching[i].npop);
		}
		
		fniche << endl << endl;
/*
		fniche << "gbest:" << endl;
		for (int i = 0; i < gbest.size(); i++)
		{
			for (int j = 0; j < nvar; j++)
				fniche << gbest[i].x[j] << "	";
			fniche << gbest[i].obj << endl;
		}

		for (int i = 0; i < niching.size(); i++)

		{
			fniche << "the nichi:" << i << endl;
			report_niching(i);
			niching[i].select_fitness(niching[i].psize, niching[i].npop);
			fniche << "the select:" << endl;
			report_niching(i);

		}
		*/
	
	}
}





void evol::evalpop(CEC2013* pfunc, vector<individual>& pop)
{
	for (int i = 0; i < pop.size(); i++)
	{
		pop[i].evaluation(pfunc);

	}
}

void evol::statisticpop(vector<individual>& pop)
{

	maxobj = pop[0].obj;
	minobj = pop[0].obj;

	for (int i = 0; i < pop.size(); i++)
	{
		//	fl << i << "	";
		if (maxobj < pop[i].obj)
			maxobj = pop[i].obj;
		if (minobj > pop[i].obj)
			minobj = pop[i].obj;
		if (bestof < pop[i].obj)
			bestof = pop[i].obj;
		if (worstof > pop[i].obj)
		{
			worstof = pop[i].obj;

		}

	}

}

int evol::calopt(double curu, CEC2013* pFunc)
{
	std::vector< std::vector< double> > pop;
	for (int i = 0; i < parentpop.size(); i++)
		pop.push_back(parentpop[i].x);
	std::vector< std::vector<double> > seeds;


	return  how_many_goptima(pop, seeds, pFunc, curu, pFunc->get_rho());

}









void evol::iteration(CEC2013* pfunc)     //迭代
{

	childpop.clear();
	double per1;


	
	{
		for (int i = 0; i < parentpop.size(); i++)
		{
			int r1, r2, r3, r4, r5;
			r1 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			r2 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			r3 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			r4 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			r5 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			while (r1 == i)
			{
				r1 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			}
			while (r2 == i && r2 == r1)
			{
				r2 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			}

			while (r3 == i && r2 == r3 && r3 == r1)
			{
				r3 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			}
			while (r4 == i && r4 == r3 && r4 == r1)
			{
				r4 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			}

			while (r5 == i && r5 == r3 && r5 == r1 && r5 == r2 && r5 == r4)
			{
				r5 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
			}
		
			individual child1, child2, child3;
			child1.compsite_de(parentpop[i], parentpop[r1], parentpop[r2], parentpop[r3], parentpop[r4], parentpop[r5], child1, child2, child3);
			childpop.push_back(child1);
			childpop.push_back(child2);
			childpop.push_back(child3);
			
			

		/*
					individual child;
			child.compsite_de1(parentpop[i], parentpop[r1], parentpop[r2], parentpop[r3], parentpop[r4], parentpop[r5], child);
			childpop.push_back(child);
			
	*/

		}



		evalpop(pfunc, childpop);

	}








}

void evol::iteration_anich(CEC2013* pfunc)
{
	fniche << "the  gen" << endl<<endl<<endl;
	for (int i = 0; i < niching.size(); i++)
	{
		//	iteration_anich(pfunc);
		fniche << "the first"  << endl;
		report_niching(i);
		niching[i].nevol(pfunc);
		fniche << "the second" << endl;
		report_niching(i);
	}
}
void evol::iteration_replace(CEC2013* pfunc)
{
	
	for (int i = 0; i < niching.size(); i++)
	{
		//	iteration_anich(pfunc);
		
		
		niching[i].revol(pfunc);
		fniche << "the second" << endl;
		
	}
}
void evol::mergepop()     //融合
{
	mixedpop.clear();



	for (int i = 0; i < childpop.size(); i++)
		mixedpop.push_back(childpop[i]);
	
	childpop.clear();

	for (int i = 0; i < parentpop.size(); i++)
		mixedpop.push_back(parentpop[i]);
	//if (step == 2)
	{
		for (int i = 0; i < niching.size(); i++)

		{
		
			for (int j = 0; j < niching[i].npop.size(); j++)
			mixedpop.push_back(niching[i].npop[j]);
			//	for (int j = 0; j < niching[i].copypop.size(); j++)
			for (int j = 0; j < niching[i].chpop.size(); j++)
						mixedpop.push_back(niching[i].chpop[j]);
		}
	}
	statisticpop(mixedpop);
	childpop.clear();
}
void  evol::generate(CEC2013* pfunc)  //生成子种群
{
	
	mergepop();
	cout << "mixsize	" << mixedpop.size()<<"	";
	vector<int> sort;
	sortfit(mixedpop, sort);

	select_elite(mixedpop, sort, pfunc);
	

	if (step == 3)
	{
		iteration_anich(pfunc);
	}

	else
		iteration(pfunc);
	/*
if (step == 1)
{
	iteration(pfunc);
}
else

{
	if(rnd_uni(&rnd_uni_init)<0.5)
		iteration_anich(pfunc);
	else
		iteration(pfunc);

}
	*/
}



void evol::run(CEC2013* pfunc)
{


	fes = 0;
	seed = (seed + 23) % 1377;
	rnd_uni_init = -(long)seed;
	evalpop(pfunc, parentpop);
	statisticpop(parentpop);


	vector<int> calo;
	int vargen = 0;
	//	fl << ID << endl << endl;
	gbest = parentpop;
	int thr;	
	step = 1;
	

	while (fes < maxfes)
	{
		cout << iter << "	the fes " << fes <<"	"<<step << endl;

		Factor = fes * 1.0 / (maxfes * 1.0);
		generate(pfunc);
		cout << gbest.size();



		//	if (step == 2)
			//	iteration_anich(pfunc);
			//	iteration_p(pfunc);
		
		iter = iter + 1;
		//	cout << iter << endl;

		if (iter % 10 == 0)
		{
			fniche << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
			fniche << nichdis << "	" << meandis << endl;;

			frank << iter << endl;
			//	frank << calopt(0.1, pfunc) << "	" << calopt(0.01, pfunc) << endl << endl;
			frank << meandis << endl;
			//cout << "parentpop" << parentpop.size() << endl;
			for (int i = 0; i < parentpop.size(); i++)
			{
				if (nvar > 1)
					frank << setiosflags(ios::fixed) << setprecision(12) << parentpop[i].x[0] << "	" << parentpop[i].x[1] << "	" << parentpop[i].obj << "	" << parentpop[i].l1 << "	" << parentpop[i].l2 << "	" << parentpop[i].lfit << "	" << parentpop[i].rank << endl;

				//	frank << setiosflags(ios::fixed) << setprecision(12) << parentpop[i].x[0] << "	" << parentpop[i].x[1] << "	" << parentpop[i].obj << "	" << endl;

			}
			frank << endl << endl;
			frank << "gebst:::" << endl;
			for (int i = 0; i < gbest.size(); i++)
			{
				if (nvar > 1)
					frank << setiosflags(ios::fixed) << setprecision(12) << gbest[i].x[0] << "	" << gbest[i].x[1] << "	" << gbest[i].obj << endl;
			}
			frank << "the dist is:	" << meandis << "	" << cmdis << endl;
		}
	}
	




}


void  evol::report_niching(int pt)
{
	if (pt > niching.size())
	{
		return;
	}
	else
	{
		fniche << endl << endl;
		fniche << niching[pt].npop.size() << " the ptth	" <<pt<< endl;
		for (int i = 0; i < niching[pt].npop.size(); i++)
		{
			for (int j = 0; j < nvar; j++)
				fniche << niching[pt].npop[i].x[j] << "	";
			fniche << niching[pt].npop[i].obj << "	";
			fniche << endl;
		}
	}

}


#endif