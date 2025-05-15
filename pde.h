#pragma once
#pragma once


#pragma once
#ifndef _cn_H_
#define _cn_H_
#include"individual.h"
#include"niches.h"
#include"CEC2013.h"
#include<iostream>
#include <iomanip> 
using namespace std;
class  problem
{
	int nx;
	
	double  evalution(vector<double> x);
	double   dychange;

	


};
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
	void localrefine(CEC2013* pfunc);
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

	void iteration1(CEC2013* pfunc);
	void iteration_pop(CEC2013* pfunc);
	void iteration_niches(CEC2013* pfunc);
	void iteration_pop_gbest(CEC2013* pfunc);
	void  generate(CEC2013* pfunc);
	void mergepop_cp();
	void mergepop_np();
	void mergepop(int st);
	void iteration_anich(CEC2013* pfunc);

	void clustering(int csize, vector<individual>& pop, vector<int>& sort);
	void find_n(int csize, vector<individual>& pop, vector<int>& sort);

	void sample(vector<individual>& lpop, CEC2013* pfunc);
	ofstream frank;
	ofstream fniche;
	ofstream fevol;
	void report(int i);
	void clustering(vector<individual>& pop, vector<int>& sort, vector<int>& dsort);
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
	sprintf_s(filename2, " .%dnch12.txt", ID);
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
	if (yita < 0.001 / nvar)
		yita = 0.001 / nvar;
	cal_distance(pop, sort);
	double av = 0;
	double maxd = 0;
	for (int i = 1; i < pop.size(); i++)
	{

		pop[sort[i]].l2 = pop[sort[i]].dis;
		av += pop[sort[i]].dis;
		pop[sort[i]].l1 = (pop[sort[i]].obj - worstof) / (bestof - worstof);

		
		if (maxd < pop[sort[i]].dis)
			maxd = pop[sort[i]].dis;
	
		for (int j = 0; j < nvar; j++)
		{
			int p = pop[sort[i]].parent;
			pop[sort[i]].grade.push_back(pop[p].x[j] - pop[sort[i]].x[j]);
		}
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

		pop[sort[i]].l2 = pop[sort[i]].dis / meandis;

	}



	for (int i = 0; i < pop.size(); i++)
	{
	
		pop[i].rank = -1;
		if (pop[i].dis > av)
			pop[i].label = 1;
		else
			pop[i].label = -1;
	}



	for (int i = 0; i < pop.size(); i++)
	{


		if (pop[i].dis == 0.0)
			pop[i].lfit = 0.0;
		else
			pop[i].lfit = pop[i].l2  *pow( pop[i].l1 ,0.1+2*Factor*Factor) ;
		

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
		if (pop[lsort[i]].dis > meandis)
			dsort.push_back(lsort[i]);
	}

	//decide the seed and the number
	gp = dsort.size();
	
	if (dsort.size() > ceil(sqrt(popsize)))

	{
		//	for (int i = 0; i < dsort.size(); i++)
				//pop[dsort[i]].lfit = pop[dsort[i]].l2 * yita + pop[dsort[i]].l1;

			//	pop[dsort[i]].lfit = pop[dsort[i]].l2;
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



	cout << gbest.size() << endl;
	
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







void evol::localrefine(CEC2013* pfunc)
{
	nich a;
	a.psize = sqrt(pops);
	int r = (rnd_uni(&rnd_uni_init) * (sqrt(pops))+(1-fes/(maxfes*1.0))*pops);
	if(r>=pops)
		r = (rnd_uni(&rnd_uni_init) * ((pops-1)) );

	r = (rnd_uni(&rnd_uni_init) * (sqrt(pops)));
	a.npop.push_back(parentpop[r]);
	for (int j = 0; j < nvar; j++)
		fniche <<parentpop[r].x[j] << "	 ";
	fniche << parentpop[r].obj << "	the s "<<endl;
	a.ievol(pfunc);
	for (int i = 0; i < a.npop.size(); i++)
	{
	//	childpop.push_back(a.npop[i]);
	}
	for (int i = 0; i < a.npop.size(); i++)
	{
		for (int j = 0; j < nvar; j++)
			fniche <<a.npop[i].x[j] << "	";
		fniche << a.npop[i].obj << "	";

		fniche << endl;
	}
	a.find_best();
	parentpop[r] = a.npop[a.best];
} 



void evol::iteration(CEC2013* pfunc)     //迭代
{

	childpop.clear();
	double per1;
	localrefine(pfunc);
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
			//compside with 3 child 

			individual child1, child2, child3;
			child1.compsite_de(parentpop[i], parentpop[r1], parentpop[r2], parentpop[r3], parentpop[r4], parentpop[r5], child1, child2, child3);
			childpop.push_back(child1);
			childpop.push_back(child2);
			childpop.push_back(child3);

			/*		individual child;
					child.compsite_de1(parentpop[i], parentpop[r1], parentpop[r2], parentpop[r3], parentpop[r4], parentpop[r5], child);
					childpop.push_back(child);
		*/
		}

		evalpop(pfunc, childpop);

	}

	




}

void evol::iteration_anich(CEC2013* pfunc)
{
	for (int i = 0; i < niching.size(); i++)
		niching[i].select_fitness(sqrt(pops), niching[i].npop);
	for (int i = 0; i < niching.size(); i++)
	{
		cout << i << endl;
		fniche << "stat" << endl;
		report_niching(i);
		niching[i].nevol(pfunc);
		report_niching(i);
		fniche << "end" << endl << endl;
	}
}
void evol::mergepop_cp()     //融合
{


	mixedpop.clear();
	for (int i = 0; i < childpop.size(); i++)
		mixedpop.push_back(childpop[i]);
	for (int i = 0; i < parentpop.size(); i++)
		mixedpop.push_back(parentpop[i]);

}
void evol::mergepop_np()
{
	mixedpop.clear();


	for (int i = 0; i < parentpop.size(); i++)
		mixedpop.push_back(parentpop[i]);
	for (int i = 0; i < niching.size(); i++)
	{
		for (int j = 0; j < niching[i].npop.size(); j++)
			mixedpop.push_back(niching[i].npop[j]);
	}

}
void  evol::generate(CEC2013* pfunc)  //生成子种群
{

	
	
		iteration(pfunc);
		statisticpop(childpop);
		mergepop_cp();
	
	vector<int> sort;
	sortfit(mixedpop, sort);
	select_elite(mixedpop, sort, pfunc);



}

void evol::run(CEC2013* pfunc)
{


	fes = 0;
	seed = (seed + 23) % 1377;
	rnd_uni_init = -(long)seed;
	evalpop(pfunc, parentpop);
	statisticpop(parentpop);

	Factor = fes * 1.0 / (maxfes * 1.0);
	vector<int> calo;
	int vargen = 0;

	/*	for (int i = 0; i < niching.size(); i++)
		{
			report_niching(i);
		}
		*/
	for (int i = 0; i < niching.size(); i++)
	{
		report_niching(i);
	}

	while (fes < maxfes)
	{
		cout << iter << "	" << fes << endl;

		Factor = fes * 1.0 / (maxfes * 1.0);

		cout << gbest.size();

		generate(pfunc);

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

	for (int i = 0; i < gbest.size(); i++)
		parentpop.push_back(gbest[i]);

	for (int i = 0; i < childpop.size(); i++)
		parentpop.push_back(childpop[i]);



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
		fniche << pt << endl;
		fniche << niching[pt].npop.size() << "	" << endl;
		if (niching[pt].l.size() != 0)
			for (int i = 0; i < nvar; i++)
			{

				fniche << niching[pt].l[i] << "	 lb " << niching[pt].b[i] << endl;
			}
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