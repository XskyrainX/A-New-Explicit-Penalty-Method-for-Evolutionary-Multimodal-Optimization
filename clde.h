#ifndef __GLOBAsL_H_
#define __GLOBAsL_H_

#pragma once
#pragma once
#pragma once
#include"individual.h"
#include"node.h"
#include"CEC2013.h"
#include<iostream>

using namespace std;
class cl
{
public:
	vector<individual> cp;
	int psize;
	cl();
	~cl();
	void select();
	void insertind(vector<int>&l, vector<individual> &mpop);
};
cl::cl()
{}
cl::~cl()
{
	cp.clear();
}
void cl::select()
{
	if (cp.size() > psize)
	{
		vector<int> sort;
		for (int i = 0; i < cp.size(); i++)
		{
			sort.push_back(i);
		}
		for (int i = 0; i < cp.size(); i++)
		{
			for (int j = i + 1; j < cp.size(); j++)
				if (0.7*cp[sort[i]].dis + 0.3*cp[sort[i]].obj < 0.7*cp[sort[i]].dis + 0.3*cp[sort[i]].obj)
				{
					int cg = sort[i];
					sort[i] = sort[j];
					sort[j] = cg;

				}
		}
		vector<individual> temp;
		for (int i = 0; i < psize; i++)
		{
			temp.push_back(cp[sort[i]]);
		}
		cp.clear();
		cp = temp;
		temp.clear();
	}
}
void cl::insertind(vector<int>&l, vector<individual> &mpop)
{
	cp.clear();
	for (int i = 0; i < l.size(); i++)
		cp.push_back(mpop[l[i]]);
}
class FDE
{
public:

	FDE();
	~FDE();


	int popsize;
	vector<individual> parentpop;
	vector<individual>  childpop;
	vector<individual>  mixedpop;
	vector<individual>  lpop;
	vector<individual>  gbest;
	vector<cl> icl;
	void sortfit(vector<individual> &pop, vector<int> &sort);
	int it;
	void cal_distance(vector<individual>& pop, vector<int>& sort);

	void select_elite(vector<individual>& pop, vector<int>& sort);

	void iteration(CEC2013* pfunc);
	void evalpop(CEC2013* pfunc, vector<individual>& pop);
	void statisticpop(vector<individual>& pop);
	int calopt(double curu, CEC2013* pFunc);
	void run(CEC2013* pfunc);

	double dist_all;
	double bestfit;
	double worstfit;
	double bestof;
	double worstof;
	ofstream frank;

	void sample(individual &a, vector<individual> &pop, CEC2013 *pfunc);
	void diff(int pp, int r1, int r2, int r3, individual& child);
	void diff2(int pp, int r1, int r2, int r3, int r4, int r5, individual& child);
	void gbestde(individual& xbest, individual& x1, individual& x2, individual& child);
	void gbestde2(individual& xbest, individual& x1, individual& x2, individual& x3, individual& x4, individual& child);
	void ls_bb(vector<individual> &cpop, vector<individual> &lp);
	void lserach(vector<individual> &lp, vector<individual> &cpop, CEC2013 *pfunc);



	void cluster();



};



FDE::FDE()
{

	stepevol = 0;
	popsize = pops;

	it = 0;
	for (int i = 0; i < popsize; i++)
	{
		individual a;
		a.rnd_init();
		parentpop.push_back(a);
	}
	/////////////////////caldistall
	dist_all = 0;
	for (int i = 0; i < nvar; i++)
		dist_all += (ubound[i] - lbound[i]) * (ubound[i] - lbound[i]);
	dist_all = sqrt(dist_all);
	bestof = -100000000000;
	worstof = 10000000000;

	char filename1[1024];


	sprintf_s(filename1, " .%dtest.txt", ID);
	frank.open(filename1, std::ios::out);
}
FDE::~FDE()
{

}


void FDE::sortfit(vector<individual>& pop, vector<int> &sort)
{

	int d = 0;
	sort.clear();
	for (int i = 0; i < pop.size(); i++)
	{
		//	int change = (int)(rnd_uni(&rnd_uni_init) * pop.size());
		//	int d = sort[i];
		//sort[i] = sort[change];
		//	sort[change] = d;
		//	sort.push_back(i);
		pop[i].grade = 0;
		pop[i].cd = 0;
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
		//	cout << pop[sort[i]].obj<<"	";
	}
	//	frank <<"the obj"<< endl;
	//	for (int i = 0; i < pop.size(); i++)
		//	frank << pop[sort[i]].obj << endl;;
	//	frank <<"###################################"<< endl<<endl;
}

void FDE::cal_distance(vector<individual>& pop, vector<int>& sort)
{
	pop[sort[0]].dis = dist_all;
	pop[sort[0]].parent = -1;
	for (int i = 1; i < sort.size(); i++)
	{
		double dis_r = 1000000000000.0;
		int pt = 0;
		for (int j = 0; j < i; j++)
		{
			double	d = pop[sort[i]].caldist(pop[sort[j]]);
			if (d <= dis_r)
			{
				pt = j;
				dis_r = d;
				pop[sort[i]].parent = sort[j];


			}


		}
		int p = pop[sort[i]].parent;
		double g = (pop[p].obj - pop[sort[i]].obj) / dis_r;
		if (g > pop[p].grade&&dis_r != 0)
		{
			pop[p].grade = g;
			pop[p].cd = dis_r;
			pop[p].dx = pop[sort[i]].x;
		}




		pop[sort[pt]].child.push_back(sort[i]);
		pop[sort[i]].dis = dis_r;
	}

	//cout<<"	"<<pop[sort[0]].dis <<"	"<< pop[sort[10]].dis<< "	" << pop[sort[20]].dis << "	" << pop[sort[30]].dis<<endl;
}

void FDE::select_elite(vector<individual>& pop, vector<int>& sort)
{

	double yita;

	yita = Factor;
	yita = 0.99 *Factor;
	cal_distance(pop, sort);
	double av = 0;

	double maxd = 0;
	for (int i = 1; i < pop.size(); i++)
	{

		//cout << avfd << "	";
		if (maxd < pop[sort[i]].dis)
			maxd = pop[sort[i]].dis;

		//	cout << av << endl;
		pop[sort[i]].dis = pop[sort[i]].dis / dist_all;
		pop[sort[i]].l2 = pop[sort[i]].dis;

		pop[sort[i]].l1 = (pop[sort[i]].obj - worstof) / (bestof - worstof);
		av += pop[sort[i]].dis;

		//	cout << worstof << " fit	" << pop[sort[i]].l1 << endl;
	}
	//cout << "i==" << maxd << "	" << av << endl;
	pop[sort[0]].grade = 0;
	pop[sort[0]].dis = maxd / dist_all;
	pop[sort[0]].l1 = (pop[sort[0]].obj - worstof) / (bestof - worstof);
	pop[sort[0]].l2 = pop[sort[0]].dis;
	av = (av + maxd / dist_all) / (pop.size());

	//	cout << av << "	the dist	"<<av<< "	"<<maxd<<"	dist_all"<< dist_all<<endl;


	int seed = 0;
	for (int i = 0; i < pop.size(); i++)
	{
	
		if (pop[i].dis > av)
		{
			pop[i].label = 1;
			seed++;
		}
				
			else
				pop[i].label = -1;
			
	}
	frank << "the seed is:" << seed<<"the mean dis:	"<<av << endl;
	vector<int> lsort;
	
	for (int i = 0; i < pop.size(); i++)
	{
		lsort.push_back(i);

		if (pop[i].dis == 0.0)
			pop[i].lfit = 0.0;
		else
			pop[i].lfit = (1 - yita)*(pop[i].l2) + yita * pop[i].l1 ;

	}




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

	gbest.clear();
	/*	for (int i = 0; i <10; i++)
		{
			if (pop[sort[i]].label == 1)
			{
				gbest.push_back(pop[sort[i]]);

			}
		}
		*/
//	frank << "the best" << endl;
	icl.clear();
	for (int i = 0; i < 10; i++)
	{

		gbest.push_back(pop[lsort[i]]);
		cl a;
		a.cp.push_back(pop[lsort[i]]);
		icl.push_back(a);
		pop[lsort[i]].label = i;
	//	frank << pop[lsort[i]].x[0] << "	" << pop[lsort[i]].x[1] << "	" << pop[lsort[i]].obj << endl;
		if (it % 10 == 0)
		frank << gbest[i].x[0] << "	" << gbest[i].x[1] << "	" << gbest[i].obj << endl;
	//	frank << sort[i] << endl;

	}

	for (int i = 0; i < pop.size(); i++)
	{
		int pt;
		double dis = 1000000;
	//	cout << "tttt" << i << endl;;
		for (int j = 0; j < icl.size(); j++)
		{
			for (int k = 0; k < icl[j].cp.size(); k++)
			{
				double d = icl[j].cp[k].caldist(pop[sort[i]]);
				if (d < dis)
				{
					dis = d;
					pt = j;
				}
			}

		}
	//	cout << pt << "	" << endl;
		if (dis != 0)
		{
			//cout << pt;
			icl[pt].cp.push_back(pop[sort[i]]);
		}

	}




	for (int i = 1; i < pop.size(); i++)
	{
		if (pop[sort[i]].label == -1)
		{
			int d = pop[sort[i]].parent;
			pop[sort[i]].label = pop[d].label;
		//	cout << "#$$$$$$$$%^^^^^^^^^^^^the lable is:	" << pop[d].label << endl;
			
	//	frank << i << "parent:	" << pop[sort[i]].parent << "	" << pop[d].label << "	" << pop[sort[i]].label << endl;
		}
	}

	for (int i = 0; i < pop.size(); i++)
	{
		gbest[pop[i].label].child.push_back(i);

		//cout << "the bable is:	" << pop[i].label << endl;
	}
	
	for (int i = 0; i < gbest.size(); i++)
	{
	//	frank << endl;
	//	for (int j = 0; j < gbest[i].child.size(); j++)
	//		frank << gbest[i].child[j] << "	";
	}
		

	cluster();
	cout << "gbest" << gbest.size() << endl;
	sort = lsort;

}
void FDE::cluster()
{
	frank << endl << endl;;
	frank << it << endl;
	int dsize = popsize / gbest.size();
	dsize = 20;
	vector<individual> cpop;
	
	for (int i = 0; i < icl.size(); i++)
	{
		
		
		if (icl[i].cp.size() < dsize)
		{
			icl[i].psize = icl[i].cp.size();
		     
			while (icl[i].cp.size() < dsize)
			{
				icl[i].cp.push_back(mixedpop[(int)(rnd_uni(&rnd_uni_init)*mixedpop.size())]);
			}
		}
		else
		{
			icl[i].psize = dsize;
			
			icl[i].select();
		}
		if(it%10==0)
		{ 
			frank << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << i << endl;;
			for (int j = 0; j < icl[i].cp.size(); j++)
				frank << icl[i].cp[j].x[0] << "	" << icl[i].cp[j].x[1] << "	" << icl[i].cp[j].obj << "	" << icl[i].cp[j].dis << endl;

		
		}
		

	}




}

void FDE::diff2(int pp, int r1, int r2, int r3, int r4, int r5, individual& child)
{
	double rate;

	double F;
	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);

	vector<double> l;
	l.push_back(0.1);
	l.push_back(0.5);
	l.push_back(0.9);
	F = l[3 * rnd_uni(&rnd_uni_init)];


	double CR;
	l.clear();

	l.push_back(0.5);
	l.push_back(0.7);
	l.push_back(0.9);

	CR = l[3 * rnd_uni(&rnd_uni_init)];
	l.clear();
	CR = 0.7;
	F = 0.5;
	if (nvar < 2)
	{
		for (int i = 0; i < nvar; i++)
		{
			child.x[i] = parentpop[r1].x[i] + F * (parentpop[r2].x[i] - parentpop[r3].x[i]) + F * (parentpop[r4].x[i] - parentpop[r5].x[i]);
			if (child.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (parentpop[r1].x[i] - lbound[i]);
			}
			if (child.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = ubound[i] - rnd * (ubound[i] - parentpop[r1].x[i]);
			}
		}
	}
	else
	{
		for (int i = 0; i < nvar; i++)
		{
			if (rnd_uni(&rnd_uni_init) <= CR || i == j)
			{
				child.x[i] = parentpop[r1].x[i] + F * (parentpop[r2].x[i] - parentpop[r3].x[i]) + F * (parentpop[r4].x[i] - parentpop[r5].x[i]);


				if (child.x[i] < lbound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
					if (rnd < 0.5)
						child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
					else
						child.x[i] = lbound[i];
				}
				if (child.x[i] > ubound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
					if (rnd < 0.5)
						child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
					else
						child.x[i] = ubound[i];
				}

			}
			else
				child.x[i] = parentpop[pp].x[i];
		}

		/*	for(int i=0;i<nvar;i++)
			{
			double rnd=rnd_uni(&rnd_uni_init);
			rnd+=fabs((rnd-0.5));
			child.x[i]=parentpop[r1].x[i]+rnd*(parentpop[r2].x[i]-parentpop[r3].x[i]);
				  if (child.x[i]<lbound[i])
					  {
				  double rnd = rnd_uni(&rnd_uni_init);
			  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);

			  }
			  if( child.x[i]>ubound[i])
			  {
				  double rnd = rnd_uni(&rnd_uni_init);
				 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
			  }
			}
			*/
	}

	//	child.obj_eval();
}
void FDE::diff(int pp, int r1, int r2, int r3, individual& child)
{
	double rate;

	double F;
	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);

	vector<double> l;
	l.push_back(0.1);
	l.push_back(0.5);
	l.push_back(0.9);
	F = l[3 * rnd_uni(&rnd_uni_init)];


	double CR;
	l.clear();

	l.push_back(0.5);
	l.push_back(0.7);
	l.push_back(0.9);

	CR = l[3 * rnd_uni(&rnd_uni_init)];
	l.clear();
	CR = 0.7;
	F = 0.5;
	if (nvar < 2)
	{
		for (int i = 0; i < nvar; i++)
		{
			child.x[i] = parentpop[r1].x[i] + F * (parentpop[r2].x[i] - parentpop[r3].x[i]);
			if (child.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (parentpop[r1].x[i] - lbound[i]);
			}
			if (child.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = ubound[i] - rnd * (ubound[i] - parentpop[r1].x[i]);
			}
		}
	}
	else
	{
		for (int i = 0; i < nvar; i++)
		{
			if (rnd_uni(&rnd_uni_init) <= CR || i == j)
			{
				child.x[i] = parentpop[r1].x[i] + F * (parentpop[r2].x[i] - parentpop[r3].x[i]);


				if (child.x[i] < lbound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
					if (rnd < 0.5)
						child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
					else
						child.x[i] = lbound[i];
				}
				if (child.x[i] > ubound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
					if (rnd < 0.5)
						child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
					else
						child.x[i] = ubound[i];
				}

			}
			else
				child.x[i] = parentpop[pp].x[i];
		}

		/*	for(int i=0;i<nvar;i++)
			{
			double rnd=rnd_uni(&rnd_uni_init);
			rnd+=fabs((rnd-0.5));
			child.x[i]=parentpop[r1].x[i]+rnd*(parentpop[r2].x[i]-parentpop[r3].x[i]);
				  if (child.x[i]<lbound[i])
					  {
				  double rnd = rnd_uni(&rnd_uni_init);
			  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);

			  }
			  if( child.x[i]>ubound[i])
			  {
				  double rnd = rnd_uni(&rnd_uni_init);
				 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
			  }
			}
			*/
	}

	//	child.obj_eval();
}



void FDE::evalpop(CEC2013* pfunc, vector<individual>& pop)
{
	for (int i = 0; i < pop.size(); i++)
	{
		pop[i].obj = pfunc->evaluate(pop[i].x);
		//	fl << i<<"	"<<pop[i].obj << "	";
		fes++;
	}
	//fl << endl << endl;;
}

void FDE::statisticpop(vector<individual>& pop)
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
			//	fl << "the obj of£º	" << pop[i].obj << "	"<<i<<"	";

			worstof = pop[i].obj;
			//	fl << "the worst of£º	" << worstof << endl;
		}

	}

}

int FDE::calopt(double curu, CEC2013* pFunc)
{
	std::vector< std::vector< double> > pop;
	for (int i = 0; i < parentpop.size(); i++)
		pop.push_back(parentpop[i].x);
	std::vector< std::vector<double> > seeds;


	return  how_many_goptima(pop, seeds, pFunc, curu, pFunc->get_rho());

}

void FDE::run(CEC2013* pfunc)
{
	fes = 0;
	seed = (seed + 23) % 1377;
	rnd_uni_init = -(long)seed;


	evalpop(pfunc, parentpop);
	statisticpop(parentpop);
	//		similar=1.0;
	int iter = 0;
	vector<int> calo;
	int vargen = 0;
	fl << ID << endl << endl;
	gbest = parentpop;
	while (fes < maxfes)

	{

		cout << it << endl;
		Factor = (fes*1.0) / maxfes;
		//	cout << bestof << "	" << worstof << endl;
		iter++;
		cout << "	iter <<:	" << iter;
		if (iter == 60)
		{
			double d = 0;
		}
		iteration(pfunc);
		evalpop(pfunc, childpop);
		if (mixedpop.size() != 0)
			mixedpop.clear();
		for (int j = 0; j < childpop.size(); j++)
			mixedpop.push_back(childpop[j]);
		for (int j = 0; j < parentpop.size(); j++)
			mixedpop.push_back(parentpop[j]);
		for (int i = 0; i < icl.size(); i++)
		{
			for (int j = 0; j < icl[i].cp.size(); j++)
				mixedpop.push_back(icl[i].cp[j]);

		}
		parentpop.clear();
		it = it+1;
		statisticpop(mixedpop);
		vector<int> sort;
		for (int i = 0; i < mixedpop.size(); i++)
			sort.push_back(i);
		sortfit(mixedpop, sort);
		cal_distance(mixedpop, sort);
		//	cal_knn(mixedpop, sort);
		select_elite(mixedpop, sort);
		for (int i = 0; i < pops; i++)
			parentpop.push_back(mixedpop[sort[i]]);
		
		sort.clear();
		if (it % 10 == 0)
		{
			frank << iter << endl;
			frank << calopt(0.1, pfunc) << "	" << calopt(0.01, pfunc) << endl << endl;
		/*	for (int i = 0; i < pops; i++)
			{
				frank << parentpop[i].grade << "	dis:" << parentpop[i].dis << "cd:	" << parentpop[i].cd << endl;
			}
			*/
			for (int i = 0; i < pops; i++)
			{

				frank << parentpop[i].x[0] << "	" << parentpop[i].x[1] << "	" << parentpop[i].obj << "	" << parentpop[i].l1 << "	" << parentpop[i].l2 << "	" << parentpop[i].lfit << "	" << parentpop[i].rank << endl;
			}
			frank << endl << endl;
			frank << "the gblest" << endl;
			for (int i = 0; i < gbest.size(); i++)
			{
				//		frank <<gbest[i].x[0] << "	" << gbest[i].x[1] << "	" << gbest[i].obj << "	" << gbest[i].l1 << "	" << gbest[i].l2 << "	" << gbest[i].lfit << "	" << gbest[i].rank << endl;
			}


			cout << iter << "		";

			cout << calopt(0.1, pfunc) << "	" << calopt(0.01, pfunc) << "	" << calopt(0.001, pfunc) << endl;
		}




		if (iter % 20 == 0)
		{
			//		fl << endl << endl;
			//		fl <<"the proble:"<<ID<<"	iter:	" <<iter << "	:" << endl;
				//	fl << ubound[0] - lbound[0] << "	" << pfunc->get_fitness_goptima() << endl;


				//	fl<<"the yita factor:"<<Factor* Factor* Factor* (ubound[0] - lbound[0]) * 100.0 * fabs((worstof - bestof)) / fabs((minobj - maxobj))<<endl;
				//	fl << "the yita factor:" << Factor * Factor * Factor * (ubound[0] - lbound[0]) * 100.0 << endl;
				//	fl<< "obj:	" << bestof			<< "	worstof:	" << worstof << "	minobj:	" << minobj << "	maxobj:	" << maxobj<<endl;
				//	fl <<  "worsof-bestof:	" << fabs((worstof - bestof)) << "	max-minobj:	" << fabs(maxobj - minobj)<<"	the rate:	"<< fabs((worstof - bestof))/ fabs(maxobj - minobj)<<endl;
				//	fl << "the opt:	" << calopt(0.1, pfunc) << "	"<<calopt(0.01, pfunc) << "	"<< calopt(0.001, pfunc) << endl;
		}

	

	}



}







void FDE::iteration(CEC2013* pfunc)
{
	if (childpop.size() != 0)
		childpop.clear();
	while (childpop.size() < popsize)
	{
		int d;
		if (it == 0)
			d = 0;
		else
		 d=(int)(rnd_uni(&rnd_uni_init) * 3);
		
		//cout << endl;
	//	cout << childpop.size() << "	d=" << d << endl;
	if (d == 0)
	{
		int r, r1, r2, r3;
		r= (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
		r1 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
		r2 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
		r3 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
		individual child;
		diff(r, r1, r2, r3, child);
		childpop.push_back(child);
	}

	if (d==1)
	{
		int c, r, r1, r2, r3;
		c= (int)(rnd_uni(&rnd_uni_init) *  icl.size());
		r = (int)(rnd_uni(&rnd_uni_init) * icl[c].cp.size());
		r1 = (int)(rnd_uni(&rnd_uni_init) * icl[c].cp.size());
		r2 = (int)(rnd_uni(&rnd_uni_init) * icl[c].cp.size());
		r3 = (int)(rnd_uni(&rnd_uni_init) * icl[c].cp.size());
	
		individual child;
		icl[c].cp[r].diff(icl[c].cp[r], icl[c].cp[r1], icl[c].cp[r2], icl[c].cp[r3], child);
		childpop.push_back(child);
	}

	if (d == 2)
	{
		int c, r, r1, r2, r3;
		c = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());

		r = (int)(rnd_uni(&rnd_uni_init) *  icl.size());
		r1 = (int)(rnd_uni(&rnd_uni_init) * icl[r].cp.size());
		r2 = (int)(rnd_uni(&rnd_uni_init) * icl[r].cp.size());
		r3 = (int)(rnd_uni(&rnd_uni_init) * icl[r].cp.size());

		//cout << c<<"	" << r << "	" << r1 << "	" << r2 << "	" << r3 << endl;
		individual child;
		parentpop[c].diff(parentpop[c], parentpop[c], icl[r].cp[r2], icl[r].cp[r3], child);
		childpop.push_back(child);
	}

	}
	




}









#endif