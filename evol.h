
#pragma once
#include"individual.h"
#include"node.h"
#include"CEC2013.h"
#include<iostream>

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
	vector<individual>  gbest;

	int evolflag;

	int num;
	int niche_size;
	vector<niche> niching;


	void sortfit(vector<individual> &pop, vector<int> &sort);

	void cal_distance(vector<individual>& pop, vector<int>& sort);

	void select_elite(vector<individual>& pop, vector<int>& sort);

	void iteration(CEC2013* pfunc);
	void evalpop(CEC2013* pfunc, vector<individual>& pop);
	void statisticpop(vector<individual>& pop);
	int calopt(double curu, CEC2013* pFunc);
	void run(CEC2013* pfunc);

	
	double bestfit;
	double worstfit;
	
	ofstream frank;

/*
	void diff(int pp, int r1, int r2, int r3, individual& child);
	void diff2(int pp, int r1, int r2, int r3, int r4, int r5, individual& child);
	void gbestde(individual& xbest, individual& x1, individual& x2, individual& child);
	void gbestde2(individual& xbest, individual& x1, individual& x2, individual& x3, individual& x4, individual& child);
	*/
	void nich_divide(CEC2013 *pfunc);
	void nich_evol(CEC2013 *pfunc);
};



evol::evol()
{

	stepevol = 0;
	popsize = pops;
	evolflag = 0;

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

	niche_size = pops / 20;
	sprintf_s(filename1, " .%dnichtest.txt", ID);
	frank.open(filename1, std::ios::out);
}
evol::~evol()
{

}


void evol::sortfit(vector<individual>& pop, vector<int> &sort)
{

	int d = 0;
	for (int i = 0; i < pop.size(); i++)
	{
		//	int change = (int)(rnd_uni(&rnd_uni_init) * pop.size());
		//	int d = sort[i];
		//sort[i] = sort[change];
		//	sort[change] = d;
		//	sort.push_back(i);

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
				pop[sort[i]].grade = (pop[sort[j]].obj - pop[sort[i]].obj) / (bestof - worstof);
			}


		}
		pop[sort[i]].dis = dis_r;
	}

	//cout<<"	"<<pop[sort[0]].dis <<"	"<< pop[sort[10]].dis<< "	" << pop[sort[20]].dis << "	" << pop[sort[30]].dis<<endl;
}

void evol::select_elite(vector<individual>& pop, vector<int>& sort)
{

	double yita;

	yita = Factor;
	yita = 0.95 *Factor;
	cal_distance(pop, sort);
	double av = 0;

	double maxd = 0;
	for (int i = 1; i < pop.size(); i++)
	{

		pop[sort[i]].l2 = pop[sort[i]].dis;
		av += pop[sort[i]].dis;
		pop[sort[i]].l1 = (pop[sort[i]].obj - worstof) / (bestof - worstof);

		//cout << avfd << "	";
		if (maxd < pop[sort[i]].dis)
			maxd = pop[sort[i]].dis;
	}
	pop[sort[0]].grade = 0;
	pop[sort[0]].dis = maxd;
	pop[sort[0]].l1 = (pop[sort[0]].obj - worstof) / (bestof - worstof);
	pop[sort[0]].l2 = pop[sort[0]].dis;
	av = (av + maxd) / (pop.size());





	for (int i = 0; i < pop.size(); i++)
	{

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
			pop[i].lfit = (1 - yita)*(pop[i].l2) / nvar + yita * pop[i].l1;

	}




	for (int i = 0; i < sort.size(); i++)
	{

		for (int j = i + 1; j < sort.size(); j++)
		{

			if (pop[sort[i]].lfit < pop[sort[j]].lfit)
			{

				int dchange = sort[i];
				sort[i] = sort[j];
				sort[j] = dchange;
			}
		}


	}

	gbest.clear();
	for (int i = 0; i < pop.size() / 2; i++)
	{
		if (pop[sort[i]].label == 1)
		{
			gbest.push_back(pop[sort[i]]);

		}
	}
	cout << "gbest" << gbest.size() << endl;


}




void evol::evalpop(CEC2013* pfunc, vector<individual>& pop)
{
	for (int i = 0; i < pop.size(); i++)
	{
		pop[i].evaluation(pfunc);
		//	fl << i<<"	"<<pop[i].obj << "	";
		
	}
	//fl << endl << endl;;
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
			//	fl << "the obj of£º	" << pop[i].obj << "	"<<i<<"	";

			worstof = pop[i].obj;
			//	fl << "the worst of£º	" << worstof << endl;
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

void evol::run(CEC2013* pfunc)
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
	int flag = -1;
	evolflag = 0;
	while (fes < maxfes)

	{
		Factor = (fes*1.0) / maxfes;
		
		iter++;
		cout << "	iter <<:	" << iter;

		if (gbest.size() < parentpop.size() / 20)
			flag = 1;
		cout << "	flag <<:	" << flag;
		if (flag == -1||evolflag==0)
			evolflag = 0;

		if (flag == 1 &&evolflag == 0)
		{
			evolflag = 1;
		
		}
		else
		{
			if (flag == 1 && evolflag == 1)
			{
			
				evolflag = 2;
			}
		}
		
		cout << "evolflag" << evolflag <<"	:"<<fes<< endl;
		
		if (evolflag == 0)
		{
			iteration(pfunc);
			evalpop(pfunc, childpop);
			if (mixedpop.size() != 0)
				mixedpop.clear();
			for (int j = 0; j < childpop.size(); j++)
				mixedpop.push_back(childpop[j]);
			for (int j = 0; j < parentpop.size(); j++)
				mixedpop.push_back(parentpop[j]);
			parentpop.clear();
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
			if (iter % 10 == 0)
			{
				frank << iter << endl;
				frank << calopt(0.1, pfunc) << "	" << calopt(0.01, pfunc) << endl << endl;
				for (int i = 0; i < pops; i++)
				{
					//	frank << parentpop[i].x[0] << "	" << parentpop[i].x[1] << "	" << parentpop[i].obj << "	" << parentpop[i].l1 << "	" << parentpop[i].l2 << "	" << parentpop[i].lfit << "	" << parentpop[i].rank << endl;
				}
				frank << endl << endl;
				frank << "the gblest" << endl;
				for (int i = 0; i < gbest.size(); i++)
				{
					frank << gbest[i].x[0] << "	" << gbest[i].x[1] << "	" << gbest[i].obj << "	" << gbest[i].l1 << "	" << gbest[i].l2 << "	" << gbest[i].lfit << "	" << gbest[i].rank << endl;
				}


				cout << iter << "		";

				cout << calopt(0.1, pfunc) << "	" << calopt(0.01, pfunc) << "	" << calopt(0.001, pfunc) << endl;
			}

		}
		
		
	
			if (evolflag == 1)
			{
				
			
				nich_divide(pfunc);
				evolflag = 2;
				system("pause");
			}
			
			
		

		 if (evolflag == 2)
		{
			cout << "thr" << endl;
			fl << "the iter:	" << iter << endl;
			nich_evol(pfunc);
		}

		

	}
	if (evolflag != 0)
	{
		parentpop.clear();
		frank << "the last niching" << endl;
		for (int i = 0; i < niching.size(); i++)
		{
			for (int j = 0; j < niching[i].ns; j++)
			{
				parentpop.push_back(niching[i].parentpop[j]);
				frank << "i	" << i << "	" << parentpop[i].x[1] << "	" << parentpop[i].x[2] << "	" << parentpop[i].obj << endl;
			}
		}
	}
	
}







void evol::iteration(CEC2013* pfunc)
{

	vector<int> l;
	if (childpop.size() != 0)
		childpop.clear();
	int pt = 0;
	double maxd = 1000000000;
	/*int count = 0;
	for (int i = 0; i < gbest.size(); i++)
	{
		if (gbest[i].obj < maxd)
		{
			maxd = gbest[i].obj;
			pt = i;
		}
	}
	cout <<"the bestworst===" <<gbest[pt].obj << endl;
	*/


	for (int i = 0; i < popsize; i++)
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

		while (r2 == i && r2 == r3 && r3 == r1)
		{
			r3 = (int)(rnd_uni(&rnd_uni_init) * parentpop.size());
		}
			   		 		
			individual child;
			diff(parentpop[i], parentpop[r1], parentpop[r2], parentpop[r3], child);
			childpop.push_back(child);
		
			}



}



void evol::nich_divide(CEC2013 *pfunc)
{
	int numc = gbest.size();
	cout << "it  is dived cal num" << endl;
	for (int i = 0; i < numc;i++)
	{
		niche a;
		a.parentpop.clear();
		a.parentpop.push_back(gbest[i]);
		double maxdis = 1000000000.0;
		int pt;
		cout << i << "	";
		for (int j = 0; j < gbest.size(); j++)
		{
			cout << j << "j	";
			if (i != j)
			{
				double d = gbest[i].caldist(gbest[j]);
				if (d < maxdis)
				{
					pt = j;
					maxdis = d;
				}
			}
		
		}
		for (int k = 0; k < nvar; k++)
		{
			a.radius.push_back(gbest[i].x[k]-gbest[pt].x[k]);
		}

		niching.push_back(a);

	}
	cout << "it  is dived cal niche" << endl;
//class		



	for (int i = 0; i < parentpop.size(); i++)
	{
		double max = 10000000.0;
		int pt = 0;
		for (int j = 0; j < niching.size(); j++)
		{
			double d = gbest[j].caldist(parentpop[i]);
			if (d < max)
			{
				max = d;
				pt = j;
			}
		}
		niching[pt].parentpop.push_back(parentpop[i]);
	}
	frank << "the divide@@@@@@@@@@@pop@@@@#£¤£¤£¤£¤£¤£¤£¤£¤£¤£¤£¤%%%%%%%%%%" << endl;
	for (int i = 0; i < parentpop.size(); i++)
	{
		

			frank << "i	" << i << "	" << mixedpop[i].x[0] << "	" << mixedpop[i].x[1] << "	" << mixedpop[i].obj << endl;
		
	}


	cout << "the divide" << endl;
	for (int i = 0; i < niching.size(); i++)
	{
		for (int j = 0; j < niching[i].parentpop.size(); j++)
		{
		
			frank << "i	" << i << "	" << niching[i].parentpop[j].x[0] << "	" << niching[i].parentpop[j].x[1] << "	" << niching[i].parentpop[j].obj << endl;
		}
	}

	//
	int n_size = parentpop.size() / niching.size();
	
	for (int i = 0; i < niching.size(); i++)
	{
		cout << "itis dived and increase" << endl;
		while (niching[i].parentpop.size() > n_size)
		{
			niching[i].parentpop.pop_back();
		}
		while (niching[i].parentpop.size() < n_size)
		{
			niching[i].sample(pfunc);
		}
		niching[i].ns = n_size;
	}
	cout << "the divide  before" << endl;
	for (int i = 0; i < niching.size(); i++)
	{
		for (int j = 0; j < niching[i].parentpop.size(); j++)
		{

			frank << "i	" << i << "	" << niching[i].parentpop[j].x[0] << "	" << niching[i].parentpop[j].x[1] << "	" << niching[i].parentpop[j].obj << endl;
		}
	}


	}

void evol::nich_evol(CEC2013 *pfunc)
{
	int num = niching.size();
	cout << "evols" << num<< endl;
	
	for (int i = 0; i < num; i++)
	{
		fl << "the iniches:	"<<i << endl;
		niching[i].devol(pfunc);
		niching[i].select();
	}



}
