#pragma once
#pragma once
#pragma once
//
/**/
/**/
#ifndef _nsE_H_
#define _nsE_H_
#include <iostream>
#include <vector>
#include"individual.h"
using namespace std;

class nich
{
public:

	nich();
	nich(individual& a, int size, double menadis);
	nich(double f);
	int psize;
	int esize;
	vector<double> l;
	vector<double> b;
	~nich();
	double mdis;
	int best;
	int st;
	double bestobj;
	double Fac;
	vector<individual> npop;
	vector<individual> cpop;
	vector<individual> chpop;
	vector<individual> elitepop;
	vector<individual>  gbest;
	vector<double> mean;
	vector<double> var;
	void  revol(CEC2013* pfunc);
	void repar(individual& a);
	void find_lb();
	//individual  best;
	double nicdis;
	vector<individual> copypop;
	vector<individual> gpop;
	void iteration_n(CEC2013* pfunc);
	//void sortfit(vector<individual> &pop, vector<int> &sort);
	//void cal_distance(vector<individual>& pop, vector<int>& sort);

	void gen_a_ind();
	//	void select_elite(vector<individual>& pop, int psize, vector<individual> &cpop);

	void find_n(double mdis, vector<individual>& bpop, individual& aind, vector<individual>& cpop);
	void operator=(const nich& a);

	void select_with_dist_f(vector<individual>& pop, int psize, vector<individual>& cpop);

	void fillup(CEC2013* pfunc);
	void nich_evol();
	void nich_select(vector<individual>& cpop, int psize);
	void cal_v_m();
	void nich_gen();
	void bestr2(individual& a);
	void bestr1(individual& a);
	void iteration_ind(CEC2013* pfunc);

	//	vector<double> var;
	//	vector<double> mean;

	void nevol(CEC2013* pfunc);


	void ievol(CEC2013* pfunc);
	void	find_best();


	void   select_fitness(int psize, vector<individual>& cpop);

};
nich::nich()
{

}
nich::nich(individual& a, int size, double meandis)
{

	psize = size;
	Fac = (1.0 * fes) / maxfes;
	mdis = meandis;
	gpop.clear();
	best = 0;
}
nich::nich(double f)
{

}
void nich::repar(individual& a)
{
	for (int i = 0; i < nvar; i++)
	{
		if (a.x[i] < l[i])
		{
			double rnd = rnd_uni(&rnd_uni_init);
			//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);

			a.x[i] = l[i] + rnd_uni(&rnd_uni_init) * (l[i] - a.x[i]);
		}
		if (a.x[i] > b[i])
		{
			double rnd = rnd_uni(&rnd_uni_init);
			//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);

			a.x[i] = b[i] - rnd_uni(&rnd_uni_init) * (a.x[i] - b[i]);
		}
	}

}
void nich::find_lb()
{
	l.clear();
	b.clear();
	if (npop.size() <= 1)
	{
		for (int i = 0; i < nvar; i++)
		{
			double max = lbound[i];
			double min = ubound[i];

			max = npop[0].x[i] + (ubound[i] - lbound[i]) / pops;
			min = npop[0].x[i] - (ubound[i] - lbound[i]) / pops;

			

			if (max > ubound[i])
				max = ubound[i];
			if (min < lbound[i])
				min = lbound[i];
			l.push_back(min);
			b.push_back(max);
		}
	}
	else
	{

		for (int i = 0; i < nvar; i++)
		{
			double max = lbound[i];
			double min = ubound[i];
			for (int j = 0; j < npop.size(); j++)
			{
				if (npop[j].x[i] > max)
					max = npop[j].x[i];
				if (npop[j].x[i] < min)
					min = npop[j].x[i];
			}

			double mid = (max - min) / 5;
			if (max - min < (ubound[i] - lbound[i]) / pops)
			{
				max = max+(ubound[i] - lbound[i]) / pops;
				min=min- (ubound[i] - lbound[i]) / pops;
			}
			max = max + mid;
			min = min - mid;
			if (max > ubound[i])
				max = ubound[i];
			if (min < lbound[i])
				min = lbound[i];
			l.push_back(min);
			b.push_back(max);
		}
	}

	for (int i = 0; i < l.size(); i++)
	{
		double ldis = (ubound[i] - lbound[i]) / pops;;
		if (fabs(b[i] - l[i]) < ldis)
		{
			double mid = (b[i] - l[i]) / 2;
			b[i] = mid + ldis;
			l[i] = mid - ldis;
		}
	}

}
nich::~nich()
{
	npop.clear();
	copypop.clear();
	gbest.clear();
	var.clear();
	mean.clear();
	gpop.clear();
}




void  nich::operator=(const nich& a)
{
	int n;

	Fac = a.Fac;

	psize = a.psize;
	npop = a.npop;
	gbest = a.gbest;
}








void nich::find_n(double mdis, vector<individual>& bpop, individual& aind, vector<individual>& cpop)
{
	cpop.clear();
	for (int i = 0; i < bpop.size(); i++)
	{
		if (aind.caldist(bpop[i]) < mdis && bpop[i].dis != 0)
		{
			cpop.push_back(bpop[i]);
		}
	}
}

void nich::nich_select(vector<individual>& cpop, int psize)
{
	copypop.clear();
	double c;
	c = fabs(cos(3.1415 * 8 * Factor));

	if (cpop.size() > psize)
	{
		double yita1 = 1 - Factor * Factor * Factor;
		for (int i = 0; i < cpop.size(); i++)
		{
			int d;


			if (cpop[i].dis == 0.0)
				cpop[i].lfit = 0.0;
			else

				//cpop[i].lfit = cpop[i].l1 + pow(cpop[i].l2, 1 + nvar * c);
				cpop[i].lfit = cpop[i].l1;

		}
		vector<int> sort;
		vector<int> lsort;
		for (int i = 0; i < cpop.size(); i++)
			lsort.push_back(i);

		for (int i = 0; i < lsort.size(); i++)
		{

			for (int j = i + 1; j < lsort.size(); j++)
			{

				if (cpop[lsort[i]].lfit < cpop[lsort[j]].lfit)
				{

					int dchange = lsort[i];
					lsort[i] = lsort[j];
					lsort[j] = dchange;
				}
			}
		}
		for (int i = 0; i < psize; i++)
		{

			copypop.push_back(cpop[lsort[i]]);
		}



		npop.clear();

		npop = copypop;
		copypop.clear();
	}


	mean.clear();
	for (int i = 0; i < nvar; i++)
	{
		double mean1 = 0;
		for (int j = 0; j < npop.size(); j++)
			mean1 += npop[j].x[i];
		mean1 = mean1 / npop.size();
		mean.push_back(mean1);
	}
	var.clear();
	for (int i = 0; i < nvar; i++)
	{
		double var1 = 0;
		for (int j = 0; j < npop.size(); j++)
			var1 += (npop[j].x[i] - mean[i]) * (npop[j].x[i] - mean[i]);
		var1 = sqrt(var1 / npop.size());
		var.push_back(var1);
	}


	/*
	copypop.clear();
	copypop.push_back(npop[0]);
	vector<individual> temp;
	for (int i = 0; i < npop.size(); i++)
		temp.push_back(npop[i]);

	vector<int> perc;
	vector<int> perc1;
	for (int i = 0; i < temp.size(); i++)
		perc.push_back(i);
	for (int i = 0; i < perc.size(); i++)
	{
		int ch1 = (int)(rnd_uni(&rnd_uni_init) * perc.size());
		int ch2 = (int)(rnd_uni(&rnd_uni_init) * perc.size());
		int ach = perc[ch1];
		perc[ch1] = perc[ch2];
		perc[ch2] = ach;
	 ch1 = (int)(rnd_uni(&rnd_uni_init) * perc.size());
		 ch2 = (int)(rnd_uni(&rnd_uni_init) * perc.size());
		int ach = perc1[ch1];
		perc1[ch1] = perc1[ch2];
		perc1[ch2] = ach;

	}
	for (int i = 0; i < perc.size(); i++)
	{
		int ch1 = (int)(rnd_uni(&rnd_uni_init) * perc.size());
		int ch2 = (int)(rnd_uni(&rnd_uni_init) * perc.size());


	}
	int c = 0;
	while (copypop.size() < psize)
	{

		copypop.push_back(temp[perc[c]]);
		c = c + 1;
	}
	npop.clear();
	npop = copypop;
	copypop.clear();
	*/
}

void nich::iteration_n(CEC2013* pfunc)
{
	copypop.clear();

	for (int i = 0; i < npop.size(); i++)
	{
		int r1, r2, r3, r4, r5;
		r1 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		r2 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		r3 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		r4 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		r5 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		while (r1 == i)
		{
			r1 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}
		while (r2 == i && r2 == r1)
		{
		r2 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}

		while (r3 == i && r2 == r3 && r3 == r1)
		{
			r3 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}
		while (r4 == i && r4 == r3 && r4 == r1 && r4 == r2)
		{
			r4 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}
		while (r5 == i && r5 == r3 && r5 == r1 && r5 == r2 && r5 == r4)
		{
			r5 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}


		
		individual child;

		individual child1, child2, child3;
		child1.compsite_de(npop[i], npop[r1], npop[r2], npop[r3], npop[r4], npop[r5], child1, child2, child3);

		//child1.repair();
	//child2.repair();
		//child3.repair(); 
		child1.repair_vec(l, b); child2.repair_vec(l, b); child3.repair_vec(l, b);
		child1.repair();
		child2.repair();
		child3.repair();
		child1.evaluation(pfunc);
		child2.evaluation(pfunc);
		child3.evaluation(pfunc);

		if (child1.obj < child2.obj)
			child = child2;
		else
			child = child1;
		if (child.obj < child3.obj)
			child = child3;
		copypop.push_back(child);
		chpop.push_back(child1);
		chpop.push_back(child2);
		chpop.push_back(child3);
	}
/*
	for (int i = 0; i < npop.size(); i++)
	{
		if (copypop[i].obj > npop[i].obj)
			npop[i] = copypop[i];
	}
*/
	copypop.clear();
}
void nich::select_fitness(int psize, vector<individual>& cpop)
{

	if (cpop.size() > psize)
	{


		vector<int> lsort;
		double yita;
		yita = (1 - Factor * Factor);
		yita = yita * yita * yita * yita * yita;
		for (int i = 0; i < cpop.size(); i++)
			lsort.push_back(i);
		for (int i = 0; i < lsort.size(); i++)
		{
			if (cpop[i].l2 == 0)

				cpop[i].lfit = -111000000000;
			cpop[i].lfit =(1-yita) *cpop[i].l1 + (yita) * cpop[i].l2 / nvar;


		}
		for (int i = 0; i < lsort.size(); i++)
		{

			for (int j = i + 1; j < lsort.size(); j++)
			{

				if (cpop[lsort[i]].lfit < cpop[lsort[j]].lfit)
				{

					int dchange = lsort[i];
					lsort[i] = lsort[j];
					lsort[j] = dchange;
				}
			}
		}
		for (int i = 0; i < psize; i++)
		{

			copypop.push_back(cpop[lsort[i]]);
		}



		npop.clear();

		npop = copypop;
		copypop.clear();

	}

}







void nich::find_best()
{
	best = 0;
	for (int i = 1; i < npop.size(); i++)
		if (npop[i].obj > npop[best].obj)
			best = i;
	bestobj = npop[best].obj;
}



void nich::fillup(CEC2013* pfunc)
{
	if (npop.size() < psize)
	{
		if (Factor < 0.001)
			Factor = 0.001;
		while (npop.size() < psize)
		{
			individual child;
			child = npop[0];
			double scal;

		

			for (int i = 0; i < nvar; i++)

			//child.x[i] = npop[0].x[i] + rnd_uni(&rnd_uni_init) * (ubound[i] - lbound[i]) / pops;
			child.x[i] =l[i]+ rnd_uni(&rnd_uni_init)*(b[i]-l[i]);
			child.repair();
			child.evaluation(pfunc);
			npop.push_back(child);
		
		}

	}
	
}



void nich::nevol(CEC2013* pfunc)
{
	int cf = 0;
	double cur_mdis;


	int evolsize;




	find_lb();


	if (npop.size() > psize)
	{
		select_fitness(psize, npop);
	}
	else
	{
		fillup(pfunc);
	}

	int iter = 0;
	fr << "~~~~~~~~2222222222222222222222222222~~~~~~~~~~the iter is generate~~22222222222222222222222222~~~~~~~~~" << endl;
	for (int i = 0; i < npop.size(); i++)
	{
		for (int j = 0; j < nvar; j++)
			fr << npop[i].x[j] << "	";
		fr << npop[i].obj << endl;
	}


	fr << endl;
	fr << gp << "	the gp " << nvar << endl;
	//while (iter <sqrt(nvar)* ceil(sqrt(pops)))
	//while (iter <= 1+ceil(sqrt(nvar)))
	//while (iter <=sqrt(nvar))
	int icount = 0;


	icount =ceil( sqrt(nvar*1.0));
	icount = 1;
	while (iter < icount)

	{
		fr << iter << endl;
		find_best();
		iteration_n(pfunc);
		fr << endl;
		fr << endl;
		for (int i = 0; i < npop.size(); i++)
		{
			for (int j = 0; j < nvar; j++)
				fr << npop[i].x[j] << "	";
			fr << npop[i].obj << endl;
		}
		iter++;
	}

}
void nich::revol(CEC2013* pfunc)
{
	int cf = 0;
	double cur_mdis;


	int evolsize;




	find_lb();


	if (npop.size() > psize)
	{
		select_fitness(psize, npop);
	}
	else
	{
		fillup(pfunc);
	}
	
	int iter = 0;
	fr << "~~~~~~~~2222222222222222222222222222~~~~~~~~~~the iter is generate~~22222222222222222222222222~~~~~~~~~" << endl;
	for (int i = 0; i < npop.size(); i++)
	{
		for (int j = 0; j < nvar; j++)
			fr << npop[i].x[j] << "	";
		fr << npop[i].obj << endl;
	}


	fr << endl;
	fr << gp << "	the gp " << nvar << endl;
	//while (iter <sqrt(nvar)* ceil(sqrt(pops)))
	//while (iter <= 1+ceil(sqrt(nvar)))
	//while (iter <=sqrt(nvar))
	int icount = 0;


	icount = sqrt(pops);
	while (iter < icount)

	{
		fr << iter << endl;
		find_best();
		iteration_n(pfunc);
		fr << endl;
		fr << endl;
		for (int i = 0; i < npop.size(); i++)
		{
			for (int j = 0; j < nvar; j++)
				fr << npop[i].x[j] << "	";
			fr << npop[i].obj << endl;
		}
		iter++;
	}

}
void nich::iteration_ind(CEC2013* pfunc)
{

	int i=(rnd_uni(&rnd_uni_init) * npop.size());
	
	{
		int r1, r2, r3, r4, r5;
		r1 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		r2 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		r3 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		r4 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		r5 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		while (r1 == i)
		{
			r1 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}
		while (r2 == i && r2 == r1)
		{
			r2 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}

		while (r3 == i && r2 == r3 && r3 == r1)
		{
			r3 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}
		while (r4 == i && r4 == r3 && r4 == r1 && r4 == r2)
		{
			r4 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}
		while (r5 == i && r5 == r3 && r5 == r1 && r5 == r2 && r5 == r4)
		{
			r5 = (int)(rnd_uni(&rnd_uni_init) * npop.size());
		}


		individual child;
		//repar(child1);
		//repar(child2);
		//repar(child3);
	/*	if (gp < sqrt(pops))
		{
			individual child1, child2, child3;
		child1.compsite_de(npop[i], npop[r1], npop[r2], npop[r3], npop[r4], npop[r5], child1,child2,child3);
			child1.evaluation(pfunc);
			child2.evaluation(pfunc);
			child3.evaluation(pfunc);

			if (child1.obj < child2.obj)
				child = child2;
			else
				child = child1;
			if (child.obj < child3.obj)
				child = child3;

		}

		else
		{

			child.compsite_de1(npop[i], npop[r1], npop[r2], npop[r3], npop[r4], npop[r5], child);
			//repar_vec(child);
			child.evaluation(pfunc);
		}
	*/
	/*	child.compsite_de1(npop[i], npop[r1], npop[r2], npop[r3], npop[r4], npop[r5], child);

		//	child.compsite_de2(npop[i], npop[r1], npop[r2], npop[r3], npop[r4], npop[r5], npop[best],child);
		child.repair_vec(l,b);
	child.evaluation(pfunc);
		copypop.push_back(child);
		*/
		//	gpop.push_back(child1);
		//	gpop.push_back(child2);
		//	gpop.push_back(child3);


		individual child1, child2, child3;
		child1.compsite_de(npop[i], npop[r1], npop[r2], npop[r3], npop[r4], npop[r5], child1, child2, child3);
	
		child1.evaluation(pfunc);
		child2.evaluation(pfunc);
		child3.evaluation(pfunc);

		if (child1.obj < child2.obj)
			child = child2;
		else
			child = child1;
		if (child.obj < child3.obj)
			child = child3;
		if (child.obj > npop[i].obj)
			npop[i] = child;

	}

	
}
#endif