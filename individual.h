#ifndef _IND_H_
#define _IND_H_

#include"global.h"
#include"CEC2013.h"
 #include <vector>

using namespace std;
class individual
{
public:
	individual();
	~individual();

	vector<  double> x;	
	vector<double> dx;
	
	int parent;
	double l1;
	double l2;
	double dis;
	double lfit;
	 double obj;
	 int kn;
	 vector<int> child;
	 vector<double> vl;
	int label;
	//void obj_eval();
	void rnd_init();
	int pdim;
	double  pdis;
    void show_objective();
	void diff(individual &pp, individual &r1, individual &r2, individual &r3, individual& child);
	void diff_r2(individual &pp, individual &r1, individual &r2, individual &r3, individual &r4, individual &r5, individual& child);
	int rank;
	double cd;
	double caldist(individual &ind1);
	 void operator=(const individual &ind2);
	 bool operator==(const individual &ind2);
	 void mutation(double yit);
	 void evaluation(CEC2013* pfunc);
	vector< double> grade;
	 void repair();
	 int islocal;
	 double ldis;
	 void diff1(individual &pp, individual &r1, individual &r2, individual &r3, individual& child);
	 void diff_r3(individual& pp, individual& r1, individual& r2, individual& r3, individual& child);
	 void diff_mut(individual &r2, individual &r3, individual& child);
	 double gaussrand();
	 void mut_gus(vector<double> &varl);
	// void evaluation(CEC2013* pfunc);
	 void compsite_de(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& r5, individual& child1, individual& child2, individual& child3);
	 void compsite_de1(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& r5, individual& child);
	 void  repair_vec(vector<double> lb, vector<double> ub);
	 double CR;
	 double F;
	 void grade_diff(individual& pp, individual& r1, individual& r2, individual& child1, individual& child2, individual& child3);
	 void diff_best_r3(individual& pp, individual& r1, individual& r2, individual& child);
	 void diff_best_r2(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& child);
	 void diff_best(individual& pp, individual& r1, individual& r2, individual& child);
	 void compsite_de_gbest(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& r5, individual& child1, individual& child2, individual& child3);
	 void de(individual& pp, individual& r1, individual& r2, individual& r3, individual& child);
	 void de2(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& r5, individual& child);
	 void de3(individual& pp, individual& r1, individual& r2, individual& r3, individual& child);
};




//void diff(individual &pp, individual &r1, individual &r2, individual &r3, individual& child);
/*
void diff2(individual &pp, individual &r1, individual &r2, individual &r3, individual &r4, individual &r5, individual& child)
{
}

void diff2_grade(individual &pp, individual &r1, individual &r2, individual &r3, individual &r4, individual &r5, individual& child)
{
}

*/

void individual::diff(individual &pp, individual &r1, individual &r2, individual &r3, individual& child)
{

	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);

	
	
	if (nvar < 2)
	{
		if(isnan(child.x[0])==1)
			child.x[0] = rnd_uni(&rnd_uni_init) * (ubound[0]-lbound[0] );
		for (int i = 0; i < nvar; i++)
		{
			child.x[i] = r1.x[i] + F * (r2.x[i] - r3.x[i]);
		//	fl << " 0	" << r1.x[i] <<"	"<< r2.x[i] <<"	"<<r3.x[i] << child.x[i];
			if (isnan(child.x[0]) == 1)
				child.x[0] = rnd_uni(&rnd_uni_init) * (ubound[0] - lbound[0]);
			if (child.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
			//	fl << "	1	" << child.x[i];
			}
			if (child.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
		//		fl << "2	" << child.x[i];
			}
		}
	}
	else
	{
		for (int i = 0; i < nvar; i++)
		{
			if (rnd_uni(&rnd_uni_init) <= CR || i == j)
			{
				child.x[i] = r1.x[i] + F * (r2.x[i] - r3.x[i]);


				if (child.x[i] < lbound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
					if (rnd < 1)
						child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
					else
						child.x[i] = lbound[i];
				}
				if (child.x[i] > ubound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
					if (rnd < 1)
						child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
					else
						child.x[i] = ubound[i];
				}

			}
			else
				child.x[i] = pp.x[i];
		}

	
	}

	
}


void individual::diff_r2(individual &pp, individual &r1, individual &r2, individual &r3, individual &r4, individual &r5, individual& child)
{
	

	
	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);
	if (nvar < 2)
	{
		for (int i = 0; i < nvar; i++)
		{
			child.x[i] = x[i] + F * (r2.x[i] - r3.x[i])+ F * (r5.x[i] - r4.x[i]);
			if (child.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) *(lbound[i] - child.x[i]);
			}
			if (child.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = ubound[i] - rnd * (ubound[i] - r1.x[i]);
			}
		}
	}
	else
	{
		for (int i = 0; i < nvar; i++)
		{
			if (rnd_uni(&rnd_uni_init) <= CR || i == j)
			{
				child.x[i] = r1.x[i] + F * (r2.x[i] - r3.x[i])+ F * (r4.x[i] - r5.x[i]);


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
				child.x[i] = pp.x[i];
		}


	}
}
void individual::diff_r3(individual& pp, individual& r1, individual& r2, individual& r3,  individual& child)
{
	double dstep = rnd_uni(&rnd_uni_init);
	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);
	for (int i = 0; i < nvar; i++)
	{

		{
			child.x[i] = pp.x[i] + dstep * (r1.x[i] - pp.x[i]) + F * (r3.x[i] - r2.x[i]);


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

				if (rnd < 0.5)
					child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
				else
					child.x[i] = ubound[i];
			}

		}
	}
}
void individual::diff1(individual &pp, individual &r1, individual &r2, individual &r3, individual& child)
{
	double rate;

	double F;
	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);


	double CR = 0.9;
	F = 0.7 ;
	if (nvar < 2)
	{
		for (int i = 0; i < nvar; i++)
		{
			child.x[i] = r1.x[i] + F * (r1.x[i] - r2.x[i])+ F * (r1.x[i] - r3.x[i]);
			if (child.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (r1.x[i] - lbound[i]);
			}
			if (child.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = ubound[i] - rnd * (ubound[i] - r1.x[i]);
			}
		}
	}
	else
	{
		for (int i = 0; i < nvar; i++)
		{
			if (rnd_uni(&rnd_uni_init) <= CR || i == j)
			{

				child.x[i] = r1.x[i] + F* (r2.x[i] - r3.x[i]) ;


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
				child.x[i] = pp.x[i];
		}


	}


}


void individual::diff_mut( individual &r2, individual &r3, individual& child)
{

	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);
	for (int i = 0; i < nvar; i++)
	{
		if (rnd_uni(&rnd_uni_init) <= 1.0 )
		{
		child.x[i] = x[i] + 0.5*(r3.x[i] - r2.x[i])*rnd_uni(&rnd_uni_init);
		//	child.x[i] = r2.x[i] + (x[i] - r2.x[i])*(2*rnd_uni(&rnd_uni_init));
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
			child.x[i] = x[i];
		
	}
}
individual::individual()
{

	for (int i = 0; i < nvar; i++)
	{
		x.push_back(0.0);
		grade.push_back(0.0);
	}
		
	
    	rank    = 0;
		
}

individual::~individual()
{
    x.clear();
	
}

void individual::rnd_init()
{
    for(int n=0;n<nvar;n++)
        x[n] = lbound[n] + rnd_uni(&rnd_uni_init)*(ubound[n] - lbound[n]);

}


void individual::show_objective()
{

	for(int n=0;n<nvar;n++)
		std::cout<<x[n]<<" ";
	std::cout<<rank<<" ";
	std::cout<<"\n";

}


bool individual::operator==(const individual& ind2)
{
	int flag = 0;
	for(int n=0;n<nvar;n++)
	{
	    if(ind2.x[n] !=x[n])
	        return false;
    }
    return true;
}
void  individual::operator=(const individual& ind2)
{   
	int n;
	
	if (ind2.x.size() != 0)
	{
		for (n = 0; n < nvar; n++)
			x[n] = ind2.x[n];
	}
	
	else
		cout<<"the error";
	obj=ind2.obj;
	dis = ind2.dis;
    rank  = ind2.rank;
	label = ind2.label;
	grade = ind2.grade;
	l1 = ind2.l1;
	l2 = ind2.l2;
	dx = ind2 .dx;
}




double individual::caldist(individual &ind1)
{
	double dist=0;
	for(int i=0;i<nvar;i++)
		dist+=(x[i]-ind1.x[i])*(x[i]-ind1.x[i]);
	return sqrt(dist);
}
void individual::mutation(double yi)
{
	
	
	for (int i = 0; i < nvar; i++)
	{
		if (rnd_uni(&rnd_uni_init)<1)
		{

			double rnd = gaussrand();
			
			x[i] = x[i] + rnd * yi;

			if (x[i] < lbound[i])
				x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(lbound[i] - x[i]);
			if (x[i] > ubound[i])
				x[i] = ubound[i] - rnd_uni(&rnd_uni_init) * (x[i] - ubound[i]);

	}
	
	}
		


}

void individual::evaluation(CEC2013* pfunc)
{
	obj = pfunc->evaluate(x);
	fes++;
}


void individual::repair()
{

	for(int i=0;i<nvar;i++)
	{
	
		while (x[i] < lbound[i])
		{
			double rnd = rnd_uni(&rnd_uni_init);
				//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);

				x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - x[i]);

		}
		while (x[i] > ubound[i])
		{
			double rnd = rnd_uni(&rnd_uni_init);
			//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);

			x[i] = ubound[i] - rnd * (x[i] - ubound[i]);

		}
	
	}
	
}
void individual::repair_vec(vector<double> lb,vector<double> ub)
{

	for (int i = 0; i < nvar; i++)
	{

		while (x[i] < lb[i])
		{
			double rnd = rnd_uni(&rnd_uni_init);
			//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);

			x[i] = lb[i] + rnd_uni(&rnd_uni_init) * (lb[i] - x[i]);

		}
		while (x[i] > ub[i])
		{
			double rnd = rnd_uni(&rnd_uni_init);
			//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);

			x[i] = ub[i] - rnd * (x[i] - ub[i]);

		}

	}

}
double individual::gaussrand()
{
	double t1, t2, a, r;
	double x;
	
	t1 = rnd_uni(&rnd_uni_init);
	t2 = rnd_uni(&rnd_uni_init);

	a = 2 * 3.14159265358979323846*t1;            //a是极坐标的角度：变成了0~2*pi的均匀分布
	r = sqrt(-2 * log(t2));   //r是极坐标的距离：变成自然对数开根号的一种分布
	/*用极坐标(a,r)转换成笛卡尔坐标(x,y)，这就是产生的高斯白噪声*/
	x = r * cos(a);
	return x;
}

void individual::mut_gus(vector<double> &varl)
{
	for (int i = 0; i < nvar; i++)
	{

		double rnd = gaussrand();
		x[i] = x[i] + rnd * varl[i];
		 rnd = rnd_uni(&rnd_uni_init);
		if (x[i] < lbound[i])
			x[i] = lbound[i] + rnd *(lbound[i] - x[i]);
		if (x[i] > ubound[i])
			x[i] = ubound[i] - rnd * (x[i] - ubound[i]);
	}

}

void individual::compsite_de(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& r5, individual& child1,individual &child2,individual &child3)
{
	vector<double> cra;
	cra.push_back(0.1);
	cra.push_back(0.5);
	cra.push_back(0.9);
	vector<double> Fa;
	Fa.push_back(0.1);
	Fa.push_back(0.5);
	Fa.push_back(0.9); 
	int sel;
	sel =(int)( rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) { CR = 0.1; F = 1.0; }
	if (sel == 1) { CR = 0.9; F = 1.0; }
	if (sel == 2) { CR = 0.2; F = 0.8; }
	sel = (int)(rnd_uni(&rnd_uni_init) * 4);	
	diff( pp,  r1, r2,  r3, child1);	
	diff_r2( pp, r1,  r2,  r3, r4,r5, child2);	
	diff_r3( pp,  r1, r2, r3, child3);

}
void individual::compsite_de1(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& r5, individual& child)
{
	int sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) { CR = 0.1; F = 1.0; }
	if (sel == 1) { CR = 0.9; F = 1.0; }
	if (sel == 2) { CR = 0.2; F = 0.8; }
	 sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel ==0) { 
		diff(pp, r1, r2, r3, child);
	}
	if (sel == 1) {
		diff_r2(pp, r1, r2, r3, r4, r5, child);
	}
	if (sel == 2) {
		diff_r3(pp, r1, r2, r3, child);
	}

}

void individual::compsite_de_gbest(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& r5, individual& child1, individual& child2,individual &child3)
{
	int sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) { CR = 0.1; }
	if (sel == 1) { CR = 0.5;  }
	if (sel == 2) { CR = 0.7; }
	sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) {  F = 0.01; }
	if (sel == 1) {  F = 0.1; }
	if (sel == 2) {  F = 1.0; }

}

void individual::diff_best(individual& pp, individual& r1, individual& r2,individual& child)
{
	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);



	{
		for (int i = 0; i < nvar; i++)
		{
			if (rnd_uni(&rnd_uni_init) <= CR || i == j)
			{
				child.x[i] = pp.x[i] + F * (r1.x[i] - r2.x[i]);


				if (child.x[i] < lbound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
					if (rnd < 1)
						child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
					else
						child.x[i] = lbound[i];
				}
				if (child.x[i] > ubound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
					if (rnd < 1)
						child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
					else
						child.x[i] = ubound[i];
				}

			}
			else
				child.x[i] = pp.x[i];
		}


	}

}


void individual::diff_best_r2(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& child)
{

	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);
	for (int i = 0; i < nvar; i++)
	{
		if (rnd_uni(&rnd_uni_init) <= CR || i == j)
		{
			child.x[i] = pp.x[i] + F * (r1.x[i] - r2.x[i])+ F * (r3.x[i] - r4.x[i]);


			if (child.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
				if (rnd < 1)
					child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
				else
					child.x[i] = lbound[i];
			}
			if (child.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
				if (rnd < 1)
					child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
				else
					child.x[i] = ubound[i];
			}

		}
		else
			child.x[i] = pp.x[i];
	}
}
void individual::diff_best_r3(individual& pp, individual& r1, individual& r2, individual& child)
{

	double d = (rnd_uni(&rnd_uni_init) );
	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);
	for (int i = 0; i < nvar; i++)
	{
		if (rnd_uni(&rnd_uni_init) <= CR || i == j)
		{
			child.x[i] = pp.x[i] + F *d* (r1.x[i] - pp.x[i])+ F * d * (r2.x[i] - pp.x[i]);


			if (child.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
				if (rnd < 1)
					child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
				else
					child.x[i] = lbound[i];
			}
			if (child.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
				if (rnd < 1)
					child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
				else
					child.x[i] = ubound[i];
			}

		}
		else
			child.x[i] = pp.x[i];
	}
}
void individual::grade_diff(individual& pp, individual& r1, individual& r2, individual& child1, individual& child2, individual& child3)
{
	int sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) { CR = 0.5; }
	if (sel == 1) { CR = 0.7; }
	if (sel == 2) { CR = 1.0; }
	sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) { F = 0.0000001; }
	if (sel == 1) { F = 0.00001; }
	if (sel == 2) { F = 0.01; }

	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);
	double d = (rnd_uni(&rnd_uni_init));
	for (int i = 0; i < nvar; i++)
	{
		if (rnd_uni(&rnd_uni_init) <= CR || i == j)
		{

			child1.x[i] = pp.x[i] + F * (0.5 - d) * r1.grade[i];


			if (child1.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
				if (rnd < 0.5)
					child1.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child1.x[i]);
				else
					child1.x[i] = lbound[i];
			}
			if (child1.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
				if (rnd < 0.5)
					child1.x[i] = ubound[i] - rnd * (child1.x[i] - ubound[i]);
				else
					child1.x[i] = ubound[i];
			}

		}
		else
			child1.x[i] = pp.x[i];
	}


	sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) { CR = 0.5; }
	if (sel == 1) { CR = 0.7; }
	if (sel == 2) { CR = 1.0; }
	sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) { F = 0.0000001; }
	if (sel == 1) { F = 0.00001; }
	if (sel == 2) { F = 0.01; }

	j = (int)(rnd_uni(&rnd_uni_init) * nvar);

	for (int i = 0; i < nvar; i++)
	{
		if (rnd_uni(&rnd_uni_init) <= CR || i == j)
		{

			child2.x[i] = pp.x[i] + F * (d - 0.5) * (r1.grade[i] - r2.grade[i]);


			if (child2.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
				if (rnd < 1)
					child2.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child2.x[i]);
				else
					child2.x[i] = lbound[i];
			}
			if (child2.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
				if (rnd < 1)
					child2.x[i] = ubound[i] - rnd * (child2.x[i] - ubound[i]);
				else
					child2.x[i] = ubound[i];
			}

		}
		else
			child2.x[i] = pp.x[i];
	}


	sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) { CR = 0.5; }
	if (sel == 1) { CR = 0.7; }
	if (sel == 2) { CR = 1.0; }
	sel = (int)(rnd_uni(&rnd_uni_init) * 3);
	if (sel == 0) { F = 0.0000001; }
	if (sel == 1) { F = 0.00001; }
	if (sel == 2) { F = 0.01; }

	for (int i = 0; i < nvar; i++)
	{
		if (rnd_uni(&rnd_uni_init) <= CR || i == j)
		{

			child3.x[i] = pp.x[i] + 0.8 * (d - 0.5) * r1.grade[i] + 0.8 * (d - 0.5) * r2.grade[i];


			if (child3.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
				if (rnd < 1)
					child3.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child3.x[i]);
				else
					child3.x[i] = lbound[i];
			}
			if (child3.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
				if (rnd < 1)
					child3.x[i] = ubound[i] - rnd * (child3.x[i] - ubound[i]);
				else
					child3.x[i] = ubound[i];
			}

		}
		else
			child2.x[i] = pp.x[i];
	}
}
void individual::de(individual& pp, individual& r1, individual& r2, individual& r3, individual& child)
{

	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);



	if (nvar < 2)
	{
		if (isnan(child.x[0]) == 1)
			child.x[0] = rnd_uni(&rnd_uni_init) * (ubound[0] - lbound[0]);
		for (int i = 0; i < nvar; i++)
		{
			child.x[i] = r1.x[i] + 0.8 * (r2.x[i] - r3.x[i]);
			//	fl << " 0	" << r1.x[i] <<"	"<< r2.x[i] <<"	"<<r3.x[i] << child.x[i];
			if (isnan(child.x[0]) == 1)
				child.x[0] = rnd_uni(&rnd_uni_init) * (ubound[0] - lbound[0]);
			if (child.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
				//	fl << "	1	" << child.x[i];
			}
			if (child.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
				//		fl << "2	" << child.x[i];
			}
		}
	}
	else
	{
		for (int i = 0; i < nvar; i++)
		{
			if (rnd_uni(&rnd_uni_init) <= 0.2 || i == j)
			{
				child.x[i] = r1.x[i] + 0.8 * (r2.x[i] - r3.x[i]);


				if (child.x[i] < lbound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//  child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init)*(parentpop[r1].x[i] - lbound[i]);
					if (rnd < 1)
						child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
					else
						child.x[i] = lbound[i];
				}
				if (child.x[i] > ubound[i])
				{
					double rnd = rnd_uni(&rnd_uni_init);
					//	 child.x[i]= ubound[i] - rnd*(ubound[i] - parentpop[r1].x[i]);
					if (rnd < 1)
						child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
					else
						child.x[i] = ubound[i];
				}

			}
			else
				child.x[i] = pp.x[i];
		}


	}


}

void individual::de2(individual& pp, individual& r1, individual& r2, individual& r3, individual& r4, individual& r5, individual& child)
{



	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);
	if (nvar < 2)
	{
		for (int i = 0; i < nvar; i++)
		{
			child.x[i] = x[i] + 0.5 * (r2.x[i] - r3.x[i]) + 0.5 * (r5.x[i] - r4.x[i]);
			if (child.x[i] < lbound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = lbound[i] + rnd_uni(&rnd_uni_init) * (lbound[i] - child.x[i]);
			}
			if (child.x[i] > ubound[i])
			{
				double rnd = rnd_uni(&rnd_uni_init);
				child.x[i] = ubound[i] - rnd * (ubound[i] - r1.x[i]);
			}
		}
	}
	else
	{
		for (int i = 0; i < nvar; i++)
		{
			if (rnd_uni(&rnd_uni_init) <= 0.5 || i == j)
			{
				child.x[i] = r1.x[i] + 0.5 * (r2.x[i] - r3.x[i]) + 0.5 * (r4.x[i] - r5.x[i]);


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
				child.x[i] = pp.x[i];
		}


	}
}

void individual::de3(individual& pp, individual& r1, individual& r2, individual& r3, individual& child)
{
	double dstep = rnd_uni(&rnd_uni_init);
	int j = (int)(rnd_uni(&rnd_uni_init) * nvar);
	for (int i = 0; i < nvar; i++)
	{

		{
			child.x[i] = pp.x[i] + dstep * (r1.x[i] - pp.x[i]) + 0.8 * (r3.x[i] - r2.x[i]);


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

				if (rnd < 0.5)
					child.x[i] = ubound[i] - rnd * (child.x[i] - ubound[i]);
				else
					child.x[i] = ubound[i];
			}

		}
	}
}

#endif