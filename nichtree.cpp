// nichtree.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<fstream>
#include<vector>
// #include "node.h"

#include"global.h"
#include"CEC2013.h"
#include"cfunction.h"
//#include"clde.h"
#include"ode.h"
#include<iostream>

using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{

	ofstream  f1("me.txt");
	ofstream  rf("result.txt");
	std::vector <int> mdtest1;
	std::vector <int> mdtest2;
	std::vector <int> mdtest3;
	std::vector <int> mdtest4;
	std::vector <int> mdtest5;
	std::vector<int> numop;
	numop.push_back(2);
	numop.push_back(5);
	numop.push_back(1);
	numop.push_back(4);
	numop.push_back(2);
	numop.push_back(18);
	numop.push_back(36);
	numop.push_back(81);
	numop.push_back(216);
	numop.push_back(12);
	numop.push_back(6);
	numop.push_back(8);
	numop.push_back(6);
	numop.push_back(6);
	numop.push_back(8);
	numop.push_back(6);
	numop.push_back(8);
	numop.push_back(6);
	numop.push_back(8);
	numop.push_back(8);
	std::vector <double> md1;
	std::vector <double> md2;
	std::vector <double> md3;
	std::vector <double> md4;
	std::vector <double> md5;
	std::vector <double> m1;
	std::vector <double> m2;
	std::vector <double> m3;
	std::vector <double> m4;
	std::vector <double> m5;
	int numiter =10;


	 for(int n=1;n<21;n++)
	{
		 ID = n;
		
				if(n==1)
		{
			pops=80;
			maxfes=5.0E+04;
		}
					if(n==2)
		{
			pops=80;
			maxfes=5.0E+04;
		}
						if(n==3)
		{
			pops=80;
			maxfes=5.0E+04;
		}
							if(n==4)
		{
			pops=80;
			maxfes=5.0E+04;
		}
								if(n==5)
		{
			pops=80;
			maxfes=5.0E+04;
		}
									if(n==6)
		{
			pops=100;
			maxfes=2.0E+05;
		}
										if(n==7)
		{
			pops=300;
		maxfes=2.0E+05;
		}
											if(n==8)
		{
			pops=300;
			maxfes=4.0E+05;
		}
												if(n==9)
		{
				pops=300;
			maxfes=4.0E+05;
		}
													if(n==10)
		{
			pops=100;
			maxfes=2.0E+05;
		}
																if(n==11)
		{
			pops=200;
			maxfes=2.0E+05;
		}
																			if(n==12)
		{
				pops=200;
			maxfes=2.0E+05;
		}
																						if(n==13)
		{
			pops=200;
			maxfes=2.0E+05;
		}
				if(n==14)
		{
				pops=200;
			maxfes=4.0E+05;
		}
																												if(n==15)
		{
			pops=200;
			maxfes=4.0E+05;
		}
															if(n==16)
		{
				pops=200;
			maxfes=4.0E+05;
		}
															if(n==17)
		{
				pops=200;
			maxfes=4.0E+05;
		}
																	if(n==18)
		{
			pops=200;
			maxfes=4.0E+05;
		}
		if(n==19)
		{
			pops=200;
			maxfes=4.0E+05;
		}
			if(n==20)
		{
			pops=200;
			maxfes=4.0E+05;
		}
																								
     for(int k=0;k<numiter;k++)
	 {
		 cout << k <<"the last is:	"<< endl;
	CEC2013 *pFunc;
	pFunc= new CEC2013(n);
	nvar=pFunc->get_dimension();
	cout<<nvar<<endl;
	cout<<"the problem is	"<<n<<endl;

	lbound=new double[nvar];
	ubound= new double[nvar];
	for(int i=0;i<pFunc->get_dimension();i++)
	{
		lbound[i]=pFunc->get_lbound(i);
		ubound[i]=pFunc->get_ubound(i);
		cout<<lbound[i]<<"	"<<ubound[i]<<endl;
	}
	worstof=-100000000000000000;
	bestof=10000000000000000000;

	evol a;
	a.run(pFunc);

	std::vector< std::vector< double> > pop;
	for(int i=0;i<a.parentpop.size();i++)
		pop.push_back(a.parentpop[i].x);
/*	cout << "-------------------------------------------------------" <<endl;
	for (std::vector< std::vector<double> >::iterator it = pop.begin(); 
			it != pop.end(); ++it) {
		cout << "Fitness: " << pFunc->evaluate(*it) << "\tGene:\t";
		for (std::vector<double>::iterator jt = it->begin();
				jt != it->end(); ++jt) {
			cout << *jt << "\t";
		}
		cout << endl;
	}
	cout << "-------------------------------------------------------" <<endl;
	*/
	/* Calculate how many global optima are in the population */
	double accuracy=0.1;
	std::vector< std::vector<double> > seeds; 
	mdtest1.push_back(how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) );
	cout << "In the current population there exist " 
		<< how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) 
		<< " global optimizers." << endl;

	 accuracy=0.01;
	//std::vector< std::vector<double> > seeds; 
	mdtest2.push_back(how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) );
	cout << "In the current population there exist " 
		<< how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) 
		<< " global optimizers." << endl;

		accuracy=0.001;
	//std::vector< std::vector<double> > seeds; 
	mdtest3.push_back(how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) );
	cout << "In the current population there exist " 
		<< how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) 
		<< " global optimizers." << endl;

	 accuracy=0.0001;
	//std::vector< std::vector<double> > seeds; 
	mdtest4.push_back(how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) );
	cout << "In the current population there exist " 
		<< how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) 
		<< " global optimizers." << endl;

		 accuracy=0.00001;
//	std::vector< std::vector<double> > seeds; 
	mdtest5.push_back(how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) );
	cout << "In the current population there exist " 
		<< how_many_goptima(pop, seeds, pFunc, accuracy, pFunc->get_rho()) 
		<< " global optimizers." << endl;
	
	/* Clean up */
	delete pFunc;
	delete []lbound;
	delete []ubound;
		
	}
///////////////////////////////////////arr1
	 f1<<"the function	"<<n<<endl;
	  int d1=0;
	  int tnum=0;
	   for(int testd=0;testd<mdtest1.size();testd++)
	   {
	  f1<<mdtest1[testd]<<endl;
	  d1+=mdtest1[testd];
	  if(mdtest1[testd]==numop[n-1])
	  {
		  tnum++;
	  }
	   }
	   md1.push_back(tnum/50.0);
	   m1.push_back((d1*1.0)/(numop[n-1]*numiter));
	   f1<<d1;
	   
	     f1<<"-------------------------------------------------------------"<<endl;

		 ///////////////////////////////////////arr2
		 int d2=0;
		 tnum=0;
	      for(int testd=0;testd<mdtest2.size();testd++)
		  {
	  f1<<mdtest2[testd]<<endl;
	  d2+=mdtest2[testd];
	    if(mdtest2[testd]==numop[n-1])
	  {
		  tnum++;
	  }
		  } 
		  md2.push_back(tnum/50.0);
	   m2.push_back((d2*1.0)/(numop[n-1]*numiter));
		  f1<<d2;
		      f1<<"-------------------------------------------------------------"<<endl;



			   ///////////////////////////////////////arr3
			  int d3=0;
			  tnum=0;
		     for(int testd=0;testd<mdtest3.size();testd++)
			 {
	  f1<<mdtest3[testd]<<endl;
			 d3+=mdtest3[testd];
			  if(mdtest3[testd]==numop[n-1])
	  {
		  tnum++;
	  }
			 }
			   md3.push_back(tnum/50.0);
	   m3.push_back((d3*1.0)/(numop[n-1]*numiter));
			  f1<<d3;
			     f1<<"-------------------------------------------------------------"<<endl;




				  ///////////////////////////////////////arr4
				 int d4=0;
				 tnum=0;
			    for(int testd=0;testd<mdtest4.size();testd++)
				{
	  f1<<mdtest4[testd]<<endl;
	  d4+=mdtest4[testd];
	   if(mdtest4[testd]==numop[n-1])
	  {
		  tnum++;
	  }
				}
				  md4.push_back(tnum/50.0);
	   m4.push_back((d4*1.0)/(numop[n-1]*numiter));
				f1<<d4;
				    f1<<"-------------------------------------------------------------"<<endl;



					 ///////////////////////////////////////arr5
					int d5=0;
					tnum=0;
				   for(int testd=0;testd<mdtest5.size();testd++)
				   {
	  f1<<mdtest5[testd]<<endl;
	  d5+=mdtest5[testd];
	   if(mdtest5[testd]==numop[n-1])
	  {
		  tnum++;
	  }
				   }
				     md5.push_back(tnum/50.0);
	   m5.push_back((d5*1.0)/(numop[n-1]*numiter));
				   	f1<<d5;
				       f1<<"-------------------------------------------------------------"<<endl;
  mdtest1.clear();
  mdtest2.clear();
  mdtest3.clear();
  mdtest4.clear();
  mdtest5.clear();

 

  }
 
  rf << "_____________________________________the accuracy=0.1" << endl;
  for (int lresult = 0; lresult < md1.size(); lresult++)
  {
	  rf << m1[lresult] << "	" << md1[lresult] << endl;
  }


  rf << "_____________________________________the accuracy=0.01" << endl;
  for (int lresult = 0; lresult < md1.size(); lresult++)
  {
	  rf << m2[lresult] << "	" << md2[lresult] << endl;
  }
  rf << "_____________________________________the accuracy=0.001" << endl;
  for (int lresult = 0; lresult < md1.size(); lresult++)
  {
	  rf << m3[lresult] << "	" << md3[lresult] << endl;
  }
  rf << "_____________________________________the accuracy=0.0001" << endl;
  for (int lresult = 0; lresult < md1.size(); lresult++)
  {
	  rf << m4[lresult] << "	" << md4[lresult] << endl;
  }
  rf << "_____________________________________the accuracy=0.00001" << endl;
  for (int lresult = 0; lresult < md1.size(); lresult++)
  {
	  rf << m5[lresult] << "	" << md5[lresult] << endl;
  }

	
	/*  pops=100;
    ID =1;
	if(ID==1)
	{
		nvar=2;
	}
	if(ID==2)
	{
		nvar=20;
	}
	if(ID==3)
	{
		nvar=2;
	}
	if(ID==4)
	{
		nvar=2;
	}
	if(ID==5)
	{
		nvar=3;
	}
	if(ID==6)
	{
		nvar=6;
	}
	if(ID==7)
	{
		nvar=20;
	}
	lbound=new double[nvar];
	ubound=new double[nvar];
	for(int i=0;i<nvar;i++)
	{
		lbound[i]=-1;
		ubound[i]=1;
	}
  Earun ea;
  ea.run();*/

    return 0;
	

	return 0;
}

