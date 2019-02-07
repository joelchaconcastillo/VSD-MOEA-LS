/*==========================================================================
// //  Author: Carlos Segura, Joel Chacón 
//     Description: 
//
// ===========================================================================*/


#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <queue>
#include <iomanip>
#include <cfloat>
#include "global.h"
#include "recomb.h"
#include "common.h"
#include "individual.h"

class MOEA
{

public:
	MOEA();
	virtual ~MOEA();

	void init_population();                  // initialize the population

	void evol_population();                                      
	void exec_emo(int run);

	void save_front(char savefilename[1024]);       // save the pareto front into files
	void save_pos(char savefilename[1024]);

        void penalize_nearest(vector<CIndividual *> &candidates, vector<CIndividual *> &penalized);
        void select_farthest_penalized(vector<CIndividual *> &candidates, vector<CIndividual *> &penalized);
        void select_best_candidate(vector<CIndividual *> &survivors, vector<CIndividual *> &candidates, vector<CIndividual *> &penalized);
        void compute_distances(vector<CIndividual *> &candidates, vector<CIndividual *> &survivors);
        void binary_tournament_selection(vector<CIndividual > &population, vector<CIndividual> &child_pop);
        void recombination(vector<CIndividual> &child_pop);
        void reproduction(vector<CIndividual> &population, vector<CIndividual> &child_pop);
        void update_diversity_factor();
        void computing_dominate_information(vector <CIndividual*> &pool);
        void select_first_survivors(vector<CIndividual*> &survivors, vector<CIndividual*> &candidates);
        void update_domianted_information(vector<CIndividual*> &survivors, vector<CIndividual*> &current);
        void update_population(vector<CIndividual*> &survivors, vector<CIndividual> &population);
	void update_archive(vector<CIndividual> &archive, vector<CIndividual> &population);

	void fast_non_dominated_sorting(vector <CIndividual*> &survivors);
	void update_archive(vector<CIndividual> &archive, vector<CIndividual> &population, vector<CIndividual> &child_pop);
	void LocalSearch1(CIndividual &old_individual);
	void LocalSearch2(CIndividual &old_individual);
	void LocalSearch3(CIndividual &old_individual);
	void improvement(vector<CIndividual> old_individuals);
	double distance( vector<double> &a, vector<double> &b);

	double distance_improvement( vector<double> &a, vector<double> &b);
	vector <CIndividual> population;
	vector<CIndividual> child_pop;	
	vector<CIndividual> archive;	
	void operator=(const MOEA &moea);

public:
//
//	// algorithm parameters
	int nfes;
//	int     nfes, max_nfes;          //  the number of function evluations

};

MOEA::MOEA()
{

}

MOEA::~MOEA()
{

}
double MOEA::distance( vector<double> &a, vector<double> &b)
{
	double dist = 0 ;
   for(int i = 0; i < a.size(); i++)
	{
	   double factor = (a[i]-b[i])/(vuppBound[i] - vlowBound[i]);
	   dist += factor*factor;
	}
   return sqrt(dist);
}
double MOEA::distance_improvement( vector<double> &a, vector<double> &b)
{
	double dist = 0 ;
	double maxd = -INFINITY;
   for(int i = 0; i < a.size(); i++)
	{
	   double factor = max(0.0,a[i]-b[i]);
	   dist += factor*factor;
	   maxd = max(maxd, max(b[i]-a[i],0.0));
	}
        if(dist == 0.0) return -maxd; //in case that this indicator is zero, this mean that it is a dominated individual...
   return dist;
   return sqrt(dist);
}
void MOEA::init_population()
{

    for(int i=0; i<pops; i++)
	{
		CIndividual indiv1, indiv2, indiv3;
		// Randomize and evaluate solution
		indiv1.rnd_init();
		indiv1.obj_eval();
		// Save in the population
		population.push_back(indiv1);

		indiv2.rnd_init();
		indiv2.obj_eval();

		child_pop.push_back(indiv2);

		indiv3.rnd_init();
		indiv3.obj_eval();

		archive.push_back(indiv3);
		nfes+=3;
	}
}
void MOEA::operator=(const MOEA &alg)
{
	//population = alg.population;
}
void MOEA::evol_population()
{
	vector<CIndividual *> penalized, survivors;
//	cout << "generation "<<endl;
	//join the offspring and parent populations
	vector<CIndividual *> candidates;
	for(int i = 0; i < pops; i++)
	{
	  candidates.push_back( &(population[i]));
	  candidates.push_back( &(child_pop[i]));
	  candidates.push_back( &(archive[i]));
	}
	computing_dominate_information(candidates); //computing the dominate count of each candidate individual...
	//select the "best" individuals that owns to candidate set and are moved in survivors set...
	select_first_survivors(survivors, candidates);
	//update the diversity-factor-parameter...	
	update_diversity_factor();
	//Pre-computing the neares distances both objective and decision spaces..
	compute_distances(candidates, survivors);
       	while( survivors.size() < pops )
	{
	  penalize_nearest(candidates, penalized);//penalize the nearest individuals.. 
	  if(candidates.empty())	  
	     select_farthest_penalized(survivors, penalized);//in case that all the individuals are penalized pick up the farstest and add it to survirvors
	  else
	    {
	     update_domianted_information(survivors, candidates); //update the rank of each candidate whitout penalized
	     select_best_candidate(survivors, candidates, penalized); // the best candidate is selected considering the improvemente distance, and the rank..
	    }
	}
	fast_non_dominated_sorting(survivors);//rank the survivors individuals..

	//this procedure is necessary since the penalized individuals
	update_population(survivors, population); //update the parent population 
	
	reproduction(population, child_pop); //generate a new population considering the survivors individuals...
	update_archive(archive, population, child_pop);
	
        //improvement of the child....
	//improvement(child_pop);
	//improvement(archive);
}

void MOEA::improvement(vector<CIndividual> child_improved) //note: is a copy of child_pop
{
//   for(int i = 0 ; i < child_improved.size(); i++)
      {
	int in =rand()%child_improved.size();
	//LocalSearch1(child_improved[i]);	
//	LocalSearch1(child_improved[in]);	
//	LocalSearch2(child_improved[in]);	
	//LocalSearch2(child_improved[rand()%child_improved.size()]);	
	LocalSearch3(child_improved[in]);	
//	LocalSearch3(child_improved[rand()%child_improved.size()]);	
      }	
	update_archive(archive, child_improved);
}

void MOEA::LocalSearch1(CIndividual &newindividual)
{
	CIndividual old_individual;
	vector<double> SearchRange(nvar);
	vector<int> Direction(nvar), order(nvar);
	for(int i = 0; i < nvar; i++)
	{
	   SearchRange[i] = (vuppBound[i]-vlowBound[i])/2.0;
	   order[i] = i;
	}
	  bool improve = false;
	  bool LocalOptima = false;
	int k=0;
	while(!LocalOptima && k < 100)
	{
	 k++;
	  if(!improve)
	  {
	   bool expand = true;
	   for(int i = 0; i < nvar; i++)  //check all variables..
	   {
		Direction[i] = (rnd_uni(&rnd_uni_init) < 0.5)?-1:1; //similar than a bernulli distrbutions, which is dedicated for the direcction..
		SearchRange[i] *=0.5;
	      if(SearchRange[i] > 1e-8) //expand again if all search range variables are lower than 1e-8, note that some variables could be zero and could be keep in this way..
	      {
		  expand = false;
	      }
	   }
	   if( expand )
	   {
	      LocalOptima = true;
	      for(int i = 0; i < nvar; i++)  
	      {
	            SearchRange[i] = (vuppBound[i]-vlowBound[i])*0.4;
	      }
	   }
	  }
	  improve = false;
	   next_permutation(order.begin(), order.end());
	   //random_shuffle(order.begin(), order.end());
	   for(int d = 0; d < nvar; d++)
	   {
		int i = order[d];
		old_individual = newindividual;
		newindividual.x_var[i] = old_individual.x_var[i] + SearchRange[i]*Direction[i];
		newindividual.x_var[i] = max(vlowBound[i], newindividual.x_var[i]);
		newindividual.x_var[i] = min(vuppBound[i], newindividual.x_var[i]);
	        newindividual.obj_eval();
	        nfes++;
		if(newindividual < old_individual)
		{
		  improve=true;
		  LocalOptima = false;
		}
	        else if( old_individual < newindividual) //weakly dominance
	        {
		   newindividual = old_individual;
	           newindividual.x_var[i] = old_individual.x_var[i] - 0.5*SearchRange[i]*Direction[i];
	           newindividual.x_var[i] = max(vlowBound[i], newindividual.x_var[i]);
	           newindividual.x_var[i] = min(vuppBound[i], newindividual.x_var[i]);
	           newindividual.obj_eval();
	           nfes++;
		   if(newindividual < old_individual)
		   {
		     improve=true;
		     LocalOptima = false;
		   }
	           else if( old_individual < newindividual)
	           {
	             newindividual = old_individual;
	           }
	         }
	     }
	}
}
void MOEA::LocalSearch2(CIndividual &newindividual)
{
	CIndividual old_individual;
	vector<double> SearchRange(nvar);
	vector<int> Direction(nvar), order(nvar);
	for(int i = 0; i < nvar; i++)
	{
	   SearchRange[i] = (vuppBound[i]-vlowBound[i])/2.0;
	   order[i] = i;
	}
	  bool improve = false;
	  bool LocalOptima = false;
	int k=0;
	while(!LocalOptima && k < 100)
	{
	 k++;
	  if(!improve)
	  {
	   bool expand = true;
	   for(int i = 0; i < nvar; i++)  //check all variables..
	   {
		Direction[i] = (rnd_uni(&rnd_uni_init) < 0.5)?-1:1; //similar than a bernulli distrbutions, which is dedicated for the direcction..
		SearchRange[i] *=0.5;
	      if(SearchRange[i] > 1e-8) //expand again if all search range variables are lower than 1e-8, note that some variables could be zero and could be keep in this way..
	      {
		  expand = false;
	      }
	   }
	   if( expand )
	   {
	      LocalOptima = true;
	      for(int i = 0; i < nvar; i++)  
	      {
	            SearchRange[i] = (vuppBound[i]-vlowBound[i])*0.4;
	      }
	   }
	  }
	  improve = false;
	   next_permutation(order.begin(), order.end());
	   //random_shuffle(order.begin(), order.end());
	   for(int d = 0; d < nvar/4; d++)
	   {
		int i = order[d];
		old_individual = newindividual;
		newindividual.x_var[i] = old_individual.x_var[i] + SearchRange[i]*Direction[i];
		newindividual.x_var[i] = max(vlowBound[i], newindividual.x_var[i]);
		newindividual.x_var[i] = min(vuppBound[i], newindividual.x_var[i]);
	  }
	  newindividual.obj_eval();
	  nfes++;
	  if(newindividual < old_individual)
	  {
	    improve=true;
	    LocalOptima = false;
	  }
	  else if( old_individual < newindividual) //weakly dominance
	  {
	   newindividual = old_individual;
	   for(int d = 0; d < nvar/4; d++)
	    {
	      int i = order[d];
	      newindividual.x_var[i] = old_individual.x_var[i] - 0.5*SearchRange[i]*Direction[i];
	      newindividual.x_var[i] = max(vlowBound[i], newindividual.x_var[i]);
	      newindividual.x_var[i] = min(vuppBound[i], newindividual.x_var[i]);
	    }
	      newindividual.obj_eval();
	      nfes++;
	      if(newindividual < old_individual)
	      {
	        improve=true;
	        LocalOptima = false;
	      }
	      else if( old_individual < newindividual)
	      {
	        newindividual = old_individual;
	      }
	   
	  }
	}
}
void MOEA::LocalSearch3(CIndividual &old_individual)
{

  CIndividual Best = old_individual;
   bool LocalOptima = false;
	int k = 0;
   while( !LocalOptima && k <10)
   {
	k++;
     vector<double> SearchU(nvar), SearchL(nvar), Disp(nvar);   
     LocalOptima = true;
     for(int i = 0; i < nvar; i++)
      {
	SearchU[i] = vuppBound[i];
	SearchL[i] = vlowBound[i];
	Disp[i] = (SearchU[i] - SearchL[i])/100.0;
      }
	vector<int> order(nvar);
	for(int i = 0; i < order.size(); i++) order[i] = i;
        next_permutation(order.begin(), order.end());
      double maxDisp = DBL_MAX;
      while(maxDisp > 1e-3)
      {
	maxDisp = -DBL_MAX;
	
	   //random_shuffle(order.begin(), order.end());for(int i = 0; i < nvar; i++) order[i] = i;
	for(int l =0; l < nvar; l++)
	{
	   int i = order[l];
	  //if(Disp[i] < 1e-3) continue;
	//  if((SearchU[i]-SearchL[i]) < DBL_EPSILON) continue;
	   for(double d = SearchL[i]+1e-5; d < SearchU[i]; d+=Disp[i])
	   { 
	      CIndividual neighbour = Best;
	      neighbour.x_var[i] = d;
	      neighbour.x_var[i] = min(neighbour.x_var[i], vuppBound[i]);
	      neighbour.x_var[i] = max(neighbour.x_var[i], vlowBound[i]);
	      neighbour.obj_eval();
	      nfes++;
	      if( neighbour < Best)
	      //if( distance_improvement(neighbour.y_obj, Best.y_obj) - distance_improvement(Best.y_obj,neighbour.y_obj ) > 1e-3 )
	      {
		bool a = (neighbour == Best), b =(neighbour< Best);
	//	neighbour.show_objective();
	//	Best.show_objective();
	//	cout <<a <<endl;
	//	cout <<b <<endl;
		  Best = neighbour;
		  LocalOptima = false;
	      }
	   }
	   SearchU[i] = min(Best.x_var[i]+2.0*Disp[i], vuppBound[i]);
	   SearchL[i] = max(Best.x_var[i]-2.0*Disp[i], vlowBound[i]);
	   Disp[i] = (SearchU[i] - SearchL[i])/10.0;
	   maxDisp = max(maxDisp, Disp[i]);
	}
      }    
      old_individual = Best;
    }

}
void MOEA::update_archive(vector<CIndividual> &archive, vector<CIndividual> &population, vector<CIndividual> &child_pop)
{
	vector<CIndividual> candidates;
	for(int i = 0; i < archive.size(); i++)
	{
	   candidates.push_back(archive[i]);
	   candidates.push_back(population[i]);
	   candidates.push_back(child_pop[i]);
	}
	archive.clear();
	vector<int> BestIndex;
	///Select the best improvement distance candidates....
	vector<bool> selected(candidates.size(), false);
	for(int m = 0; m < nobj; m++)
	{
		int indxmaxim;
		double bestvector = INFINITY;
		for(int i = 0; i <  candidates.size(); i++)
		 {	
			// if(candidates[i]->times_dominated != 0) continue; //just consider the first front
		 	  if(selected[i])continue;
		        double s = 0.0;	
		        double maxv = -INFINITY;
		        for(int k = 0; k < nobj; k++)
		        {
		      	   double fi = fabs(candidates[i].y_obj[k]);
		      	   s += fi;
		      	   double ti = (k==m)?fi:1e5*fi;
			    if(ti > maxv)   maxv=ti;
		        }
		         maxv = maxv + 0.0001*s;
		        if(bestvector > maxv)
		        { indxmaxim = i; bestvector = maxv;}
		  }
		selected[indxmaxim] = true;
		archive.push_back(candidates[indxmaxim]);	
	}
	vector< double > min_improvement_dist(candidates.size(), DBL_MAX);
	for(int i = 0; i < candidates.size(); i++)
	{
	   if(selected[i]) continue; 
	   for(int j = 0; j < archive.size(); j++)
	   {
		min_improvement_dist[i] = min(min_improvement_dist[i], distance_improvement(archive[j].y_obj, candidates[i].y_obj));
	   } 
	}
	while(archive.size() < population.size())
	{
	   double maxd = -DBL_MAX;
	   int maxi =-1;
	   for(int i = 0; i < min_improvement_dist.size(); i++)
	   {
		if( selected[i] )continue;
	        if( maxd < min_improvement_dist[i])
		{
		   maxd = min_improvement_dist[i];
		   maxi = i;
		}
	   }
	   selected[maxi] = true;
	   archive.push_back(candidates[maxi]);
	   for(int i = 0; i < min_improvement_dist.size(); i++)
	   {
	     if(selected[i])continue;
		min_improvement_dist[i] = min(min_improvement_dist[i], distance_improvement(candidates[maxi].y_obj, candidates[i].y_obj));
	   }
	}

}
void MOEA::update_archive(vector<CIndividual> &archive, vector<CIndividual> &population)
{
	vector<CIndividual> candidates;
	for(int i = 0; i < archive.size(); i++)
	{
	   candidates.push_back(archive[i]);
	   candidates.push_back(population[i]);
	}
	archive.clear();
	vector<int> BestIndex;
	///Select the best improvement distance candidates....
	vector<bool> selected(candidates.size(), false);
	for(int m = 0; m < nobj; m++)
	{
		int indxmaxim;
		double bestvector = INFINITY;
		for(int i = 0; i <  candidates.size(); i++)
		 {	
			// if(candidates[i]->times_dominated != 0) continue; //just consider the first front
		 	  if(selected[i])continue;
		        double s = 0.0;	
		        double maxv = -INFINITY;
		        for(int k = 0; k < nobj; k++)
		        {
		      	   double fi = fabs(candidates[i].y_obj[k]);
		      	   s += fi;
		      	   double ti = (k==m)?fi:1e5*fi;
			    if(ti > maxv)   maxv=ti;
		        }
		         maxv = maxv + 0.0001*s;
		        if(bestvector > maxv)
		        { indxmaxim = i; bestvector = maxv;}
		  }
		selected[indxmaxim] = true;
		archive.push_back(candidates[indxmaxim]);	
	}
	vector< double > min_improvement_dist(candidates.size(), DBL_MAX);
	for(int i = 0; i < candidates.size(); i++)
	{
	   if(selected[i]) continue; 
	   for(int j = 0; j < archive.size(); j++)
	   {
		min_improvement_dist[i] = min(min_improvement_dist[i], distance_improvement(archive[j].y_obj, candidates[i].y_obj));
	   } 
	}
	while(archive.size() < population.size())
	{
	   double maxd = -DBL_MAX;
	   int maxi =-1;
	   for(int i = 0; i < min_improvement_dist.size(); i++)
	   {
		if( selected[i] )continue;
	        if( maxd < min_improvement_dist[i])
		{
		   maxd = min_improvement_dist[i];
		   maxi = i;
		}
	   }
	   selected[maxi] = true;
	   archive.push_back(candidates[maxi]);
	   for(int i = 0; i < min_improvement_dist.size(); i++)
	   {
	     if(selected[i])continue;
		min_improvement_dist[i] = min(min_improvement_dist[i], distance_improvement(candidates[maxi].y_obj, candidates[i].y_obj));
	   }
	}

}
void MOEA::fast_non_dominated_sorting(vector <CIndividual*> &survivors)
{
   vector< vector < int > > dominate_list(survivors.size()); //in the worst case the number of fronts is the same as the survivors size
   vector< int > dominated_count (survivors.size(), 0), currentfront;
   for(int i = 0; i < survivors.size(); i++)
   {
	   for(int j = 0; j < survivors.size(); j++)
	  {
		if(i==j) continue;
	       if( *(survivors[i]) < *(survivors[j]))
	   	    dominate_list[i].push_back(j);
		else if (*(survivors[j]) < *(survivors[i]))
		   dominated_count[i]++;
 	  }
	if(dominated_count[i] == 0 ) currentfront.push_back(i);// get the first front
   }
   int rank = 0;
   while(!dominate_list[rank].empty())
   {
	vector<int> nextFront;
	for(int i = 0; i < currentfront.size(); i++)
	{
	   survivors[currentfront[i]]->rank = rank;
	   for(int j = 0; j < dominate_list[currentfront[i]].size(); j++)
	   {
		dominated_count[dominate_list[currentfront[i]][j]]--;
		if( dominated_count[dominate_list[currentfront[i]][j]] == 0) nextFront.push_back(dominate_list[currentfront[i]][j]);
		
	   }
	}	
	rank++;
	currentfront = nextFront;
   }
}
void MOEA::update_population(vector<CIndividual*> &survivors, vector<CIndividual> &population)
{
	vector<CIndividual> pool;
   for(int i = 0; i < survivors.size(); i++) pool.push_back(*(survivors[i]));
   for(int i = 0; i < population.size(); i++) population[i] = pool[i];
}
void MOEA::update_domianted_information(vector<CIndividual*> &survivors, vector<CIndividual*> &candidates)
{
     bool firstfrontcurrent = false; 
     for(int i = 0; i < candidates.size(); i++) if(candidates[i]->times_dominated==0) firstfrontcurrent = true; //check if there exists at least one candidate in the lowest current front
     
     if( !firstfrontcurrent) //this indicates that there is not able a current in the lowest front, so the next front is to be considered
	{	
	  
	   for(int i = 0; i < survivors.size(); i++)
	   {
		if(survivors[i]->times_dominated == 0)
		{
		      for(int j = 0; j < survivors[i]->ptr_dominate.size(); j++)
		  	   {
		  		survivors[i]->ptr_dominate[j]->times_dominated--;
		   	   }
		  	   survivors[i]->times_dominated--;
		}
	   }
	}
}
void MOEA::select_first_survivors(vector<CIndividual*> &survivors, vector<CIndividual*> &candidates)
{
	vector<int> BestIndex;
	///Select the best improvement distance candidates....
	vector<bool> grid(candidates.size(), false);
	for(int m = 0; m < nobj; m++)
	{
		int indxmaxim;
		double bestvector = INFINITY;
		for(int i = 0; i <  candidates.size(); i++)
		 {	
			 if(candidates[i]->times_dominated != 0) continue; //just consider the first front
		 	  if(grid[i])continue;
		        double s = 0.0;	
		        double maxv = -INFINITY;
		        for(int k = 0; k < nobj; k++)
		        {
		      	   double fi = fabs(candidates[i]->y_obj[k]);
		      	   s += fi;
		      	   double ti = (k==m)?fi:1e5*fi;
			    if(ti > maxv)   maxv=ti;
		        }
		         maxv = maxv + 0.0001*s;
		        if(bestvector > maxv)
		        { indxmaxim = i; bestvector = maxv;}
		  }
		grid[indxmaxim] = true;
		BestIndex.push_back(indxmaxim);
	}
	
		sort(BestIndex.begin(), BestIndex.end()); //sort the indexes and remove from candidates

		for(int i = BestIndex.size()-1; i >= 0 ; i--)
		{
			int index = BestIndex[i];
			survivors.push_back( candidates[index]);
			iter_swap(candidates.begin()+index, candidates.end()-1);
			candidates.pop_back();
		}
}
//get the rank of each individual...
void MOEA::computing_dominate_information(vector <CIndividual*> &pool)
{
    for(int i = 0; i < pool.size(); i++)
    {
	pool[i]->times_dominated = 0;
	pool[i]->ptr_dominate.clear();
	for(int j = 0; j < pool.size(); j++)
	{
	    if(i == j) continue;
	    if( *(pool[i]) < *(pool[j]) ) //the check if pop[i] dominates pop[j], tht '<' is overloaded
	    {
		pool[i]->ptr_dominate.push_back(pool[j]);
	    }
	    else if( *(pool[j]) < *(pool[i]) )
	   {
		pool[i]->times_dominated++;	
	   }
	}
    }
}
//updates the lowest distance factor of the diversity explicitly promoted
void MOEA::update_diversity_factor()
{
	double ratio = ((double) nfes)/max_nfes;
	lowestDistanceFactor = Initial_lowest_distance_factor - Initial_lowest_distance_factor*(ratio/0.9);
//	double frequency = 100;
//        lowestDistanceFactor = Initial_lowest_distance_factor*( 0.5- 0.5*cos(M_PI*cos( frequency*(ratio)*M_PI )*(1.0-ratio)));

}
void MOEA::reproduction(vector<CIndividual> &population, vector<CIndividual> &child_pop)
{
   //binary tournament selction procedure
   binary_tournament_selection(population, child_pop);
   //recombination of the individuals, through the SBX code (taken from the nsga-ii code), also the evaluation of the population is performed
   recombination(child_pop); 
}
void MOEA::recombination(vector<CIndividual> &child_pop)
{
   vector<CIndividual> child_pop2 = child_pop;
	
   for(int i = 0; i < child_pop.size(); i+=2)
    {
       int indexa = int(rnd_uni(&rnd_uni_init)*pops);
       int indexb = int(rnd_uni(&rnd_uni_init)*pops);	
       real_sbx_xoverA( child_pop2[indexa], child_pop2[indexb], child_pop[i], child_pop[i+1]);//the crossover probability and index distribution eta are configured in the global.h file
       realmutation(child_pop[i]); //the index distribution (eta) and  mutation probability are configured in the global.h file
       realmutation(child_pop[i+1]);
       child_pop[i].obj_eval();
       child_pop[i+1].obj_eval();
    }
}
void MOEA::binary_tournament_selection(vector<CIndividual > &population, vector<CIndividual> &child_pop)
{
   for(int i = 0; i < population.size(); i++)
	{
	   int indexa = int(rnd_uni(&rnd_uni_init)*pops);
	   int indexb = int(rnd_uni(&rnd_uni_init)*pops);
	   if(population[indexa].rank < population[indexb].rank)
	      child_pop[i] = population[indexa];
	   else if(population[indexa].rank > population[indexb].rank)
	      child_pop[i] = population[indexb];
	   else 
	   {
	      child_pop[i] = (rnd_uni(&rnd_uni_init) < 0.5  )? population[indexa] : population[indexb];
	   }	
	}
}
void MOEA::compute_distances(vector<CIndividual *> &candidates, vector<CIndividual *> &survivors)
{
	for(int i = 0; i < candidates.size(); i++)
	{
	    candidates[i]->nearest_variable_distance = INFINITY;
	    candidates[i]->neares_objective_distance = INFINITY;
	   for(int j = 0; j < survivors.size(); j++)
	   {
		candidates[i]->nearest_variable_distance = min( candidates[i]->nearest_variable_distance, distance(candidates[i]->x_var, survivors[j]->x_var));
		candidates[i]->neares_objective_distance = min( candidates[i]->neares_objective_distance, distance_improvement(survivors[j]->y_obj, candidates[i]->y_obj));
	   }
	}	

}
void MOEA::select_best_candidate(vector<CIndividual *> &survivors, vector<CIndividual *> &candidates, vector<CIndividual *> &penalized)
{
	int best_index_lastfront = -1;//the index of current with the farthes improvement distance
	double max_improvement = -INFINITY;
	  for(int i = 0 ; i < candidates.size(); i++)
	    {
		   if(candidates[i]->times_dominated != 0) continue;
			if(  max_improvement < candidates[i]->neares_objective_distance  )
			{
				max_improvement = candidates[i]->neares_objective_distance;
				best_index_lastfront= i;
			}
	    }
	 if(best_index_lastfront == -1) return; //this occurs when the first m-survirvors are dominated bewteen them, thus there are not candidates availables to pick, therefore this iteration is skiped, so in the next iteration will be available some candidates...

	//update distances of Current and penalized
	  for(int i = 0 ; i < candidates.size(); i++)
	   {
		if( i != best_index_lastfront) // Avoid to be updated by itself..
	        {
		 candidates[i]->nearest_variable_distance = min( candidates[i]->nearest_variable_distance, distance(candidates[i]->x_var, candidates[best_index_lastfront]->x_var ) );
		 candidates[i]->neares_objective_distance = min( candidates[i]->neares_objective_distance, distance_improvement(candidates[best_index_lastfront]->y_obj, candidates[i]->y_obj));
		}
	   }
	  for(int i = 0 ; i < penalized.size(); i++)
	  {
		penalized[i]->nearest_variable_distance = min( penalized[i]->nearest_variable_distance, distance( penalized[i]->x_var, candidates[best_index_lastfront]->x_var ) )  ;
		penalized[i]->neares_objective_distance  =  min( penalized[i]->neares_objective_distance, distance_improvement(candidates[best_index_lastfront]->y_obj, penalized[i]->y_obj));
	  }
	  survivors.push_back(candidates[best_index_lastfront]);
	  iter_swap(candidates.begin()+best_index_lastfront, candidates.end()-1);
	  candidates.pop_back();
}
void MOEA::select_farthest_penalized(vector<CIndividual *> &survivors, vector<CIndividual *> &penalized)
{
    	double largestDCN = -INFINITY;
	int index_largestDCN=0;
	for(int i = 0; i < (int)penalized.size(); i++) // get the index of penalized with larges DCN
	{
		if(penalized[i]->nearest_variable_distance >  largestDCN )
		{
			index_largestDCN = i;
			largestDCN = penalized[i]->nearest_variable_distance;
		}
	}

	for(int i = 0 ; i < (int)penalized.size(); i++) //update the nearest distance once that the penalized is moved to candidate (thereafter to survivors)
	{
		if( i != index_largestDCN )
		penalized[i]->nearest_variable_distance = min( penalized[i]->nearest_variable_distance, distance( penalized[i]->x_var, penalized[index_largestDCN]->x_var));
	}	

	survivors.push_back(penalized[index_largestDCN]);
	iter_swap(penalized.begin()+index_largestDCN, penalized.end()-1);
	penalized.pop_back();
}
void MOEA::penalize_nearest(vector<CIndividual *> &candidates, vector<CIndividual *> &penalized)
{
   	for(int i = candidates.size()-1; i >=0; i--)
	{	
		if( candidates[i]->nearest_variable_distance < lowestDistanceFactor )
		{
			penalized.push_back(candidates[i]);
			for(int j = 0; j < candidates[i]->ptr_dominate.size(); j++)
			{
				candidates[i]->ptr_dominate[j]->times_dominated--; //decreasing the times in which survivors is dominated, this since penalized individuals are not considered..
			}
			//remove the candidate with index "i"
			iter_swap(candidates.begin()+i, candidates.end()-1);
			candidates.pop_back();
		}
	}
}
void MOEA::exec_emo(int run)
{
       char filename1[5024];
       char filename2[5024];
		seed = run;
	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;

	nfes      = 0;
	init_population(); //Initialize individuals...

	sprintf(filename1,"%s/POS/POS_MOEAD_%s_RUN%d_seed_%d_nobj_%d.dat_bounded",currentPATH, strTestInstance,run, seed, nobj);
	sprintf(filename2,"%s/POF/POF_MOEAD_%s_RUN%d_seed_%d_nobj_%d.dat_bounded",currentPATH, strTestInstance,run, seed, nobj);
	while(nfes < max_nfes )
	{
		evol_population();
		nfes += pops;
//	    if( !(nfes % (max_nfes/10)  ))
	    {
	      cout << "nfes... "<< nfes <<endl;
		cout << lowestDistanceFactor <<endl;
              save_front(filename2); //save the objective space information
	    }
	}
	save_pos(filename1); //save the decision variable space information
        save_front(filename2); //save the objective space information
	population.clear();
}
void MOEA::save_front(char saveFilename[1024])
{

    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename,fstream::app|fstream::out );
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<archive[n].y_obj[k]<<"  ";
	for(int k=0;k<nobj;k++)
			fout<<child_pop[n].y_obj[k]<<"  ";
	for(int k=0;k<nobj;k++)
			fout<<population[n].y_obj[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void MOEA::save_pos(char saveFilename[1024])
{
    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename, fstream::app|fstream::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<population[n].x_var[k] << "  ";
			//fout<<population[n].indiv.x_var[k]<< fixed << setprecision(30) << "  ";
//	  for(int k=0;k<nvar;k++)
//			fout<<child_pop[n].x_var[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}



#endif
