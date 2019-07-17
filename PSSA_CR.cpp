/* 
 The population vector is augmented with two components representing the number of sources and the volume.
 Accordingly, the stoichiometric matrix is augmented with two columns representing the variation of the sources and volume.
 Since the sources never change (I always have one source), its stoichiometry is always 0; similarly, since the volume 
increases deterministically, its stoichiometry is always 0.
In this way I can model reactions whose rates are proportional to the volume or that comes from a source (e.g. 0->A)

I suppose that the reactions are proportional to the number of reactants and not concentration. This amounts to say that the
concentration of ribosomes and polymerases are constant, i.e. I'm modelling an "intrinsic noise" model.

Right after a cell division, the next reaction time is recomputed

ToDo
Remove the time interval skip errror message
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

using namespace std;

//Create the random distributions
std::random_device rd; 
std::mt19937 rnd_gen (rd ());
std::uniform_real_distribution<> dist_unif (0, 1);
std::exponential_distribution<> dist_exp (1);
std::normal_distribution<> dist_lambda;
std::normal_distribution<> dist_DV;
std::binomial_distribution<> dist_binom;

int parse_file(string filename, vector<double>& population, vector<double>& rates, vector<int>& nu, vector<vector<int>>& reactants,
               vector<double>& growth, double& time_step, vector<double>& time_range,
               vector<string>& species_names, vector<string>& reactions_names)
{
  ifstream file(filename.c_str());
  if(!file.is_open()){
    cerr << "Error: couldn't open the reaction file" << endl;
    return(-1);
  }
  
  string line;
  int tmp_int;
  double tmp_double;
  string tmp_string;
  
  //These are used to check that I've read all the required info
  bool got_population=false, got_matrix=false, got_rates=false, got_reactants=false, got_timestep=false;
  //These are optional arguments
  bool got_growth=false, got_times=false, got_species_names=false, got_reactions_names=false;
  
  while(getline(file, line)){
    if(line[0] == '#')
      continue;
    if(line.compare("!population")==0){
      getline(file, line);
      istringstream iss(line);
      while(iss >> tmp_double) population.push_back(tmp_double);
      got_population = true;
    }
    if(line.compare("!rates")==0){
      getline(file, line);
      istringstream iss(line);
      while(iss >> tmp_double) rates.push_back(tmp_double);
      got_rates = true;
    }
    if(line.compare("!matrix")==0){
      getline(file, line, ';');
      istringstream iss(line);
      while(iss >> tmp_int) nu.push_back(tmp_int);
      got_matrix = true;
    }
    if(line.compare("!reactants")==0){
      getline(file, line);
      while(!line.empty()){
        istringstream iss(line);
        reactants.push_back(vector<int>());
        while(iss >> tmp_int) reactants.back().push_back(tmp_int);
        getline(file, line);
      }
      got_reactants = true;
    }
    if(line.compare("!cell_growth")==0){
      getline(file, line);
      istringstream iss(line);
      while(iss >> tmp_double) growth.push_back(tmp_double);
      got_growth = true;
    }
    if(line.compare("!time_step")==0){
      file >> time_step;
      got_timestep = true;
    }
    if(line.compare("!time_range")==0){
      getline(file, line);
      istringstream iss(line);
      while(iss >> tmp_double) time_range.push_back(tmp_double);
      got_times = true;
    }
    if(line.compare("!species_names")==0){
      getline(file, line);
      istringstream iss(line);
      while(iss >> tmp_string) species_names.push_back(tmp_string);
      got_species_names = true;
    }
    if(line.compare("!reactions_names")==0){
      getline(file, line);
      istringstream iss(line);
      while(iss >> tmp_string) reactions_names.push_back(tmp_string);
      got_reactions_names = true;
    }      
  }
  
  //Check that the parameters are given correctly
  if(got_population && got_rates && got_matrix && got_reactants && got_timestep == false){
    cerr << "Error, some parameters are missing!" << endl;
    return(-1);
  }
  if(population.size()*rates.size()!=nu.size()){
    cerr << "The dimension of the stoichiometric matrix is not correct" << endl;
    return(-1);
  }
  if(got_growth && growth.size()!=5){
    cerr << "The growth of the cell is not correctly specified" << endl;
    return(-1);
  }
  if(got_species_names && species_names.size()!=population.size()){
    cerr << "The number of species is different from the number of names" << endl;
    return(-1);
  }
  if(got_reactions_names && reactions_names.size()!=rates.size()){
    cerr << "The number of reactions is different from the number of names" << endl;
    return(-1);
  }
  
  
  //Based on the reactant list I have to compute the  propensities.
  //Reactants marked as -1 are the sources, but since I cannot have negative index in the partial propensities array,
  //I rename -1 to max_species (i.e. the source corresponds to species number max_species)
  int num_species = population.size();
  for(auto it=reactants.begin(); it!=reactants.end(); it++){
    transform(it->begin(), it->end(), it->begin(), [num_species](int x){
      if(x==-1) return(num_species);
      if(x==-2) return(num_species+1);
      return(x);
    });
  }
  
  //If the time range is not specified, I set it to a default value of 0 to 1000
  if(!got_times){
    time_range.push_back(0);
    time_range.push_back(1000);
  }
  
  //If the species and/or reactions names are not given, use default ones
  if(!got_species_names){
    for(size_t i=0; i<population.size(); i++){
      std::stringstream ss;
      ss << "Species_" << i;
      species_names.push_back(ss.str());
    }
  }
  if(!got_reactions_names){
    for(size_t i=0; i<rates.size(); i++){
      std::stringstream ss;
      ss << "Reaction_" << i;
      reactions_names.push_back(ss.str());
    }
  }
  
  //Now I add two rows to the stoichiometric matrix which are made of zeros and represent the source and volume term.
  for(size_t j=population.size(); j<=rates.size()*(population.size()+2); j+=(population.size()+2))  nu.insert(nu.begin()+j, 2, 0);
  
  //I add the source term and the volume term to the population
  double V0 = (growth.size()==0) ? 1 : growth[0];
  population.push_back(1);
  population.push_back(V0);
  
  //Set the parameters of the distribution for the growth rate and added volume, if the cell is growing
  if(growth.size()!=0){
    dist_lambda.param(std::normal_distribution<>::param_type(growth[1], growth[2]));
    dist_DV.param(std::normal_distribution<>::param_type(growth[3], growth[4]));
  }
  
  return(0);
}

template <class T>
inline void print_matrix(T matrix, int n_cols){
  int n_rows = matrix.size()/n_cols;
  for(int row = 0; row < n_rows; row++){
    for(auto it = matrix.begin()+row*n_cols; it<matrix.begin() + (row+1)*n_cols && it<matrix.end(); it++)
      cerr << *it << "\t";
    cerr << endl;
  }
}

void print_params(vector<double> population, vector<double> rates, vector<int> nu, vector<vector<int>> reactants_list,
                  vector<double> growth, double time_step, vector<double> time_range)
{
  cerr << "*----- Parameters -----*" << endl;
  cerr << "Extended initial population:  ";
  for(auto it : population) cerr << it << " ";
  cerr << "\nSpecific rates: ";
  for(auto it : rates) cerr << it << " ";
  cerr << "\nExtended stoichiometric matrix:\n";
  print_matrix(nu, population.size());
  cerr << "\nReactants list:\n";
  for(auto it : reactants_list){
    for(auto it2 : it) cerr << it2 << " ";
    cerr << endl;
  }
  if(growth.size()!=0){
    cerr << "\nGrowth paramaters (initial volume, mean and standard deviation of growth rate and added volume)\n";
    for(auto it : growth) cerr << it << " ";
  }
  cerr << "\nTime step for printing: " << time_step;
  cerr << "\nTime range to report: " << time_range[0] << " - " << time_range[1] << endl;
}

inline void print_results(vector<double> population, vector<double> rates, double time, int cell, int parent){
  cout << time << " ";
  for(auto it=population.begin(); it<population.end(); it++) cout << *it << "\t";
  for(auto it : rates) cout << it << "\t";
  cout << cell << "\t" << parent << endl;
}

void run_simulation(vector<double>& population, vector<double> rates, vector<int> nu, vector<vector<int>> reactants_list,
                    vector<double> growth, double time_step, vector<double> time_range){
  size_t N = population.size();
  size_t M = rates.size();
  
  double next_t, cum_time=0, time_next_reaction; //Simulation times
  int current_interval = 0;
  
  vector<int> reactants_molecules(M, 1); //Number of molecules associated to each reaction
  vector<double> propensities(M, 0.); // Propensities = rate*number of molecules of reactants
  double prop_tot; //Sum of all the propensities
  vector<double> propensity_cum_sum(M, 0); // Cumulative sum of the propensity vector
  
  double r; //Uniform random number to select next reaction
  int next_reaction;
  
  //Growth parameters
  double lambda;
  double &V = population.back(); //Everytime I modify V, I modify the last element of the population
  double DV;
  double Vdiv;
  int cell;
  int parent;
  if(growth.size()!=0){
    lambda = dist_lambda (rnd_gen);
    V = growth[0]; DV=dist_DV (rnd_gen);
  }
  else{
    lambda = 0;
    V=1; DV=2;
  }
  Vdiv = V+DV;
  cell=0; parent=-1;
  
  
  //PERFORM THE REACTION STEPS
  while(cum_time < time_range[1]){
    //Find the number of molecules associated to each reaction
    reactants_molecules.assign(M, 1);
    for(size_t i=0; i<M; i++){
      for(size_t j=0; j<reactants_list[i].size(); j++) reactants_molecules[i]*=population[reactants_list[i][j]];
    }
    
    //Find propensities
    transform(rates.begin(), rates.end(), reactants_molecules.begin(), propensities.begin(), multiplies<double>());
    
    //Find next reaction time
    prop_tot = accumulate(propensities.begin(), propensities.end(), 0.);
    dist_exp.param(std::exponential_distribution<>::param_type(prop_tot)); //Update rate of exponential
    next_t = dist_exp (rnd_gen);
    time_next_reaction = cum_time + next_t;
    
    //From current time until division or next reaction time, nothing changes, except the volume
    if(V*exp(lambda*next_t) < Vdiv){
      V *= exp(lambda*next_t);  //Remember: V is a reference to the last element of population
      cum_time = time_next_reaction;
      
      //Find next reaction
      partial_sum (propensities.begin(), propensities.end(), propensity_cum_sum.begin());
      
      r = dist_unif(rnd_gen)*prop_tot;
      next_reaction = distance(propensity_cum_sum.begin(),
                               find_if(propensity_cum_sum.begin(), propensity_cum_sum.end(), [r](double &x){return r<x;}));
      
      //Update the population
      transform(nu.begin()+next_reaction*N, nu.begin()+(next_reaction+1)*N, population.begin(), population.begin(), plus<double>());
    }
    else{   //If a division has occured before I reached the next reaction, I divide the cell
      cum_time += log(Vdiv/V)/lambda;  //Go to the time when the cell reaches Vdiv
      //Generate new growth values
      lambda = dist_lambda (rnd_gen);
      DV = dist_DV (rnd_gen);
      V = Vdiv/2; //It's been shown that the volume of the daughter is half the one of the mother with very low noise (2%)
      Vdiv = V + DV;
      
      //Assign new cell number
      parent = cell;
      cell += 1;  
      
      //Binomial partition of the species
      for(auto it=population.begin(); it<population.end()-2; it++){
        dist_binom.param(std::binomial_distribution<>::param_type(*it, 0.5)); 
        *it = dist_binom(rnd_gen);
      }
    }
    
    //Print the results if I am in a new time interval
    if(time_step<=0){
      print_results(population, rates, cum_time, cell, parent);
    }else{
      if(floor(cum_time/time_step) > current_interval){
        if(floor(cum_time/time_step) != current_interval+1) cerr << "Error! One time interval has been skipped" << endl;
        current_interval = floor(cum_time/time_step);
        if(cum_time>time_range[0])  print_results(population, rates, current_interval*time_step, cell, parent);
      }
    }
    
  }  
}

void print_header(vector<string> species_names, vector<string> reactions_names){
  cout << "time" << "\t";
  for(auto it : species_names) cout << it << "\t";
  cout << "source" << "\t" <<  "length" << "\t";
  for(auto it : reactions_names) cout << it << "\t";
  cout << "cell_id" << "\t" << "parent_id" << endl;
}

int main(int argc, char**argv){
  if(argc!=2){
    cerr << "Number of arguments not correct" << endl;
    return(-1);
  }
  
  string filename = argv[1];
  
  //Store the initial population, reaction rates and stoichiometric matrix
  vector<double> population;
  vector<double> rates;
  vector<int> nu;
  vector<vector<int>> reactants_list;
  vector<double> growth;
  vector<double> time_range;
  vector<string> species_names;
  vector<string> reactions_names;
  
  //Store some parameters for the simulation
  double time_step;
  
  if(parse_file(filename, population, rates, nu, reactants_list, growth, time_step, time_range, species_names, reactions_names) != 0)
    return(-1);
  print_header(species_names, reactions_names);
  print_params(population, rates, nu, reactants_list, growth, time_step, time_range);
  run_simulation(population, rates, nu, reactants_list, growth, time_step, time_range);
  
  return(0);
}
