/*
 * Uses Rcpp and BH for ode integration
 * Defines the dynamical system
 * Borrows concepts from Boost odeint examples
 * Borrows concepts and code snippets from StackOverflow mostly answered by
 * Dirk Eddelbuttel
 */

#include "header.h"
#include <Rcpp.h>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/lexical_cast.hpp>

// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;
using boost::lexical_cast;
using boost::bad_lexical_cast;
// forward declartion -- defined in header.cpp
size_t convertAdjMatToVector(
    IntegerMatrix gene_interaction, std::vector<size_t>& tgtGene,
    std::vector<std::pair<size_t,size_t> >& intSrcType);

void generateThresholds(IntegerMatrix gene_interaction,
                        NumericVector threshold_gene,
                        Rcpp::List simulationData);

extern double Hs_Racipe(double A, double AB0, size_t n_ab, double lambda_ab);

// Define container type for the dynamical system
typedef std::vector< double > container_type;

// Define Euler Maruyama stepper for stochastic integration
template<size_t N> class stochastic_euler
{
public:

  typedef container_type state_type;
  typedef container_type deriv_type;
  typedef double value_type;
  typedef double time_type;
  typedef unsigned short order_type;
  typedef boost::numeric::odeint::stepper_tag stepper_category;

  static order_type order( void ) { return 1; }

  template< class System >
  void do_step(
      System system , container_type &x , time_type t , time_type dt ) const
  {
    deriv_type det , stoch ;
    system.first( x , det );
    system.second( x , stoch );
    for( size_t i=0 ; i<x.size() ; ++i )
      x[i] += dt * det[i] + sqrt( dt ) * stoch[i];
  }
};

struct stoch
{
  boost::mt19937 &m_rng;
  boost::random::normal_distribution<> m_dist;

  stoch( boost::mt19937 &rng , double sigma ) : m_rng( rng ) ,
  m_dist( 0.0 , sigma ) { }

  void operator()( const container_type &x , container_type &dxdt )
  {
    dxdt[0] = m_dist( m_rng );
  }
};

// Define the dynamical system
class genRegNet
{
  container_type params;
  //  container_type state;
  std::vector<size_t> tgtGene;
  //vector containing source and type of nth interaction
  std::vector<std::pair<size_t,size_t> > intSrcType;
  //number of interactions

  // variables to shift the pointer location to the desired category of params
  size_t nInteractions;
  size_t number_gene;
  size_t kStart;
  size_t nStart;
  size_t lambdaStart;
  size_t threshStart;


public:
  //default constructor
  // genRegNet() = default;

  // Initialize the dynamical system with the network
  void setNetwork(const IntegerMatrix gene_interaction){
    //call the function to initialize tftGene and intSrcType using the
    //adjacency matrix gene_interaction
    nInteractions = convertAdjMatToVector(gene_interaction, tgtGene, intSrcType);
    number_gene = gene_interaction.ncol();
    //  state.resize(number_gene);
    params.resize(2*number_gene + 3*nInteractions);
    kStart = number_gene;
    nStart = 2*number_gene;
    lambdaStart = 2*number_gene + nInteractions;
    threshStart = 2*number_gene + 2*nInteractions;
  }


  // set the parameters and initial conditions of the dynamical system
  void setParamsIc(const NumericVector& ic, const NumericVector modelParams )
  {
    //   if(state.size() == ic.size()){
    //    for (size_t i=0; i<ic.size(); i++){
    //    state[i] = ic[i];
    //    }
    //    }
    if(params.size() == modelParams.size()){
      for (size_t i=0; i<modelParams.size(); i++){
        params[i] = modelParams[i];
      }
    }

    //    if(params.size() != modelParams.size() || state.size() != ic.size())
    if(params.size() != modelParams.size())

    {stop("Size mismatch between initial conditions and state.");}
  }


  //  void setParams(const NumericVector modelParams ) {
  //    for (size_t i=0; i<modelParams.size(); i++){
  //    params.push_back(modelParams[i]);
  //      }
  //  }

// parameter accessor
  container_type getParams( void ) const { return params; }

  // A print fuction -- used in debugging
  void print(){
    Rcout<<"Number of interactions:\t"<<nInteractions
    <<"\t Number of genes \t"<<number_gene<<"\n";
    //--Debugging mode
    //Rcout<<kStart<<"\t"<<nStart<<"\t"<<lambdaStart<<"\t"<<threshStart<<"\n";
    Rcout<<"Number of parameters:\t"<<params.size()<<"\n";

    //--Debugging mode
    //Rcout<<state.size()<<"\n";
    Rcout<<"Parameters:\t";
    for(size_t i=0;i<params.size();i++)
    {Rcout<<params[i]<<"\t";}
    Rcout<<"\n";
    //--Debugging mode
    //  for(size_t i=0;i<state.size();i++)
    //  {Rcout<<state[i]<<"\t";}
  }

  // Operator overloading for use in boost odeint integrators
  void operator()(
      const container_type& state, container_type& dxdt,
      double /* t */ ) const
  {

    size_t intCounter = 0;
    for( size_t i = 0 ; i < state.size(); i++ ){
      double hillMultiplier = 1.0;
      // Check if the ith gene is target of intCounter interaction
      // If yes, then include its effect
      while(tgtGene[intCounter] <= i){
        if(tgtGene[intCounter] == i){
          double gene_lambda = params[lambdaStart + intCounter];
          if(intSrcType[intCounter].second == 2){
            // For an inhibitory link invert the fold change
            gene_lambda = 1./gene_lambda;
          }
          // Multiply the hillMultiplier with the evaluated hill function
          hillMultiplier = hillMultiplier*Hs_Racipe(
            state[intSrcType[intCounter].first],
            params[threshStart + intCounter],
            static_cast<size_t>(params[nStart + intCounter]), gene_lambda);
        }
        intCounter++;
      }
      // Define dxdt = g*H - kx
      dxdt[i] = hillMultiplier*params[i] - params[kStart + i]*state[i];
    }
  }
};


// function to set initial conditions
// [[Rcpp::export]]
void selectIc(NumericMatrix& ic, const IntegerMatrix& gene_interaction,
                     const NumericMatrix& modelParameters)
{
  // Rcout<<"genModelParamsSelectIcTrue"<<"\n";
  size_t number_gene = gene_interaction.ncol();
  size_t model_count_max = modelParameters.nrow();
  size_t nIc = std::round(ic.ncol()/number_gene);

  //vector containing tgtGene of nth interaction
  std::vector<size_t> tgtGene;
  //vector containing source and type of nth interaction
  std::vector<std::pair<size_t,size_t> > intSrcType;
  //number of interactions
  size_t nInteractions =0;
  //call the function to initialize tftGene and intSrcType using the
  //adjacency matrix gene_interaction
  nInteractions = convertAdjMatToVector(gene_interaction, tgtGene, intSrcType);

  for(size_t model_counter=0; model_counter<model_count_max;model_counter++){
    // vectors to hold the range of initial conditions for a specific model
    // the initial conditions will be chosen from this range (log distribution)
    std::vector<double> min_gene(number_gene);
    std::vector<double> max_gene(number_gene);

    for(size_t gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      double value = modelParameters(model_counter, gene_count1)/
        modelParameters(model_counter, (number_gene + gene_count1));
      max_gene[gene_count1] = value;
      min_gene[gene_count1] = value;
    }
    for(size_t intCounter = 0; intCounter < nInteractions; intCounter++)
      {
      // Divide the minimum value by the fold change of that interaction
      min_gene[tgtGene[intCounter]] =
        min_gene[tgtGene[intCounter]]/modelParameters(
            model_counter,(2*number_gene + intCounter + nInteractions));
      }
    // Now select and store initial conditions
    for(size_t gene_count1=0; gene_count1<number_gene; gene_count1++)
    {
    for(size_t icCounter=0;icCounter< nIc; icCounter++){
      ic(model_counter, (number_gene*icCounter + gene_count1)) =
        std::exp(std::log(min_gene[gene_count1])+(log(max_gene[gene_count1]) -
        log(min_gene[gene_count1]))*u_distribution(u_generator));
      }
    }
  }
}


// continer for storing gene expression values at each time point
// Borrowed from odeint examples
// struct push_back_state_and_time
// {
//   std::vector< container_type >& m_states;
//   std::vector< double >& m_times;
//   size_t m_steps;
//  size_t& max_steps;
//
//   push_back_state_and_time(
//     std::vector< container_type > &states , std::vector< double > &times, size_t& max_steps )
//     : m_states( states ) , m_times( times ), max_steps( max_steps ), m_steps(0) { }
//
//   void operator()( const container_type &x , double t )
//   {
//     m_states.push_back( x );
//     m_times.push_back( t );
//     ++m_steps;
//     if( m_steps > max_steps ) throw std::runtime_error( "Simulation Complete" );
//   }
// };

struct push_back_state_and_time
{
  std::vector< container_type >& allStates;
  std::vector< double >& allTimes;

  push_back_state_and_time( std::vector< container_type > &states , std::vector< double > &times )
    : allStates( states ) , allTimes( times ) { }

  void operator()( const container_type &x , double t )
  {
    allStates.push_back( x );
    allTimes.push_back( t );
    Rcout<<allTimes.size()<<"\n";
  }
};

struct expObserver
{
  NumericMatrix geneExp;
  size_t timeCounter;
  expObserver(NumericMatrix geneExpression) :
    geneExp(geneExpression), timeCounter(0){};

  //push_back_state_and_time( std::vector< container_type > &states , std::vector< double > &times )
  //  : allStates( states ) , allTimes( times ) { }

  void operator()( const container_type &x , double t )
  {
    for(size_t geneCounter =0; geneCounter < geneExp.ncol(); geneCounter++)
    {
      geneExp(timeCounter, geneCounter) = x[geneCounter];
    }

  timeCounter++;
  }
};



void write_cout( const container_type &x , double t, size_t step )
{
  Rcout << step<<"\t"<<t << '\t' << x[1] << "\n";
}

///////////////////////////////////////////////////////////////////////////////
// A do-all function. bool paramters to control calculations.
// [[Rcpp::export]]
void simulateNetwork(NumericMatrix geneExpression,
    const IntegerMatrix gene_interaction, NumericMatrix& modelParameters,
    NumericMatrix& ic, Rcpp::List simulationData,
    const bool genThresh = true,
    const bool genModelParams = true, const bool genIc = true,
    const bool integrate = true,
    const bool timeSeries = false,
    const bool annealing = false)
  {
  // Initialize the network
  size_t number_gene = gene_interaction.ncol();
  //vector containing tgtGene of nth interaction
  std::vector<size_t> tgtGene;
  //vector containing source and type of nth interaction
  std::vector<std::pair<size_t,size_t> > intSrcType;
  //number of interactions
  size_t nInteractions =0;
  //call the function to initialize tftGene and intSrcType using the
  //adjacency matrix gene_interaction
  nInteractions = convertAdjMatToVector(gene_interaction, tgtGene, intSrcType);
  // numer of genes in the network

// Initialize the parameters from R objects

  NumericVector simulationParameters =
    as<NumericVector>(simulationData["simParams"]);
  NumericVector stochasticParameters =
    as<NumericVector>(simulationData["stochParams"]);
  NumericVector hyperParameters =
    as<NumericVector>(simulationData["hyperParams"]);

  NumericVector threshold_gene;
  if(genThresh){
    generateThresholds( gene_interaction,
                        threshold_gene,
                       simulationData);
  }
  else{
    threshold_gene = as<NumericVector>(simulationData["thresholds"]);
  }


  size_t model_count_max = static_cast<size_t>(simulationParameters[0]);
  double parameter_range = simulationParameters[1];
  double dt = simulationParameters[2];
  double tot_time = simulationParameters[3];
  size_t INITIAL_CONDITIONS = static_cast<size_t> (simulationParameters[4]);
  size_t output_precision = static_cast<size_t> (simulationParameters[5]);

  size_t D_levels = static_cast<size_t>(stochasticParameters[0]);
  double D_scaling = stochasticParameters[1];
  double D_max = stochasticParameters[2];
  bool gene_noise_scaling = static_cast<bool>(stochasticParameters[3]);
  double D_shot = static_cast<bool>(stochasticParameters[4]);

  double g_min = hyperParameters[0];
  double g_max = hyperParameters[1];
  double k_min = hyperParameters[2];
  double k_max = hyperParameters[3];
  double lambda_min = hyperParameters[4];
  double lambda_max = hyperParameters[5];
  size_t n_min = static_cast<size_t>(hyperParameters[6]);
  size_t n_max = static_cast<size_t>(hyperParameters[7]);
  size_t possible_interactions = static_cast<size_t>(hyperParameters[8]);
  size_t threshold_max = static_cast<size_t>(hyperParameters[9]);
  double standard_deviation_factor = hyperParameters[10];

  size_t max_steps = static_cast<size_t>(tot_time/dt);
  applyParameterRange(parameter_range, g_min, g_max, k_min, k_max,
                      lambda_min, lambda_max);

 // If genModelParams is true, generate model paramters
 // and store in modelParameters
if(genModelParams){

  size_t kStart = number_gene;
  size_t nStart = 2*number_gene;
  size_t lambdaStart = 2*number_gene + nInteractions;
  size_t threshStart = 2*number_gene + 2*nInteractions;

  for(size_t model_counter=0; model_counter < model_count_max; model_counter++){
    for(size_t gene_count1=0;gene_count1<number_gene;gene_count1++){
      //Initialize production rate of genes
      modelParameters(model_counter, gene_count1) = g_min +
        (g_max-g_min)*u_distribution(u_generator);
   //   Rcout<<modelParameters(model_counter, gene_count1)<<"\n";

      //Initialize degradation rate of genes
      modelParameters(model_counter, (kStart + gene_count1)) =
        k_min+(k_max-k_min)*u_distribution(u_generator);
    }

      for(size_t int_counter=0;int_counter<nInteractions;int_counter++)
      {

        // n for the inward links for genes starting from gene 1
        modelParameters(model_counter, (nStart + int_counter)) =
          static_cast<double> ( round(n_min +
          (n_max-n_min)*u_distribution(u_generator)));

        // lambda for the inward links for genes starting from gene 1
        modelParameters(model_counter, (lambdaStart + int_counter)) =
                          (lambda_max-lambda_min)*u_distribution(u_generator) +
                          lambda_min;

        //Initialize threshold for each interaction
        // Thresholds for inward links for genes starting from gene 1
        //(thresholds for all inward links of gene 1 and so on)
        modelParameters(model_counter, (threshStart + int_counter)) =
          (1-standard_deviation_factor*sqrt(3))*
          threshold_gene[intSrcType[int_counter].first] +
          (2*sqrt(3)*standard_deviation_factor*
          threshold_gene[intSrcType[int_counter].first])*
          u_distribution(u_generator);
      }
  }
}

if(genIc){
  selectIc( ic, gene_interaction, modelParameters);
}
if(integrate){
  genRegNet model;
  model.setNetwork(gene_interaction);
  for(size_t modelCounter = 0; modelCounter < model_count_max; modelCounter++)
  {
    NumericVector modelIc = ic(modelCounter,_);
    NumericVector modelParams = modelParameters(modelCounter,_);
    model.setParamsIc(modelIc, modelParams);
    //--Debug Mode
    // model.print();
    container_type x;
    for(size_t i=0; i < number_gene; i++){
      x.push_back(modelIc[i]);
    }


    if(!timeSeries)
    {
      size_t  steps = boost::numeric::odeint::integrate(
      model , x, 0.0 , tot_time , dt);

    for(size_t i=0; i < number_gene; i++){
      geneExpression(modelCounter, i) = x[i];
    }
    }
    else {

      size_t steps = 0;
      // x_vec and times used only for time series analysis
  //    std::vector<container_type> x_vec;
      std::vector<double> times;
      expObserver obs(geneExpression);

      // calculate some transients steps
      //   size_t  steps = boost::numeric::odeint::integrate(testCase , x, 0.0 , 100.0 , dt,
      //                                                   push_back_state_and_time( x_vec , times ) );

 //     try {
//        steps = boost::numeric::odeint::integrate( model, x, 0.0, 50000, dt, push_back_state_and_time( x_vec, times, max_steps) );
        //Rcout<<max_steps<<"\n";

   //if(integrate)
     {
        boost::numeric::odeint::runge_kutta4< container_type > stepper;


//     boost::numeric::odeint::integrate_n_steps(
//       stepper , model , x , 0.0 , dt , max_steps-1,
//       write_cout(x,t,steps) );

//        boost::numeric::odeint::integrate_n_steps(
//          stepper , model , x , 0.0 , dt , max_steps-1,
//          push_back_state_and_time( x_vec , times) );

        boost::numeric::odeint::integrate_n_steps(
          stepper ,boost::ref( model) , x , 0.0 , dt , max_steps-1,
          boost::ref(obs));

        std::vector<std::string> name;
        double tempT=0;
        for(size_t counter=0; counter<geneExpression.nrow();counter++)
        {
          name.push_back(std::to_string(tempT));
          tempT+=dt;
        }
        rownames(geneExpression) = Rcpp::wrap(name);
   }
  // return;
  //      steps = boost::numeric::odeint::integrate(
  //       model , x, 0.0 , 10000.0 , dt,
  //       push_back_state_and_time( x_vec , times, max_steps ) );
  //    } catch( const std::runtime_error& error ) { steps = max_steps;

//        return;}
    //  Rcout<<"Size"<<x_vec.size()<<"\t"<<times.size()<<"\n";
    //  Rcout<<"gESize"<<geneExpression.nrow()<<"\n";

    //  for(size_t stepCount = 0; stepCount<geneExpression.nrow(); stepCount++){
    //  if(stepCount < x_vec.size()){
    //  name.push_back(std::to_string(times[stepCount]));
    //  for(size_t i=0; i < number_gene; i++){
    //    geneExpression(stepCount, i) = x_vec[stepCount][i];
    //    }
    //  }
    //  }



    }
    // Rcout<<"genModelParamsTrue"<<"\n";
  }

}

 return;
}
