/*
 * Uses Rcpp and BH for ode integration
 * Defines the dynamical system
 * Borrows concepts from Boost odeint examples
 * Borrows concepts and code snippets from StackOverflow mostly answered by
 * Dirk Eddelbuttel
 */

#include"header.h"
#include <Rcpp.h>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/lexical_cast.hpp>

// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;

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

// // Define Euler Maruyama stepper for stochastic integration
// template<size_t N> class stochastic_euler
// {
// public:
//
//   typedef container_type state_type;
//   typedef container_type deriv_type;
//   typedef double value_type;
//   typedef double time_type;
//   typedef unsigned short order_type;
//   typedef boost::numeric::odeint::stepper_tag stepper_category;
//
//   static order_type order( void ) { return 1; }
//
//   template< class System >
//   void do_step(
//       System system , container_type &x , time_type t , time_type dt ) const
//   {
//     deriv_type det , stoch ;
//     system.first( x , det );
//     system.second( x , stoch );
//     for( size_t i=0 ; i<x.size() ; ++i )
//       x[i] += dt * det[i] + sqrt( dt ) * stoch[i];
//   }
// };
//
// struct stoch
// {
//   boost::mt19937 &m_rng;
//   boost::random::normal_distribution<> m_dist;
//
//   stoch( boost::mt19937 &rng , double sigma ) : m_rng( rng ) ,
//   m_dist( 0.0 , sigma ) { }
//
//   void operator()( const container_type &x , container_type &dxdt )
//   {
//     dxdt[0] = m_dist( m_rng );
//   }
// };

// Define the dynamical system
class genRegNet
{
  std::vector<double> params;
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
  void setParams(const std::vector<double>& modelParams )
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

size_t getParamsSize(void){return params.size();}
// parameter accessor
std::vector<double> getParams( void ) const { return params; }

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
      const std::vector<double>& state, std::vector<double>& dxdt,
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
void selectIc(const IntegerMatrix& gene_interaction, const size_t& nIc,
              const size_t& model_count_max, const size_t& outputPrecision,
              std::string fileName)
{
  std::fstream outIc("./tmp/ic" + fileName + ".txt",std::ios::out);

   std::ifstream inParameters("./tmp/parameters" + fileName + ".txt",std::ifstream::in);

  // Rcout<<"genModelParamsSelectIcTrue"<<"\n";
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
  size_t paramsSize = 2*number_gene + 3*nInteractions;

  for(size_t model_counter=0; model_counter<model_count_max;model_counter++){
    std::vector<double> modelParams;
    modelParams.resize(paramsSize);
    double tmpParameter = 0.0;
    for(size_t tmpParCounter=0; tmpParCounter < paramsSize; tmpParCounter++){
      inParameters >> tmpParameter;
      modelParams[tmpParCounter] = tmpParameter;
    }
    // vectors to hold the range of initial conditions for a specific model
    // the initial conditions will be chosen from this range (log distribution)
    std::vector<double> min_gene(number_gene);
    std::vector<double> max_gene(number_gene);

    for(size_t gene_count1=0;gene_count1<number_gene;gene_count1++)
    {
      double value = modelParams[gene_count1]/
        modelParams[number_gene + gene_count1];
      max_gene[gene_count1] = value;
      min_gene[gene_count1] = value;
    }
    for(size_t intCounter = 0; intCounter < nInteractions; intCounter++)
      {
      // Divide the minimum value by the fold change of that interaction
      min_gene[tgtGene[intCounter]] =
        min_gene[tgtGene[intCounter]]/modelParams[(2*number_gene + intCounter + nInteractions)];
      }
    // Now select and store initial conditions

    for(size_t icCounter=0;icCounter< nIc; icCounter++){
      for(size_t gene_count1=0; gene_count1<number_gene; gene_count1++)
      {
      outIc<<std::setprecision(outputPrecision)<<
        std::exp(std::log(min_gene[gene_count1])+(log(max_gene[gene_count1]) -
        log(min_gene[gene_count1]))*u_distribution(u_generator))<<"\t";
      }
      outIc<<"\n";
    }
  }

  outIc.close();
  inParameters.close();
}


struct expObserver
{


  size_t outputPrecision;
  std::fstream& outGE;

  expObserver(size_t& outputPrec, std::fstream& outFile) :
    outputPrecision(outputPrec ),
    outGE(outFile) {  };

  void operator()( const std::vector<double> &x , double t )
  {
    if(outGE.is_open()){
    outGE<<std::setprecision(outputPrecision)<<t<<"\t";
    for(size_t geneCounter =0; geneCounter < x.size(); geneCounter++)
    {
      outGE<<std::setprecision(outputPrecision)<<x[geneCounter]<<"\t";
    }
    outGE<<"\n";
    }
    else{(stop("Error in opening the output file"));}
  }
};

//
//
// void write_cout( const std::vector<double> &x , double t, size_t step )
// {
//   Rcout << step<<"\t"<<t << '\t' << x[1] << "\n";
// }

///////////////////////////////////////////////////////////////////////////////
// A do-all function. bool paramters to control calculations.
// [[Rcpp::export]]
void simulateNetwork(
    const IntegerMatrix gene_interaction,
    Rcpp::List simulationData, String tmpFileName,
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
  // if(genThresh){
  //   generateThresholds( gene_interaction,
  //                       threshold_gene,
  //                      simulationData);
  // }
 //  else{
    threshold_gene = as<NumericVector>(simulationData["thresholds"]);
 // }

  size_t model_count_max = static_cast<size_t>(simulationParameters[0]);
  double parameter_range = simulationParameters[1];
  double dt = simulationParameters[2];
  double tot_time = simulationParameters[3];
  size_t nIc = static_cast<size_t> (simulationParameters[4]);
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

  std::string fileName= tmpFileName;


 // If genModelParams is true, generate model paramters
 // and store in modelParameters
if(genModelParams){

  std::fstream outParameters("./tmp/parameters" + fileName + ".txt",std::ios::out);

  size_t kStart = number_gene;
  size_t nStart = 2*number_gene;
  size_t lambdaStart = 2*number_gene + nInteractions;
  size_t threshStart = 2*number_gene + 2*nInteractions;

  for(size_t model_counter=0; model_counter < model_count_max; model_counter++){
    for(size_t gene_count1=0;gene_count1<number_gene;gene_count1++){
      //Initialize production rate of genes
      outParameters<<std::setprecision(output_precision)<<
        (g_min + (g_max-g_min)*u_distribution(u_generator))<<"\t";
    }
    for(size_t gene_count1=0;gene_count1<number_gene;gene_count1++){
      //Initialize degradation rate of genes
      outParameters<<std::setprecision(output_precision)<<
        k_min+(k_max-k_min)*u_distribution(u_generator)<<"\t";
    }

      for(size_t int_counter=0;int_counter<nInteractions;int_counter++)
      {

        // n for the inward links for genes starting from gene 1

        outParameters<<std::setprecision(1)<<
          static_cast<double> ( round(n_min +
          (n_max-n_min)*u_distribution(u_generator)))<<"\t";

        // lambda for the inward links for genes starting from gene 1
        outParameters<<std::setprecision(output_precision)<<
                          (lambda_max-lambda_min)*u_distribution(u_generator) +
                          lambda_min<<"\t";

        //Initialize threshold for each interaction
        // Thresholds for inward links for genes starting from gene 1
        //(thresholds for all inward links of gene 1 and so on)
        outParameters<<std::setprecision(output_precision)<<
          (1-standard_deviation_factor*sqrt(3))*
          threshold_gene[intSrcType[int_counter].first] +
          (2*sqrt(3)*standard_deviation_factor*
          threshold_gene[intSrcType[int_counter].first])*
          u_distribution(u_generator)<<"\t";
      }
      outParameters<<"\n";
  }
  outParameters.close();
}

if(genIc){
  selectIc( gene_interaction, nIc, model_count_max, output_precision, fileName);
}

if(integrate){
  std::ifstream inIc("./tmp/ic" + fileName + ".txt",std::ifstream::in);

  std::ifstream inParameters("./tmp/parameters" + fileName + ".txt",std::ifstream::in);

  std::fstream outGE("./tmp/GE" + fileName + ".txt",std::ios::out);

  genRegNet model;
  model.setNetwork(gene_interaction);
  size_t paramsSize = model.getParamsSize();
  for(size_t modelCounter = 0; modelCounter < model_count_max; modelCounter++)
  {
    std::vector<double> modelParams;
    modelParams.resize(paramsSize);
    double tmpParameter = 0.0;
    for(size_t tmpParCounter=0; tmpParCounter < paramsSize; tmpParCounter++){
      if(inParameters >> tmpParameter){
      modelParams[tmpParCounter] = tmpParameter;}
      else {stop("Error in reading parameter file!");}
    }

    std::vector<double> modelIc;
      modelIc.resize(number_gene);
    double tmpIc = 0.0;
    for(size_t tmpParCounter=0; tmpParCounter < number_gene; tmpParCounter++){
      if(inIc >> tmpIc){
        modelIc[tmpParCounter] = tmpIc;}
      else {stop("Error in reading initial condition file!");}
    }

    model.setParams(modelParams);
    //--Debug Mode
    // model.print();
    std::vector<double> x;
    for(size_t i=0; i < number_gene; i++){
      x.push_back(modelIc[i]);
    }


    if(!timeSeries)
    {
   //   std::fstream outGE("./tmp/GE" + fileName + ".txt",std::ios::out);
      size_t  steps = boost::numeric::odeint::integrate(
        boost::ref( model) , x, 0.0 , tot_time , dt);
      for(size_t i=0; i < number_gene; i++){
      outGE<<std::setprecision(output_precision)<<x[i]<<"\t";
        }
      outGE<<"\n";
      }
    else if(timeSeries){
      expObserver obs(output_precision, outGE);
      size_t  steps = boost::numeric::odeint::integrate(
          boost::ref( model) , x, 0.0 , tot_time , dt, boost::ref( obs));

        // boost::numeric::odeint::runge_kutta4< std::vector<double> > stepper;
        // boost::numeric::odeint::integrate_n_steps(
        //   stepper ,boost::ref( model) , x , 0.0 , dt , max_steps-1,
        //   boost::ref(obs));

    }
  }
  outGE<<"\n";

}

 return;
}
