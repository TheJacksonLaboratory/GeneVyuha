// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
using namespace Rcpp;

// Defines the gene regulatory network as a dynamical system

#include <iostream>

#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
// forward declartion -- defined in header.cpp
size_t convertAdjMatToVector(
    IntegerMatrix gene_interaction, std::vector<size_t>& tgtGene,
    std::vector<std::pair<size_t,size_t> >& intSrcType);
extern double Hs_Racipe(double A, double AB0, size_t n_ab, double lambda_ab);

typedef std::vector< double > container_type;

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
  void do_step( System system , container_type &x , time_type t , time_type dt ) const
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

  stoch( boost::mt19937 &rng , double sigma ) : m_rng( rng ) , m_dist( 0.0 , sigma ) { }

  void operator()( const container_type &x , container_type &dxdt )
  {
    dxdt[0] = m_dist( m_rng );
  }
};

class genRegNet
{
  container_type params;
//  container_type state;
  std::vector<size_t> tgtGene;
  //vector containing source and type of nth interaction
  std::vector<std::pair<size_t,size_t> > intSrcType;
  //number of interactions
  size_t nInteractions;
  size_t number_gene;
  size_t kStart;
  size_t nStart;
  size_t lambdaStart;
  size_t threshStart;


public:
  genRegNet() = default;
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

//  void setParams(const NumericVector modelParams ) {
//    for (size_t i=0; i<modelParams.size(); i++){
//    params.push_back(modelParams[i]);
//      }
//  }


  container_type getParams( void ) const { return params; }

void print(){
  Rcout<<nInteractions<<"\t"<<number_gene<<"\t"<<kStart<<"\t"<<nStart<<"\t"
  <<lambdaStart<<"\t"<<threshStart<<"\n";
//  Rcout<<params.size()<<"\t"<<state.size()<<"\n";
  for(size_t i=0;i<params.size();i++)
  {Rcout<<params[i]<<"\t";}
  Rcout<<"\n";
//  for(size_t i=0;i<state.size();i++)
//  {Rcout<<state[i]<<"\t";}
  Rcout<<"\n";
}

  void operator()(
      const container_type& state, container_type& dxdt,
                double /* t */ ) const
  {
    size_t intCounter = 0;
    for( size_t i = 0 ; i < state.size(); i++ ){
      double hillMultiplier = 1.0;
      while(tgtGene[intCounter] <= i){
        if(tgtGene[intCounter] == i){
          double gene_lambda = params[lambdaStart + intCounter];
          if(intSrcType[intCounter].second == 2){
            gene_lambda = 1./gene_lambda;
          }
          hillMultiplier = hillMultiplier*Hs_Racipe(
            state[intSrcType[intCounter].first], params[threshStart + intCounter],
            static_cast<size_t>(params[nStart + intCounter]), gene_lambda);
        }
        intCounter++;
      }
      dxdt[i] = hillMultiplier*params[i] - params[kStart + i]*state[i];
    }
  }
};

struct push_back_state_and_time
{
  std::vector< container_type >& m_states;
  std::vector< double >& m_times;

  push_back_state_and_time( std::vector< container_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

  void operator()( const container_type &x , double t )
  {
    m_states.push_back( x );
    m_times.push_back( t );
  }
};

// Used in debugging mode
// Check the Class functionality
// [[Rcpp::export]]
void testDS(IntegerMatrix gene_interaction, NumericVector params, NumericVector ic){
  genRegNet testCase;
  testCase.setNetwork(gene_interaction);
  testCase.setParamsIc(ic,params);
  testCase.print();
  container_type x;
for(size_t i=0; i<ic.size(); i++){
  x.push_back(ic[i]);
}
  const double dt = 0.1;
std::vector<container_type> x_vec;
std::vector<double> times;
  // calculate some transients steps
size_t  steps = boost::numeric::odeint::integrate(testCase , x, 0.0 , 100.0 , dt,
                                            push_back_state_and_time( x_vec , times ) );
//size_t  steps = boost::numeric::odeint::integrate(testCase , x, 0.0 , 100.0 , dt);

/* output */
 for( size_t i=0; i<=steps; i++ )
 {
   Rcout << times[i] << '\t' << x_vec[i][0] << '\t' << x_vec[i][1] << '\n';
 }
 Rcout <<"\n"<<"\n"<< x[0]<<'\t'<<x[1] << '\n';
//
 }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#test(gene_interaction)
*/
