#include"header.h"
using namespace Rcpp;
size_t convertAdjMatToVector(
    IntegerMatrix gene_interaction, std::vector<size_t>& tgtGene,
    std::vector<std::pair<size_t,size_t> >& intSrcType);


// [[Rcpp::export]]
void generateThresholds(IntegerMatrix gene_interaction,
                                 NumericVector threshold_gene,
                                 Rcpp::List simulationData)
{
  NumericVector simulationParameters = as<NumericVector>(simulationData["simParams"]);
  NumericVector stochasticParameters = as<NumericVector>(simulationData["stochParams"]);
  NumericVector hyperParameters = as<NumericVector>(simulationData["hyperParams"]);

    size_t model_count_max = static_cast<size_t>(simulationParameters[0]);
    double parameter_range = simulationParameters[1];

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

//Rcout<<g_min<<"\t"<<g_max<<"\t"<<possible_interactions<<"\t"<<n_max<<"\n";
//possible_interactions = 5;
//Rcout<<hyperParameters[6]<<"\t"<<possible_interactions<<"\n";
//return 0;
  size_t number_gene=gene_interaction.nrow();

  // Adjust the parameters minimum and maximum values according to
  // parameter range
  applyParameterRange(parameter_range, g_min, g_max, k_min, k_max,
                      lambda_min, lambda_max);

// First step : Find the median expression of isolated gene
    double gene_isolated_median=(g_min+g_max)/(k_min+k_max);
  for(size_t i=0;i<number_gene;i++){threshold_gene[i]= gene_isolated_median;}

  // Use the median expressions of isolated genes to incorporate the effects
  //of all incoming interactions on a gene

  //an array of interaction type and their number
  size_t interactions_number[number_gene][possible_interactions];
  for(size_t i=0;i<number_gene;i++)
    for(size_t j=0;j<possible_interactions;j++)
      interactions_number[i][j]=0;
  for(size_t it1=0;it1<number_gene;it1++)
    {
    for(size_t it2=0;it2<number_gene;it2++)
      {
      size_t g_1=gene_interaction(it1,it2);
      interactions_number[it1][g_1]+=1;
      }
    }

  const size_t model_count_max2=std::max(model_count_max,threshold_max);
  for(size_t gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    std::vector<double> Af;
    for(size_t model_count=0;model_count<model_count_max2;model_count++)
    {
      double g_a = g_min+(g_max-g_min)*u_distribution(u_generator);
      double k_a = k_min+(k_max-k_min)*u_distribution(u_generator);
      Af.push_back(g_a/k_a);
    // Effect of excitatory interactions
      for(size_t counter_interaction1=0;
          counter_interaction1<interactions_number[gene_count1][1];
          counter_interaction1++)
      {
        double g_b=g_min+(g_max-g_min)*u_distribution(u_generator);
        double k_b=k_min+(k_max-k_min)*u_distribution(u_generator);
        double BA0=gene_isolated_median*(1-standard_deviation_factor*sqrt(3)) +
          (2*sqrt(3)*standard_deviation_factor*gene_isolated_median)*
          u_distribution(u_generator);

        size_t n_ba =
          static_cast<size_t>(u_distribution(u_generator)*(n_max-n_min))+n_min;
        double gene_lambda=(lambda_max-lambda_min)*u_distribution(u_generator)
          + lambda_min;

        Af[model_count] = Af[model_count]*Hs_Racipe(g_b/k_b, BA0, n_ba,
                                                    gene_lambda)/gene_lambda;
      }

// Effect of inhibitory interactions
      for(size_t counter_interaction1=0;counter_interaction1 <
        interactions_number[gene_count1][2];counter_interaction1++)
      {
        double g_b=g_min+(g_max-g_min)*u_distribution(u_generator);
        double k_b=k_min+(k_max-k_min)*u_distribution(u_generator);
        double BA0=gene_isolated_median*(1-standard_deviation_factor*sqrt(3)) +
          (2*sqrt(3)*standard_deviation_factor*gene_isolated_median)*
          u_distribution(u_generator);

        size_t n_ba =
          static_cast<size_t>(u_distribution(u_generator)*(n_max-n_min))+n_min;
      double gene_lambda=1./((lambda_max-lambda_min)*(n_max-n_min)+lambda_min);
        double hill_eval=Hs_Racipe(g_b/k_b, BA0, n_ba, gene_lambda);
        Af[model_count] = Af[model_count]*hill_eval;
      }
    }
// Find median
    std::sort(Af.begin(), Af.end());
    threshold_gene[gene_count1] =
      0.5*(Af[static_cast<size_t>((model_count_max2)*0.5)] +
      Af[static_cast<size_t>((model_count_max2+2)*0.5)]);
    Af.clear();

  }
return;
}
