#include "header.h"
#include <Rcpp.h>


using namespace Rcpp;
void write_parameters(int number_gene, int output_precision, std::vector<double> &g_gene,
                      std::vector<double> &k_gene, std::vector<std::vector<int> > &n_gene,
                      std::vector<std::vector<double> > &lambda_gene,
                      std::vector<std::vector<double> > &threshold_gene_log,
                      std::fstream &out_0)
{
  // out_0<<"Parameters"<<"\n";
  ///////////////////////////////////////////////////////////////////////////////////////

  //Writing parameters to file
  ///////////////////////////////////////////////////////////////////////////////////////
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {out_0<<std::setprecision(output_precision)<<g_gene[gene_count1]<<"\t";} //production rate of each gene

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {out_0<<std::setprecision(output_precision)<<k_gene[gene_count1]<<"\t";} // degradation rate of each gene

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {if(threshold_gene_log[gene_count1][gene_count2]>0)
    out_0<<std::setprecision(output_precision)<<threshold_gene_log[gene_count1][gene_count2]<<"\t";}
  // above--thresholds for the inteaction links, thresholds for inward links for genes are written, starting from gene 1 (thresholds for all inward links of gene 1 and so on)
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {if(n_gene[gene_count1][gene_count2]>0)
    out_0<<std::setprecision(1)<<n_gene[gene_count1][gene_count2]<<"\t";}
  // above--n for the inteaction links, n for inward links for genes are written, starting from gene 1
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {if(lambda_gene[gene_count1][gene_count2]>0)
    out_0<<std::setprecision(output_precision)<<lambda_gene[gene_count1][gene_count2]<<"\t";}
  // above--lambda for the inteaction links, lambda for inward links for genes are written, starting from gene 1
  out_0<<"\n";

}

void read_parameters(IntegerMatrix gene_interaction, int number_gene, std::vector<double> &g_gene,
                     std::vector<double> &k_gene, std::vector<std::vector<int> > &n_gene,
                     std::vector<std::vector<double> > &lambda_gene,
                     std::vector<std::vector<double> > &threshold_gene_log,
                     std::fstream &in_parameters)
{
  ///////////////////////////////////////////////////////////////////////////////////////

  //Reading parameters from file
  ///////////////////////////////////////////////////////////////////////////////////////
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {    in_parameters >>  g_gene[gene_count1];}

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  { in_parameters >>  k_gene[gene_count1]; }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      if(gene_interaction(gene_count1,gene_count2)!=0)
      {
        in_parameters >>  threshold_gene_log[gene_count1][gene_count2];
        //Rcout<<"Here"<< threshold_gene_log[gene_count1][gene_count2]<<"\n";
      }
    }
  }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {double test;
    if(gene_interaction(gene_count1,gene_count2)!=0) { in_parameters >>  test; n_gene[gene_count1][gene_count2] = std::round(test);}}

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
  {if(gene_interaction(gene_count1,gene_count2)!=0) in_parameters >>  lambda_gene[gene_count1][gene_count2];}

}

void select_parameters(IntegerMatrix gene_interaction, NumericVector threshold_gene,
                       double g_min, double g_max,
                       double k_min, double k_max, int possible_interactions,
                       long model_count_max,long threshold_max,
                       double h, double lambda_min,
                       double lambda_max, int n_min, int n_max,
                       double tot_time,
                       double standard_deviation_factor, int number_gene,
                       double D_max,  double D_shot_scaling,
                       int GENE_NOISE_SCALING, int file_writing_interval,
                       int D_levels, double D_scaling,
                       int output_precision, int ANNEALING,
                       int INITIAL_CONDITIONS, String filename, std::vector<double> &g_gene,
                       std::vector<double> &k_gene, std::vector<std::vector<int> > &n_gene,
                       std::vector<std::vector<double> > &lambda_gene,
                       std::vector<std::vector<double> > &threshold_gene_log,
                       std::fstream &out_0)
{
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++){g_gene[gene_count1]=g_min+(g_max-g_min)*u_distribution(u_generator);}

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++){k_gene[gene_count1]=k_min+(k_max-k_min)*u_distribution(u_generator);}

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      if(gene_interaction(gene_count1,gene_count2)==0)
      {
        n_gene[gene_count1][gene_count2]=0;
      }
      else n_gene[gene_count1][gene_count2]=int((n_max-n_min)*u_distribution(u_generator))+n_min;
    }
  }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      if(gene_interaction(gene_count1,gene_count2)==0)
      {
        lambda_gene[gene_count1][gene_count2]=0;
      }
      else lambda_gene[gene_count1][gene_count2]=(lambda_max-lambda_min)*u_distribution(u_generator)+lambda_min;
    }
  }

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      if(gene_interaction(gene_count1,gene_count2)==0)
      {
        threshold_gene_log[gene_count1][gene_count2]=0;
      }
      else
      {

        threshold_gene_log[gene_count1][gene_count2]=(1-standard_deviation_factor*sqrt(3))*threshold_gene[gene_count2]+(2*sqrt(3)*standard_deviation_factor*threshold_gene[gene_count2])*u_distribution(u_generator);
      }
    }
  }


}


void select_IC_range(int number_gene, IntegerMatrix gene_interaction,
                     std::vector<double> &g_gene,
                     std::vector<double> &k_gene,
                     std::vector<std::vector<double> > &lambda_gene,
                     std::vector<double> &max_gene,
                     std::vector<double> &min_gene)
{
  for(int gene_count1=0;gene_count1<number_gene;gene_count1++){max_gene[gene_count1]=g_gene[gene_count1]/k_gene[gene_count1];}

  for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
  {
    double min_gene_multiplier_final=1;

    for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
    {
      double gene_min_multiplier=1;
      double gene_lambda=lambda_gene[gene_count1][gene_count2];
      switch(gene_interaction(gene_count1,gene_count2))
      {
      case 0:
        gene_min_multiplier=1.0;
        break;

      case 2:
        gene_lambda=1./gene_lambda;
        gene_min_multiplier=gene_lambda;
        break;

      case 1:
        gene_min_multiplier=1./gene_lambda;
        break;
      default :
        Rcout << "Invalid Interation code for Gene"<<gene_count1<<" and gene"<<gene_count2<<" interaction"<<"\n";
      }

      min_gene_multiplier_final=min_gene_multiplier_final*gene_min_multiplier;
    }

    min_gene[gene_count1]=g_gene[gene_count1]/k_gene[gene_count1]*min_gene_multiplier_final;

  }

}

void select_IC(int number_gene, IntegerMatrix gene_interaction,
               std::vector<double> &g_gene,
               std::vector<double> &k_gene,
               std::vector<std::vector<int> > &lambda_gene,
               std::vector<double> &max_gene,
               std::vector<double> &min_gene,
               std::fstream &out_ic)
{

}
// [[Rcpp::export]]

int simulate_GRN(IntegerMatrix gene_interaction, NumericVector threshold_gene,
                 double g_min, double g_max,
                 double k_min, double k_max, int possible_interactions,
                 long model_count_max,long threshold_max,
                 double h, double lambda_min,
                 double lambda_max, int n_min, int n_max,
                 double tot_time,
                 double standard_deviation_factor, int number_gene,
                 double D_max,  double D_shot_scaling,
                 int GENE_NOISE_SCALING,
                 int D_levels, double D_scaling,
                 int output_precision, int ANNEALING,
                 int INITIAL_CONDITIONS, String filename, bool parameters_file,
                 bool readIC, bool timeSeries)

{
  try {
    std::string file_name= filename;
    Rcout<<"Running time evolution simulations for "<<std::to_string(number_gene)<<" genes..."<<"\n";
    //  if(INITIAL_CONDITIONS>1)
    int file_writing_interval=1;
    D_max=D_max/sqrt(number_gene); //Scale the maximum noise by the number of genes.
    double D=D_max; // setting noise to maximum noise level
    double Darray[number_gene]; //array to scale the noise level in each gene
    //Rcout<<"parameter_file"<<parameters_file<<"\n";
    // if(parameters_file) Rcout<<"If true"<<"\n";
    // Scale the noise level in each gene if GENE_NOISE_SCALING=1
    if(GENE_NOISE_SCALING==1){

      for(int i=0; i<number_gene;i++)
      {
        Darray[i]=threshold_gene[i];}
      Rcout<<"Using a noise level that is proportional to median expression of the gene"<<"\n";

    }
    else{
      for(int i=0; i<number_gene;i++)
      {
        Darray[i]=1.0;
      }
      Rcout<<"Using same noise level for each gene"<<"\n";

    }

    //Create output file if not there already
    std::fstream out_0;
    if(parameters_file)
    { out_0.open("./tmp/parameters" + file_name +".txt",std::ios::in);}
    else {
  out_0.open("./tmp/parameters" + file_name +".txt",std::ios::out);
    }
    if(!out_0) {     Rcout <<
      "Cannot open output file for writing parameters.\n";  return 1;}


    std::fstream out_ic;
    if(readIC)
    {  out_ic.open("./tmp/ic"+file_name+".txt",std::ios::in);
    }
    else{
      out_ic.open("./tmp/ic"+file_name+".txt",std::ios::out);
    }
    if(!out_ic) {Rcout <<
      "Cannot open output file for writing initial conditions.\n";  return 1;}

    std::fstream out_1("./tmp/GE"+file_name+".txt",std::ios::out);
    if(!out_1) {     Rcout << "Cannot open output file.\n";  return 1;}


    // std::ifstream in_parameters("./tmp/parameters"+file_name+".txt",
    //                             std::ifstream::in);
    // std::ifstream in_ic("./tmp/ic"+file_name+".txt",std::ifstream::in);

    // Rcout<<in_ic.is_open()<<"in_ic open true"<<"\n";

    //File writing is minimized to speed up the program.
    //Output is written to file after every FILE_WRITING_INTERVAL times.
    // Not needed...will be removed later
    int super_count_max=int(model_count_max/file_writing_interval);
    //  Rcout<<"Super Count="<<super_count_max<<"\n";
    model_count_max=file_writing_interval;
    for(int super_count=0;super_count<super_count_max;super_count++)
    {

      //A temporary array to store gene expression values
     double expression_gene_final[file_writing_interval][D_levels][number_gene];
      for(int fwi_count=0;fwi_count<file_writing_interval;fwi_count++)
        for(int fn_count=0;fn_count<D_levels;fn_count++)
          for(int gf=0;gf<number_gene;gf++)
            expression_gene_final[fwi_count][fn_count][gf]=0;

      for(long model_count=0;model_count<model_count_max;model_count++)
      {


        //Initialize production rate of genes
        std::vector<double> g_gene(number_gene);

        //Initialize degradation rate of genes
        std::vector<double> k_gene(number_gene);

        //Initialize hill coefficient for each interaction
        std::vector<std::vector<int> > n_gene(number_gene,
                                              std::vector<int>(number_gene));

        //Initialize fold change for each interaction
        std::vector<std::vector<double> > lambda_gene(
            number_gene, std::vector<double>(number_gene));

        //Initialize threshold for each interaction
        std::vector<std::vector<double> > threshold_gene_log(
            number_gene, std::vector<double>(number_gene));

          if(parameters_file){

            if(!out_0) {Rcout <<
              "Cannot open input file for reading parameters.\n";  return 1;}

            read_parameters( gene_interaction, number_gene, g_gene,
                             k_gene, n_gene,
                             lambda_gene,
                             threshold_gene_log, out_0);

            //    Rcout<<in_ic.is_open()<<"in_ic open true"<<"\n";

          }
          else
          {

            select_parameters( gene_interaction,  threshold_gene,
                               g_min,  g_max,
                               k_min,  k_max,  possible_interactions,
                               model_count_max, threshold_max,
                               h,  lambda_min,
                               lambda_max,  n_min,  n_max,
                               tot_time,
                               standard_deviation_factor,  number_gene,
                               D_max,   D_shot_scaling,
                               GENE_NOISE_SCALING,  file_writing_interval,
                               D_levels,  D_scaling,
                               output_precision,  ANNEALING,
                               INITIAL_CONDITIONS, filename, g_gene,
                               k_gene, n_gene,
                               lambda_gene,
                               threshold_gene_log, out_0);

            write_parameters( number_gene, output_precision, g_gene,
                              k_gene, n_gene,
                              lambda_gene,
                              threshold_gene_log, out_0);

          }


           ///////////////////////////////////////////////////////////////////////////////////////

        //Initial condition range selection
        ///////////////////////////////////////////////////////////////////////////////////////
        std::vector<double> max_gene(number_gene);
        std::vector<double> min_gene(number_gene);
        if(!readIC)
        {
          select_IC_range( number_gene,  gene_interaction, g_gene, k_gene,
                           lambda_gene, max_gene, min_gene);
        }

        ///////////////////////////////////////////////////////////////////////////////////////

        //Initial condition  selection
        ///////////////////////////////////////////////////////////////////////////////////////
        for(int ic_count=0;ic_count<INITIAL_CONDITIONS;ic_count++)
        {
          std::vector <double> expression_gene(number_gene);
          //array for current gene expression
          std::vector <double> expression_gene0(number_gene);
          //array for initial gene expression
          if(readIC)
          {
            if(!out_ic) {Rcout <<
              "Cannot open input file for reading initial conditions.\n";
              return 1;}
            for(size_t ic_counter=0;ic_counter <number_gene; ic_counter++)
            {
              out_ic >> expression_gene0[ic_counter];
            }
          }


          else
          {
            for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
            {
              expression_gene0[gene_count1]=exp(log(min_gene[gene_count1]) +
                (log(max_gene[gene_count1])-
                log(min_gene[gene_count1]))*u_distribution(u_generator));
              expression_gene[gene_count1]=expression_gene0[gene_count1];
            }

          }



          ///////////////////////////////////////////////////////////////////////////////////////

          //Writing initial condition to file
          ///////////////////////////////////////////////////////////////////////////////////////
  if(!readIC){
 //   out_ic.clear();
   // out_ic.seekg(0, std::ios::beg);
          for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
          {
            out_ic<<std::setprecision(output_precision)<<
              expression_gene0[gene_count1]<<"\t";
          }
          out_ic<<"\n";
}


          ///////////////////////////////////////////////////////////////////////////////////////

          //Time Evolution
          ///////////////////////////////////////////////////////////////////////////////////////
          //        for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
          //      {out_ic<<std::setprecision(output_precision)<<expression_gene0[gene_count1]<<"\t";} //initial condition of each gene
          D=D_max; //Start with maximum noise level for each model
          for(int file_count=0;file_count<D_levels;file_count++)
          {

            double expression_gene_h[number_gene]; //array for temp gene expression values
            if(ANNEALING==1) {}
            else {
              for(int gene_count_temp=0;gene_count_temp<number_gene;gene_count_temp++)
              { expression_gene[gene_count_temp]=expression_gene0[gene_count_temp];

                expression_gene_h[gene_count_temp]=expression_gene0[gene_count_temp];
                // out_ic<<std::setprecision(output_precision)<<expression_gene0[gene_count_temp]<<"\t";
              }
            }
            //Rcout<<"D="<<D<<"\n";
            double i=0.0;
            do
            {
              for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
              {
                double final_multiplier=1;

                for(int gene_count2=0;gene_count2<number_gene;gene_count2++)
                {
                  double gene_value=expression_gene[gene_count2];
                  double gene_threshold=threshold_gene_log[gene_count1][gene_count2];
                  int gene_n=n_gene[gene_count1][gene_count2];
                  double gene_lambda=lambda_gene[gene_count1][gene_count2];
                  double gene_activation_multiplier=1;

                  switch(gene_interaction(gene_count1,gene_count2))
                  {
                  case 0:
                    gene_activation_multiplier=1.0;
                    break;

                  case 2:
                    gene_lambda=1./gene_lambda;
                    gene_activation_multiplier=gene_lambda+(1.-gene_lambda)*1./(1.+pow((gene_value/gene_threshold),gene_n));
                    break;

                  case 1:
                    gene_activation_multiplier=(gene_lambda+(1.-gene_lambda)*1./(1.+pow((gene_value/gene_threshold),gene_n)))/gene_lambda;
                    break;

                  default :
                    Rcout << "Invalid Interation code for Gene"<<gene_count1<<" and gene"<<gene_count2<<" interaction"<<"\n";
                  }

                  final_multiplier=final_multiplier*gene_activation_multiplier;
                }


                expression_gene_h[gene_count1]=expression_gene[gene_count1]+h*(g_gene[gene_count1]*final_multiplier-k_gene[gene_count1]*expression_gene[gene_count1])+
                  D*sqrt(h)*g_distribution(g_generator)*Darray[gene_count1]+
                  D_shot_scaling*D*sqrt(h)*g_distribution(g_generator)*Darray[gene_count1]*expression_gene[gene_count1];
                if(expression_gene_h[gene_count1]<0) expression_gene_h[gene_count1]=0;
              }

              for(int gene_count1=0;gene_count1<number_gene;gene_count1++){
                expression_gene[gene_count1]=expression_gene_h[gene_count1];}




              if(timeSeries) {
                out_1<<std::setprecision(output_precision)<<i<<"\t";
                for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
                  {
                  out_1<<std::setprecision(output_precision)<<
                    expression_gene[gene_count1]<<"\t";
                  }
                out_1<<"\n";
                }
              i+=h;
            }while(i<tot_time);
            // Rcout<< "Noise Level" << file_count<<"\t"<<D<<"\n";
            if(file_count==D_levels-2){D=0;}
            else {D=D*D_scaling;}
            // if(file_writing_interval==1)
            if(!timeSeries)
            {
              //Rcout<<"IC";
              for(int gene_count1=0;gene_count1<number_gene;gene_count1++)
              {
                out_1<<std::setprecision(output_precision)<<expression_gene[gene_count1]<<"\t";
              }
              // out_1<<"\n";
            }


          }


          out_1<<"\n";
        }

      }

    }
    out_1.close();
    out_0.close();
    out_ic.close();
//    if(in_ic.is_open()) in_ic.close();
//    if(in_parameters.is_open()) in_parameters.close();

    Rcout<<"Simulations completed successfully. Data files are in results folder.\n";

  }
  catch (Rcpp::internal::InterruptedException& e)
  {
    Rcout << "Caught an interrupt!" << std::endl;
  }
  return 0;
}


