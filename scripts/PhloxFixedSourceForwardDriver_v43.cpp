// ANTO-authored diagnostic driver for original 2014 PhloxSelection forward
// C++ code ONLY. No copied original author source or full MCMC.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/random.hpp>
#include "World.h"
#include "likelihood.h"
#include "Params.h"

// Required external RNG symbol used by original code. The fixed forward
// simulation must not sample it; deterministic replay checked in v43 CI.
boost::mt19937 gen;

static std::vector<double> original_start() {
    std::ifstream fin("seedSet.params");
    if (!fin.is_open()) throw std::runtime_error("ORIGINAL_SEEDSET_FILE_MISSING");
    std::vector<double> x;
    double v;
    while (fin >> v) x.push_back(v);
    const double expected[] = {
      .279, .362, .612, .639, .362, 0, .72, .88, .91, .90, 3.7, 6.3, 5.5
    };
    if (x.size()!=13) throw std::runtime_error("ORIGINAL_SEEDSET_DIMENSION_DRIFT");
    for (int i=0;i<13;++i) {
        if (!std::isfinite(x[i]) || std::fabs(x[i]-expected[i])>1e-14)
            throw std::runtime_error("ORIGINAL_SEEDSET_VALUES_DRIFT");
    }
    return x;
}

static void good_prob_vector(const std::vector<double>& a, const char* which,
                             double& max_row_sum_error) {
    if (a.size()!=4) throw std::runtime_error(std::string("MODEL_BAD_4CLASS_DIMENSION_")+which);
    double sum=0;
    for (double x:a) {
        if (!std::isfinite(x) || x<-1e-10 || x>1.0+1e-10)
            throw std::runtime_error(std::string("MODEL_INVALID_PROBABILITY_")+which);
        sum+=x;
    }
    double delta=std::fabs(sum-1);
    max_row_sum_error=std::max(max_row_sum_error,delta);
    if (delta>1e-6)
        throw std::runtime_error(std::string("MODEL_NOT_NORMALIZED_")+which);
}

int main() {
 try {
    std::cout<<std::setprecision(17);
    Params params;  // original author inPhlox_nuissance.params parser
    if (params.nWorlds!=5 || params.wName.size()!=5 || params.totDemes!=36 ||
        params.nZones!=2 || params.migDis!=4 || params.maxGens!=100000 ||
        std::fabs(params.maxChange-1e-8)>1e-17 ||
        std::fabs(params.startPop-1.5)>1e-10 ||
        std::fabs(params.deltaDeme-.5)>1e-10) {
        throw std::runtime_error("ORIGINAL_NUISANCE_PARAMETERS_DRIFT");
    }
    for (int i=0;i<5;i++) {
        if(params.wName[i]!=i+1)
            throw std::runtime_error("ORIGINAL_TRANSECT_ID_DRIFT");
    }
    std::vector<double> start=original_start();
    params.addUntransformedPars(start);
    double grand_like=0, grand_a=0,grand_f=0,max_delta=0;
    int grand_pop=0,grand_moms=0,grand_adult=0,grand_children=0;
    std::cout<<"ANTO_V43_DRIVER_MODEL_SCOPE ORIGINAL_FORWARD_ONLY_NOT_MCMC"<<std::endl;
    for (int i=0;i<params.nWorlds;i++) {
        World world(&params,i);
        world.simulate(&params,i);
        const std::vector<DemeOut> predicted=world.MakesimDemeOut();
        std::vector<ObsData> observed=readObsFiles(params.wName[i]);
        if (predicted.size()!=36)
            throw std::runtime_error("ORIGINAL_DEME_COUNT_DRIFT");
        double l_source=calcTotLike(const_cast<std::vector<DemeOut>&>(predicted),observed,&params);
        double l_adult=0, l_family=0;
        double sum_row_error=0;
        double adults_pred[4]={0,0,0,0}, kids_pred[4]={0,0,0,0};
        int n_adult=0,n_family=0,n_child=0;
        for (unsigned j=0;j<observed.size();j++) {
            ObsData& ob=observed[j];
            DemeOut interpolated=calcBar(const_cast<std::vector<DemeOut>&>(predicted),
                                          ob.position,&params);
            good_prob_vector(interpolated.simPheno,"ADULT",sum_row_error);
            if (ob.PhenoCounts.size()!=4)
                throw std::runtime_error("SOURCE_POP_PHENOTYPE_DIMENSION_DRIFT");
            int count=0;
            for (int k=0;k<4;k++) {
                count+=ob.PhenoCounts[k];
                if (ob.PhenoCounts[k]<0)
                    throw std::runtime_error("SOURCE_NEGATIVE_ADULT_COUNT");
            }
            if (count<=0)
                throw std::runtime_error("SOURCE_EMPTY_ADULT_POPULATION");
            l_adult+=MakeLikelihood(interpolated.simPheno,ob.PhenoCounts);
            for (int k=0;k<4;k++)
                adults_pred[k]+=count*interpolated.simPheno[k];
            n_adult+=count;
            for (auto& fam:ob.momOff) {
                if (fam.momPheno<0 || fam.momPheno>=4 || fam.offCount.size()!=4)
                    throw std::runtime_error("MALFORMED_SOURCE_MOTHER_FAMILY_ROW");
                int children=0;
                for (int k: fam.offCount) {
                    if (k<0)throw std::runtime_error("NEGATIVE_CHILD_COUNT");
                    children+=k;
                }
                if (children<=0)
                    throw std::runtime_error("SOURCE_ALL_ZERO_OR_PADDED_MOTHER_ROW");
                std::vector<double> f_expected=interpolated.ktable[fam.momPheno];
                good_prob_vector(f_expected,"MOTHER_CONDITIONAL_CHILD",sum_row_error);
                l_family+=MakeLikelihood(f_expected,fam.offCount);
                for(int k=0;k<4;k++)
                    kids_pred[k]+=children*f_expected[k];
                n_family++; n_child+=children;
            }
        }
        if (std::fabs(l_source-(l_adult+l_family))>1e-6)
            throw std::runtime_error("SOURCE_ORIGINAL_LIKELIHOOD_COMPONENTS_DISAGREE");
        if (!std::isfinite(l_source) || l_source>0)
            throw std::runtime_error("ORIGINAL_SOURCE_NONFINITE_OR_POSITIVE_LOGLIKE");
        int generations=(int)world.output_Gens();
        if (generations<=0 || generations>params.maxGens)
            throw std::runtime_error("MODEL_GENERATION_COUNT_BAD");
        if (generations==params.maxGens)
            throw std::runtime_error("FIXED_SOURCE_WORLD_HIT_MAX_GENERATIONS_NOT_CONFIRMED_EQUILIBRIUM");
        if (observed.empty() || n_family==0 || n_adult==0 || n_child==0)
            throw std::runtime_error("EMPTY_ORIGINAL_TRANSECT_DATA");
        std::cout<<"ANTO_V43_ORIGINAL_SOURCE_TRANSECT "
                 <<params.wName[i]<<" "<<observed.size()
                 <<" "<<n_adult<<" "<<n_family<<" "<<n_child
                 <<" "<<generations<<" "<<l_source<<" "<<l_adult
                 <<" "<<l_family<<" "<<sum_row_error;
        for(int k=0;k<4;k++)std::cout<<" "<<adults_pred[k]/n_adult;
        for(int k=0;k<4;k++)std::cout<<" "<<kids_pred[k]/n_child;
        std::cout<<std::endl;
        grand_like+=l_source;grand_a+=l_adult;grand_f+=l_family;
        grand_pop+=observed.size();grand_adult+=n_adult;
        grand_moms+=n_family;grand_children+=n_child;
        max_delta=std::max(max_delta,sum_row_error);
    }
    if (grand_pop!=32 || grand_adult!=11379 ||
        grand_moms!=371 || grand_children!=2496)
        throw std::runtime_error("ORIGINAL_AUTHOR_EXAMPLE_DENOMINATOR_DRIFT");
    if (std::fabs(grand_like-grand_a-grand_f)>1e-6)
        throw std::runtime_error("ORIGINAL_MULTITRANSECT_LIKELIHOOD_SUM_DRIFT");
    std::cout<<"ANTO_V43_ORIGINAL_SOURCE_TOTAL "
             <<grand_pop<<" "<<grand_adult<<" "<<grand_moms<<" "<<grand_children
             <<" "<<grand_like<<" "<<grand_a<<" "<<grand_f
             <<" "<<max_delta<<std::endl;
    return 0;
 } catch (const std::exception& exc) {
    std::cerr<<"ANTO_V43_ORIGINAL_FORWARD_TECHNICAL_FAILURE "<<exc.what()<<std::endl;
    return 3;
 }
}
