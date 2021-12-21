#include "functions.hpp"
using namespace std;

int main() {
    
    //SIR MODEL FILE NAMES FOR OUTPUT
    string sirClassFiles[numberOfFiles] = {"Class_SIR_RK1.csv",
                                           "Class_SIR_RK2.csv",
                                           "Class_SIR_RK4.csv"};
    string sirCovidFiles[numberOfFiles] = {"Covid_SIR_RK1.csv",
                                           "Covid_SIR_RK2.csv",
                                           "Covid_SIR_RK4.csv"};
    string sirDeltaFiles[numberOfFiles] = {"Delta_SIR_RK1.csv",
                                           "Delta_SIR_RK2.csv",
                                           "Delta_SIR_RK4.csv"};
    //SEIR MODEL FILE NAMES FOR OUTPUT
    string seirClassFiles[numberOfFiles] = {"Class_SEIR_RK1.csv",
                                            "Class_SEIR_RK2.csv",
                                            "Class_SEIR_RK4.csv"};
    string seirCovidFiles[numberOfFiles] = {"Covid_SEIR_RK1.csv",
                                            "Covid_SEIR_RK2.csv",
                                            "Covid_SEIR_RK4.csv"};
    string seirDeltaFiles[numberOfFiles] = {"Delta_SEIR_RK1.csv",
                                            "Delta_SEIR_RK2.csv",
                                            "Delta_SEIR_RK4.csv"};

    //SERS CLASS/COVID/DELTA MODEL
    
    //R_0=Reproductive rate of the virus (Secondary infections per infection)
    //f=fraction of those who recover from the disease
    //gama=recovery rate from infected (reciprocal is infectious period)
    
    
    //CLASS PARAMETERS
    float sers_class_R_0 = 1.5;
    float sers_class_f = 1;
    float sers_class_gama = 0.2;
    float sers_class_population = 100000;
    float sersClassParams[4] = {sers_class_R_0,
                                sers_class_f,
                                sers_class_gama,
                                sers_class_population};
    generateSERData(sirClassFiles, sersClassParams);

    //COVID PARAMETERS
    float sers_covid_R_0 = 2.2; //New England Journal of Mecidine
    float sers_covid_f = 0.9963; // National Library of medicine
    float sers_covid_gama = (1.0/10.0); //CDC JAMA
    float sers_covid_population = 320000000.0;//US Population
    float sersCovidParams[4] = {sers_covid_R_0,
                                sers_covid_f,
                                sers_covid_gama,
                                sers_covid_population};
    generateSERData(sirCovidFiles, sersCovidParams);

    //DELTA PARAMETERS
    float sers_delta_R_0 = 5.08; //Journal of Travel Medicine
    float sers_delta_f = 0.9963; //Unknown at this time
    float sers_delta_gama = (1.0/10.0); //Unknown at this time
    float sers_delta_population = 320000000.0;
    float sersDeltaParams[4] = {sers_delta_R_0,
                                sers_delta_f,
                                sers_delta_gama,
                                sers_delta_population};
    generateSERData(sirDeltaFiles,sersDeltaParams);

    
    
    //SEIRS COVID/DELTA MODEL
    
    
    //R_0=Reproductive rate of the virus (Secondary infections per infection)
    //A=Per capita birth rate
    //u=Per capita death rate
    //p=virus induced fatality rate
    //e=rate of progression from exposed to infected (reciprocal is the incubation period)
    //y=recovery rate from infected (reciprocal is infectious period)
    //n = Population
    
    
    //CONVERGENCE PARAMETERS
    float seirs_convergence_R_0 = 5.72;
    float seirs_convergence_A = (10267.0/320000000.0);
    float seirs_convergence_u = (8500.0/320000000.0);
    float seirs_convergence_p = 0.006;
    float seirs_convergence_e = (1.0/3.0);
    float seirs_convergence_y = (1.0/8.0);
    float seirs_convergence_n = 10000000;
    float seirConvergenceParams[7] = {seirs_convergence_R_0,
                           seirs_convergence_A,
                           seirs_convergence_u,
                           seirs_convergence_p,
                           seirs_convergence_e,
                           seirs_convergence_y,
                           seirs_convergence_n};
    generateSEIRData(seirClassFiles,seirConvergenceParams);
    
    //COVID PARAMETERS
    float seirs_covid_R_0 = 2.2; //New England Journal of Medicine
    float seirs_covid_A = (10267.0/320000000.0); //Per capita birth rate per day (US) CDC
    float seirs_covid_u = (8500.0/320000000.0);//Per capita death rate per day (US) CDC
    float seirs_covid_p = 0.0037/10; //MedRxiv
    float seirs_covid_e = (1.0/5.1); //Annals of Internal Medicine
    float seirs_covid_y = (1.0/10.0); // CDC Jama
    float seirs_covid_n = 320000000.0;
    float seirCovidParams[7] = {seirs_covid_R_0,
                           seirs_covid_A,
                           seirs_covid_u,
                           seirs_covid_p,
                           seirs_covid_e,
                           seirs_covid_y,
                           seirs_covid_n};
    generateSEIRData(seirCovidFiles,seirCovidParams);

    //DELTA PARAMETERS
    float seirs_delta_R_0 = 5.08; //Journal of Travel Medicine
    float seirs_delta_A = 0.001;//Per capita birth rate per day (US) CDC
    float seirs_delta_u = 0.00001;//Per capita death rate per day (US) CDC
    float seirs_delta_p = 0.0037/10;//MedRxiv
    float seirs_delta_e = (1.0/4.0);//CCDC
    float seirs_delta_y = (1.0/10.0);
    float seirs_delta_n = 320000000.0;
    float seirDeltaParams[7] = {seirs_delta_R_0,
                                seirs_delta_A,
                                seirs_delta_u,
                                seirs_delta_p,
                                seirs_delta_e,
                                seirs_delta_y,
                                seirs_delta_n};
    generateSEIRData(seirDeltaFiles, seirDeltaParams);

}
