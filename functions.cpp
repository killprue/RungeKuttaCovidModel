#include "functions.hpp"

void rungeKutta1(vector <float> &S, vector <float> &I, vector <float> &R, float Params[4] ){
    float R_0 = Params[0], f = Params[1], gama = Params[2], beta;
    float DS1 = 0, DI1 = 0, DR1 = 0;
    for (int i = 0; i < N-1; i++){
        beta = R_0*gama/S[0];
    
        DS1 = -beta * S[i] * I[i];
        S[i + 1] = S[i] + h * DS1;
        DI1 = beta * S[i] * I[i] - gama*I[i];
        I[i + 1] = I[i] + h*DI1;
        DR1 = f * gama * I[i];
        R[i + 1] = R[i] + h*DR1;
    };
}

void rungeKutta2(vector <float> &S, vector <float> &I, vector <float> &R,float Params[4]) {
    float R_0 = Params[0], f = Params[1], gama = Params[2], beta;
    float DS1 = 0, DI1 = 0, DR1 = 0;
    float DS2 = 0, DI2 = 0, DR2 = 0;
    for (int i = 0; i < N-1; i++){
        beta = R_0*gama/S[0];

        DS1 = -beta * S[i] * I[i];
        S[i + 1] = S[i] + h * DS1/2.0;
        DI1 = beta * S[i] * I[i] - gama*I[i];
        I[i + 1] = I[i] + h*DI1/2.0;
        DR1 = f * gama * I[i];
        R[i + 1] = R[i] + h*DR1/2.0;


        DS2 = -beta*S[i+1]*I[i+1];
        DI2 = beta * S[i+1] * I[i+1] - gama*I[i+1];
        DR2 = f*gama*I[i+1];
        S[i+1] = S[i] + h*DS2;
        I[i+1] = I[i] + h*DI2;
        R[i+1] = R[i] + h*DR2;

    };
}

void rungeKutta4(vector <float> &S, vector <float> &I, vector <float> &R,float Params[4]){
    float R_0 = Params[0], f = Params[1], gama = Params[2], beta;
    float DS1 = 0, DI1 = 0, DR1 = 0;
    float DS2 = 0, DI2 = 0, DR2 = 0;
    float DS3 = 0, DI3 = 0, DR3 = 0;
    float DS4 = 0, DI4 = 0, DR4 = 0;
    for (int i = 0; i < N-1; i++){
        beta = R_0*gama/S[0];

        DS1=-beta*S[i]*I[i];
        S[i+1] = S[i] + (h*DS1/2.0);
        DI1=beta*S[i]*I[i]-(gama*I[i]);
        I[i+1]=I[i]+h*DI1/2.0;
        DR1=f*gama*I[i];
        R[i+1]=R[i] + h*DR1/2.0;

        DS2=-beta*S[i+1]*I[i+1];
        DI2=beta*S[i+1]*I[i+1]-gama*I[i+1];
        DR2=f*gama*I[i+1];
        S[i+1]=S[i]+h*DS2/2.0;
        I[i+1]=I[i]+h*DI2/2.0;
        R[i+1]=R[i]+h*DR2/2.0;

        DS3=-beta*S[i+1]*I[i+1];
        DI3=beta*S[i+1]*I[i+1]-gama*I[i+1];
        DR3=f*gama*I[i+1];
        S[i+1]=S[i]+h*DS3;
        I[i+1]=I[i]+h*DI3;
        R[i+1]=R[i]+h*DR3;

        DS4=-beta*S[i+1]*I[i+1];
        DI4=beta*S[i+1]*I[i+1]-gama*I[i+1];
        DR4=f*gama*I[i+1];
        S[i+1] = S[i] + ((1.0/6.0)*(DS1+(2*DS2)+(2*DS3)+DS4)*h) ;
        I[i+1] = I[i]+(1.0/6.0)*(DI1+2*DI2+2*DI3+DI4)*h;
        R[i+1] = R[i]+(1.0/6.0)*(DR1+2*DR2+2*DR3+DR4)*h;
        
    };
}

void rungeKutta1(vector <float> &S, vector <float> &I, vector <float> &R, float Params[6], vector <float> &E, vector <float> &D, vector<float> &Q) {
    float R_0 = Params[0], A = Params[1], u = Params[2], p = Params[3];
    float e = Params[4], y = Params[5], n = Params[6], beta;
    float DS1 = 0, DE1 = 0, DI1 = 0, DR1 = 0, DQ1 = 0;
    //p_m = percent masked, p_um = percent unmasked
    float p_m = 0.0, p_um = 1-p_m;
    //i_l = 1/(time to lose immunity in days), beta_coeff = reduced transmisibility after maskeing
    float i_l = 0.0055, beta_coeff = 0.21;
    //p_r = positivty rate of tests
    float p_r = 0, tests = 10;
    A = 0;
    u = 0;
    for(int i = 0; i < N-1; i++){
        beta = (1+abs((cos(0.0086*i*h))/2))*(1+abs(cos(0.0349*i*h)/2))*((R_0*((e+y)*(y+p+u)))/e);
        
        DS1 = A - u*S[i] + i_l*R[i] - beta*(beta_coeff*I[i]*p_m*S[i]*p_m + I[i]*p_um*S[i]*p_um + 2*beta_coeff*I[i]*S[i]*p_um*p_m + 2*beta_coeff*I[i]*S[i]*p_m*p_um)/n;
        S[i + 1] = S[i] + h*DS1;
        DQ1 = p_r*tests - y*Q[i];
        Q[i+1] = Q[i] + h*DQ1;
        DE1 = beta*(beta_coeff*I[i]*p_m*S[i]*p_m + I[i]*p_um*S[i]*p_um + 2*beta_coeff*I[i]*S[i]*p_um*p_m + 2*beta_coeff*I[i]*S[i]*p_m*p_um)/n - (u + e)*E[i];
        E[i + 1] = E[i] + h*DE1;
        DI1 = e*E[i] - p_r*tests - (y + u + p) * I[i];
        I[i + 1] = I[i] + h*DI1;
        DR1 = y*I[i] - u*R[i] - i_l*R[i] + y*Q[i];
        R[i + 1] = R[i] + h*DR1;
        
        D[i+1] = n - (S[i+1] + I[i+1] + R[i+1] + E[i+1] + Q[i+1]);
        
        //Arbitrary positivity rate
        p_r = (I[i]/n)*4;
        
    };
}

void rungeKutta2(vector <float> &S, vector <float> &I, vector <float> &R, float Params[6], vector <float> &E, vector <float> &D, vector<float> &Q) {
    float R_0 = Params[0], A = Params[1], u = Params[2], p = Params[3];
    float e = Params[4], y = Params[5], n = Params[6], beta;
    float DS1 = 0, DE1 = 0, DI1 = 0, DR1 = 0, DQ1 = 0;
    float DS2 = 0, DE2 = 0, DI2 = 0, DR2 = 0, DQ2 = 0;
    //p_m = percent masked, p_um = percent unmasked
    float p_m = 0.0, p_um = 1-p_m;
    //i_l = 1/(time to lose immunity in days), beta_coeff = reduced transmisibility after maskeing
    float i_l = 0.0055, beta_coeff = 0.21;
    //p_r = positivty rate of tests
    float p_r = 0, tests = 10;
    A = 0;
    u = 0;
    for(int i = 0; i < N-1; i++){
        beta = (1+abs((cos(0.0086*i*h))/2))*(1+abs(cos(0.0349*i*h)/2))*((R_0*((e+y)*(y+p+u)))/e);
        
        DS1 = A - u*S[i] + i_l*R[i] - beta*(beta_coeff*I[i]*p_m*S[i]*p_m + I[i]*p_um*S[i]*p_um + 2*beta_coeff*I[i]*S[i]*p_um*p_m + 2*beta_coeff*I[i]*S[i]*p_m*p_um)/n;
        S[i + 1] = S[i] + h*DS1/2.0;
        DQ1 = p_r*tests - y*Q[i];
        Q[i+1] = Q[i] + h*DQ1/2.0;
        DE1 = beta*(beta_coeff*I[i]*p_m*S[i]*p_m + I[i]*p_um*S[i]*p_um + 2*beta_coeff*I[i]*S[i]*p_um*p_m + 2*beta_coeff*I[i]*S[i]*p_m*p_um)/n - (u + e)*E[i];
        E[i + 1] = E[i] + h*DE1/2.0;
        DI1 = e*E[i] - p_r*tests - (y + u + p) * I[i];
        I[i + 1] = I[i] + h*DI1/2.0;
        DR1 = y*I[i] - u*R[i] - i_l*R[i] + y*Q[i];
        R[i + 1] = R[i] + h*DR1/2.0;

        DS2 = A - u*S[i+1] + i_l*R[i+1] - beta*(beta_coeff*I[i+1]*p_m*S[i+1]*p_m + I[i+1]*p_um*S[i+1]*p_um + 2*beta_coeff*I[i+1]*S[i+1]*p_um*p_m + 2*beta_coeff*I[i+1]*S[i+1]*p_m*p_um)/n;
        DQ2 = p_r*tests - y*Q[i+1];
        DE2 = beta*(beta_coeff*I[i+1]*p_m*S[i+1]*p_m + I[i+1]*p_um*S[i+1]*p_um + 2*beta_coeff*I[i+1]*S[i+1]*p_um*p_m + 2*beta_coeff*I[i+1]*S[i+1]*p_m*p_um)/n - (u + e)*E[i+1];
        DI2 = e*E[i+1] - p_r*tests - (y + u + p) * I[i+1];
        DR2 = y*I[i+1] - u*R[i+1] - i_l*R[i+1] + y*Q[i+1];
        
        S[i + 1] = S[i] + h * DS2;
        Q[i + 1] = Q[i] + h * DQ2;
        E[i + 1] = E[i] + h * DE2;
        I[i + 1] = I[i] + h * DI2;
        R[i + 1] = R[i] + h * DR2;
        D[i+1] = n - (S[i+1] + I[i+1] + R[i+1] + E[i+1] + Q[i+1]);
        
        //Arbitrary positivity rate
        p_r = (I[i]/n)*4;
        //50 days after the start of the infection. Approximately mid-march.
        //Major lockdowns take effect, arbitrary reduction of R_0
        if(i == 5000){
            R_0 = R_0/6;
        }
        //CDC recommends Masks as an effective way of fighting covid 19
        if(i == 7200){
            p_m = 0.75;
            p_um = 1 - p_m;
            tests = 10000;
        }
        //Approximately May 1st Lockdowns are lifted in much of the country. Guard is let down and mask
        //Compliance slackens
        if(i == 10000){
            R_0 = (R_0 * 6);
            tests = 100000;
        }
        //Beginning of August new wave of lockdowns take effect after resurgence likely
        //Due to July 4th and wave of protests
        if(i == 16000){
            R_0 = (R_0)/6;
            tests = 700000;
        }
        //September 1st, Classes begin for much of the youth across the United States and a weird form
        //Of normalcy takes hold
        if(i == 19000){
            R_0 = R_0 * 6;
        }

        //Between September first and December 1st testing per day tripled later within the month of December vaccines were beginning to be distributed
        if(i == 26000){
            tests = 2100000;
        }
        if(i == 32000){
            tests = 1200000;
        }
        //Delta becomes the dominant form of covid in early july of 2021, but its R_0 is somewhat diminished due to a large portion of the population being vaccinated
        if(i == 53200 && R_0 != 5.08){
            R_0 = 5.08;
        }
    };
}

void rungeKutta4(vector <float> &S, vector <float> &I, vector <float> &R, float Params[6], vector <float> &E, vector <float> &D, vector<float> &Q) {
    float R_0 = Params[0], A = Params[1], u = Params[2], p = Params[3];
    float e = Params[4], y = Params[5], n = Params[6], beta;
    float DS1 = 0, DE1 = 0, DI1 = 0, DR1 = 0, DQ1 = 0;
    float DS2 = 0, DE2 = 0, DI2 = 0, DR2 = 0, DQ2 = 0;
    float DS3 = 0, DE3 = 0, DI3 = 0, DR3 = 0, DQ3 = 0;
    float DS4 = 0, DE4 = 0, DI4 = 0, DR4 = 0, DQ4 = 0;
    //p_m = percent masked, p_um = percent unmasked
    float p_m = 0.0, p_um = 1-p_m;
    //i_l = 1/(time to lose immunity in days), beta_coeff = reduced transmisibility after maskeing
    float i_l = 0.0055, beta_coeff = 0.21;
    //p_r = positivty rate of tests
    float p_r = 0, tests = 10;
    A = 0;
    u = 0;
    for(int i = 0; i < N-1; i++){
        beta = (1+abs((cos(0.0086*i*h))/2))*(1+abs(cos(0.0349*i*h)/2))*((R_0*((e+y)*(y+p+u)))/e);
        
        DS1 = A - u*S[i] + i_l*R[i] - beta*(beta_coeff*I[i]*p_m*S[i]*p_m + I[i]*p_um*S[i]*p_um + 2*beta_coeff*I[i]*S[i]*p_um*p_m + 2*beta_coeff*I[i]*S[i]*p_m*p_um)/n;
        S[i + 1] = S[i] + h*DS1/2.0;
        DQ1 = p_r*tests - y*Q[i];
        Q[i+1] = Q[i] + h*DQ1/2.0;
        DE1 = beta*(beta_coeff*I[i]*p_m*S[i]*p_m + I[i]*p_um*S[i]*p_um + 2*beta_coeff*I[i]*S[i]*p_um*p_m + 2*beta_coeff*I[i]*S[i]*p_m*p_um)/n - (u + e)*E[i];
        E[i + 1] = E[i] + h*DE1/2.0;
        DI1 = e*E[i] - p_r*tests - (y + u + p) * I[i];
        I[i + 1] = I[i] + h*DI1/2.0;
        DR1 = y*I[i] - u*R[i] - i_l*R[i] + y*Q[i];
        R[i + 1] = R[i] + h*DR1/2.0;

        DS2 = A - u*S[i+1] + i_l*R[i+1] - beta*(beta_coeff*I[i+1]*p_m*S[i+1]*p_m + I[i+1]*p_um*S[i+1]*p_um + 2*beta_coeff*I[i+1]*S[i+1]*p_um*p_m + 2*beta_coeff*I[i+1]*S[i+1]*p_m*p_um)/n;
        DQ2 = p_r*tests - y*Q[i+1];
        DE2 = beta*(beta_coeff*I[i+1]*p_m*S[i+1]*p_m + I[i+1]*p_um*S[i+1]*p_um + 2*beta_coeff*I[i+1]*S[i+1]*p_um*p_m + 2*beta_coeff*I[i+1]*S[i+1]*p_m*p_um)/n - (u + e)*E[i+1];
        DI2 = e*E[i+1] - p_r*tests - (y + u + p) * I[i+1];
        DR2 = y*I[i+1] - u*R[i+1] - i_l*R[i+1] + y*Q[i+1];
        
        S[i + 1] = S[i] + h * DS2/2.0;
        Q[i + 1] = Q[i] + h * DQ2/2.0;
        E[i + 1] = E[i] + h * DE2/2.0;
        I[i + 1] = I[i] + h * DI2/2.0;
        R[i + 1] = R[i] + h * DR2/2.0;

        DS3 = A - u*S[i+1] + i_l*R[i+1] - beta*(beta_coeff*I[i+1]*p_m*S[i+1]*p_m + I[i+1]*p_um*S[i+1]*p_um + 2*beta_coeff*I[i+1]*S[i+1]*p_um*p_m + 2*beta_coeff*I[i+1]*S[i+1]*p_m*p_um)/n;
        DQ3 = p_r*tests - y*Q[i+1];
        DE3 = beta*(beta_coeff*I[i+1]*p_m*S[i+1]*p_m + I[i+1]*p_um*S[i+1]*p_um + 2*beta_coeff*I[i+1]*S[i+1]*p_um*p_m + 2*beta_coeff*I[i+1]*S[i+1]*p_m*p_um)/n - (u + e)*E[i+1];
        DI3 = e*E[i+1] - p_r*tests - (y + u + p) * I[i+1];
        DR3 = y*I[i+1] - u*R[i+1] - i_l*R[i+1] + y*Q[i+1];
        
        S[i + 1] = S[i] + h * DS3;
        Q[i + 1] = Q[i] + h * DQ3;
        E[i + 1] = E[i] + h * DE3;
        I[i + 1] = I[i] + h * DI3;
        R[i + 1] = R[i] + h * DR3;

        DS4 = A - u*S[i+1] + i_l*R[i+1] - beta*(beta_coeff*I[i+1]*p_m*S[i+1]*p_m + I[i+1]*p_um*S[i+1]*p_um + 2*beta_coeff*I[i+1]*S[i+1]*p_um*p_m + 2*beta_coeff*I[i+1]*S[i+1]*p_m*p_um)/n;
        DQ4 = p_r*tests - y*Q[i+1];
        DE4 = beta*(beta_coeff*I[i+1]*p_m*S[i+1]*p_m + I[i+1]*p_um*S[i+1]*p_um + 2*beta_coeff*I[i+1]*S[i+1]*p_um*p_m + 2*beta_coeff*I[i+1]*S[i+1]*p_m*p_um)/n - (u + e)*E[i+1];
        DI4 = e*E[i+1] - p_r*tests - (y + u + p) * I[i+1];
        DR4 = y*I[i+1] - u*R[i+1] - i_l*R[i+1] + y*Q[i+1];
        
        S[i+1] = S[i] + ((1.0/6.0)*(DS1+(2*DS2)+(2*DS3)+DS4)*h);
        Q[i+1] = Q[i] + ((1.0/6.0)*(DQ1+(2*DQ2)+(2*DQ3)+DQ4)*h);
        E[i+1] = E[i] + ((1.0/6.0)*(DE1+(2*DE2)+(2*DE3)+DE4)*h);
        I[i+1] = I[i] + ((1.0/6.0)*(DI1+2*DI2+2*DI3+DI4)*h);
        R[i+1] = R[i] + ((1.0/6.0)*(DR1+2*DR2+2*DR3+DR4)*h);
        
        D[i+1] = n - (S[i+1] + I[i+1] + R[i+1] + E[i+1] + Q[i+1]);
        
        //Arbitrary positivity rate
        p_r = (I[i]/n)*4;
        //50 days after the start of the infection. Approximately mid-march.
        //Major lockdowns take effect, arbitrary reduction of R_0
        if(i == 5000){
            R_0 = R_0/6;
        }
        //CDC recommends Masks as an effective way of fighting covid 19
        if(i == 7200){
            p_m = 0.75;
            p_um = 1 - p_m;
            tests = 10000;
        }
        //Approximately May 1st Lockdowns are lifted in much of the country. Guard is let down and mask
        //Compliance slackens
        if(i == 10000){
            R_0 = (R_0 * 6);
            tests = 100000;
        }
        //Beginning of August new wave of lockdowns take effect after resurgence likely
        //Due to July 4th and wave of protests
        if(i == 16000){
            R_0 = (R_0)/6;
            tests = 700000;
        }
        //September 1st, Classes begin for much of the youth across the United States and a weird form
        //Of normalcy takes hold
        if(i == 19000){
            R_0 = R_0 * 6;
        }

        //Between September first and December 1st testing per day tripled later within the month of December vaccines were beginning to be distributed
        if(i == 26000){
            tests = 2100000;
        }
        if(i == 32000){
            tests = 1200000;
        }
        //Delta becomes the dominant form of covid in early july of 2021, but its R_0 is somewhat diminished due to a large portion of the population being vaccinated
        if(i == 53200 && R_0 != 5.08){
            R_0 = 5.08;
        }

    };
}

void resetArrays(vector <float> &S, vector <float> &I, vector <float> &R) {
    float S_0 = S[0], I_0 = I[0];
    fill(S.begin(),S.end(),0);
    fill(I.begin(),I.end(),0);
    fill(R.begin(),R.end(),0);
    S[0] = S_0;
    I[0] = I_0;
}
void resetArrays(vector <float> &S, vector <float> &I, vector <float> &R, vector <float> &E, vector <float> &D, vector <float> &Q) {
    float S_0 = S[0], E_0 = E[0], I_0 = I[0];
    fill(S.begin(),S.end(),0);
    fill(I.begin(),I.end(),0);
    fill(R.begin(),R.end(),0);
    fill(E.begin(),E.end(),0);
    fill(D.begin(),D.end(),0);
    fill(Q.begin(),Q.end(),0);
    S[0] = S_0;
    E[0] = E_0;
    I[0] = I_0;
}

void write(string location,string fileName,vector <float> S,vector <float> I,vector <float> R){
    vector <float> X(N,0);
    for (int i = 0; i < N; i++){ X[i] = i+1; };
    ofstream output;
    output.open(location + fileName);
    output << "X," << "Infected," << "Recovered," << "Susceptible," << endl;
    for (int i = 0; i < N; i++){
        output << X[i]*h << "," << I[i] << "," << R[i] << "," << S[i] <<  endl;
    }
    output.close();
}
void write(string location,string fileName,vector <float> S,vector <float> I,vector <float> R, vector <float> E, vector <float> &D, vector <float> &Q){
    float X[N] = {};
    for (int i = 0; i < N; i++){ X[i] = i+1; };
    ofstream output;
    output.open(location + fileName);
    output << "X," << "Infected," << "Recovered," << "Susceptible," << "Exposed," << "Dead,"<< "Quarantined"<<endl;
    for (int i = 0; i < N; i++){
        output << X[i]*h << "," << I[i] << "," << R[i] << "," << S[i] << "," << E[i] << "," << D[i] << "," << Q[i] << endl;
    }
    output.close();
}

void generateSERData(string fileList[numberOfFiles], float Params[4]) {

    vector <float> S(N,0), I(N,0), R(N,0), D(N,0);

    
    S[0] = Params[3] - 1;
    I[0] = 1;

    rungeKutta1(S, I, R, Params);
    write(saveLocation,fileList[0],S,I,R);
    resetArrays(S, I, R);


    rungeKutta2(S, I, R, Params);
    write(saveLocation,fileList[1],S,I,R);
    resetArrays(S, I, R);

    rungeKutta4(S, I, R, Params);
    write(saveLocation,fileList[2],S,I,R);

}

void generateSEIRData(string fileList[numberOfFiles], float Params[7]) {
    vector <float> S(N,0), I(N,0), R(N,0), E(N,0), D(N,0), Q(N,0);

    S[0] = Params[6] - (Params[6]*0.0000002) - 1;
    E[0] = Params[6]*0.0000002;
    I[0] = 1;


    rungeKutta1(S, I, R, Params, E,D,Q);
    write(saveLocation,fileList[0],S,I,R,E,D,Q);
    resetArrays(S, I, R, E,D,Q);

    rungeKutta2(S, I, R, Params, E,D,Q);
    write(saveLocation,fileList[1],S ,I,R,E,D,Q);
    resetArrays(S, I, R, E,D,Q);

    rungeKutta4(S, I, R, Params, E,D, Q);
    write(saveLocation,fileList[2],S,I,R,E,D, Q);
}


