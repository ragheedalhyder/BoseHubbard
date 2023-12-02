//
//  main.cpp
//  MISFtransitionBH
//
//  Created by Ragheed Alhyder on 24/10/2022.
//

#include "FunctionsSF.h"

int main(int argc, const char * argv[]) {
    
    string Dir = "/Users/ragheed/Library/Mobile Documents/com~apple~CloudDocs/Aarhus/Projects/MottSuperfluidBosePolaron/MISFtransitionBH/DataPoints/FinalData/Fixed01/";

    
    fstream imported("/Users/ragheed/Library/Mobile Documents/com~apple~CloudDocs/Aarhus/Projects/MottSuperfluidBosePolaron/MISFtransitionBH/DataPoints/NewData/nbar01.txt", std::ios_base::in);
    double a, b;
    int dim = 100;
    double dJUs[dim];
    double muUs[dim];
    int i = 0;
    while ( imported >> a >> b)
    {
        dJUs[i] = a;
        muUs[i] = b;
//        cout << dJUs[i] << '\t' << muUs[i] << endl;
        i++;
    }
    
    
    
    double dkx = 2 * pi / L;
    double dky = 2 * pi / L;
    double KXs[M], KYs[M], dkxs[M], dkys[M];
    
    for (int i = 0; i < M ; i++){
        KXs[i] = -pi + i * dkx;
        KYs[i] = -pi + i * dky;
        dkxs[i] = dkx;
        dkys[i] = dky;
        if (i == 0 || i == M - 1){
            dkxs[i] = dkx / 2.0;
            dkys[i] = dky / 2.0;
        }
        cout << KXs[i] << endl;
    }
    
    double muU =  sqrt(2) - 1;
    double UIB = 0.2;
    double dJU = 0.39; // 2 d J / U with d the dimension
    
    int Max = 100;
    //    for(int count = 0; count < 50 ; count ++){
    //
    //        cout << (double) (count + 1) / Max + 0.168 << endl;
    //    }
    mat n0s = mat(Max, Max, fill::zeros);
    vec sigma0s, sigma1s, sigma2s, sigmapols, T11s, T12s, T21s, T2200s, T22s;
    sigma0s = sigma1s = sigma2s = sigmapols = T11s =  T12s =  T21s = T2200s =  T22s = vec(Max, fill::zeros);
    
    double Epol = 0.14;
    
    
    for (int cutoff = 2; cutoff <= 2 ;cutoff ++)
    {
        std::string s1 = std::to_string(cutoff);
        std::string s2 = std::to_string(M);
        string ext = ".txt";
        string dash = "_";
        
        
        ostringstream oss;
        oss << Dir << "Pert" << ext;
        string var = oss.str();
        ofstream myfile1;
        myfile1.open(var); // Change this too
       
        for(int count = 1; count < 69 ; count ++){

            dJU = (double) (count) / Max ;
            dJU = dJUs[count];
            muU = muUs[count];

            cout << count << "\t dJU = " << dJU << "\t muU = " << muU << endl;
            
            mat Mat = mat(N, N, fill::zeros);
            mat  A  = mat(N, N, fill::zeros);
            mat  B  = mat(N, N, fill::zeros);
            double uks[N][M][M][N], vks[N][M][M][N], omegaklambda[N][M][M];
            double psi0 = 1;
            
            cx_vec Eigvals;
            vec cns;
            cx_mat Eigvecs;
            
            for (int k = 0 ; k < 2; k++){      // Solving the equation of motion self-consistently we find the coefficients cn
                for (int i = 0 ; i < N; i++){
                    for (int j = 0 ; j < N; j++){
                        Mat(i,j) = (0.5 * i * (i - 1) - muU * i) * delta( i , j ) -  dJU * sqrt( i ) * psi0 * delta( i , j + 1) - dJU * sqrt( i + 1 ) * psi0 * delta( i + 1, j );
//                        cout << dJU * sqrt( i + 1 ) * psi0 * delta( i + 1, j ) << "\t" << Mat(i,j) << endl;
                    }
                }
                
                eig_gen(Eigvals, Eigvecs, Mat);
                uvec ind = find( real(Eigvals) == min(real(Eigvals)) );
                cns = abs(real( Eigvecs.col(ind(0)) ));
                
                psi0 = Psi0(cns);
            }
            
            
            double n0 = 0;
            double omega0U = 0;
            
            for (int n = 0 ; n < N; n++){
                n0 = n0 + n * real(cns(n)) * real(cns(n));
                omega0U = omega0U + ( 0.5 * n * (n - 1) - muU * n)  * real(cns(n)) * real(cns(n));
            }
            
            n0s(count) = n0;
            omega0U = - 2 * dJU * psi0 * psi0 + omega0U;
            cout << cns << endl;
            cout << "n0 = " << n0 << endl;
            cout << "Omega0U =" << omega0U << endl;
            
            
            for (int kx = 0; kx < M ; kx ++){
                for (int ky = 0; ky < M ; ky ++){
                    double x = epsI( KXs[kx], KYs[ky] );
                    for (int n = 0; n < N ; n++){
                        for (int m = 0 ; m < N ; m++){
                            A ( n , m ) = ( 0.5 * n * (n - 1) - muU * n - omega0U) * delta( n , m ) -  JkU(dJU, 0) * psi0 * ( sqrt( n ) * delta( n , m + 1) + sqrt( m ) * delta( n + 1 , m ) )- JkU(dJU, x) * ( sqrt( n ) * sqrt( m ) * cn( cns , m - 1 ) * cn( cns , n - 1 ) + sqrt( n + 1 ) * sqrt( m + 1 ) * cn( cns , m + 1 ) * cn( cns , n + 1 ) );
                            B ( n , m ) = - JkU(dJU, x) * ( sqrt( n ) * sqrt( m + 1 ) * cn( cns , m + 1 ) * cn( cns , n - 1 ) + sqrt( n + 1 ) * sqrt( m ) * cn( cns , m - 1 ) * cn( cns , n + 1 ) );
                        }
                    }
                    
                    mat AB = join_rows(A,B);
                    mat BA = join_rows(-B,-A);
                    mat MatAB = join_cols(AB,BA);
                    eig_gen(Eigvals, Eigvecs, MatAB);
                    vec omega00 = sort(real(Eigvals));
                    vec omega0 = omega00.subvec( N , 2 * N - 1);
                    
                    for (int lambda = 0; lambda < N ; lambda ++){
                        omegaklambda[lambda][kx][ky]= omega0(lambda);
                        bool doublezero = false;
                        if (lambda == 0){
                            uvec ind = find(floor(10 * abs(real(Eigvals))) == floor(10 * abs(omega0(lambda))) );
                            
                            if (ind.size() > 2){
                                doublezero = true;
                            }
                            
                            for (int n = 0; n < N; n++){
                                uks[n][kx][ky][lambda] = cns(n);
                                vks[n][kx][ky][lambda] = 0;
                            }
                            omegaklambda[lambda][kx][ky] = 0;
                            //                        cout << cns << endl;
                            
                        }
                        else{
                                
                            uvec ind = find(real(Eigvals) == omega0(lambda) );
                            
                            uword ind1 = ind(0);
                            
                            if (ind.size() > 1){
                                ind1 = ind(1);
                            }
                            double Norm = 0;
                            
                            for (int n = 0; n < N; n++){
                                //                            cout << real(Eigvecs(n, ind1))  << '\t' << real(Eigvecs(N + n, ind1))  << endl;
                                Norm = Norm + pow(  real(Eigvecs(n, ind1)), 2 )  - pow( real(Eigvecs(N + n, ind1)) , 2 );
                            }
                            if (round(10E6 * Norm) <= 0)
                                Norm = 1;
                            for (int n = 0; n < N; n++){
                                uks[n][kx][ky][lambda] = real(Eigvecs(n, ind1)) / sqrt(Norm);
                                vks[n][kx][ky][lambda] = real(Eigvecs(N + n, ind1)) / sqrt(Norm);
//                                if(kx == 5 && ky == 5 && lambda == 1)
//                                    cout << real(Eigvecs(n, ind1)) << "\t" << sqrt(Norm) << endl;
                            }
                        }
                    }
                }
            }

            //     Sigma0 zeroth term of self energy
            
            double deltan2 = 0;
            for (int kx = 0 ; kx < M; kx++){
                for (int ky = 0 ; ky < M; ky++){
                    for (int lambda = 1 ; lambda <= cutoff; lambda++){
                        deltan2 =  deltan2 + 1.0 / pow( 2 * pi, 2 ) * dkxs[kx] * dkys[ky] * V(kx, ky, lambda, kx, ky, lambda, vks, KXs, KYs, n0);
                    }
                }
            }
            
            sigma0s(count) = UIB * (n0 +  deltan2);
            cout << "sig0 = " << sigma0s(count) << endl;

            double Sigma1 = 0;
            double res = 0;
            
            for (int kx = 0 ; kx < M; kx++){
                for (int ky = 0 ; ky < M; ky++){
                    for (int lambda = 1 ; lambda <= cutoff; lambda++){
                        double Nkres = Nk(kx, ky, lambda, cns, uks, vks, KXs, KYs);
                        res = real(pow(UIB / (2 * pi),2) * dkxs[kx] * dkys[ky] * abs( pow( Nkres , 2 )) / ( -  omegaklambda[ lambda][ kx][ ky ]  - dJU * epsI( KXs[kx], KYs[ky] ) + 1i * 0.0001));
                        Sigma1 = Sigma1 + res;
                        if (isnan(res))
                        {
                            cout <<  Nkres;
                        }
                    }
                }
            }
            
            sigma1s (count) = Sigma1;
            cout << "sig1 = " << sigma1s(count) << endl;
            double Sigma2 = 0;
            
            for (int lambda = 1 ; lambda <= cutoff; lambda++){
                for (int lambda1 = 1 ; lambda1 <= cutoff; lambda1++){
                    for (int kx = 0 ; kx < M; kx++){
                        for (int ky = 0 ; ky < M; ky++){

                            for (int px = 0 ; px < M; px++){
                                for (int py = 0 ; py < M; py++){

                                    double W1 = W(kx, ky, lambda, px, py, lambda1, uks, vks, KXs, KYs, n0);
                                    double W2 = W(px, py, lambda1, kx, ky, lambda, uks, vks, KXs, KYs, n0);
                                    Sigma2 = Sigma2 + real(pow(UIB,2) / (2 * pow(2 * pi,4)) * (dkxs[kx] * dkys[ky] * dkxs[px] * dkys[py]) *  abs( pow(W1 + W2 ,2)) / ( -  omegaklambda[lambda][kx][ky] -  omegaklambda[lambda1][px][py] - dJU * epsI( KXs[kx] + KXs[px], KYs[ky] + KYs[py])+ 1i * 0.0001)) ;
                                }
                            }
                        }
                    }
                }
            }
            
            sigma2s (count) = Sigma2;
            
//            myfile << dJU << "\t" << sigma0s(count) + sigma1s (count) + sigma2s (count) << endl;
//            myfile.close();
            cout << "sig2 = " << Sigma2 << endl;
//            cout << cutoff << "\t" << dJU << "\t" << sigma0s(count) << "\t" << sigma1s(count) << "\t" << sigma2s(count) << endl;
            myfile1 << dJU << "\t" << sigma0s(count) << "\t" << sigma1s(count) << "\t" << sigma2s(count) << endl; // Write out perturbative result


            for (Epol = 0.01; Epol <= 0.30; Epol += 0.00005){
                cx_vec sig1vec = SigmaPolaron (Epol , KXs,  KYs,  dkxs, dkys, uks, vks, omegaklambda, dJU, n0, UIB, cutoff);
//                double sig1 = real(sig1vec(0));
//                cout << sig1 << endl;
                //cout << dJU << "\t" << cutoff << "\t" << Epol << "\t" << sig1vec(0) << "\t" << sig1vec(1) << "\t" << sig1vec(2) << "\t" << sig1vec(3) << "\t" << sig1vec(4) << "\t" << sig1vec(5) << "\t" << sig1vec(6) << endl;
//                myfile << dJU << "\t" << Epol << "\t" << real(sig1vec(0)) << endl;
                // return {real(Epol - sigpol), T1100(60,0), T12(60,0), T2100(60,0), T2200(60,0) ,Sigma22, sigpol};
                myfile << dJU << "\t" << Epol << "\t" << real(sig1vec(0)) << "\t" << real(sig1vec(1)) << "\t" << real(sig1vec(2)) << "\t" << real(sig1vec(3)) << "\t" << real(sig1vec(4)) << "\t" << real(sig1vec(5)) << "\t" << real(sig1vec(6)) << "\t" << imag(sig1vec(6)) << endl; // Write out output files for T matrix calculation
            }

            
        }
        myfile1.close();
    }
    return 0;
}
