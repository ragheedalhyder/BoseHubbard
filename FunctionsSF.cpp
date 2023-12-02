#include "FunctionsSF.h"

//int L = 9; // Number of lattice sites in the X and Y directions (Lx = Ly = L)
//int M = L + 1;

int delta(int x, int y) {
    if (x == y)
        return 1;
    else
        return 0;
}
///
///
double cn(vec cns, int n){
    if(n == - 1)
        return 0;
    if(n == N)
        return 0;
    else
        return real(cns(n));
        
}
///
///
double Psi0(vec cns){
    double Psi0 = 0;
    for(long i = 1; i < N; i++){
        Psi0 = Psi0 + sqrt(i) * cns( i - 1 ) * cns( i );
    }
    return Psi0;
}
///
///
double JkU(double dJU, double x){
    return dJU - dJU * x;
}
///
///
double epsI(double kx, double ky){
    return pow(sin( kx / 2 ),2) + pow(sin( ky / 2 ),2);
}
///
///
double U(int kx, int ky, int lambda, int qx, int qy, int lambda1, double uks[N][M][M][N], double KXs[M], double KYs[M], double n0 ){
    double Ures = 0;
//    if(KXs[kx] == 0 && KYs[ky] == 0 && lambda == 1)
//        return 0;
//    if(KXs[qx] == 0 && KYs[qy] == 0 && lambda1 == 1)
//        return 0;
    for(long i = 0; i < N; i++){
        Ures = Ures + ( i - n0 * (1 - delta(lambda, lambda1) * delta(lambda, 0))) * uks[i][kx][ky][lambda] * uks[i][qx][qy][lambda1];
    }
    return Ures;
}
///
///
double V(int kx, int ky, int lambda, int qx, int qy, int lambda1, double vks[N][M][M][N], double KXs[M], double KYs[M], double n0 ){
    double Vres = 0;
//    if(KXs[kx] == 0 && KYs[ky] == 0 && lambda == 1)
//        return 0;
//    if(KXs[qx] == 0 && KYs[qy] == 0 && lambda1 == 1)
//        return 0;
    for(long i = 0; i < N; i++){
        Vres = Vres + ( i - n0 ) * vks[i][kx][ky][lambda] * vks[i][qx][qy][lambda1];
    }
    return Vres;
}
///
///
double W(int kx, int ky, int lambda, int qx, int qy, int lambda1, double uks[N][M][M][N], double vks[N][M][M][N], double KXs[M], double KYs[M], double n0 ){
    double Wres = 0;
//    if(KXs[kx] == 0 && KYs[ky] == 0 && lambda == 1)
//        return 0;
//    if(KXs[qx] == 0 && KYs[qy] == 0 && lambda1 == 1)
//        return 0;
    for(long i = 0; i < N; i++){
        Wres = Wres + ( i - n0 ) * uks[i][kx][ky][lambda] * vks[i][qx][qy][lambda1];
//
//        int indk = lambda * M * M + kx * M + ky;
//        int indq = lambda1 * M * M + qx * M + qy;
//        if(indk == 303 && indq == 303)
//            cout << uks[i][kx][ky][lambda] << "\t" << vks[i][qx][qy][lambda1] << endl;
    }
    return Wres;
}
///
///
double Nk(int kx, int ky, int lambda, vec cns, double uks[N][M][M][N], double vks[N][M][M][N], double KXs[M], double KYs[M]){
    double Nkres = 0;
//    if(KXs[kx] == 0 && KYs[ky] == 0 && lambda == 1)
//        return 0;
//    cout << kx << "\t" << ky << "\t" << lambda << "\t" << Nkres << "\t";
    for(long i = 0; i < N; i++){
        Nkres = Nkres + i * cns(i) * ( uks[i][kx][ky][lambda] + vks[i][kx][ky][lambda] );
//        cout << cns(i) << "\t" << uks[i][kx][ky][lambda] << "\t" << vks[i][kx][ky][lambda] << "\t";
    }
//    cout << "\t" << Nkres << endl;
    return Nkres;
}
///
///


cx_vec SigmaPolaron (double Epol, double KXs[M], double KYs[M], double dkxs[M], double dkys[M], double uks[N][M][M][N], double vks[N][M][M][N], double omegaklambda[N][M][M], double dJU, double n0, double UIB, int cutoff){
//    cutoff = 5;
    int dim =  M * M * (cutoff + 1) ;
//    double eta = 0.00000001;
    double eta = 0.003;
    eta = 0.0001;
    mat Us = mat(dim, dim, fill::zeros);
    mat Vs = mat(dim, dim, fill::zeros);
    mat Ws = mat(dim, dim, fill::zeros);
    
    cx_mat PI1100, PI1121, PI21, PI12, PI2200, PI22, inv11, inv22, inv12;
    
    cx_mat T1100, T1200, T1121, T1222, T2100, T2200, T12, T22 ;
    
    PI1100 = PI1121 = PI21 = PI12 = PI2200 = PI22 = inv11 =inv12 = inv22 = T1100 = T1200 = T1121 = T1222 = T2100 = T2200 = T12 = T22 = cx_mat(dim, dim, fill::zeros);
    //cx_vec T22diag = cx_vec(dim, fill::zeros);
    
    int kx, ky, qx, qy, indk, indq, indk0, indq0 ;
    indk0 = 0;
    indq0 = 0;
    double Uelem, Velem, Welem, coeff;
    complex<double> Den11, Den22, Den;
//    int ZeroInd = 181; // Zero momentum in the Goldstone mode diverges so we remove it by hand.
    
    double epsplus = 0;
    double epsminus = 0;
    bool ZeroIndx1 = true;
    bool ZeroIndx2 = true;
    for (int lambda = 0; lambda <= cutoff; lambda++){
        for (int lambda1 = 0 ; lambda1 <= cutoff; lambda1++){
            for (int kx1 = 0 ; kx1 < M; kx1++){
                for (int ky1 = 0 ; ky1 < M; ky1++){
                    for (int qx1 = 0 ; qx1 < M; qx1++){
                        for (int qy1 = 0 ; qy1 < M; qy1++){
                            
                            kx = kx1;
                            ky = ky1;
                            qx = qx1;
                            qy = qy1;
                            
                            indk = lambda * M * M + kx * M + ky;
                            
                            indq = lambda1 * M * M + qx * M + qy;

                            if (lambda == 0 && KXs[kx] == 0 && KYs[ky] == 0 && ZeroIndx1 == true){
                                indk0 = indk;
                                ZeroIndx1 = false;
                                //cout << "indk0 = " << indk0 << endl;
                            }
                            
                            if (lambda1 == 0 && KXs[qx] == 0 && KYs[qy] == 0 && ZeroIndx1 == true){
                                indq0 = indq;
                                ZeroIndx2 = false;

                            }

                            
//                            if (lambda == 0)
//                                kx = ky = 5;

                            Uelem = UIB * U(kx, ky, lambda, qx, qy, lambda1, uks, KXs, KYs, n0);

                            Velem = UIB * V(kx, ky, lambda, qx, qy, lambda1, vks, KXs, KYs, n0);

                            Welem = UIB * ( W(kx, ky, lambda, qx, qy, lambda1, uks, vks, KXs, KYs, n0) + W(qx, qy, lambda1, kx, ky, lambda, uks, vks, KXs, KYs, n0));
                            
                            epsplus = epsI(  KXs[kx] + KXs[qx],  KXs[ky] + KXs[qy]);
                            epsminus = epsI(  KXs[kx] - KXs[qx],  KXs[ky] - KXs[qy]);
                            
                            if (lambda1 == 0)
                            {
                                epsplus = epsI(  KXs[kx], KXs[ky]);
                                epsminus = epsI(  KXs[kx], KXs[ky]);
                            }
                            if (lambda == 0)
                            {
                                epsplus = epsI(  KXs[qx],  KXs[qy]);
                                epsminus = epsI( - KXs[qx],  - KXs[qy]);
                            }
                            
                            Den11 = Epol - omegaklambda[lambda1][qx][qy] - dJU * epsI(  KXs[qx],  KXs[qy] ) + eta * 1i;
                            
                            Den22 = Epol - omegaklambda[lambda1][qx][qy] - dJU * epsI(  KXs[qx],  KXs[qy] ) + eta * 1i;
                            
                            Den   = Epol - omegaklambda[lambda][kx][ky] - omegaklambda[lambda1][qx][qy] - dJU * epsplus + eta * 1i;
//                            if (indk == indq){
//                                cout << real(Den11) << endl;
////                                cout << "1." << real(Den11) << endl;
////                                cout << "2." << real(Den22) << endl;
////                                cout << "3." << real(Den) << endl;
//                            }
                            Us(indk,indq) =  Uelem;
                            Vs(indk,indq) =  Velem;
                            Ws(indk,indq) =  Welem;
                            

                            //cout << dkxs[kx] << '\t' << dkys[ky] << endl;
                            coeff = 1 / (4 * pi * pi) * dkxs[qx] * dkys[qy];
                            
                            PI1100(indk,indq) = coeff * (Uelem / Den11); // Pair propagator for T1100 and T1200
                            
                            PI1121(indk,indq) = coeff *  (Uelem / Den22); // Pair propagator for T11 to calculate T21
                            
                            PI21(indk,indq) = coeff * (Welem / Den22);   // Pair propagator to calculate T21 from T11
                            
                            PI2200(indk,indq) = coeff * (Welem / Den22);
                            
                            
                            PI12(indk,indq) = coeff * (Uelem / Den);     // Pair propagator for T12 to calculate T22
                            
                            PI22(indk,indq) = coeff * (Welem / Den);     // Pair propagator for T22 as a function of T12
                            
                            
//                            if(indk == 180 ){
//                                if (real(Den) == 0.00540844)
//                                    cout << omegaklambda[lambda][kx][ky] << "\t" << omegaklambda[lambda1][qx][qy] << "\t" <<  dJU * epsI( KXs[kx] +  KXs[qx], KXs[ky] + KXs[qy]) ;
////                                cout << lambda << kx << ky << endl;
////                                cout << PI12(indk,indq) << endl;
//                                cout << "Den = " << Den << endl;
//                            }
                            
//                            cout << PI11(indk,indq) << "\t" << PI22(indk,indq) << "\t" << omegaklambda[lambda][kx][ky] << "\t" << omegaklambda[lambda1][qx][qy] << endl;
                            if (lambda1 == 0){
                                PI1100(indk,indq) = 0;
                                PI1121(indk,indq) = 0;
                                PI21(indk,indq) = 0;
                                PI12(indk,indq) = 0;
                                PI2200(indk,indq) = 0;
                                PI22(indk,indq) = 0;
                            }
                            
                        }
                    }
                }
            }
        }
    }
    
    mat IMat = mat(dim,dim,fill::eye);
    
    inv11 = inv ( IMat  - PI1100 );
    inv22 = inv ( IMat  - PI1121 );
//    inv11.print();
//    cout << inv22.max() << endl;
    
    T1100 = inv11 * Us;
    T1200 = inv11 * Ws;
    
    T1121 = inv22 * Us;
    T2100 = Ws + PI21 * T1121;
    
    T1222 = inv22 * Ws;
    T2200 = Vs + PI21 * T1222;
    
    inv12 =  inv(IMat  - PI12);
    T12 = inv12 * Ws;
    T22 = Vs + 0.5 * PI22 * T12;
    
    cx_vec T22diag = diagvec(T22) ;
//    cout << T22diag << endl;
    complex<double> Sigma22 = 0;
    
    for (int lambda = 1 ; lambda <= cutoff; lambda++){
        for (int kx = 0 ; kx < M; kx++){
            for (int ky = 0 ; ky < M; ky++){
                indk = lambda * M * M + kx * M + ky;
//                if (indk == 170)
//                    cout << lambda;
                Sigma22 = Sigma22 + 1.0 / (4 * pi * pi) * dkxs[kx] * dkxs[ky] * T22diag(indk);
//                cout << T22diag(indk) << endl;
            }
        }
    }
//    cout << T1100(60,0) << endl;
//    cout << real(T11(0,0)) << "\t" << real(T12(0,0)) << "\t" << real(T21(0,0)) << "\t" << real(T2200(0,0)) << "\t" << Sigma22 << endl;
    complex<double> sigpol = T1100(indk0,indq0) + T1200(indk0,indq0) + T2100(indk0,indq0) + T2200(indk0,indq0) + Sigma22 ;
//    sigpol = Sigma22 ;
//    return Sigma22;
    return {real(Epol - sigpol), T1100(indk0,indq0), T12(indk0,indq0), T2100(indk0,indq0), T2200(indk0,indq0) ,Sigma22, sigpol};
}


cx_vec FixedPoint(double sig1, double sig2, double Epol, double dEpol, double KXs[M], double KYs[M],  double dkxs[M], double dkys[M], double uks[N][M][M][N], double vks[N][M][M][N], double omegaklambda[N][M][M], double dJU, double n0, double UIB, int cutoff){
    cout << "FP \t " << Epol << endl;
    cx_vec sigmidvec = SigmaPolaron(Epol + dEpol / 2, KXs, KYs, dkxs, dkys, uks, vks, omegaklambda, dJU, n0, UIB, cutoff);

    if (abs(dEpol) < epsilon){
        cout << "Epol = "<< Epol + dEpol / 2.0 << endl;
        return {Epol + dEpol / 2.0, sigmidvec(1), sigmidvec(2), sigmidvec(3), sigmidvec(4) ,sigmidvec(5)};
    }

    double sigmid = real(sigmidvec(0));
    
    if (sig1 * sig2 > 0){
        cout << " No Trimer Energy Values in the defined interval" << endl;
        return {Epol + dEpol / 2.0, sigmidvec(1), sigmidvec(2), sigmidvec(3), sigmidvec(4) ,sigmidvec(5)};
    }
    if (sig1 * sigmid < 0 ){
        return FixedPoint(sig1, sigmid, Epol, dEpol / 2.0, KXs,  KYs,  dkxs, dkys, uks, vks, omegaklambda, dJU, n0, UIB, cutoff);
    }
    else{
        return FixedPoint(sigmid, sig2, Epol + dEpol / 2.0, dEpol / 2.0, KXs,  KYs,  dkxs, dkys, uks, vks, omegaklambda, dJU, n0, UIB, cutoff);
    }
}


cx_vec LimitFinder(double sig1, double Epol, double dEpol, double KXs[M], double KYs[M],  double dkxs[M], double dkys[M], double uks[N][M][M][N], double vks[N][M][M][N], double omegaklambda[N][M][M], double dJU, double n0, double UIB, int cutoff){
    cout << Epol << "\t" << sig1 << endl;
    double sig2;//, sig3, sig4, sig5, sig6= 0.;
    
    cx_vec sig2vec = SigmaPolaron ( Epol + dEpol , KXs,  KYs,  dkxs, dkys, uks, vks, omegaklambda, dJU, n0, UIB, cutoff);

    sig2 = real(sig2vec(0));

    if (sig1 > sig2){
        return FixedPoint(sig1, sig2, Epol, dEpol, KXs,  KYs,  dkxs, dkys, uks, vks, omegaklambda, dJU, n0, UIB, cutoff);;
    }

    if (sig1 * sig2 < 0 ){
        return FixedPoint(sig1, sig2, Epol, dEpol, KXs,  KYs,  dkxs, dkys, uks, vks, omegaklambda, dJU, n0, UIB, cutoff);;
    }
    else if (sig1 * sig2 > 0 ){
        return LimitFinder(sig2, Epol + dEpol, dEpol, KXs,  KYs,  dkxs, dkys, uks, vks, omegaklambda, dJU, n0, UIB, cutoff);
    }
    return 0;
}
