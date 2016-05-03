//
//  Population Synthesis of Planet Formation
//
//  gimpact.c
//
//
//  Created by Takanori Sasaki
//  2016/05/03
//
//

#include "header.h"

void gimpact(double e[], double a[], double Mp[], double peri[], double apo[], double *M, double *t, int *n)
{
    int i, j, k, l, tmp_i, kk, ll;
    double e_h, A, B, rH, Tk[100][100], mu, e0, aa, bb, theta, theta_i, theta_j, omega_i, omega_j, dt_system;
    double tau_short, e_esc, rho, s, E_0, E, Sig_a, tmp, ran, prob, lambda, epsi, lambda2, epsi2;
    int em_i, em_j, e_enc_in[100], e_enc_out[100], em_max_in, em_max_out, em_max_in2, em_max_out2, g1, g2, gtmp;
    double a_i0, a_j0, Wj, fj, Cj, Min, Mout, a_tmp_in[100], a_tmp_out[100], tau_cross[100][100];

    //int nnn = 15;
    rho = 3.0;
    Cj = 2.0/3.0;

    //printf ("%f ", a[nnn]);
    //printf ("%f ", *t);
    //printf ("\n");

    for (i=0; i<*n; i++) {
    if (Mp[i] > 1e-2*ME) {
        e_h = pow(Mp[i]/(3.0*(*M)), 1.0/3.0);
        e[i] = e_h*sqrt(-1.0*log((double)rand()/((double)RAND_MAX+1)));
        peri[i] = (1.0 - e[i])*a[i];
        apo[i] = (1.0 + e[i])*a[i];
    }
    }

    do {

        tau_short = 1.0e9;
    for (i=0; i<(*n)-1; i++) {
        if (Mp[i] > 1e-2*ME) {
        for (j=i+1; j<(*n); j++) {
            if (Mp[j] > 1e-2*ME) {
            aa = sqrt(a[i]*a[j]);
            bb = fabs(a[i]-a[j]);
            mu = 0.5*(Mp[i]+Mp[j])/(*M);
            e0 = 0.5*(e[i]+e[j])*aa/bb;
            A = -2.0 + e0 - 0.27*log10(mu);
            B = 18.7 + 1.1*log10(mu) - (16.8+1.2*log10(mu))*e0;
            Tk[i][j] = (sqrt(4.0*M_PI*M_PI*aa*AU*aa*AU*aa*AU/(G*(*M))))/YEAR;
            if (a[i] < a[j]) rH = pow((Mp[i]+Mp[j])/(3.0*(*M)), 1.0/3.0)*a[i];
            else rH = pow((Mp[i]+Mp[j])/(3.0*(*M)), 1.0/3.0)*a[j];
            tau_cross[i][j] = pow(10.0, A + B*log10(bb/(2.3*rH)))*Tk[i][j];
            if (tau_cross[i][j] < tau_short) {
                tau_short = tau_cross[i][j];
                em_i = i;
                em_j = j;

            }
        }
    }
        }
    }
        if (tau_short > pow(10.0, 5.5)*Tk[em_i][em_j]) dt_system = tau_short;
        else dt_system = pow(10.0, 5.5)*Tk[em_i][em_j] + (pow(10.0, 6.5)-pow(10.0, 5.5))*((double)rand()/((double)RAND_MAX+1))*Tk[em_i][em_j];
        //printf ("%f %f\n", *t-2.0e7, dt_system);
        //printf ("%d %d %f %f\n", em_i, em_j, e[em_i], e[em_j]);
        //printf ("%f %f %f %f %f %e %e\n", tau_short, dt_system, Tk[em_i][em_j], a[em_i], a[em_j], Mp[em_i]/ME, Mp[em_j]/ME);
        //printf("%f %e %f %f\n", sqrt(a[em_i]*a[em_j]), 0.5*(Mp[em_i]+Mp[em_j])/(*M), 0.5*(e[em_i]+e[em_j])*sqrt(a[em_i]*a[em_j])/fabs(a[em_i]-a[em_j]), pow((Mp[em_i]+Mp[em_j])/(3.0*(*M)), 1.0/3.0)*sqrt(a[em_i]*a[em_j]));

        //printf ("%e %e\n", Mp[em_i]/ME, Mp[em_j]/ME);

    if (*t + dt_system > 1.0e9) return;

        if (a[em_i] > a[em_j]) {
            tmp_i = em_i;
            em_i = em_j;
            em_j = em_i;
        }

    e_esc = 0.28*pow((Mp[em_i]+Mp[em_j])/ME, 1.0/3.0)*pow(rho/3.0, 1.0/6.0)*pow(sqrt(a[em_i]*a[em_j]), 0.5);
    e[em_i] = Mp[em_j]*e_esc*sqrt(-1.0*log((double)rand()/((double)RAND_MAX+1)))/(Mp[em_i]+Mp[em_j]);
    e[em_j] = Mp[em_i]*e_esc*sqrt(-1.0*log((double)rand()/((double)RAND_MAX+1)))/(Mp[em_i]+Mp[em_j]);

    a_i0 = a[em_i];
    a_j0 = a[em_j];

        a_tmp_in[0] = a_i0;
        a_tmp_out[0] = a_j0;
        a[em_i] -= e[em_i]*a[em_i];
        a[em_j] += e[em_j]*a[em_j];

    peri[em_i] = (1.0 - e[em_i])*a[em_i];
    apo[em_i] = (1.0 + e[em_i])*a[em_i];
    peri[em_j] = (1.0 - e[em_j])*a[em_j];
    apo[em_j] = (1.0 + e[em_j])*a[em_j];

    k = 0;
    l = 0;
    for (i=0; i<*n; i++) {
        if (Mp[i] > 1.0e-2*ME) {
        if (i!=em_i && i!=em_j) {


                if ((apo[em_j]>peri[i]) && (peri[em_j]<apo[i])) {
                    e_enc_out[k] = i;
                    k++;
                }
                if ((apo[em_i]>peri[i]) && (peri[em_i]<apo[i])) {
                    e_enc_in[l] = i;
                    l++;
                }
                e_enc_out[k] = em_j;
                e_enc_in[l] = em_i;
                    }
    }
    }


    //printf("%d %d\n", em_i, em_j);

    if (k!=0) {
    em_max_out = k;
    for (i=0; i<k; i++) {
        if (Mp[e_enc_out[i]] > Mp[e_enc_out[em_max_out]]) em_max_out = i;
    }
    if (k == em_max_out) {
        em_max_out2 = 0;
        for (i=1; i<k; i++) {
            if (Mp[e_enc_out[i]] > Mp[e_enc_out[em_max_out2]]) em_max_out2 = i;
        }
    }
    else {
        em_max_out2 = k;
        for (i=0; i<k; i++) {
            if (i!=em_max_out) {
                if (Mp[e_enc_out[i]] > Mp[e_enc_out[em_max_out2]]) em_max_out2 = i;
            }
        }
    }

    for (i=0; i<k; i++) {
        if (i!=em_max_out) {
            e_esc = 0.28*pow((Mp[e_enc_out[i]]+Mp[e_enc_out[em_max_out]])/ME, 1.0/3.0)*pow(rho/3.0, 1.0/6.0)*pow(sqrt(a[e_enc_out[i]]*a[e_enc_out[em_max_out]]), 0.5);
            e[e_enc_out[i]] = Mp[e_enc_out[em_max_out]]*e_esc*sqrt(-1.0*log((double)rand()/((double)RAND_MAX+1)))/(Mp[e_enc_out[i]]+Mp[e_enc_out[em_max_out]]);
        }
        else {
            e_esc = 0.28*pow((Mp[e_enc_out[i]]+Mp[e_enc_out[em_max_out2]])/ME, 1.0/3.0)*pow(rho/3.0, 1.0/6.0)*pow(sqrt(a[e_enc_out[i]]*a[e_enc_out[em_max_out2]]), 0.5);
            e[e_enc_out[i]] = Mp[e_enc_out[em_max_out2]]*e_esc*sqrt(-1.0*log((double)rand()/((double)RAND_MAX+1)))/(Mp[e_enc_out[i]]+Mp[e_enc_out[em_max_out2]]);
        }
    }


    for (i=0; i<k+1; i++) {
        Min = 0.0;
        Mout = 0.0;
        for (j=0; j<k+1; j++) {
            if (i!=j) {
            if (a[e_enc_out[i]] < a[e_enc_out[j]]) Mout += Mp[e_enc_out[j]];
            else Min += Mp[e_enc_out[j]];
            }
        }
        fj = Cj*(Min - Mout)/(Min + Mout);
        Wj = fj + (1.0 - fabs(fj))*(1.0 - 2.0*(double)rand()/((double)RAND_MAX+1));
        a_tmp_out[i] = a[e_enc_out[i]];
        a[e_enc_out[i]] += Wj*e[e_enc_out[i]]*a[e_enc_out[i]];
        //printf ("%f\n", Wj*e[e_enc_out[i]]*a[e_enc_out[i]]);

    }

    /*
    E_0 = 0.0;
    E = 0.0;
    for (i=0; i<k+1; i++) {
        E_0 += 1.0/a_tmp_out[i];
        E += 1.0/a[e_enc_out[i]];
    }
    s = E/E_0;

     for (i=0; i<k+1; i++) {
        a[e_enc_out[i]] *= s;
        peri[e_enc_out[i]] = (1.0 - e[e_enc_out[i]])*a[e_enc_out[i]];
        apo[e_enc_out[i]] = (1.0 + e[e_enc_out[i]])*a[e_enc_out[i]];
    }*/
    }

    if (l!=0) {
    em_max_in = l;
    for (i=0; i<l; i++) {
        if (Mp[e_enc_in[i]] > Mp[e_enc_in[em_max_in]]) em_max_in = i;
    }
    if (l == em_max_in) {
        em_max_in2 = 0;
        for (i=1; i<l; i++) {
                if (Mp[e_enc_in[i]] > Mp[e_enc_in[em_max_in2]]) em_max_in2 = i;
        }
    }
    else {
        em_max_in2 = l;
        for (i=0; i<l; i++) {
            if (i!=em_max_in) {
                if (Mp[e_enc_in[i]] > Mp[e_enc_in[em_max_in2]]) em_max_in2 = i;
            }
        }
    }

    for (i=0; i<l; i++) {
        if (i!=em_max_in) {
            e_esc = 0.28*pow((Mp[e_enc_in[i]]+Mp[e_enc_in[em_max_in]])/ME, 1.0/3.0)*pow(rho/3.0, 1.0/6.0)*pow(sqrt(a[e_enc_in[i]]*a[e_enc_in[em_max_in]]), 0.5);
            e[e_enc_in[i]] = Mp[e_enc_in[em_max_in]]*e_esc*sqrt(-1.0*log((double)rand()/((double)RAND_MAX+1)))/(Mp[e_enc_in[i]]+Mp[e_enc_in[em_max_in]]);
        }
        else {
            e_esc = 0.28*pow((Mp[e_enc_in[i]]+Mp[em_max_in2])/ME, 1.0/3.0)*pow(rho/3.0, 1.0/6.0)*pow(sqrt(a[e_enc_in[i]]*a[em_max_in2]), 0.5);
            e[e_enc_in[i]] = Mp[em_max_in2]*e_esc*sqrt(-1.0*log((double)rand()/((double)RAND_MAX+1)))/(Mp[e_enc_in[i]]+Mp[em_max_in2]);
        }
    }

        for (i=0; i<l+1; i++) {
            Min = 0.0;
            Mout = 0.0;
            for (j=0; j<l+1; j++) {
                if (i!=j) {
                    if (a[e_enc_in[i]] < a[e_enc_in[j]]) Mout += Mp[e_enc_in[j]];
                    else Min += Mp[e_enc_in[j]];
                }
            }
            fj = Cj*(Min - Mout)/(Min + Mout);
            Wj = fj + (1.0 - fabs(fj))*(1.0 - 2.0*(double)rand()/((double)RAND_MAX+1));
            a_tmp_in[i] = a[e_enc_in[i]];
            a[e_enc_in[i]] += Wj*e[e_enc_in[i]]*a[e_enc_in[i]];
            //printf ("%f\n", Wj*e[e_enc_in[i]]*a[e_enc_in[i]]);
        }

/*
    E_0 = 0.0;
    E = 0.0;
    for (i=0; i<l+1; i++) {
        E_0 += 1.0/a_tmp_in[i];
        E += 1.0/a[e_enc_in[i]];
    }
    s = E/E_0;

    for (i=0; i<l+1; i++) {
        a[e_enc_in[i]] *= s;
        peri[e_enc_in[i]] = (1.0 - e[e_enc_in[i]])*a[e_enc_in[i]];
        apo[e_enc_in[i]] = (1.0 + e[e_enc_in[i]])*a[e_enc_in[i]];
    }*/
    }

        kk = k;
        ll = l;
        //printf ("%d %d\n", em_i, em_j);
        do {
    k = 0;
    l = 0;
    for (i=0; i<(*n)+1; i++) {
        if (Mp[i] > 1.0e-2*ME) {

        if (i!=em_i) {
            if ((apo[em_i]>peri[i]) && (peri[em_i]<apo[i])) {
                e_enc_in[l] = i;
                l++;
            }
        }

        if (i!=em_j) {
            if ((apo[em_j]>peri[i]) && (peri[em_j]<apo[i])) {
                e_enc_out[k] = i;
                k++;
            }
        }
        }
    }


        if (kk!=0) {
    g1 = em_j;
    tmp = 0.0;
    for (i=0; i<k; i++) {
        tmp += pow(a[e_enc_out[i]], -3.0);
    }
    tmp = 1.0/tmp;
    ran = (double)rand()/((double)RAND_MAX+1);
    prob = 0.0;
    for (i=0; i<k; i++) {
        if (ran > prob) g2 = e_enc_out[i];
        prob += tmp*pow(a[e_enc_out[i]], -3.0);
    }

    if (a[g1] > a[g2]) {
        gtmp = g1;
        g1 = g2;
        g2 = gtmp;
    }


        //printf("%d %d\n", g1, g2);


            theta = 0.0;
            do {
                theta += 0.01*M_PI;
                lambda = a[g2]*pow(1.0-e[g2], 2.0)/(1.0+e[g2]);
                epsi = 2.0*sqrt(e[g2])/(1.0+e[g2]);
                lambda2 = a[g1]*pow(1.0-e[g1], 2.0)/(1.0+e[g1]);
                epsi2 = 2.0*sqrt(e[g1])/(1.0+e[g1]);
                theta_i = lambda*sin(theta)/(1.0+epsi*cos(theta));
                theta_j = lambda2*sin(M_PI-theta)/(1.0-epsi2);
            } while (theta_i < theta_j);

            theta *= 3.0;
            theta = M_PI - theta;
            theta = (M_PI - theta) + 2.0*theta*(double)rand()/((double)RAND_MAX+1);

    omega_i = 2.0*M_PI*(double)rand()/((double)RAND_MAX+1);
    omega_j = omega_i + theta;

            //printf ("%f %f\n", omega_i, omega_j);

    e[g1] = sqrt((pow(Mp[g2]*e[g2]*cos(omega_i)+Mp[g1]*e[g1]*cos(omega_j), 2.0) + pow(Mp[g2]*e[g2]*sin(omega_i)+Mp[g1]*e[g1]*sin(omega_j), 2.0))/pow(Mp[g2]+Mp[g1], 2.0));
    e[g2] = 0.0;

    a[g1] = (Mp[g1] + Mp[g2])/((Mp[g1]/a[g1])+(Mp[g2]/a[g2]));
    peri[g1] = (1.0 - e[g1])*a[g1];
    apo[g1] = (1.0 + e[g1])*a[g1];
    a[g2] = 50.0 + 50.0*(double)rand()/((double)RAND_MAX+1);
    peri[g2] = (1.0 - e[g2])*a[g2];
    apo[g2] = (1.0 + e[g2])*a[g2];

    Mp[g1] = Mp[g1] + Mp[g2];
    Mp[g2] = 1.0e20;
        }

        if (ll!=0) {
    g1 = em_i;
    tmp = 0.0;
    for (i=0; i<l; i++) {
        tmp += pow(a[e_enc_in[i]], -3.0);
    }
    tmp = 1.0/tmp;
    ran = (double)rand()/((double)RAND_MAX+1);
    prob = 0.0;
    for (i=0; i<l; i++) {
        if (ran > prob) g2 = e_enc_in[i];
        prob += tmp*pow(a[e_enc_in[i]], -3.0);
    }

    if (a[g1] > a[g2]) {
        gtmp = g1;
        g1 = g2;
        g2 = gtmp;
    }
        //printf("%d %d\n", g1, g2);

                theta = 0.0;
            do {
                theta += 0.01*M_PI;
                lambda = a[g2]*pow(1.0-e[g2], 2.0)/(1.0+e[g2]);
                epsi = 2.0*sqrt(e[g2])/(1.0+e[g2]);
                lambda2 = a[g1]*pow(1.0-e[g1], 2.0)/(1.0+e[g1]);
                epsi2 = 2.0*sqrt(e[g1])/(1.0+e[g1]);
                theta_i = lambda*sin(theta)/(1.0+epsi*cos(theta));
                theta_j = lambda2*sin(M_PI-theta)/(1.0-epsi2);
            } while (theta_i < theta_j);

            theta *= 3.0;
            theta = M_PI - theta;
            theta = (M_PI - theta) + 2.0*theta*(double)rand()/((double)RAND_MAX+1);


    omega_i = 2.0*M_PI*(double)rand()/((double)RAND_MAX+1);
    omega_j = omega_i + theta;

    e[g1] = sqrt((pow(Mp[g2]*e[g2]*cos(omega_i)+Mp[g1]*e[g1]*cos(omega_j), 2.0) + pow(Mp[g2]*e[g2]*sin(omega_i)+Mp[g1]*e[g1]*sin(omega_j), 2.0))/pow(Mp[g2]+Mp[g1], 2.0));
    e[g2] = 0.0;

    a[g1] = (Mp[g1] + Mp[g2])/((Mp[g1]/a[g1])+(Mp[g2]/a[g2]));
    peri[g1] = (1.0 - e[g1])*a[g1];
    apo[g1] = (1.0 + e[g1])*a[g1];

            Mp[g2] = 1.0e-20*ME;
    /*
    a[g2] = a[g2] = 50.0 + 50.0*(double)rand()/((double)RAND_MAX+1);
    peri[g2] = (1.0 - e[g2])*a[g2];
    apo[g2] = (1.0 + e[g2])*a[g2];
      */
        //    printf ("Impact!\n");

    Mp[g1] = Mp[g1] + Mp[g2];
            Mp[g2] = 1.0e20;
        }
            //printf("%d %d\n", k, l);
        }while (k!=0 && l!=0);

        //printf("a\n");
        //for (i=0; i<*n; i++) printf("%f %f\n", a[i], Mp[i]/ME);

        *t += dt_system;
            //if (a[nnn] < 3.0) {
                //printf ("%f ", a[nnn]);
                //printf ("%e\n", *t);
            //}

    //printf ("\n");
    /*
        printf("%f ", *t);
        for (i=0; i<(*n); i++) printf("%f ", a[i]);
        printf("\n");
*/
        //printf ("%f %f %f %f %f %f %f %f\n", *t-2.0e7, a[0], a[5], a[7], a[9], a[10], a[16], a[17]);

    } while (*t < 1.0e9);

}
