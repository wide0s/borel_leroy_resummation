void series_z_zz(double zz[], double z[], int n)
{
    double ETA[3];
    double R[3];

    double ETA_amplitude, R_amplitude, amplitude;

    // Helpers:
    // n^2, n^3, n^4, (n+2), (n+8), 1/(n+8), (1/(n+8))^2
    double dn = (double) n;
    double dnn = dn * dn;
    double dnnn = dnn * dn;
    double dnnnn = dnnn * dn;
    double S = dn + 2.0;
    double L = dn + 8.0;
    double reverse_l = 1.0 / L;
    double LL = reverse_l * reverse_l;

    // Eta = Amplitude * (eta0 + eta1 * e + eta2 * e^2)
    ETA_amplitude = 0.5 * S * LL;
    ETA[0] = 1.0;
    ETA[1] = 0.25 * LL * (56.0 * dn + 272.0 - dnn);
    ETA[2] = 0.0625 * LL * LL * 
             (1124.0 * dnn + 17920.0 * dn + 46144.0 -
              5.0 * dnnnn - 230.0 * dnnn -
              461.5898508 * L * (5.0 * dn + 22.0));

//    printf ("ETA: %.10f * (%.10f + %.10f*e + %.10f*e^2)\n", ETA_amplitude, ETA[0], ETA[1], ETA[2]);

    // R = Amplitude * (r0 + r1 * e + r2 * e^2)
    R_amplitude = 6.0 * log (4.0 / 3.0) - 1.0;
    R[0] = 1.0;
    R[1] = -0.188483;
    R[2] = -0.0999529 + LL * (4.772294768 * dn + 21.49187274);

//    printf ("R: %.10f * (%.10f + %.10f*e + %.10f*e^2)\n", R_amplitude, R[0], R[1], R[2]);

    // z = 2.0 + Amplitude * e^2 * zz(e), where zz
    // zz = zz0 + zz1 * e + zz2 * e^2
    zz[0] = ETA[0] * R[0];
    zz[1] = ETA[1] + R[1];
    zz[2] = ETA[2] + ETA[1] * R[1] + R[2];

//    printf ("zz: %.10f + %.10f*e + %.10f*e^2\n", zz[0], zz[1], zz[2]);

    amplitude = ETA_amplitude * R_amplitude;

    z[0] = 2.0;
    z[1] = 0.0;
    z[2] = amplitude;
    z[3] = amplitude * zz[1];
    z[4] = amplitude * zz[2];

//    printf ("z: %.10f + %.10f*e^1 + %.10f*e^2 + %.10f*e^3 + %.10f*e^4\n", z[0], z[1], z[2], z[3], z[4]);

//    zzz[0] = z[2];
//    zzz[1] = z[3];
//    zzz[2] = z[4];
}
