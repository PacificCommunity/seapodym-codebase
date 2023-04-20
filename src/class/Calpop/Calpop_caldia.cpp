#include "calpop.h"

/////////////////////////////////////////////////////////////////////////////////
//   -------------------------------------------------------------
//     sous programme 'caldia' : calcul des coefficients diagonaux
//     -----------------------------------------------------------
//     de t=t a t=t+1/2   (a - b - c)
//     de t=t+1/2 a t=t+1 (d - e - f)

//     voir Fournier et Sibert (1994) pour le calcul détaillé

//	si type_diff=0 la diffusion est constante pour source et forage
//	si type_diff=1 la diffusion est fonction de l'indice d'habitat pour les
//thons 	l'argument sp permet de definir l'espece :ex: 0 skj, 1, yft...
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// AVERTISSEMENT !! Les fonction caldia et calrec ont ete developpees sous
// fortran et utilisent l'ordre d'incrementation i pour colonnes (longitudes) et
// j pour lignes (latitudes). Les matrices sont donc lues dans le sens fortran
// et ecrites dans dans le sens C
////////////////////////////////////////////////////////////////////////////////

void CCalpop::caldia(
    const PMap& map, const CParam& param, const DMATRIX& diffusion_x,
    const DMATRIX& advection_x, const DMATRIX& diffusion_y,
    const DMATRIX& advection_y) {
    if (map.global == 1) {
        caldia_GO(
            map, param, diffusion_x, advection_x, diffusion_y, advection_y);
        return;
    }

    a.initialize();
    b.initialize();
    c.initialize();
    d.initialize();
    e.initialize();
    f.initialize();

    CBord bord;
    int dt = 0;  // 2*iterationNumber;
    double dx = param.deltaX;
    double dxdx = dx * dx;
    double twodxdx = 2 * dxdx;

    for (int j = map.jmin; j <= map.jmax; j++) {
        for (int i = map.iinf[j]; i <= map.isup[j]; i++) {
            bord.b = map.bord_cell[i][j];
            char pos = bord.cotex();
            double sigma = diffusion_x[i][j];
            double sigmam = diffusion_x[i - 1][j];
            double sigmap = diffusion_x[i + 1][j];
            double uinf = advection_x[i - 1][j];
            double u = advection_x[i][j];
            double usup = advection_x[i + 1][j];

            a[j][i] = d1(pos, sigmam, sigma, uinf, twodxdx, dxdx, dx);
            b[j][i] = d2(pos, sigmam, sigma, sigmap, u, twodxdx, dxdx, dx, dt);
            c[j][i] = d3(pos, sigma, sigmap, usup, twodxdx, dxdx, dx);
        }
    }
    double dy = param.deltaY;
    double dydy = dy * dy;
    double twodydy = 2 * dydy;

    for (int i = map.imin; i <= map.imax; i++) {
        for (int j = map.jinf[i]; j <= map.jsup[i]; j++) {
            bord.b = map.bord_cell[i][j];
            char pos = bord.cotey();
            double sigma = diffusion_y[i][j];
            double sigmam = diffusion_y[i][j - 1];
            double sigmap = diffusion_y[i][j + 1];
            double vinf = advection_y[i][j - 1];
            double v = advection_y[i][j];
            double vsup = advection_y[i][j + 1];

            d[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
            e[i][j] = d2(pos, sigmam, sigma, sigmap, v, twodydy, dydy, dy, dt);
            f[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
        }
    }
    ybet_comp(map);
    // cout << norm(a) << " " <<  norm(c) << " " << norm(d) << " " << norm(e) <<
    // " " << norm(f) <<" " << norm(ybet) << endl;
}

void CCalpop::caldia_GO(
    const PMap& map, const CParam& param, const DMATRIX& diffusion_x,
    const DMATRIX& advection_x, const DMATRIX& diffusion_y,
    const DMATRIX& advection_y) {
    // cout << __LINE__ << endl;
    a.initialize();
    b.initialize();
    c.initialize();
    d.initialize();
    e.initialize();
    f.initialize();

    CBord bord;
    int dt = 0;  // 2*iterationNumber;
    double dx = param.deltaX;
    double dxdx = dx * dx;
    double twodxdx = 2 * dxdx;
    int nti = param.get_nbi();

    for (int j = map.jmin; j <= map.jmax; j++) {
        for (int i = map.iinf[j]; i <= map.isup[j]; i++) {
            bord.b = map.bord_cell[i][j];
            char pos = bord.cotex();
            double sigma = diffusion_x[i][j];
            double sigmam = diffusion_x[i - 1][j];
            double sigmap = diffusion_x[i + 1][j];
            double uinf = advection_x[i - 1][j];
            double u = advection_x[i][j];
            double usup = advection_x[i + 1][j];

            if ((i == 1) && (map.carte[nti - 2][j])) {
                // cout << i << " " << j << " " << std::tolower(pos) << " " <<
                // map.bord_cell[i][j] << endl;
                sigmam = diffusion_x[nti - 2][j];
                uinf = advection_x[nti - 2][j];
            }
            if ((i == nti - 2) && (map.carte[1][j])) {
                // cout << i << " " << j << " " << std::tolower(pos) << " " <<
                // map.bord_cell[i][j] << endl;
                sigmap = diffusion_x[1][j];
                usup = advection_x[1][j];
            }

            a[j][i] = d1(pos, sigmam, sigma, uinf, twodxdx, dxdx, dx);
            b[j][i] = d2(pos, sigmam, sigma, sigmap, u, twodxdx, dxdx, dx, dt);
            c[j][i] = d3(pos, sigma, sigmap, usup, twodxdx, dxdx, dx);
        }
    }
    // exit(1);
    double dy = param.deltaY;
    double dydy = dy * dy;
    double twodydy = 2 * dydy;

    for (int i = map.imin; i <= map.imax; i++) {
        for (int j = map.jinf[i]; j <= map.jsup[i]; j++) {
            bord.b = map.bord_cell[i][j];
            char pos = bord.cotey();
            double sigma = diffusion_y[i][j];
            double sigmam = diffusion_y[i][j - 1];
            double sigmap = diffusion_y[i][j + 1];
            double vinf = advection_y[i][j - 1];
            double v = advection_y[i][j];
            double vsup = advection_y[i][j + 1];

            d[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
            e[i][j] = d2(pos, sigmam, sigma, sigmap, v, twodydy, dydy, dy, dt);
            f[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
        }
    }
    ybet_comp(map);
}

double CCalpop::d1(
    const char pos, const double sigmam, const double sigma, const double uinf,
    const double twodd, const double dd, const double d) {
    double value = 0;
    switch (pos) {
        case SANS:     // si frontieres ouvertes
        case D_FERME:  // ou fermee a droite seulement
            if (uinf > 0)
                value = -((sigmam + sigma) / twodd) - (uinf / d);
            else
                value = -((sigmam + sigma) / twodd);
            break;
    }
    return value;
}

double CCalpop::d2(
    const char pos, const double sigmam, const double sigma,
    const double sigmap, const double u, const double twodd, const double dd,
    const double d, const int dt) {
    double value = 0;
    switch (pos) {
        case SANS:
            if (u > 0)
                value = ((sigmam + 2 * sigma + sigmap) / twodd) + (u / d) + dt;
            else
                value = ((sigmam + 2 * sigma + sigmap) / twodd) - (u / d) + dt;
            break;

        case G_FERME:
            if (u > 0)
                value = ((sigma + sigmap) / twodd) + (u / d) + dt;
            else
                value = ((sigma + sigmap) / twodd) + dt;
            break;

        case D_FERME:
            if (u > 0)
                value = ((sigmam + sigma) / twodd) + dt;
            else
                value = ((sigmam + sigma) / twodd) - (u / d) + dt;
            break;

        case TERRE:
            value = dt;
            break;
    }
    return value;
}

double CCalpop::d3(
    const char pos, const double sigma, const double sigmap, const double usup,
    const double twodd, const double dd, const double d) {
    double value = 0;
    switch (pos) {
        case SANS:
        case G_FERME:
            if (usup > 0)
                value = -((sigma + sigmap) / twodd);
            else
                value = -((sigma + sigmap) / twodd) + (usup / d);
            break;
    }
    return value;
}

void CCalpop::xbet_comp(const PMap& map) {
    int dt = 2 * iterationNumber;
    for (int j = map.jmin; j <= map.jmax; j++) {
        int inf = map.iinf[j];
        xbet[j][inf] = 1 / (bm[j][inf] + dt);
        for (int i = inf + 1; i <= map.isup[j]; i++) {
            xbet[j][i] = bm[j][i] + dt - c[j][i - 1] * a[j][i] * xbet[j][i - 1];
            xbet[j][i] = 1 / xbet[j][i];
        }
    }
}

void CCalpop::ybet_comp(const PMap& map) {
    int dt = 2 * iterationNumber;
    for (int i = map.imin; i <= map.imax; i++) {
        int inf = map.jinf[i];
        ybet[i][inf] = 1 / (e[i][inf] + dt);
        for (int j = inf + 1; j <= map.jsup[i]; j++) {
            ybet[i][j] = e[i][j] + dt - f[i][j - 1] * d[i][j] * ybet[i][j - 1];
            ybet[i][j] = 1 / ybet[i][j];
        }
    }
}

// FIN DU SOUS PROGRAMME CALDIA
