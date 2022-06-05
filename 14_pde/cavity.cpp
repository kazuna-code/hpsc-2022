#include <vector> // vector
#include <math.h> // pow

using namespace std;
typedef vector<vector<float>> matrix;

int main() {
    const int nx = 41;
    const int ny = 41;
    const int nt = 500;
    const int nit = 50;
    float dx = 2 / (nx - 1);
    float dy = 3 / (ny - 1);
    const float dt = 0.01;
    const int rho = 1;
    const float nu = 0.02;

    matrix u(ny, vector<float>(nx));
    matrix v(ny, vector<float>(nx));
    matrix p(ny, vector<float>(nx));
    matrix b(ny, vector<float>(nx));

    matrix un(ny, vector<float>(nx));
    matrix vn(ny, vector<float>(nx));
    matrix pn(ny, vector<float>(nx));

    for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
            u[j][i] = 0.0;
            v[j][i] = 0.0;
            p[j][i] = 0.0;
            b[j][i] = 0.0;
        }
    }

    for (int n=0; n<nt; n++) {
        for (int j=1; j<ny-1; j++) {
            for (int i=1; i<nx-1; i++) {
                b[j][i] = rho * (1 / dt *\
                          ((u[j][i+1] - u[j][i-1]) / (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -\
                          pow(((u[j][i+1] - u[j][i-1]) / (2 * dx)), 2.0) - 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy) *\
                          (v[j][i+1] - v[j][i-1]) / (2 * dx)) - pow(((v[j+1][i] - v[j-1][i]) / (2 * dy)), 2.0));
            }
        }
        for (int it=0; it<nit; it++) {
            // pn = p.copy()
            for (int j=0; j<ny; j++) {
                for (int i=0; i<nx; i++) {
                    pn[j][i] = p[j][i];
                }
            }
            for (int j=1; j<ny-1; j++) {
                for (int i=1; i<nx-1; i++) {
                    p[j][i] = (pow(dy, 2.0) * (pn[j][i+1] + pn[j][i-1]) +\
                               pow(dx, 2.0) * (pn[j+1][i] + pn[j-1][i]) -\
                               b[j][i] * pow(dx, 2.0) * pow(dy, 2.0))\
                               / (2 * (pow(dx, 2.0) + pow(dy, 2.0)));
                }
            }
            // p[:, -1] = p[:, -2]
            // p[:, 0] = p[:, 1]
            for (int j=0; j<ny; j++) {
                p[j][nx-1] = p[j][nx-2];
                p[j][0] = p[j][1];
            }
            // p[0, :] = p[1, :]
            // p[-1, :] = 0
            for (int i=0; i<nx; i++) {
                p[0][i] = p[1][i];
                p[ny-1][i] = 0.0;
            }
        }
        // un = u.copy()
        // vn = v.copy()
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                un[j][i] = u[j][i];
                vn[j][i] = v[j][i];
            }
        }
        for (int j=1; j<ny-1; j++) {
            for (int i=1; i<nx-1; i++) {
                u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i - 1])\
                                   - un[j][i] * dt / dy * (un[j][i] - un[j - 1][i])\
                                   - dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1])\
                                   + nu * dt / pow(dx, 2.0) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1])\
                                   + nu * dt / pow(dy, 2.0) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]);
                v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1])\
                                   - vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i])\
                                   - dt / (2 * rho * dx) * (p[j+1][i] - p[j-1][i])\
                                   + nu * dt / pow(dx, 2.0) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1])\
                                   + nu * dt / pow(dy, 2.0) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]);
            }
        }
        // u[:, 0] = 0
        // u[:, -1] = 0
        // v[:, 0] = 0
        // v[:, -1] = 0
        for (int j=0; j<ny; j++) {
            u[j][0] = u[j][nx-1] = v[j][0] = v[j][nx-1] = 0.0;
        }
        // u[0, :] = 0
        // u[-1, :] = 1
        // v[0, :] = 0
        // v[-1, :] = 0
        for (int i=0; i<nx; i++) {
            u[0][i] = v[0][i] = v[ny-1][i] = 0.0;
            u[ny-1][i] = 1.0;
        }
    }

    return 0;
}