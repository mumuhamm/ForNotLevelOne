#include <iostream>
#include <cmath>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>

const double mu = 1.0; // Van der Pol parameter

// Van der Pol oscillator equations
void vanderpol(double t, const double y[], double dydt[]) {
    dydt[0] = y[1];
    dydt[1] = mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
}

// Fourth-order Runge-Kutta method for numerical integration
void RungeKutta(double t0, double y0[], double h, int steps, TGraph* graph) {
    double t = t0;
    double y[2];
    double k1[2], k2[2], k3[2], k4[2];
    graph->SetPoint(0, t, y0[0]);

    for (int i = 1; i <= steps; ++i) {
        vanderpol(t, y0, k1);
        for (int j = 0; j < 2; ++j)
            y[j] = y0[j] + 0.5 * h * k1[j];
        
        vanderpol(t + 0.5 * h, y, k2);
        for (int j = 0; j < 2; ++j)
            y[j] = y0[j] + 0.5 * h * k2[j];

        vanderpol(t + 0.5 * h, y, k3);
        for (int j = 0; j < 2; ++j)
            y[j] = y0[j] + h * k3[j];

        vanderpol(t + h, y, k4);
        for (int j = 0; j < 2; ++j)
            y0[j] += (h / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);

        t += h;
        graph->SetPoint(i, t, y0[0]);
    }
}

int main() {
    const double t0 = 0.0;         // Initial time
    const double y0[2] = {1.0, 0}; // Initial conditions: position and velocity
    const double h = 0.01;         // Step size
    const int steps = 1000;        // Number of integration steps

    TCanvas *canvas = new TCanvas("canvas", "Van der Pol Oscillator", 800, 600);
    TGraph *graph = new TGraph();
    graph->SetTitle("Van der Pol Oscillator;Time;t");
    graph->SetMarkerStyle(8);

    RungeKutta(t0, const_cast<double*>(y0), h, steps, graph);

    graph->Draw("APL");
    graph->GetXaxis()->SetTitle("Time");
    graph->GetYaxis()->SetTitle("Position");

    canvas->Draw();
    canvas->SaveAs("myplot.png");
    
    return 0;
}

