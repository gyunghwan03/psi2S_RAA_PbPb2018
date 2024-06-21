#include "TMath.h"

double cal_p_value(double nll_low_dim, double nll_high_dim, int diff_dim)
{
    double nll_diff = 2*(nll_low_dim - nll_high_dim);
    double p_value = 100*TMath::Prob(nll_diff, diff_dim);
    return p_value;
}