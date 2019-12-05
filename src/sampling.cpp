#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// Get probability given variable
vec pnormVec(const vec& a, double mu, double sigma, int lt, int lg){
    return as<vec>(wrap(pnorm(as<NumericVector>(wrap(a)), mu, sigma, lt, lg)));
}

// Get normal random variables
vec rnormVec(const vec& mu, const vec& sigma){
    return as<vec>(wrap(mapply(as<NumericVector>(wrap(mu)), as<NumericVector>(wrap(sigma)), R::rnorm)));
}

// Get n bernoulli trials
vec rbinomVec(const vec& n, const vec& p){
    return as<vec>(wrap(mapply(as<NumericVector>(wrap(n)), as<NumericVector>(wrap(p)), R::rbinom)));
}

// Get multinomial random variables
IntegerVector rmultinomC(int n, NumericVector prob){
    //n is size instead of total number of random samples
    //stackoverflow.com/questions/24618370/using-rmultinom-with-rcpp
    int k = prob.size();
    IntegerVector ans(k);
    R::rmultinom(n, prob.begin(), k, ans.begin());
    return ans;
}


// Sample g
void gSample(vec& g, int ng, mat& x, mat& dz, mat& sig2, rowvec& s, vec& gammam, vec& denom, uvec& idx_na, double nu_g, double tausq_g, vec& upd, const vec& step_g){
    vec g_temp = rnormVec(g, step_g);
    vec g_diff = g_temp - g;
    vec g2_diff = g_diff % (g_temp + g);
    mat msz = dz.each_row() + s;
    mat ldiff = x - msz;
    ldiff = ((ldiff.each_col() % g_diff).each_col() - g2_diff / 2) / sig2;
    ldiff(idx_na) = pnormVec(((msz.each_col() + g_temp)(idx_na) * gammam(0) + gammam(1)) / denom, 0, 1, 1, 1) - pnormVec(((msz.each_col() + g)(idx_na) * gammam(0) + gammam(1)) / denom, 0, 1, 1, 1);
    uvec idx = find(sum(ldiff, 1) + (nu_g * g_diff - g2_diff / 2) / tausq_g - as<vec>(wrap(log(runif(ng)))) > 0);
    if (idx.n_elem > 0){
        g(idx) = g_temp(idx);
        upd(idx) += 1;
    }
}


// Sample g hyperparameter νg from
double nu_gSample(double m_g, double dsq_g, double tausq_g, vec& g, int ng){
    double v = 1 / (1 / dsq_g + ng / tausq_g);
    return R::rnorm((m_g / dsq_g + sum(g) / tausq_g) * v, sqrt(v));
}


// Sample g hyperparameter τg2 from
double tausq_gSample(double a_g, double b_g, double nu_g, vec& g, int ng){
    return 1 / R::rgamma(a_g + ng / 2, 1 / (b_g + sum(pow(g - nu_g, 2)) / 2));
}


// Sample Gene-specific expression increment δj
void deltaSample(vec& delta, int ng, int ns, mat& mis, mat& z, mat& psi, vec& sigma1sq, mat& mus, vec& gammam, double alpha_delta, double beta_delta, vec& upd, const vec& step_delta){
    vec delta_temp = rnormVec(delta, step_delta);
    uvec idx_pos = find(delta_temp > 0);
    int npos = idx_pos.n_elem;
    if (npos > 0){
        mat psi_p = psi.rows(idx_pos), mus_p = mus.rows(idx_pos);
        uvec idx_pos_na = find(mis.rows(idx_pos) % z.rows(idx_pos) == 1);
        vec delta_diff = delta_temp(idx_pos) - delta(idx_pos);
        vec delta2_diff = delta_diff % (delta_temp(idx_pos) + delta(idx_pos));
        mat ldiff = z.rows(idx_pos) % (((psi_p.each_col() % delta_diff).each_col() - delta2_diff / 2).each_col() / sigma1sq(idx_pos));
        ldiff.replace(datum::nan, 0);
        if (idx_pos_na.n_elem > 0){
            //mat temp = repmat(sqrt(1 + pow(gammam(0), 2) * sigma1sq(idx_pos)), 1, ns);
            //vec denom = temp(idx_pos_na);
            umat sub_na = ind2sub(size(psi_p), idx_pos_na);
            vec denom = sqrt(1 + pow(gammam(0), 2) * sigma1sq(idx_pos(sub_na.row(0).t())));
            ldiff(idx_pos_na) = pnormVec((gammam(0) * (mus_p.each_col() + delta_temp(idx_pos))(idx_pos_na) + gammam(1)) / denom, 0, 1, 1, 1) - pnormVec((gammam(0) * (mus_p.each_col() + delta(idx_pos))(idx_pos_na) + gammam(1)) / denom, 0, 1, 1, 1);
        }
        uvec idx = idx_pos(find(sum(ldiff, 1) + (alpha_delta - 1) * (log(delta_temp(idx_pos)) - log(delta(idx_pos))) - beta_delta * delta_diff - as<vec>(wrap(log(runif(npos)))) > 0));
        if (idx.n_elem > 0){
            delta(idx) = delta_temp(idx);
            upd(idx) += 1;
        }
    }
}


// Sample δ hyperparameter βδ from
double beta_deltaSample(double a_l, double b_l, double alpha_delta, vec& delta, int ng){
    return R::rgamma(a_l + ng * alpha_delta, 1 / (b_l + sum(delta)));
}


// Sample σ0 hyperparameter β0 from
void sigma0sqSample(vec& sigma0sq, int ns, mat& mis, mat& z, mat& psi, mat& mus, vec& gammam, double alpha0, double beta0, vec& upd, const vec& step_s0){
    vec sigma0sq_temp = rnormVec(sigma0sq, step_s0);
    uvec idx_pos = find(sigma0sq_temp > 0);
    int npos = idx_pos.n_elem;
    if (npos > 0){
        uvec idx_pos_na = find(mis.rows(idx_pos) % (1 - z.rows(idx_pos)) == 1);
        vec a0t = alpha0 + sum((1 - mis.rows(idx_pos)) % (1 - z.rows(idx_pos)), 1) / 2;
        mat temp = (1 - z.rows(idx_pos)) % pow(psi.rows(idx_pos), 2);
        temp.replace(datum::nan, 0);
        vec b0t = beta0 + sum(temp, 1) / 2;
        vec ldiff = - (a0t + 1) % (log(sigma0sq_temp(idx_pos)) - log(sigma0sq(idx_pos))) - b0t / sigma0sq_temp(idx_pos) + b0t / sigma0sq(idx_pos);
        if (idx_pos_na.n_elem > 0){
            mat temp = mus.rows(idx_pos);
            vec num = gammam(0) * temp(idx_pos_na) + gammam(1);
            mat lmiss = zeros<mat>(npos, ns);
            //temp = ones<mat>(npos, ns);
            //lmiss(idx_pos_na) = pnormVec(num / sqrt((temp.each_col() +  pow(gammam(0), 2) * sigma0sq_temp(idx_pos))(idx_pos_na)), 0, 1, 1, 1) - pnormVec(num / sqrt((temp.each_col() + pow(gammam(0), 2) * sigma0sq(idx_pos))(idx_pos_na)), 0, 1, 1, 1);
            umat sub_na = ind2sub(size(temp), idx_pos_na);
            lmiss(idx_pos_na) = pnormVec(num / sqrt(1 + pow(gammam(0), 2) * sigma0sq_temp(idx_pos(sub_na.row(0).t()))), 0, 1, 1, 1) - pnormVec(num / sqrt(1 + pow(gammam(0), 2) * sigma0sq(idx_pos(sub_na.row(0).t()))), 0, 1, 1, 1);
            ldiff += sum(lmiss, 1);
        }
        uvec idx = idx_pos(find(ldiff - as<vec>(wrap(log(runif(npos)))) > 0));
        if (idx.n_elem > 0){
            sigma0sq(idx) = sigma0sq_temp(idx);
            upd(idx) += 1;
        }
    }
}


// Sample beta0
double beta0Sample(double a_l, double b_l, double alpha0, vec& sigma0sq, int ng){
    return R::rgamma(alpha0 * ng + a_l, 1 / (sum(1 / sigma0sq) + b_l));
}


// Sample σ1 hyperparameter β1 from
void sigma1sqSample(vec& sigma1sq, int ns, mat& mis, mat& z, mat& psi, mat& mus, vec& delta, vec& gammam, double alpha1, double beta1, vec& upd, const vec& step_s1){
    vec sigma1sq_temp = rnormVec(sigma1sq, step_s1);
    uvec idx_pos = find(sigma1sq_temp > 0);
    int npos = idx_pos.n_elem;
    if (npos > 0){
        uvec idx_pos_na = find(mis.rows(idx_pos) % z.rows(idx_pos) == 1);
        vec a1t = alpha1 + sum((1 - mis.rows(idx_pos)) % z.rows(idx_pos), 1) / 2;
        mat psi_p = psi.rows(idx_pos);
        mat temp = z.rows(idx_pos) % pow(psi_p.each_col() - delta(idx_pos), 2);
        temp.replace(datum::nan, 0);
        vec b1t = beta1 + sum(temp, 1) / 2;
        vec ldiff = - (a1t + 1) % (log(sigma1sq_temp(idx_pos)) - log(sigma1sq(idx_pos))) - b1t / sigma1sq_temp(idx_pos) + b1t / sigma1sq(idx_pos);
        if (idx_pos_na.n_elem > 0){
            temp = mus.rows(idx_pos);
            vec num = gammam(0) * (temp.each_col() + delta(idx_pos))(idx_pos_na) + gammam(1);
            mat lmiss = zeros<mat>(npos, ns);
            //temp = ones<mat>(npos, ns);
            //lmiss(idx_pos_na) = pnormVec(num / sqrt((temp.each_col() + pow(gammam(0), 2) * sigma1sq_temp(idx_pos))(idx_pos_na)), 0, 1, 1, 1) - pnormVec(num / sqrt((temp.each_col() + pow(gammam(0), 2) * sigma1sq(idx_pos))(idx_pos_na)), 0, 1, 1, 1);
            umat sub_na = ind2sub(size(temp), idx_pos_na);
            lmiss(idx_pos_na) = pnormVec(num / sqrt(1 + pow(gammam(0), 2) * sigma1sq_temp(idx_pos(sub_na.row(0).t()))), 0, 1, 1, 1) - pnormVec(num / sqrt(1 + pow(gammam(0), 2) * sigma1sq(idx_pos(sub_na.row(0).t()))), 0, 1, 1, 1);
            ldiff += sum(lmiss, 1);
        }
        uvec idx = idx_pos(find(ldiff - as<vec>(wrap(log(runif(npos)))) > 0));
        sigma1sq(idx) = sigma1sq_temp(idx);
        upd(idx) += 1;
    }
}


// Sample beta1
double beta1Sample(double a_l, double b_l, double alpha1, vec& sigma1sq, int ng){
    return R::rgamma(alpha1 * ng + a_l, 1 / (sum(1 / sigma1sq) + b_l));
}


// Sample binary expression status zij
void zSample(int nuc, int ng, List& idx_c, mat& psi, mat& mus, mat& mis, mat& z, vec& delta, vec& sigma0sq, vec& sigma1sq, vec& gammam, double beta_pi, double alpha_pi){
    uvec idx_k, j;
    int nk;
    vec nzk, nzki, prob, ratio, mkj;
    for (int k = 0; k < nuc; k++){
        idx_k = as<uvec>(wrap(idx_c[k]));
        nk = idx_k.n_elem;
        if (nk == 1){
            prob = 1 / (1 + exp(pow(psi.col(idx_k(0)) - delta, 2) / (2 * sigma1sq) - pow(psi.col(idx_k(0)), 2) / (2 * sigma0sq)) * beta_pi / alpha_pi % sqrt(sigma1sq / sigma0sq));
            j = find(mis.col(idx_k(0)) == 1);
            if (j.n_elem > 0){
                mkj = mus.col(idx_k(0));
                mkj = mkj(j);
                prob(j) = 1 / (1 + exp(pnormVec((gammam(0) * mkj + gammam(1)) / sqrt(1 + pow(gammam(0), 2) * sigma0sq(j)), 0, 1, 1, 1) - pnormVec((gammam(0) * (mkj + delta(j)) + gammam(1)) / sqrt(1 + pow(gammam(0), 2) * sigma1sq(j)), 0, 1, 1, 1)) * beta_pi / alpha_pi);
                z.col(idx_k(0)) = rbinomVec(ones<vec>(ng), prob);
            }
        } else{
            nzk = sum(z.cols(idx_k), 1);
            for (int i = 0; i < nk; i++){
                nzki = nzk - z.col(idx_k(i));
                ratio = sqrt(sigma1sq / sigma0sq) % exp(pow(psi.col(idx_k(i)) - delta, 2) / (2 * sigma1sq) - pow(psi.col(idx_k(i)), 2) / (2 * sigma0sq));
                j = find(mis.col(idx_k(i)) == 1);
                if (j.n_elem > 0){
                    mkj = mus.col(idx_k(i));
                    mkj = mkj(j);
                    ratio(j) = exp(pnormVec((gammam(0) * mkj + gammam(1)) / sqrt(1 + pow(gammam(0), 2) * sigma0sq(j)), 0, 1, 1, 1) - pnormVec((gammam(0) * (mkj + delta(j)) + gammam(1)) / sqrt(1 + pow(gammam(0), 2) * sigma1sq(j)), 0, 1, 1, 1));
                }
                prob = 1 / (1 + ratio % (beta_pi + nk - 1 - nzki) / (alpha_pi + nzki));
                z.col(idx_k(i)) = rbinomVec(ones<vec>(ng), prob);
                nzk = nzki + z.col(idx_k(i));
            }
        }
    }
}


// Sample dropout probability parameters γ = (γ1,γ2)T
void gammamSample(vec& gammam, double nu1, double tau1sq, double nu2, double tau2sq, uvec& idx_na, uvec& idx_nna, mat& sig2, mat& x, mat& mut, double& upd, double step_gamma){
    vec gammam_temp = rnormVec(gammam, ones<vec>(2) * step_gamma);
    if ((pow(gammam(0) - nu1, 2) - pow(gammam_temp(0) - nu1, 2)) / (2 * tau1sq) + (pow(gammam(1) - nu2, 2) - pow(gammam_temp(1) - nu2, 2)) / (2 * tau2sq) + sum(pnormVec((gammam_temp(0) * mut(idx_na) + gammam_temp(1)) / sqrt(1 + pow(gammam_temp(0), 2) * sig2(idx_na)), 0, 1, 1, 1) - pnormVec((gammam(0) * mut(idx_na) + gammam(1)) / sqrt(1 + pow(gammam(0), 2) * sig2(idx_na)), 0, 1, 1, 1)) + sum(pnormVec(gammam_temp(0) * x(idx_nna) + gammam_temp(1), 0, 1, 0, 1) - pnormVec(gammam(0) * x(idx_nna) + gammam(1), 0, 1, 0, 1)) > log(R::runif(0, 1))){
        gammam = gammam_temp;
        upd += 1;
    }
}


// Sample component label ci
void cluSample(Row<int>& clu, int ns, int ng, double alpha_pi, double beta_pi, mat& z, double alpha_lambda){
    uvec idx, idx_k;
    Col<int> ci, uci;
    int nci, nk;
    vec prob, lprob, zab;
    mat zi;
    for (int i = 0; i < ns; i++){
        if (i == 0){
            idx = regspace<uvec>(1, 1, ns - 1);
        } else if(i == ns - 1){
            idx = regspace<uvec>(0, 1, ns - 2);
        } else{
            idx = join_cols(regspace<uvec>(0, 1, i - 1), regspace<uvec>(i + 1, 1, ns - 1));
        }
        ci = clu(idx);
        uci = unique(ci);
        nci = uci.n_elem;
        zab = alpha_pi * z.col(i) + beta_pi * (1 - z.col(i));
        prob = zeros<vec>(nci + 1);
        lprob = zeros<vec>(nci + 1);
        zi = z.cols(idx);
        zi = 1 - abs(zi.each_col() - z.col(i));
        for (int k = 0; k < nci; k++){
            idx_k = find(ci == uci(k));
            nk = idx_k.n_elem;
            lprob(k) = log(nk) + sum(log(sum(zi.cols(idx_k), 1) + zab)) - ng * log(nk + alpha_pi + beta_pi);
        }
        lprob(nci) = log(alpha_lambda) + sum(log(zab)) - ng * log(alpha_pi + beta_pi);
        for (int k = 0; k < nci + 1; k++){
            prob(k) = 1 / sum(exp(lprob - lprob(k)));
        }
        idx = find(as<vec>(wrap(rmultinomC(1, as<NumericVector>(wrap(prob))))) == 1);
        if (idx(0) < nci){
            clu(i) = uci(idx(0));
        } else{
            clu(i) = max(uci) + 1;
        }
    }
}


// Sample expression probability π hyperparameter βπ
List beta_piSample(double beta_pi, double alpha_pi, double a_l, double b_l, int ng, Row<int>& uc, List idx_c, mat& z, double& upd, double step_bp){
    double beta_pi_temp = R::rnorm(beta_pi, step_bp);
    double lprob = ng * uc.n_elem * (lgamma(alpha_pi + beta_pi) - lgamma(alpha_pi) - lgamma(beta_pi)) + (a_l - 1) * log(beta_pi) - b_l * beta_pi;
    double lprob_temp = R_NegInf;
    vec la, zsk;
    uvec idx_k;
    int nk;
    if (beta_pi_temp > 0){
        lprob_temp = ng * uc.n_elem * (lgamma(alpha_pi + beta_pi_temp) - lgamma(alpha_pi) - lgamma(beta_pi_temp)) + (a_l - 1) * log(beta_pi_temp) - b_l * beta_pi_temp;
    }
    for (int k = 0; k < idx_c.size(); k++){
        idx_k = as<uvec>(wrap(idx_c[k]));
        nk = idx_k.n_elem;
        zsk = sum(z.cols(idx_k), 1);
        la = lgamma(alpha_pi + zsk);
        lprob += sum(la + lgamma(beta_pi + nk - zsk)) - ng * lgamma(alpha_pi + beta_pi + nk);
        if (beta_pi_temp > 0){
            lprob_temp += sum(la + lgamma(beta_pi_temp + nk - zsk)) - ng * lgamma(alpha_pi + beta_pi_temp + nk);
        }
    }
    if (lprob_temp - lprob > log(R::runif(0, 1))){
        upd += 1;
        return List::create(_["beta_pi"] = beta_pi_temp, _["lprob"] = lprob_temp);
    } else{
        return List::create(_["beta_pi"] = beta_pi, _["lprob"] = lprob);
    }
}


// Sample Drichlet-process concentration parameter αλ
List alpha_lambdaSample(double alpha_lambda, double a_l, double b_l, int ns, Row<int>& uc, double& upd, double step_al){
    double alpha_lambda_temp = R::rnorm(alpha_lambda, step_al);
    double lprob = (a_l + uc.n_elem - 1) * log(alpha_lambda) - b_l * alpha_lambda + lgamma(alpha_lambda) - lgamma(alpha_lambda + ns);
    double lprob_temp = R_NegInf;
    if (alpha_lambda_temp > 0){
        lprob_temp = (a_l + uc.n_elem - 1) * log(alpha_lambda_temp) - b_l * alpha_lambda_temp + lgamma(alpha_lambda_temp) - lgamma(alpha_lambda_temp + ns);
    }
    if (lprob_temp - lprob > log(R::runif(0, 1))){
        upd += 1;
        return List::create(_["alpha_lambda"] = alpha_lambda_temp, _["lprob"] = lprob_temp);
    } else{
        return List::create(_["alpha_lambda"] = alpha_lambda, _["lprob"] = lprob);
    }
}


// [[Rcpp::export]]
void sampling(mat& x, rowvec& s, int NT, int nthin, int nupd, int Nburn,
              vec& g, double nu_g, double tausq_g, vec& delta, double beta_delta,
              vec& sigma0sq, double beta0, vec& sigma1sq, double beta1, mat& z,
              vec& gammam, Row<int>& clu, double beta_pi, double alpha_lambda,
              mat& g_seq, mat& delta_seq, mat& sigma0sq_seq, mat& sigma1sq_seq,
              mat& z_seq, Mat<int>& c_seq, mat& gamma_seq, vec& beta_pi_seq,
              vec& alpha_lambda_seq, vec& nu_g_seq, vec& tausq_g_seq,
              vec& beta_delta_seq, vec& beta0_seq, vec& beta1_seq, vec& llambda, vec& lpi,
              vec& lobs, vec& lother, vec& lpost, double m_g, double dsq_g, double a_g,
              double b_g, double alpha_delta, double alpha0, double alpha1, double a_l, double b_l,
              double nu1, double tau1sq, double nu2, double tau2sq, double alpha_pi) {

    // Check finite/non-finite values (NaN, Inf) etc using Armadillo
    uvec idx_na = find_nonfinite(x), idx_nna = find_finite(x), idx;
    int ng = x.n_rows, ns = x.n_cols, l;

    mat mis(ng, ns, fill::zeros);
    mis(idx_na).fill(1.0);

    mat dz = z.each_col() % delta;
    mat sig2 = as<mat>(wrap(1 - z)).each_col() % sigma0sq + z.each_col() % sigma1sq;
    vec denom = sqrt(1 + pow(gammam(0), 2) * sig2(idx_na));
    Row<int> uc = unique(clu);
    int nuc = uc.n_elem;
    List idx_c = List(nuc);
    for (int k = 0; k < nuc; k++){
        idx_c[k] = find(clu == uc(k));
    }

    mat psi, mut;
    List temp;
    vec step(4 * ng + 3, fill::ones);
    step(span(0, 4 * ng - 1)) *= 0.5;
    step(4 * ng) = 0.01;
    step(span(4 * ng + 1, 4 * ng + 2)) *= 0.05;
    vec upd(4 * ng + 3, fill::zeros), upd_g(ng, fill::zeros), upd_d(ng, fill::zeros), upd_s0(ng, fill::zeros), upd_s1(ng, fill::zeros);
    double lprob_alpha_lambda, lprob_beta_pi;
    mat mus = zeros<mat>(ng, ns);

    // Do MCMC iterations; NT --> number of iterations
    for (int nt = 0; nt < NT; nt++){
        gSample(g, ng, x, dz, sig2, s, gammam, denom, idx_na, nu_g, tausq_g, upd_g, step.rows(0, ng - 1));

        mus.zeros();
        mus = (mus.each_col() + g).each_row() + s;
        psi = x - mus;

        nu_g = nu_gSample(m_g, dsq_g, tausq_g, g, ng);
        tausq_g = tausq_gSample(a_g, b_g, nu_g, g, ng);

        deltaSample(delta, ng, ns, mis, z, psi, sigma1sq, mus, gammam, alpha_delta, beta_delta, upd_d, step.rows(ng, 2 * ng - 1));

        beta_delta = beta_deltaSample(a_l, b_l, alpha_delta, delta, ng);

        sigma0sqSample(sigma0sq, ns, mis, z, psi, mus, gammam, alpha0, beta0, upd_s0, step.rows(2 * ng, 3 * ng - 1));

        beta0 = beta0Sample(a_l, b_l, alpha0, sigma0sq, ng);

        sigma1sqSample(sigma1sq, ns, mis, z, psi, mus, delta, gammam, alpha1, beta1, upd_s1, step.rows(3 * ng, 4 * ng - 1));

        beta1 = beta1Sample(a_l, b_l, alpha1, sigma1sq, ng);

        zSample(nuc, ng, idx_c, psi, mus, mis, z, delta, sigma0sq, sigma1sq, gammam, beta_pi, alpha_pi);

        dz = z.each_col() % delta;
        mut = mus + dz;
        sig2 = as<mat>(wrap(1 - z)).each_col() % sigma0sq + z.each_col() % sigma1sq;

        gammamSample(gammam, nu1, tau1sq, nu2, tau2sq, idx_na, idx_nna, sig2, x, mut, upd(4 * ng), step(4 * ng));
        denom = sqrt(1 + pow(gammam(0), 2) * sig2(idx_na));

        cluSample(clu, ns, ng, alpha_pi, beta_pi, z, alpha_lambda);
        uc = unique(clu);
        nuc = uc.n_elem;
        idx_c = List(nuc);
        for (int k = 0; k < nuc; k++){
            idx_c[k] = find(clu == uc(k));
        }
        uc = linspace<Row<int> >(0, nuc - 1, nuc);
        for (int k = 0; k < nuc; k++){
            clu(as<uvec>(wrap(idx_c[k]))).fill(k);
        }

        temp = beta_piSample(beta_pi, alpha_pi, a_l, b_l, ng, uc, idx_c, z, upd(4 * ng + 1), step(4 * ng + 1));
        beta_pi = temp["beta_pi"];
        lprob_beta_pi = temp["lprob"];

        temp = alpha_lambdaSample(alpha_lambda, a_l, b_l, ns, uc, upd(4 * ng + 2), step(4 * ng + 2));
        alpha_lambda = temp["alpha_lambda"];
        lprob_alpha_lambda = temp["lprob"];

        if (nt < Nburn * 0.8 && nt - nupd * floor(nt / nupd) == 0){
            upd.rows(0, ng - 1) = upd_g;
            upd.rows(ng, 2 * ng - 1) = upd_d;
            upd.rows(2 * ng, 3 * ng - 1) = upd_s0;
            upd.rows(3 * ng, 4 * ng - 1) = upd_s1;
            upd /= nupd;
            idx = find(upd < 0.15);
            step(idx) /= 1.5;
            idx = find(upd > 0.25);
            step(idx) *= 1.5;
            upd.zeros();
            upd_g.zeros();
            upd_d.zeros();
            upd_s0.zeros();
            upd_s1.zeros();
        }

        l = floor((nt - Nburn) / nthin);
        if (nt >= Nburn && nt - Nburn - nthin * l == 0){
            delta_seq.col(l) = delta;
            sigma0sq_seq.col(l) = sigma0sq;
            sigma1sq_seq.col(l) = sigma1sq;
            g_seq.col(l) = g;
            z_seq.col(l) = as<vec>(wrap(z));
            c_seq.col(l) = clu.t();
            gamma_seq.col(l) = gammam;
            beta_pi_seq(l) = beta_pi;
            alpha_lambda_seq(l) = alpha_lambda;
            nu_g_seq(l) = nu_g;
            tausq_g_seq(l) = tausq_g;
            beta_delta_seq(l) = beta_delta;
            beta0_seq(l) = beta0;
            beta1_seq(l) = beta1;

            llambda(l) = lprob_alpha_lambda;
            for (int k = 0; k < idx_c.size(); k++){
                llambda(l) += Rcpp::internal::lfactorial(as<uvec>(wrap(idx_c[k])).n_elem - 1);
            }
            lpi(l) = lprob_beta_pi;
            lobs(l) = sum(pnormVec(gammam(0) * x(idx_nna) + gammam(1), 0, 1, 0, 1) - pow(x(idx_nna) - mut(idx_nna), 2) / (2 * sig2(idx_nna)) - log(sig2(idx_nna)) / 2) + sum(pnormVec((gammam(0) * mut(idx_na) + gammam(1)) / denom, 0, 1, 1, 1));
            lother(l) = -sum(pow(g - nu_g, 2)) / (2 * tausq_g) - log(tausq_g) * (ng / 2 + a_g + 1) - b_g / tausq_g - pow(nu_g - m_g, 2) / (2 * dsq_g) + alpha_delta * log(beta_delta) + (alpha_delta - 1) * sum(log(delta)) - beta_delta * sum(delta) + (a_l - 1) * log(beta_delta) - b_l * beta_delta + alpha0 * log(beta0) - (alpha0 + 1) * sum(log(sigma0sq)) - beta0 * sum(1 / sigma0sq) + alpha1 * log(beta1) - (alpha1 + 1) * sum(log(sigma1sq)) - beta1 * sum(1 / sigma1sq) + (a_l - 1) * (log(beta0) + log(beta1)) - b_l * (beta0 + beta1) - pow(gammam(0) - nu1, 2) / (2 * tau1sq) - pow(gammam(1) - nu2, 2) / (2 * tau2sq);
        }
    }
    lpost = llambda + lpi + lobs + lother;
}

