// Branch-length (node-height) inference for the WF hidden Markov tree.
//
// The existing scorers (mitodrift.cpp) use one transition matrix for every edge. Here every edge e has its own
// matrix A(t_e), t_e = h(parent) - h(child) in generations, interpolated linearly (probability space) between the
// two bracketing matrices of a precomputed grid A(g_1), ..., A(g_G) (as interpolate_matrices() did in the original
// branch-length work). Layouts follow mitodrift.cpp:
//   E      : edge vector c(parents, children), 1-indexed, postorder (reorder_phylo), root = parent of the last edge
//   logP_l : per locus, state-major (logP[c * n + node]); the root row fixes the root state
//   A_g    : row-major probabilities, A[r * C + j] = P(child state j | parent state r)
// Node heights: tips 0, root = htotal (fixed here); an internal node lies between its highest child and its parent.
//
// Sampler: one internal non-root node v per iteration. Its position f = (h_p - h_v) / (h_p - lb_v) in (0, 1)
// (lb_v = highest child) gets a Gaussian random-walk step on logit(f); prior uniform on valid heights, so the
// Hastings ratio is f'(1 - f') / (f (1 - f)). Only the edges at v (its child edges and the edge above it) change,
// and only messages on the path from v to the root: cost (deg(v) + depth(v)) * C^2 per locus, not a full pass.
// Loci are scored in parallel (RcppParallel parallelFor; per-locus results summed serially, so runs are
// reproducible for a given seed whatever the thread count).

#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;

namespace {

struct BrlTree {
    int n = 0, m = 0, C = 0, L = 0, root = -1, nt = 0;
    std::vector<int> par, ch;                 // per edge (0-indexed)
    std::vector<int> up_edge;                 // per node: edge whose child is the node (-1 at root)
    std::vector<std::vector<int>> child_edges; // per node
    std::vector<int> parent_node;             // per node (-1 at root)
};

BrlTree build_tree(arma::Col<int> E, int C, int L, int n) {
    BrlTree T; T.C = C; T.L = L; T.n = n; T.m = E.n_elem / 2;
    E -= 1;
    T.par.resize(T.m); T.ch.resize(T.m);
    T.up_edge.assign(n, -1); T.parent_node.assign(n, -1); T.child_edges.assign(n, {});
    for (int i = 0; i < T.m; ++i) {
        T.par[i] = E(i); T.ch[i] = E(T.m + i);
        T.up_edge[T.ch[i]] = i; T.parent_node[T.ch[i]] = T.par[i];
        T.child_edges[T.par[i]].push_back(i);
    }
    T.root = E(T.m - 1);
    T.nt = 0; for (int v = 0; v < n; ++v) if (T.child_edges[v].empty()) ++T.nt;
    return T;
}

// Interpolated edge matrix in "exp-shifted" form: out_exp[r*C+j] = A[r][j] / max_j A[r][j], out_rowmax[r] = log max.
inline void edge_matrix(double t, const std::vector<double>& grid, const std::vector<std::vector<double>>& Ag,
                        int C, double* out_exp, double* out_rowmax) {
    const int G = static_cast<int>(grid.size());
    int lo = 0; double w = 0.0;
    if (t <= grid[0]) { lo = 0; w = 0.0; }
    else if (t >= grid[G - 1]) { lo = G - 1; w = 0.0; }
    else {
        lo = static_cast<int>(std::upper_bound(grid.begin(), grid.end(), t) - grid.begin()) - 1;
        w = (t - grid[lo]) / (grid[lo + 1] - grid[lo]);
    }
    const double* A0 = Ag[lo].data();
    const double* A1 = (w > 0.0) ? Ag[lo + 1].data() : A0;
    for (int r = 0; r < C; ++r) {
        double mx = 0.0;
        for (int j = 0; j < C; ++j) {
            const double a = (1.0 - w) * A0[r * C + j] + w * A1[r * C + j];
            out_exp[r * C + j] = a; if (a > mx) mx = a;
        }
        if (mx <= 0.0) mx = 1e-300;
        for (int j = 0; j < C; ++j) out_exp[r * C + j] /= mx;
        out_rowmax[r] = std::log(mx);
    }
}

// Message a child sends to its parent: out[c] = log sum_j A[c][j] exp(logP[j*n + child] + msg_child[j]).
inline void child_message(const double* logP, int n, int child, const double* msg_child, const double* expA,
                          const double* rowmax, int C, double* out, double* temp, double* u) {
    double mt = -std::numeric_limits<double>::infinity();
    for (int j = 0, off = child; j < C; ++j, off += n) {
        const double v = logP[off] + msg_child[j]; temp[j] = v; if (v > mt) mt = v;
    }
    if (!std::isfinite(mt)) { for (int c = 0; c < C; ++c) out[c] = -std::numeric_limits<double>::infinity(); return; }
    for (int j = 0; j < C; ++j) u[j] = std::exp(temp[j] - mt);
    for (int c = 0; c < C; ++c) {
        const double* row = expA + c * C; double s = 0.0;
        for (int j = 0; j < C; ++j) s += row[j] * u[j];
        out[c] = rowmax[c] + mt + std::log(s);
    }
}

inline double root_logz(const double* logP, int n, int root, const double* msg_root, int C) {
    double mx = -std::numeric_limits<double>::infinity(); std::vector<double> t(C);
    for (int c = 0; c < C; ++c) { t[c] = logP[c * n + root] + msg_root[c]; if (t[c] > mx) mx = t[c]; }
    if (!std::isfinite(mx)) return mx;
    double s = 0.0; for (int c = 0; c < C; ++c) s += std::exp(t[c] - mx);
    return mx + std::log(s);
}

// Full pass for one locus: fills msg (n*C, node-major) and contrib (m*C), returns logZ.
double full_pass(const BrlTree& T, const double* logP, const std::vector<double>& expA, const std::vector<double>& rowmax,
                 std::vector<double>& msg, std::vector<double>& contrib) {
    const int C = T.C;
    std::fill(msg.begin(), msg.end(), 0.0);
    std::vector<double> temp(C), u(C);
    for (int e = 0; e < T.m; ++e) {      // postorder: children before parents
        const int c = T.ch[e], p = T.par[e];
        child_message(logP, T.n, c, &msg[c * C], &expA[(size_t)e * C * C], &rowmax[(size_t)e * C], C, &contrib[(size_t)e * C], temp.data(), u.data());
        for (int s = 0; s < C; ++s) msg[p * C + s] += contrib[(size_t)e * C + s];
    }
    return root_logz(logP, T.n, T.root, &msg[T.root * C], C);
}

struct FullPassWorker : public Worker {
    const BrlTree& T; const std::vector<std::vector<double>>& logP;
    const std::vector<double>& expA; const std::vector<double>& rowmax;
    std::vector<std::vector<double>>& msg; std::vector<std::vector<double>>& contrib; std::vector<double>& logz;
    FullPassWorker(const BrlTree& T, const std::vector<std::vector<double>>& logP, const std::vector<double>& expA,
                   const std::vector<double>& rowmax, std::vector<std::vector<double>>& msg,
                   std::vector<std::vector<double>>& contrib, std::vector<double>& logz)
        : T(T), logP(logP), expA(expA), rowmax(rowmax), msg(msg), contrib(contrib), logz(logz) {}
    void operator()(std::size_t b, std::size_t e) {
        for (std::size_t l = b; l < e; ++l) logz[l] = full_pass(T, logP[l].data(), expA, rowmax, msg[l], contrib[l]);
    }
};

// Proposal evaluation for one node move, per locus. New values go to per-locus scratch; nothing is committed.
struct MoveWorker : public Worker {
    const BrlTree& T; const std::vector<std::vector<double>>& logP;
    const std::vector<std::vector<double>>& msg; const std::vector<std::vector<double>>& contrib;
    const std::vector<double>& expA; const std::vector<double>& rowmax;      // current edge matrices
    const std::vector<int>& cedges;            // child edges of v
    const std::vector<double>& cexp; const std::vector<double>& crow;        // new matrices for child edges
    const std::vector<double>& uexp; const std::vector<double>& urow;        // new matrix for v's up edge
    const std::vector<int>& path;              // v, parent(v), ..., root
    std::vector<std::vector<double>>& scratch; // per locus: [child contribs | path msgs | path up-contribs]
    std::vector<double>& logz_new;
    MoveWorker(const BrlTree& T, const std::vector<std::vector<double>>& logP, const std::vector<std::vector<double>>& msg,
               const std::vector<std::vector<double>>& contrib, const std::vector<double>& expA, const std::vector<double>& rowmax,
               const std::vector<int>& cedges, const std::vector<double>& cexp, const std::vector<double>& crow,
               const std::vector<double>& uexp, const std::vector<double>& urow, const std::vector<int>& path,
               std::vector<std::vector<double>>& scratch, std::vector<double>& logz_new)
        : T(T), logP(logP), msg(msg), contrib(contrib), expA(expA), rowmax(rowmax), cedges(cedges), cexp(cexp), crow(crow),
          uexp(uexp), urow(urow), path(path), scratch(scratch), logz_new(logz_new) {}
    void operator()(std::size_t b, std::size_t e) {
        const int C = T.C, n = T.n, nd = cedges.size(), np = path.size();
        std::vector<double> temp(C), u(C);
        for (std::size_t l = b; l < e; ++l) {
            const double* lp = logP[l].data(); const std::vector<double>& M = msg[l]; const std::vector<double>& K = contrib[l];
            std::vector<double>& S = scratch[l];
            S.resize((size_t)(nd + 2 * np) * C);
            double* sc = S.data();                       // nd * C new child contributions
            double* sm = sc + (size_t)nd * C;            // np * C new messages on the path
            double* su = sm + (size_t)np * C;            // np * C new up-contributions on the path (last unused)
            // v: new child contributions, new msg
            double* mv = sm;
            std::fill(mv, mv + C, 0.0);
            for (int i = 0; i < nd; ++i) {
                const int ce = cedges[i], c = T.ch[ce];
                child_message(lp, n, c, &M[(size_t)c * C], &cexp[(size_t)i * C * C], &crow[(size_t)i * C], C, sc + (size_t)i * C, temp.data(), u.data());
                for (int s = 0; s < C; ++s) mv[s] += sc[(size_t)i * C + s];
            }
            // walk up
            for (int k = 0; k + 1 < np; ++k) {
                const int node = path[k], pnode = path[k + 1], ue = T.up_edge[node];
                const double* ex = (k == 0) ? uexp.data() : &expA[(size_t)ue * C * C];
                const double* rm = (k == 0) ? urow.data() : &rowmax[(size_t)ue * C];
                child_message(lp, n, node, sm + (size_t)k * C, ex, rm, C, su + (size_t)k * C, temp.data(), u.data());
                double* mp = sm + (size_t)(k + 1) * C;
                for (int s = 0; s < C; ++s) mp[s] = M[(size_t)pnode * C + s] - K[(size_t)ue * C + s] + su[(size_t)k * C + s];
            }
            logz_new[l] = root_logz(lp, n, T.root, sm + (size_t)(np - 1) * C, C);
        }
    }
};

void all_edge_matrices(const BrlTree& T, const std::vector<double>& h, const std::vector<double>& grid,
                       const std::vector<std::vector<double>>& Ag, std::vector<double>& expA, std::vector<double>& rowmax) {
    const int C = T.C;
    expA.assign((size_t)T.m * C * C, 0.0); rowmax.assign((size_t)T.m * C, 0.0);
    for (int e = 0; e < T.m; ++e) edge_matrix(h[T.par[e]] - h[T.ch[e]], grid, Ag, C, &expA[(size_t)e * C * C], &rowmax[(size_t)e * C]);
}

} // namespace

//' Per-locus log-likelihood of a tree with per-edge branch lengths (internal)
//' @param E edge vector c(parents, children), 1-indexed postorder (reorder_phylo)
//' @param logP_list per-locus log-probability vectors (state-major), as convert_logliks_to_logP_list()
//' @param A_grid list of row-major transition matrices (probabilities), one per grid point
//' @param grid increasing generation values of A_grid
//' @param heights node heights (length n_tips + n_nodes; tips 0)
//' @keywords internal
// [[Rcpp::export]]
NumericVector score_tree_brl_cpp(arma::Col<int> E, const std::vector<std::vector<double>>& logP_list,
                                 const std::vector<std::vector<double>>& A_grid, const std::vector<double>& grid,
                                 const std::vector<double>& heights) {
    const int L = logP_list.size(); const int C = std::lround(std::sqrt((double)A_grid[0].size()));
    const int n = logP_list[0].size() / C;
    BrlTree T = build_tree(E, C, L, n);
    std::vector<double> expA, rowmax; all_edge_matrices(T, heights, grid, A_grid, expA, rowmax);
    std::vector<std::vector<double>> msg(L, std::vector<double>((size_t)n * C)), contrib(L, std::vector<double>((size_t)T.m * C));
    std::vector<double> logz(L);
    FullPassWorker w(T, logP_list, expA, rowmax, msg, contrib, logz);
    parallelFor(0, L, w);
    return wrap(logz);
}

//' Random-walk MCMC over node heights on a fixed topology (internal)
//' @param heights_init node heights (tips 0, root = total height, which stays fixed)
//' @param n_iter iterations (one node move each); thin: keep every thin-th state
//' @param step_sd sd of the logit-scale random-walk step
//' @param refresh full likelihood recomputation every refresh iterations (removes rounding drift)
//' @param p_scale probability of a root-height (global scaling) move per iteration; 0 keeps the root fixed
//' @param scale_sd sd of log(c) for the scaling move
//' @keywords internal
// [[Rcpp::export]]
List tree_brl_mcmc_cpp(arma::Col<int> E, const std::vector<std::vector<double>>& logP_list,
                       const std::vector<std::vector<double>>& A_grid, const std::vector<double>& grid,
                       std::vector<double> heights_init, int n_iter, double step_sd = 1.0, int thin = 1,
                       int seed = 1, int refresh = 1000, double p_scale = 0.0, double scale_sd = 0.1) {
    const int L = logP_list.size(); const int C = std::lround(std::sqrt((double)A_grid[0].size()));
    const int n = logP_list[0].size() / C;
    BrlTree T = build_tree(E, C, L, n);
    std::vector<double> h = heights_init;
    // validity
    for (int e = 0; e < T.m; ++e) if (h[T.par[e]] < h[T.ch[e]]) stop("heights_init: a parent is below its child");
    std::vector<int> movable;
    for (int v = 0; v < n; ++v) if (!T.child_edges[v].empty() && v != T.root) movable.push_back(v);
    if (movable.empty()) stop("no movable internal node");

    std::vector<double> expA, rowmax; all_edge_matrices(T, h, grid, A_grid, expA, rowmax);
    std::vector<std::vector<double>> msg(L, std::vector<double>((size_t)n * C)), contrib(L, std::vector<double>((size_t)T.m * C));
    std::vector<double> logz(L);
    { FullPassWorker w(T, logP_list, expA, rowmax, msg, contrib, logz); parallelFor(0, L, w); }
    double cur = 0.0; for (double z : logz) cur += z;

    std::mt19937_64 rng(seed); std::normal_distribution<double> nrm(0.0, step_sd); std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::uniform_int_distribution<int> pick(0, (int)movable.size() - 1);
    const int n_save = n_iter / thin;
    NumericMatrix trace(n_save, (int)movable.size());
    NumericVector ll_trace(n_save), root_trace(n_save); IntegerVector acc_trace(n_save);
    std::vector<std::vector<double>> scratch(L); std::vector<double> logz_new(L);
    long accepted = 0, n_scale = 0, acc_scale = 0; int s = 0;
    std::normal_distribution<double> nrm_s(0.0, scale_sd); const double hmax = grid.back();

    for (int it = 1; it <= n_iter; ++it) {
        bool acc = false;
        if (p_scale > 0 && unif(rng) < p_scale) {
            // root-height move: scale every height by c (relative positions kept); with a 1/H prior on the root height
            // and a log-scale random walk the acceptance ratio is the likelihood ratio
            ++n_scale;
            const double c = std::exp(nrm_s(rng));
            if (h[T.root] * c <= hmax) {
                std::vector<double> h2(h); for (double& x : h2) x *= c;
                std::vector<double> expA2, rowmax2; all_edge_matrices(T, h2, grid, A_grid, expA2, rowmax2);
                std::vector<std::vector<double>> msg2(L, std::vector<double>((size_t)n * C)), contrib2(L, std::vector<double>((size_t)T.m * C));
                std::vector<double> logz2(L);
                { FullPassWorker w(T, logP_list, expA2, rowmax2, msg2, contrib2, logz2); parallelFor(0, L, w); }
                double prop = 0.0; for (double z : logz2) prop += z;
                if (std::isfinite(prop) && std::log(unif(rng)) < prop - cur) {
                    acc = true; ++accepted; ++acc_scale; cur = prop;
                    h.swap(h2); expA.swap(expA2); rowmax.swap(rowmax2); msg.swap(msg2); contrib.swap(contrib2); logz.swap(logz2);
                }
            }
        } else {
            const int v = movable[pick(rng)];
            const int p = T.parent_node[v];
            double lb = 0.0; for (int ce : T.child_edges[v]) lb = std::max(lb, h[T.ch[ce]]);
            const double span = h[p] - lb;
            if (span > 1e-12) {
                double f = (h[p] - h[v]) / span; f = std::min(std::max(f, 1e-9), 1 - 1e-9);
                const double z = std::log(f / (1 - f)) + nrm(rng);
                const double f2 = 1.0 / (1.0 + std::exp(-z));
                const double hv2 = h[p] - f2 * span;
                // new matrices for v's child edges and up edge
                const std::vector<int>& ced = T.child_edges[v]; const int nd = ced.size();
                std::vector<double> cexp((size_t)nd * C * C), crow((size_t)nd * C), uexp((size_t)C * C), urow(C);
                for (int i = 0; i < nd; ++i) edge_matrix(hv2 - h[T.ch[ced[i]]], grid, A_grid, C, &cexp[(size_t)i * C * C], &crow[(size_t)i * C]);
                edge_matrix(h[p] - hv2, grid, A_grid, C, uexp.data(), urow.data());
                std::vector<int> path; for (int x = v; x != -1; x = T.parent_node[x]) path.push_back(x);
                MoveWorker w(T, logP_list, msg, contrib, expA, rowmax, ced, cexp, crow, uexp, urow, path, scratch, logz_new);
                parallelFor(0, L, w);
                double prop = 0.0; for (double zz : logz_new) prop += zz;
                const double log_r = prop - cur + std::log(f2 * (1 - f2)) - std::log(f * (1 - f));
                if (std::isfinite(prop) && std::log(unif(rng)) < log_r) {
                    acc = true; ++accepted; cur = prop; h[v] = hv2;
                    for (int i = 0; i < nd; ++i) {
                        std::copy(&cexp[(size_t)i * C * C], &cexp[(size_t)(i + 1) * C * C], &expA[(size_t)ced[i] * C * C]);
                        std::copy(&crow[(size_t)i * C], &crow[(size_t)(i + 1) * C], &rowmax[(size_t)ced[i] * C]);
                    }
                    const int ue = T.up_edge[v];
                    std::copy(uexp.begin(), uexp.end(), &expA[(size_t)ue * C * C]); std::copy(urow.begin(), urow.end(), &rowmax[(size_t)ue * C]);
                    const int np = path.size();
                    for (int l = 0; l < L; ++l) {
                        const double* sc = scratch[l].data(); const double* sm = sc + (size_t)nd * C; const double* su = sm + (size_t)np * C;
                        for (int i = 0; i < nd; ++i) std::copy(sc + (size_t)i * C, sc + (size_t)(i + 1) * C, &contrib[l][(size_t)ced[i] * C]);
                        for (int k = 0; k < np; ++k) std::copy(sm + (size_t)k * C, sm + (size_t)(k + 1) * C, &msg[l][(size_t)path[k] * C]);
                        for (int k = 0; k + 1 < np; ++k) std::copy(su + (size_t)k * C, su + (size_t)(k + 1) * C, &contrib[l][(size_t)T.up_edge[path[k]] * C]);
                        logz[l] = logz_new[l];
                    }
                }
            }
        }
        if (refresh > 0 && it % refresh == 0) {
            FullPassWorker w(T, logP_list, expA, rowmax, msg, contrib, logz); parallelFor(0, L, w);
            cur = 0.0; for (double z : logz) cur += z;
        }
        if (it % thin == 0 && s < n_save) {
            for (std::size_t j = 0; j < movable.size(); ++j) trace(s, j) = h[movable[j]];
            ll_trace[s] = cur; acc_trace[s] = acc; root_trace[s] = h[T.root]; ++s;
        }
        if (it % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    // final full recomputation for checking the incremental bookkeeping
    std::vector<double> logz_full(L);
    { FullPassWorker w(T, logP_list, expA, rowmax, msg, contrib, logz_full); parallelFor(0, L, w); }
    double full = 0.0; for (double z : logz_full) full += z;
    IntegerVector mv(movable.size()); for (std::size_t j = 0; j < movable.size(); ++j) mv[j] = movable[j] + 1;
    return List::create(_["heights"] = trace, _["logL"] = ll_trace, _["accepted"] = acc_trace, _["nodes"] = mv,
                        _["heights_final"] = wrap(h), _["acceptance"] = (double)accepted / n_iter,
                        _["logL_incremental"] = cur, _["logL_full"] = full, _["root_height"] = root_trace,
                        _["scale_moves"] = (double)n_scale, _["scale_acceptance"] = n_scale > 0 ? (double)acc_scale / n_scale : NA_REAL);
}

namespace {

// One node's height change, prepared and scored incrementally (same bookkeeping as the sampler's node move).
struct NodeMove {
    int v = -1; double hv2 = 0.0; int nd = 0;
    std::vector<double> cexp, crow, uexp, urow; std::vector<int> path;
};

void prepare_move(const BrlTree& T, const std::vector<double>& h, int v, double hv2, const std::vector<double>& grid,
                  const std::vector<std::vector<double>>& Ag, NodeMove& M) {
    const int C = T.C; M.v = v; M.hv2 = hv2;
    const std::vector<int>& ced = T.child_edges[v]; M.nd = ced.size();
    M.cexp.resize((size_t)M.nd * C * C); M.crow.resize((size_t)M.nd * C); M.uexp.resize((size_t)C * C); M.urow.resize(C);
    for (int i = 0; i < M.nd; ++i) edge_matrix(hv2 - h[T.ch[ced[i]]], grid, Ag, C, &M.cexp[(size_t)i * C * C], &M.crow[(size_t)i * C]);
    edge_matrix(h[T.parent_node[v]] - hv2, grid, Ag, C, M.uexp.data(), M.urow.data());
    M.path.clear(); for (int x = v; x != -1; x = T.parent_node[x]) M.path.push_back(x);
}

double score_move(const BrlTree& T, const std::vector<std::vector<double>>& logP, const std::vector<std::vector<double>>& msg,
                  const std::vector<std::vector<double>>& contrib, const std::vector<double>& expA, const std::vector<double>& rowmax,
                  const NodeMove& M, std::vector<std::vector<double>>& scratch, std::vector<double>& logz_new) {
    MoveWorker w(T, logP, msg, contrib, expA, rowmax, T.child_edges[M.v], M.cexp, M.crow, M.uexp, M.urow, M.path, scratch, logz_new);
    parallelFor(0, logP.size(), w);
    double s = 0.0; for (double z : logz_new) s += z; return s;
}

void commit_move(const BrlTree& T, std::vector<double>& h, std::vector<std::vector<double>>& msg, std::vector<std::vector<double>>& contrib,
                 std::vector<double>& expA, std::vector<double>& rowmax, std::vector<double>& logz, const NodeMove& M,
                 const std::vector<std::vector<double>>& scratch, const std::vector<double>& logz_new) {
    const int C = T.C, L = msg.size(), nd = M.nd, np = M.path.size(); const std::vector<int>& ced = T.child_edges[M.v];
    h[M.v] = M.hv2;
    for (int i = 0; i < nd; ++i) {
        std::copy(&M.cexp[(size_t)i * C * C], &M.cexp[(size_t)(i + 1) * C * C], &expA[(size_t)ced[i] * C * C]);
        std::copy(&M.crow[(size_t)i * C], &M.crow[(size_t)(i + 1) * C], &rowmax[(size_t)ced[i] * C]);
    }
    const int ue = T.up_edge[M.v];
    std::copy(M.uexp.begin(), M.uexp.end(), &expA[(size_t)ue * C * C]); std::copy(M.urow.begin(), M.urow.end(), &rowmax[(size_t)ue * C]);
    for (int l = 0; l < L; ++l) {
        const double* sc = scratch[l].data(); const double* sm = sc + (size_t)nd * C; const double* su = sm + (size_t)np * C;
        for (int i = 0; i < nd; ++i) std::copy(sc + (size_t)i * C, sc + (size_t)(i + 1) * C, &contrib[l][(size_t)ced[i] * C]);
        for (int k = 0; k < np; ++k) std::copy(sm + (size_t)k * C, sm + (size_t)(k + 1) * C, &msg[l][(size_t)M.path[k] * C]);
        for (int k = 0; k + 1 < np; ++k) std::copy(su + (size_t)k * C, su + (size_t)(k + 1) * C, &contrib[l][(size_t)T.up_edge[M.path[k]] * C]);
        logz[l] = logz_new[l];
    }
}

} // namespace

//' Maximum-likelihood node heights by coordinate ascent (internal)
//'
//' Sweeps over the internal non-root nodes; each node's position f between its highest child and its parent is
//' maximised by golden-section search on logit(f) (incremental path-to-root scoring per evaluation), as the original
//' estimate_brl_iter() did per node with Brent and full rescoring. After each sweep, if optimize_root, the global
//' scale (root height) is maximised by golden-section search on log(c) within [h_root / max_factor, grid max].
//' Stops when a sweep gains less than tol_ll log-likelihood units, or after max_sweeps.
//' @keywords internal
// [[Rcpp::export]]
List tree_brl_optimize_cpp(arma::Col<int> E, const std::vector<std::vector<double>>& logP_list,
                           const std::vector<std::vector<double>>& A_grid, const std::vector<double>& grid,
                           std::vector<double> heights_init, int max_sweeps = 50, double tol_ll = 0.01,
                           bool optimize_root = true, double max_factor = 50.0, int golden_iter = 30) {
    const int L = logP_list.size(); const int C = std::lround(std::sqrt((double)A_grid[0].size()));
    const int n = logP_list[0].size() / C;
    BrlTree T = build_tree(E, C, L, n);
    std::vector<double> h = heights_init;
    std::vector<int> movable;
    for (int v = 0; v < n; ++v) if (!T.child_edges[v].empty() && v != T.root) movable.push_back(v);
    std::vector<double> expA, rowmax; all_edge_matrices(T, h, grid, A_grid, expA, rowmax);
    std::vector<std::vector<double>> msg(L, std::vector<double>((size_t)n * C)), contrib(L, std::vector<double>((size_t)T.m * C));
    std::vector<double> logz(L);
    auto full = [&](const std::vector<double>& hh, std::vector<double>& ex, std::vector<double>& rm,
                    std::vector<std::vector<double>>& ms, std::vector<std::vector<double>>& ct, std::vector<double>& lz) {
        all_edge_matrices(T, hh, grid, A_grid, ex, rm);
        FullPassWorker w(T, logP_list, ex, rm, ms, ct, lz); parallelFor(0, L, w);
        double s = 0.0; for (double z : lz) s += z; return s;
    };
    double cur = full(h, expA, rowmax, msg, contrib, logz);
    std::vector<std::vector<double>> scratch(L); std::vector<double> logz_new(L);
    NodeMove M;
    const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
    std::vector<double> ll_sweep; ll_sweep.push_back(cur); std::vector<double> root_sweep; root_sweep.push_back(h[T.root]);
    int sweep = 0;
    for (sweep = 1; sweep <= max_sweeps; ++sweep) {
        const double start = cur;
        // postorder (children before parents) node sweep
        for (int v : movable) {
            const int p = T.parent_node[v];
            double lb = 0.0; for (int ce : T.child_edges[v]) lb = std::max(lb, h[T.ch[ce]]);
            const double span = h[p] - lb; if (span <= 1e-12) continue;
            auto eval = [&](double z) { const double f = 1.0 / (1.0 + std::exp(-z)); prepare_move(T, h, v, h[p] - f * span, grid, A_grid, M);
                                        return score_move(T, logP_list, msg, contrib, expA, rowmax, M, scratch, logz_new); };
            double a = -12.0, b = 12.0, x1 = b - gr * (b - a), x2 = a + gr * (b - a), f1 = eval(x1), f2 = eval(x2);
            for (int g = 0; g < golden_iter; ++g) {
                if (f1 > f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = eval(x1); }
                else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = eval(x2); }
            }
            const double zb = (f1 > f2) ? x1 : x2; const double fb = std::max(f1, f2);
            if (fb > cur) { eval(zb); commit_move(T, h, msg, contrib, expA, rowmax, logz, M, scratch, logz_new); cur = fb; }
        }
        if (optimize_root) {
            const double H0 = h[T.root], hmax = grid.back();
            auto evalc = [&](double lc) { std::vector<double> h2(h); for (double& x : h2) x *= std::exp(lc);
                std::vector<double> ex, rm, lz(L); std::vector<std::vector<double>> ms(L, std::vector<double>((size_t)n * C)), ct(L, std::vector<double>((size_t)T.m * C));
                return full(h2, ex, rm, ms, ct, lz); };
            double a = std::log(1.0 / max_factor), b = std::log(hmax / H0);
            if (b > a) {
                double x1 = b - gr * (b - a), x2 = a + gr * (b - a), f1 = evalc(x1), f2 = evalc(x2);
                for (int g = 0; g < golden_iter; ++g) {
                    if (f1 > f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = evalc(x1); }
                    else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = evalc(x2); }
                }
                const double lc = (f1 > f2) ? x1 : x2; const double fb = std::max(f1, f2);
                if (fb > cur) { for (double& x : h) x *= std::exp(lc); cur = full(h, expA, rowmax, msg, contrib, logz); }
            }
        }
        ll_sweep.push_back(cur); root_sweep.push_back(h[T.root]);
        Rcpp::checkUserInterrupt();
        if (cur - start < tol_ll) break;
    }
    std::vector<double> lz2(L); std::vector<double> ex2, rm2; std::vector<std::vector<double>> ms2(L, std::vector<double>((size_t)n * C)), ct2(L, std::vector<double>((size_t)T.m * C));
    const double check = full(h, ex2, rm2, ms2, ct2, lz2);
    return List::create(_["heights"] = wrap(h), _["logL"] = cur, _["logL_full"] = check, _["logL_sweeps"] = wrap(ll_sweep),
                        _["root_sweeps"] = wrap(root_sweep), _["sweeps"] = std::min(sweep, max_sweeps));
}

//' Log-likelihood and its gradient with respect to internal node heights (internal)
//'
//' dlogL/dh_v for every internal non-root node by central differences, each evaluated incrementally (only v's edges
//' and its path to the root are rescored; nothing is committed). Used by optimize_brl(method = "joint"), which chains
//' these into the gradient of the top-down (x_to_length) parameterisation.
//' @param rel_step step as a fraction of the node's admissible span (between its highest child and its parent)
//' @return list(logL, grad (named by 1-indexed node), nodes)
//' @keywords internal
// [[Rcpp::export]]
List brl_grad_heights_cpp(arma::Col<int> E, const std::vector<std::vector<double>>& logP_list,
                          const std::vector<std::vector<double>>& A_grid, const std::vector<double>& grid,
                          std::vector<double> heights, double rel_step = 1e-3) {
    const int L = logP_list.size(); const int C = std::lround(std::sqrt((double)A_grid[0].size()));
    const int n = logP_list[0].size() / C;
    BrlTree T = build_tree(E, C, L, n);
    std::vector<double> expA, rowmax; all_edge_matrices(T, heights, grid, A_grid, expA, rowmax);
    std::vector<std::vector<double>> msg(L, std::vector<double>((size_t)n * C)), contrib(L, std::vector<double>((size_t)T.m * C));
    std::vector<double> logz(L);
    { FullPassWorker w(T, logP_list, expA, rowmax, msg, contrib, logz); parallelFor(0, L, w); }
    double cur = 0.0; for (double z : logz) cur += z;
    std::vector<std::vector<double>> scratch(L); std::vector<double> logz_new(L);
    NodeMove M; std::vector<int> nodes; std::vector<double> grad;
    for (int v = 0; v < n; ++v) {
        if (T.child_edges[v].empty() || v == T.root) continue;
        const int p = T.parent_node[v];
        double lb = 0.0; for (int ce : T.child_edges[v]) lb = std::max(lb, heights[T.ch[ce]]);
        const double span = heights[p] - lb;
        double g = 0.0;
        if (span > 1e-12) {
            const double d = std::max(rel_step * span, 1e-9);
            const double hp = std::min(heights[v] + d, heights[p]), hm = std::max(heights[v] - d, lb);
            prepare_move(T, heights, v, hp, grid, A_grid, M); const double lp = score_move(T, logP_list, msg, contrib, expA, rowmax, M, scratch, logz_new);
            prepare_move(T, heights, v, hm, grid, A_grid, M); const double lm = score_move(T, logP_list, msg, contrib, expA, rowmax, M, scratch, logz_new);
            if (hp > hm) g = (lp - lm) / (hp - hm);
        }
        nodes.push_back(v + 1); grad.push_back(g);
    }
    NumericVector gr = wrap(grad); gr.attr("names") = wrap(nodes);
    return List::create(_["logL"] = cur, _["grad"] = gr, _["nodes"] = wrap(nodes));
}
