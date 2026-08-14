#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <limits>
#include <algorithm>
#include <cmath>
#include <RcppArmadillo.h>
#include <array>
#include <string>
#include <vector>
#include <sstream>
#include <numeric>
#include <mutex>
#include <memory>
#include <random>
#include <chrono>
#include <cstdint>
#include <cstdlib>

using namespace Rcpp;
using namespace RcppParallel;

namespace {

struct ScratchBuffers {
    double* temp;
    double* u;
    double* F_p2;
    double* F_p1;
    double* childF;
    double* prev_childF;
    double* F_v;
    double* s_full;
    double* max_t;
};

struct ThreadScratchBuffers {
    std::vector<double> temp;
    std::vector<double> u;
    std::vector<double> F_p2;
    std::vector<double> F_p1;
    std::vector<double> childF;
    std::vector<double> prev_childF;
    std::vector<double> F_v;
    std::vector<double> s_full;
    std::vector<double> max_t;

    void ensure(std::size_t CL, int L) {
        auto ensure_vec = [](std::vector<double>& vec, std::size_t target) {
            if (vec.size() < target) vec.resize(target);
        };
        ensure_vec(temp, CL);
        ensure_vec(u, CL);
        ensure_vec(F_p2, CL);
        ensure_vec(F_p1, CL);
        ensure_vec(childF, CL);
        ensure_vec(prev_childF, CL);
        ensure_vec(F_v, CL);
        ensure_vec(s_full, CL);
        ensure_vec(max_t, static_cast<std::size_t>(L));
    }

    ScratchBuffers view() {
        return ScratchBuffers{
            temp.data(), u.data(), F_p2.data(), F_p1.data(),
            childF.data(), prev_childF.data(), F_v.data(), s_full.data(), max_t.data()
        };
    }
};

inline ThreadScratchBuffers& thread_scratch_buffers() {
    thread_local ThreadScratchBuffers scratch;
    return scratch;
}

inline bool mitodrift_profile_enabled() {
    static bool enabled = []() {
        const char* env = std::getenv("MITODRIFT_PROFILE");
        if (!env) return false;
        if (env[0] == '\0') return false;
        return !(env[0] == '0' && env[1] == '\0');
    }();
    return enabled;
}

struct OptionalTimer {
    std::chrono::steady_clock::time_point start;
    long long* accum;

    explicit OptionalTimer(long long* target) : accum(target) {
        if (accum) start = std::chrono::steady_clock::now();
    }

    ~OptionalTimer() {
        if (!accum) return;
        const auto end = std::chrono::steady_clock::now();
        const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        *accum += static_cast<long long>(ns);
    }
};

} // namespace

/////////////////////////////////////// NNI ////////////////////////////////////////

// Inline helper performing the recursive postorder traversal.
// Parameters:
//   node, nTips: as before.
//   e1_ptr and e2_ptr: raw pointers to the parent's and child's arrays.
//   neworder: vector (of length nEdges) that will be filled with the new row ordering (stored 1-indexed).
//   L_ptr, xi_ptr, xj_ptr: auxiliary arrays computed for the internal nodes.
//   iii: current index to fill in neworder (initialized to nEdges-1).
inline void bar_reorderRcpp_inline(int node, int nTips,
    const int* e1_ptr, const int* e2_ptr,
    std::vector<int>& neworder,
    const int* L_ptr, const int* xi_ptr, const int* xj_ptr,
    int &iii) {
    
    int i = node - nTips - 1;
    // Output the current node's block (edges for this internal node) in reverse order.
    for (int j = xj_ptr[i] - 1; j >= 0; j--) {
        neworder[iii--] = L_ptr[xi_ptr[i] + j] + 1; // +1 converts to R's 1-indexing
    }
    // Recursively process each child.
    for (int j = 0; j < xj_ptr[i]; j++) {
        int k = e2_ptr[ L_ptr[xi_ptr[i] + j] ];
        if (k > nTips)
            bar_reorderRcpp_inline(k, nTips, e1_ptr, e2_ptr, neworder, L_ptr, xi_ptr, xj_ptr, iii);
    }
}

// Helper function to reorder the rows of a column-major edge matrix.
// E is a one-dimensional arma::Col<int> where the first nEdges entries are the parent column
// and the next nEdges entries are the child column.
// 'order' is a vector of 1-indexed row numbers in the desired order.
// The function returns a new arma::Col<int> with rows rearranged (still column-major).
arma::Col<int> reorder_rows(const arma::Col<int>& E, const std::vector<int>& order) {
    int nEdges = E.n_elem / 2;
    arma::Col<int> newE(E.n_elem);
    for (int i = 0; i < nEdges; i++) {
        int orig = order[i] - 1; // convert from 1-indexed to 0-indexed row number
        newE[i] = E[orig];             // parent's value (first column)
        newE[nEdges + i] = E[nEdges + orig];  // child's value (second column)
    }
    return newE;
}

// E is a one-dimensional vector representing an edge matrix in column-major order.
// [[Rcpp::export]]
arma::Col<int> reorderRcpp(const arma::Col<int>& E) {
    int nEdges = E.n_elem / 2;
    int nTips = nEdges / 2 + 1;
    int root = nTips + 1;
    
    arma::Col<int> parent = E.rows(0, nEdges - 1);
    arma::Col<int> child = E.rows(nEdges, E.n_elem - 1);
    
    // The maximum parent label tells us the total number of nodes.
    int m_val = parent.max();
    int nnode = m_val - nTips; // number of internal nodes

    // Allocate working arrays.
    std::vector<int> L(nEdges);            // Will store edge indices for each internal node.
    std::vector<int> neworder(nEdges);       // The final reordering (stored as 1-indexed row numbers).
    std::vector<int> pos(nnode, 0);          // Current fill position for each internal node.
    std::vector<int> xi(nnode, 0);           // Starting index for each internal node in L.
    std::vector<int> xj(nnode, 0);           // Count of children per internal node.

    // First pass: count children per internal node using parent's values.
    for (int i = 0; i < nEdges; i++) {
        int idx = parent[i] - nTips - 1;
        xj[idx]++;
    }
    
    // Compute starting positions xi as cumulative sums.
    for (int i = 1; i < nnode; i++) {
        xi[i] = xi[i - 1] + xj[i - 1];
    }
    
    // Fill L: For each edge, assign its row index to the appropriate block.
    for (int i = 0; i < nEdges; i++) {
        int k = parent[i] - nTips - 1;
        int j = pos[k];
        L[xi[k] + j] = i;
        pos[k]++;
    }
    
    // Reset the new order index.
    int iii = nEdges - 1;
    
    // Get raw pointers for fast access.
    const int* e1_ptr = parent.memptr();
    const int* e2_ptr = child.memptr();
    const int* L_ptr   = L.data();
    const int* xi_ptr  = xi.data();
    const int* xj_ptr  = xj.data();
    
    // Run the recursive postorder traversal.
    bar_reorderRcpp_inline(root, nTips, e1_ptr, e2_ptr, neworder, L_ptr, xi_ptr, xj_ptr, iii);
    
    // Use the computed new order to reorder the rows of E.
    arma::Col<int> newE = reorder_rows(E, neworder);
    
    return newE;
}

// [[Rcpp::export]]
std::vector<arma::Col<int>> nnin_cpp(const arma::Col<int>& E, const int n) {

    const int numEdges = E.n_elem / 2;
    const int nTips = numEdges / 2 + 1;
    
    // Create subviews for parent's and child's columns without copying.
    arma::subview_col<int> parent = E.subvec(0, numEdges - 1);
    arma::subview_col<int> child  = E.subvec(numEdges, E.n_elem - 1);
    
    // Find internal edges (child > nTips) using vectorized operations.
    arma::uvec internalEdges = arma::find(child > nTips);
    if (internalEdges.n_elem < (unsigned int)n)
        stop("n is larger than the number of valid internal edges.");
    
    // Select the nth internal edge (0-indexed)
    const int ind = internalEdges(n - 1);
    
    // Retrieve parent's value (p1) and child's value (p2) for the chosen edge.
    const int p1 = parent(ind);
    const int p2 = child(ind);
    
    // Find indices where parent equals p1.
    arma::uvec indices_p1 = arma::find(parent == p1);
    // Choose the index that is not equal to 'ind'
    const int ind1 = (indices_p1(0) == (unsigned int) ind) ? indices_p1(1) : indices_p1(0);
    
    // Find indices where parent equals p2.
    arma::uvec indices_p2 = arma::find(parent == p2);
    const int ind2_0 = indices_p2(0);
    const int ind2_1 = indices_p2(1);
    
    // Retrieve the child values for these swap candidates.
    const int e1_val = child(ind1);
    const int e2_val = child(ind2_0);
    const int e3_val = child(ind2_1);
    
    // Create copies for the two alternative topologies.
    arma::Col<int> E1 = E;  // copy of E
    arma::Col<int> E2 = E;  // another copy
    
    // Topology 1: swap child at ind1 with that at ind2_0.
    E1(numEdges + ind1) = e2_val;
    E1(numEdges + ind2_0) = e1_val;
    
    // Topology 2: swap child at ind1 with that at ind2_1.
    E2(numEdges + ind1) = e3_val;
    E2(numEdges + ind2_1) = e1_val;
    
    // Reorder each topology (e.g., into postorder) using the helper reorderRcpp.    
    std::vector<arma::Col<int>> res(2);
    res[0] = reorderRcpp(E1);
    res[1] = reorderRcpp(E2);
    return res;
}

/////////////////////////////////////// ML Tree Search ////////////////////////////////////////

//' definitions for logSumExp function
// https://github.com/helske/seqHMM/blob/master/src/logSumExp.cpp
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

//' logSumExp function for a vector
//'
//' @param x NumericVector
//' @return double logSumExp of x
//' @keywords internal
// [[Rcpp::export]]
double logSumExp(const arma::vec& x) {
    unsigned int maxi = x.index_max();
    LDOUBLE maxv = x(maxi);
    if (!(maxv > -arma::datum::inf)) {
        return -arma::datum::inf;
    }
    LDOUBLE cumsum = 1.0; // Include the max element's contribution (exp(0)=1)
    for (unsigned int i = 0; i < maxi; i++) {
        cumsum += EXPL(x(i) - maxv);
    }
    for (unsigned int i = maxi + 1; i < x.n_elem; i++) {
        cumsum += EXPL(x(i) - maxv);
    }
    return maxv + std::log(cumsum);
}


static inline double logsumexp_array(const double* x, int len) {
	double maxv = -std::numeric_limits<double>::infinity();
	for (int i = 0; i < len; ++i) {
		if (x[i] > maxv) maxv = x[i];
	}
	if (!(maxv > -std::numeric_limits<double>::infinity())) {
		return -std::numeric_limits<double>::infinity();
	}
	double sum = 0.0;
	for (int i = 0; i < len; ++i) {
		sum += std::exp(x[i] - maxv);
	}
	return maxv + std::log(sum);
}

// bp: Belief-propagation function.
// logP is a flattened likelihood matrix (row-major; dimensions: C x n)
// logA is a flattened transition matrix (dimensions: C x C)
// E is a flattened edge list in postorder (each edge: parent, child) stored in column-major order.
// n: number of nodes, C: number of states, m: number of edges, root: index of the root node.
// [[Rcpp::export]]
double score_tree_bp(const arma::Col<int>& E,
                     const std::vector<double>& logP, 
                     const std::vector<double>& logA,
                     const int n, const int C, const int m, const int root) {
    // Allocate memory for messages and temporary state values.
    std::vector<double> log_messages(C * n, 0.0);
    std::vector<double> state_log_values(C);
    std::vector<double> temp(C);
    int idx;

    for (int i = 0; i < m; i++) {
        int par  = E(i);         // parent's value for edge i
        int node = E(m + i);     // child's value for edge i

        // Precompute common expression for each c_child.
        for (int c_child = 0; c_child < C; c_child++) {
            idx = c_child * n + node;
            temp[c_child] = logP[idx] + log_messages[idx];
        }
        
        for (int c = 0; c < C; c++) {
            for (int c_child = 0; c_child < C; c_child++) {
                state_log_values[c_child] = logA[c * C + c_child] + temp[c_child];
            }
            log_messages[c * n + par] += logSumExp(state_log_values);
        }
    }
  
    for (int c = 0; c < C; c++) {
        idx = c * n + root;
        state_log_values[c] = logP[idx] + log_messages[idx];
    }
    return logSumExp(state_log_values);
}

// Core routine for belief propagation with precomputed exp-shifted A and row maxes.
// Messages are stored in node-major layout: msg[node * C + c]
static inline double score_tree_bp2_core(const int* e_ptr, int m, int n, int C, int root,
                                         const double* logP_ptr,
                                         const double* expA_shifted,   // size C*C, row-major
                                         const double* row_maxA,       // size C
                                         double* msg_nm,               // size n*C, node-major
                                         double* temp,                 // size C (scratch)
                                         double* u) {                  // size C (scratch)
	for (int i = 0; i < m; ++i) {
		const int par  = e_ptr[i];
		const int node = e_ptr[m + i];

		// Build temp[c] = logP[c*n + node] + msg_nm[node*C + c], and find its max for stability.
		double max_t = -std::numeric_limits<double>::infinity();
		int offP = node;
		const int base_child = node * C;
		for (int c = 0; c < C; ++c, offP += n) {
			const double v = logP_ptr[offP] + msg_nm[base_child + c];
			temp[c] = v;
			if (v > max_t) max_t = v;
		}
		// u[c] = exp(temp[c] - max_t)  (shared across all rows of A)
		for (int c = 0; c < C; ++c) u[c] = std::exp(temp[c] - max_t);

		// For each parent state c, accumulate log-sum-exp as:
		// row_maxA[c] + max_t + log( dot( expA_shifted_row_c , u ) )
		const int base_par = par * C;
		const double* rowA = expA_shifted;
		for (int c = 0; c < C; ++c, rowA += C) {
			double s = 0.0;
			// dot product between rowA (length C) and u (length C)
			for (int j = 0; j < C; ++j) s += rowA[j] * u[j];
			msg_nm[base_par + c] += row_maxA[c] + max_t + std::log(s);
		}
	}

	// Aggregate at root over states: logsumexp_c( logP[c*n + root] + msg_nm[root*C + c] )
	int offP = root;
	const int base_root = root * C;
	for (int c = 0; c < C; ++c, offP += n) {
		temp[c] = logP_ptr[offP] + msg_nm[base_root + c];
	}
	return logsumexp_array(temp, C);
}

// bp: Belief-propagation function.
// bp2: Belief-propagation function.
// logP_list: List of flattened likelihood matrices (each in row-major order)
// logA: Flattened transition matrix (row-major order)
// [[Rcpp::export]]
double score_tree_bp_wrapper2(arma::Col<int> E,
                              const std::vector< std::vector<double> >& logP_list,
                              const std::vector<double>& logA) {
	const int L = static_cast<int>(logP_list.size());
	const int C = static_cast<int>(std::sqrt(logA.size()));
	const int n = static_cast<int>(logP_list[0].size() / C);
	const int m = static_cast<int>(E.n_elem / 2);

	// 1-indexed (R) -> 0-indexed (C++)
	E -= 1;
	const int root = E(m - 1);
	const int* e_ptr = E.memptr();

	// Precompute exp-shifted rows of A and per-row max once (shared across loci).
	std::vector<double> row_maxA(C);
	std::vector<double> expA_shifted(static_cast<size_t>(C) * C);
	for (int r = 0; r < C; ++r) {
		const double* Arow = logA.data() + static_cast<size_t>(r) * C;
		double mr = Arow[0];
		for (int j = 1; j < C; ++j) if (Arow[j] > mr) mr = Arow[j];
		row_maxA[r] = mr;
		double* out = expA_shifted.data() + static_cast<size_t>(r) * C;
		for (int j = 0; j < C; ++j) out[j] = std::exp(Arow[j] - mr);
	}

	// Reusable buffers
	std::vector<double> msg_nm(static_cast<size_t>(n) * C);
	std::vector<double> temp(C);
	std::vector<double> u(C);

	double logZ = 0.0;
	for (int l = 0; l < L; ++l) {
		std::fill(msg_nm.begin(), msg_nm.end(), 0.0);
		logZ += score_tree_bp2_core(e_ptr, m, n, C, root,
		                            logP_list[l].data(),
		                            expA_shifted.data(), row_maxA.data(),
		                            msg_nm.data(), temp.data(), u.data());
	}
	return logZ;
}


/* ===========================
 * Message-caching NNI engine
 * ===========================
 *
 * We cache, for every locus ℓ and node v, the length‑C vector F_ℓ[v,·] which is the
 * child→parent contribution added to the parent's message (i.e., the quantity we
 * accumulate into msg_nm[parent,·] inside score_tree_bp2_core). For a parent p with two
 * children a and b, msg_nm[p,·] = F[a,·] + F[b,·]. An NNI around edge (p1→p2) only changes
 * the sets of children of p1 and p2, so we can recompute F at p2, then p1, and then walk
 * upward to the root, updating F along that single path. This yields O(height × C × L)
 * updates per proposal instead of O(n × C × L).
 */

struct NNICache {
	// static tree info
	int n;					// #nodes
	int m;					// #edges
    int C;				// #states
    int L;				// #loci
    std::size_t CL;			// states × loci per node
    int nTips;			// #tips
	int root;				// root node id (0-indexed)
	arma::Col<int> E;		// edge list, 0-indexed, postorder (parent col [0..m-1], child col [m..2m-1])


	// topology helpers
	std::vector<int> parent_of;					// size n
	std::vector< std::array<int,2> > children_of;		// for internal nodes: exactly two children
    std::vector<int> internal_edge_indices;       // internal-edge ordinal -> position in E
	// transition precompute
	std::vector<double> row_maxA;				// size C
    arma::Mat<double> expA_shifted_t;    // stores exp(A) rows as columns (C x C)

	// data likelihoods
	std::vector< std::vector<double> > logP_list;	// L × (C*n), row-major per locus

	// per-locus cached child→parent contributions F (n*C) and root logZ
    // F layout: [node][c][l] -> node * (C*L) + c*L + l
	std::vector<double> F;		// n * C * L
	std::vector<double> logZ;	// L
    std::vector<double> logP_storage; // n * C * L
    double current_total_loglik_val;

    // Scratch space for calculations to avoid reallocations
    // Size C * L
    mutable std::vector<double> scratch_temp;
    mutable std::vector<double> scratch_u;
    mutable std::vector<double> scratch_F_p2;
    mutable std::vector<double> scratch_F_p1;
    mutable std::vector<double> scratch_childF;
    mutable std::vector<double> scratch_prev_childF;
    mutable std::vector<double> scratch_F_v;
    mutable std::vector<double> scratch_s_full;

    // Size L
    mutable std::vector<double> scratch_max_t;
    mutable std::mutex scratch_mutex;

    inline ScratchBuffers shared_scratch() const {
        return ScratchBuffers{
            scratch_temp.data(),
            scratch_u.data(),
            scratch_F_p2.data(),
            scratch_F_p1.data(),
            scratch_childF.data(),
            scratch_prev_childF.data(),
            scratch_F_v.data(),
            scratch_s_full.data(),
            scratch_max_t.data()
        };
    }

    inline ScratchBuffers acquire_scratch(bool stage_mode) const {
        if (stage_mode) {
            return shared_scratch();
        }
        auto& local = thread_scratch_buffers();
        local.ensure(CL, L);
        return local.view();
    }

    mutable std::vector<int> staged_nodes;
    mutable std::vector<double> staged_F;
    mutable bool staged_ready = false;
    mutable int staged_edge_n = -1;
    mutable int staged_edge_index = -1;
    mutable int staged_which = -1;
    mutable int staged_p1 = -1;
    mutable int staged_p2 = -1;
    mutable int staged_c1 = -1;
    mutable int staged_cX = -1;
    mutable int staged_cStay = -1;
    mutable double staged_total_loglik = -std::numeric_limits<double>::infinity();

    inline void reset_staged_state_unlocked() const {
        staged_ready = false;
        staged_nodes.clear();
        staged_F.clear();
        staged_edge_n = -1;
        staged_edge_index = -1;
        staged_which = -1;
        staged_p1 = staged_p2 = staged_c1 = staged_cX = staged_cStay = -1;
        staged_total_loglik = -std::numeric_limits<double>::infinity();
    }

    inline void stage_node_F(int node_id, const double* src) const {
        staged_nodes.push_back(node_id);
        staged_F.insert(staged_F.end(), src, src + CL);
    }

	NNICache(arma::Col<int> E_in,
		const std::vector< std::vector<double> >& logP_in,
		const std::vector<double>& logA_in,
        bool reorder = true) {

        L = static_cast<int>(logP_in.size());
        C = static_cast<int>(std::sqrt(logA_in.size()));
        n = static_cast<int>(logP_in[0].size() / C);
        CL = static_cast<std::size_t>(C) * static_cast<std::size_t>(L);

        // Initialize scratch space
        scratch_temp.resize(CL);
        scratch_u.resize(CL);
        scratch_F_p2.resize(CL);
        scratch_F_p1.resize(CL);
        scratch_childF.resize(CL);
        scratch_prev_childF.resize(CL);
        scratch_F_v.resize(CL);
        scratch_s_full.resize(CL);
		
        scratch_max_t.resize(L);
        reset_staged_state_unlocked();

		// Bring E to postorder and 0-indexed
        if (reorder) {
		    E_in = reorderRcpp(E_in);
        }
		E = E_in - 1;

		m = static_cast<int>(E.n_elem / 2);
		nTips = m / 2 + 1;
		root = E(m - 1);

		// parent/children
		parent_of.assign(n, -1);
		children_of.assign(n, std::array<int,2>{-1,-1});
		for (int i = 0; i < m; ++i) {
			const int p = E[i];
			const int c = E[m + i];
			parent_of[c] = p;
			// fill two slots
			if (children_of[p][0] == -1) children_of[p][0] = c;
			else children_of[p][1] = c;
		}
        rebuild_internal_edge_indices();

		// A precompute
        row_maxA.resize(C);
        expA_shifted_t.set_size(C, C);
        for (int r = 0; r < C; ++r) {
            const double* Arow = logA_in.data() + static_cast<size_t>(r) * C;
            double mr = Arow[0];
            for (int j = 1; j < C; ++j) if (Arow[j] > mr) mr = Arow[j];
            row_maxA[r] = mr;
            for (int j = 0; j < C; ++j) {
                expA_shifted_t(j, r) = std::exp(Arow[j] - mr);
            }
        }

		// Transpose logP
        logP_storage.resize(n * C * L);
        for (int l = 0; l < L; ++l) {
            const auto& P_l = logP_in[l];
            for (int node = 0; node < n; ++node) {
                for (int c = 0; c < C; ++c) {
                    // Old: P_l[c * n + node]
                    // New: logP_storage[node * C * L + c * L + l]
                    logP_storage[node * C * L + c * L + l] = P_l[c * n + node];
                }
            }
        }

		// allocate caches
		F.assign(n * C * L, 0.0);
		logZ.assign(L, 0.0);

        std::vector<double> msg_nm(n * C * L, 0.0); // Accumulator for children messages
        
        // Scratch buffers for initialization
        std::vector<double> init_temp(C * L);
        std::vector<double> init_u(C * L);
        std::vector<double> init_max_t(L);
        std::vector<double> init_s_full(C * L);
        
        for (int i = 0; i < m; ++i) {
            const int par  = E[i];
            const int node = E[m + i];
            
            // Compute F[node]
            // temp = P[node] + msg_nm[node]
            
            const double* P_ptr = logP_storage.data() + node * C * L;
            const double* msg_ptr = msg_nm.data() + node * C * L;
            double* F_ptr = F.data() + node * C * L;
            double* msg_par_ptr = msg_nm.data() + par * C * L;
            
            std::fill(init_max_t.begin(), init_max_t.end(), -std::numeric_limits<double>::infinity());
            
            for (int c = 0; c < C; ++c) {
                const double* p = P_ptr + c * L;
                const double* m_ = msg_ptr + c * L;
                double* t = init_temp.data() + c * L;
                for (int l = 0; l < L; ++l) {
                    double val = p[l] + m_[l];
                    t[l] = val;
                    if (val > init_max_t[l]) init_max_t[l] = val;
                }
            }
            
            for (int c = 0; c < C; ++c) {
                double* t = init_temp.data() + c * L;
                double* u = init_u.data() + c * L;
                for (int l = 0; l < L; ++l) {
                    if (init_max_t[l] == -std::numeric_limits<double>::infinity()) {
                        u[l] = 0.0;
                    } else {
                        u[l] = std::exp(t[l] - init_max_t[l]);
                    }
                }
            }
            
            arma::Mat<double> init_U(init_u.data(), L, C, false, true);
            arma::Mat<double> init_S(init_s_full.data(), L, C, false, true);
            init_S = init_U * expA_shifted_t;
            const double* s_ptr = init_S.memptr();
            for (int r = 0; r < C; ++r) {
                double row_max = row_maxA[r];
                const double* s_col = s_ptr + static_cast<size_t>(r) * L;
                double* f = F_ptr + r * L;
                double* mp = msg_par_ptr + r * L;
                for (int l = 0; l < L; ++l) {
                    if (init_max_t[l] == -std::numeric_limits<double>::infinity() || s_col[l] <= 0.0) {
                        f[l] = -std::numeric_limits<double>::infinity();
                        mp[l] += -std::numeric_limits<double>::infinity();
                    } else {
                        double val = row_max + init_max_t[l] + std::log(s_col[l]);
                        f[l] = val;
                        mp[l] += val;
                    }
                }
            }
        }
        
        // Root marginal
        const double* P_root = logP_storage.data() + root * C * L;
        const double* msg_root = msg_nm.data() + root * C * L;
        
        std::fill(init_max_t.begin(), init_max_t.end(), -std::numeric_limits<double>::infinity());
        
        for (int c = 0; c < C; ++c) {
            const double* p = P_root + c * L;
            const double* m_ = msg_root + c * L;
            double* t = init_temp.data() + c * L;
            for (int l = 0; l < L; ++l) {
                double val = p[l] + m_[l];
                t[l] = val;
                if (val > init_max_t[l]) init_max_t[l] = val;
            }
        }
        
        current_total_loglik_val = 0.0;
        for (int l = 0; l < L; ++l) {
            if (init_max_t[l] == -std::numeric_limits<double>::infinity()) {
                logZ[l] = -std::numeric_limits<double>::infinity();
            } else {
                double sum_exp = 0.0;
                for (int c = 0; c < C; ++c) {
                    sum_exp += std::exp(init_temp[c * L + l] - init_max_t[l]);
                }
                logZ[l] = init_max_t[l] + std::log(sum_exp);
            }
            current_total_loglik_val += logZ[l];
        }
	}

	inline int other_child(int parent, int child) const {
		const auto &ch = children_of[parent];
		return (ch[0] == child) ? ch[1] : ch[0];
	}

    inline void rebuild_internal_edge_indices() {
        internal_edge_indices.clear();
        internal_edge_indices.reserve(static_cast<std::size_t>(nTips - 2));
        for (int i = 0; i < m; ++i) {
            if (E[m + i] >= nTips) internal_edge_indices.push_back(i);
        }
    }

    inline int locate_internal_edge_index(int edge_n) const {
        if (edge_n < 1 || static_cast<std::size_t>(edge_n) > internal_edge_indices.size()) return -1;
        return internal_edge_indices[static_cast<std::size_t>(edge_n - 1)];
    }

    inline void compute_F_vectorized(
        const double* F_c1, const double* F_c2,
        const double* P_base,
        double* outF,
        const ScratchBuffers& scratch
    ) const {
        double* u = scratch.u;
        double* temp = scratch.temp;
        double* max_t = scratch.max_t;
        double* s_full = scratch.s_full;
        const double neg_inf = -std::numeric_limits<double>::infinity();
        std::fill(max_t, max_t + L, neg_inf);
		
        for (int c = 0; c < C; ++c) {
            const double* p_ptr = P_base + c * L;
            const double* f1_ptr = F_c1 + c * L;
            const double* f2_ptr = F_c2 + c * L;
            double* t_ptr = temp + c * L;
            for (int l = 0; l < L; ++l) {
                double val = p_ptr[l] + f1_ptr[l] + f2_ptr[l];
                t_ptr[l] = val;
                if (val > max_t[l]) max_t[l] = val;
            }
        }
		
        for (int c = 0; c < C; ++c) {
            double* t_ptr = temp + c * L;
            double* u_ptr = u + c * L;
            for (int l = 0; l < L; ++l) {
                u_ptr[l] = (max_t[l] == neg_inf) ? 0.0 : std::exp(t_ptr[l] - max_t[l]);
            }
        }
		
        arma::Mat<double> U_view(u, L, C, false, true);
        arma::Mat<double> S_view(s_full, L, C, false, true);
        S_view = U_view * expA_shifted_t;
        const double* s_ptr = S_view.memptr();
		
        for (int r = 0; r < C; ++r) {
            double row_max = row_maxA[r];
            const double* s_col = s_ptr + static_cast<size_t>(r) * L;
            double* out_ptr = outF + r * L;
            for (int l = 0; l < L; ++l) {
                if (max_t[l] == neg_inf || s_col[l] <= 0.0) {
                    out_ptr[l] = neg_inf;
                } else {
                    out_ptr[l] = row_max + max_t[l] + std::log(s_col[l]);
                }
            }
        }
    }

    inline double compute_root_logZ_vectorized(
        const double* F_c1, const double* F_c2,
        const double* P_base,
        const ScratchBuffers& scratch
    ) const {
        double* temp = scratch.temp;
        double* max_t = scratch.max_t;
        std::fill(max_t, max_t + L, -std::numeric_limits<double>::infinity());
        
        for (int c = 0; c < C; ++c) {
            const double* p_ptr = P_base + c * L;
            const double* f1_ptr = F_c1 + c * L;
            const double* f2_ptr = F_c2 + c * L;
            double* t_ptr = temp + c * L;
            
            for (int l = 0; l < L; ++l) {
                double val = p_ptr[l] + f1_ptr[l] + f2_ptr[l];
                t_ptr[l] = val;
                if (val > max_t[l]) max_t[l] = val;
            }
        }
        
        double total_logZ = 0.0;
        for (int l = 0; l < L; ++l) {
            if (max_t[l] == -std::numeric_limits<double>::infinity()) {
                total_logZ += -std::numeric_limits<double>::infinity();
            } else {
                double sum_exp = 0.0;
                for (int c = 0; c < C; ++c) {
                    sum_exp += std::exp(temp[c * L + l] - max_t[l]);
                }
                total_logZ += max_t[l] + std::log(sum_exp);
            }
        }
        if (std::isnan(total_logZ)) return -std::numeric_limits<double>::infinity();
        return total_logZ;
    }

    double compute_new_loglik_impl(int edge_n, int which, bool stage_results, const ScratchBuffers& scratch) const {
        const int ind = locate_internal_edge_index(edge_n);
        if (ind < 0) stop("edge_n out of range in compute_new_loglik");

        const int p1 = E[ind];
        const int p2 = E[m + ind];
        const int c1 = other_child(p1, p2);
        const int c2 = children_of[p2][0];
        const int c3 = children_of[p2][1];

        const int cX    = (which == 0 ? c2 : c3);
        const int cStay = (which == 0 ? c3 : c2);

        if (stage_results) {
            staged_edge_n = edge_n;
            staged_edge_index = ind;
            staged_which = which;
            staged_p1 = p1;
            staged_p2 = p2;
            staged_c1 = c1;
            staged_cX = cX;
            staged_cStay = cStay;
        }

        double* scratch_F_p2_ptr = scratch.F_p2;
        double* scratch_F_p1_ptr = scratch.F_p1;

        compute_F_vectorized(
            F.data() + c1 * C * L,
            F.data() + cStay * C * L,
            logP_storage.data() + p2 * C * L,
            scratch_F_p2_ptr,
            scratch
        );
        if (stage_results) stage_node_F(p2, scratch_F_p2_ptr);

        compute_F_vectorized(
            scratch_F_p2_ptr,
            F.data() + cX * C * L,
            logP_storage.data() + p1 * C * L,
            scratch_F_p1_ptr,
            scratch
        );
        if (stage_results) stage_node_F(p1, scratch_F_p1_ptr);

        if (parent_of[p1] == -1) {
            double new_total = compute_root_logZ_vectorized(
                scratch_F_p2_ptr,
                F.data() + cX * C * L,
                logP_storage.data() + p1 * C * L,
                scratch
            );
            if (stage_results) {
                staged_total_loglik = new_total;
                staged_ready = true;
            }
            return new_total;
        }

        int child = p1;
        double* p_childF = scratch.childF;
        double* p_prev_childF = scratch.prev_childF;
        double* p_F_v = scratch.F_v;

        std::copy(scratch_F_p1_ptr, scratch_F_p1_ptr + CL, p_childF);

        int prev_child = -1;

        while (true) {
            const int v = parent_of[child];
            if (v == -1) {
                const int root_node = child;
                const int sib = other_child(root_node, prev_child);
                double new_total = compute_root_logZ_vectorized(
                    p_prev_childF,
                    F.data() + sib * C * L,
                    logP_storage.data() + root_node * C * L,
                    scratch
                );
                if (stage_results) {
                    staged_total_loglik = new_total;
                    staged_ready = true;
                }
                return new_total;
            }

            const int sib = other_child(v, child);
            compute_F_vectorized(
                p_childF,
                F.data() + sib * C * L,
                logP_storage.data() + v * C * L,
                p_F_v,
                scratch
            );
            if (stage_results) stage_node_F(v, p_F_v);

            prev_child = child;
            double* temp = p_prev_childF;
            p_prev_childF = p_childF;
            p_childF = p_F_v;
            p_F_v = temp;
            child = v;
        }
    }

    // Compute new total loglik if we perform the NNI "which" (0 or 1) at the nth internal edge.
    double compute_new_loglik(int edge_n, int which, bool stage_results = false) const {
        if (stage_results) {
            std::lock_guard<std::mutex> guard(scratch_mutex);
            reset_staged_state_unlocked();
            return compute_new_loglik_impl(edge_n, which, true, shared_scratch());
        }
        return compute_new_loglik_impl(edge_n, which, false, acquire_scratch(false));
    }

    void commit_staged_nni() {
        std::lock_guard<std::mutex> guard(scratch_mutex);
        if (!staged_ready) stop("No staged NNI proposal to apply");
        const std::size_t block = CL;
        double* F_ptr = F.data();
        const double* staged_ptr = staged_F.data();
        for (std::size_t idx = 0; idx < staged_nodes.size(); ++idx) {
            const int node_id = staged_nodes[idx];
            double* dest = F_ptr + static_cast<std::size_t>(node_id) * block;
            std::copy(staged_ptr, staged_ptr + block, dest);
            staged_ptr += block;
        }
        current_total_loglik_val = staged_total_loglik;

        const int p1 = staged_p1;
        const int p2 = staged_p2;
        const int c1 = staged_c1;
        const int cX = staged_cX;
        const int cStay = staged_cStay;

        auto &ch1 = children_of[p1];
        if (ch1[0] == c1) ch1[0] = cX; else ch1[1] = cX;
        auto &ch2 = children_of[p2];
        if (ch2[0] == cX) ch2[0] = c1; else ch2[1] = c1;
        parent_of[c1] = p2;
        parent_of[cX] = p1;

        for (int i = 0; i < m; ++i) {
            if (E[i] == p1 && E[m + i] == c1) { E[m + i] = cX; break; }
        }
        for (int i = 0; i < m; ++i) {
            if (E[i] == p2 && E[m + i] == cX) { E[m + i] = c1; break; }
        }

        arma::Col<int> E1 = E + 1;
        E1 = reorderRcpp(E1);
        E = E1 - 1;
        root = E(m - 1);
        rebuild_internal_edge_indices();

        reset_staged_state_unlocked();
    }

    void discard_staged_nni() const {
        std::lock_guard<std::mutex> guard(scratch_mutex);
        reset_staged_state_unlocked();
    }

    // Commit the NNI: update topology and cached F/logZ across loci.
    void apply_nni(int edge_n, int which) {
        compute_new_loglik(edge_n, which, true);
        commit_staged_nni();
    }

	double total_loglik() const {
		return current_total_loglik_val;
	}
};

// ---- Rcpp exported wrappers around NNICache ----

SEXP nni_cache_create(arma::Col<int> E,
	const std::vector< std::vector<double> >& logP,
	const std::vector<double>& logA) {
	Rcpp::XPtr<NNICache> ptr(new NNICache(E, logP, logA), true);
	return ptr;
}

double nni_cache_loglik(SEXP xp) {
	Rcpp::XPtr<NNICache> ptr(xp);
	return ptr->total_loglik();
}

double nni_cache_delta(SEXP xp, int edge_n, int which) {
	Rcpp::XPtr<NNICache> ptr(xp);
    double new_ll = ptr->compute_new_loglik(edge_n, which);
    double cur_ll = ptr->total_loglik();
    if (cur_ll == -std::numeric_limits<double>::infinity()) {
        if (new_ll > -std::numeric_limits<double>::infinity()) return std::numeric_limits<double>::infinity();
        else return 0.0;
    }
	return new_ll - cur_ll;
}

void nni_cache_apply(SEXP xp, int edge_n, int which) {
	Rcpp::XPtr<NNICache> ptr(xp);
	ptr->apply_nni(edge_n, which);
}

arma::Col<int> nni_cache_current_E(SEXP xp) {
	Rcpp::XPtr<NNICache> ptr(xp);
	return ptr->E + 1; // back to 1-indexed
}

// --------------------------------------------------------------------------
// Worker struct that scores multiple trees in parallel.
// Each worker processes a subset of trees from the input list.
struct ScoreTreesWorker : public Worker {
    // Inputs are passed by const reference.
    const std::vector<arma::Col<int>>& trees;
    const std::vector< std::vector<double> >& logP;
    const std::vector<double>& logA;
    
    // Output: scores for each tree
    RVector<double> scores;
    
    // Constructor.
    ScoreTreesWorker(const std::vector<arma::Col<int>>& trees,
                     const std::vector< std::vector<double> >& logP,
                     const std::vector<double>& logA,
                     NumericVector scores)
      : trees(trees), logP(logP), logA(logA), scores(scores) {}
    
    // Operator() for processing tree indices [begin, end)
    void operator()(std::size_t begin, std::size_t end) {
        // Each iteration scores one tree.
        for (std::size_t i = begin; i < end; i++) {
            scores[i] = score_tree_bp_wrapper2(trees[i], logP, logA);
        }
    }
};

// --------------------------------------------------------------------------
// Parallel wrapper: scores multiple trees in parallel.
// Parameters:
//   trees: vector of edge matrices (each arma::Col<int>) in column-major order.
//   logP: list of flattened likelihood matrices (each in row-major order).
//   logA: flattened transition matrix (row-major order).
// Returns a vector of scores (one per tree).
// [[Rcpp::export]]
NumericVector score_trees_parallel(const std::vector<arma::Col<int>>& trees,
                                   const std::vector< std::vector<double> >& logP,
                                   const std::vector<double>& logA) {
    
    int n_trees = trees.size();
    
    // Prepare the output vector to hold scores for all trees.
    NumericVector scores(n_trees);
    
    // Create the worker that will score trees in parallel.
    ScoreTreesWorker worker(trees, logP, logA, scores);
    
    // Launch parallelFor over tree indices [0, n_trees).
    parallelFor(0, n_trees, worker);
    
    return scores;
}


// E: Edge matrix (each row is a (parent, child) pair, 1-indexed from R) provided as an arma::Col<int> in column-major order.
// logP_list: List of flattened likelihood matrices (each in row-major order)
// logA: Flattened transition matrix (row-major order)
// Computes:
//   - L: number of loci (length of logP_list),
//   - C: number of states (inferred from logA),
//   - n: number of nodes (inferred from the first likelihood matrix),
//   - m: number of edges (number of rows in the original edge matrix),
//   - root: parent's value from the last edge.
// [[Rcpp::export]]
double score_tree_bp_wrapper_multi(arma::Col<int> E,
                             const std::vector< std::vector<double> >& logP_list,
                             const std::vector< std::vector<double> >& logA_list) {
    // Compute number of loci.
    int L = logP_list.size();
    // Infer number of states from the transition matrix.
    int C = std::sqrt(logA_list[0].size());
    // Infer number of nodes from the first likelihood matrix.
    int n = logP_list[0].size() / C;
    // Compute m: number of edges.
    int m = E.n_elem / 2;

    // Adjust the edge matrix from 1-indexing (R) to 0-indexing (C++).
    E = E - 1;
    // Find the root
    int root = E(m - 1);

    double logZ = 0.0;
    for (int l = 0; l < L; l++) {
        logZ += score_tree_bp(E, logP_list[l], logA_list[l], n, C, m, root);
    }
    return logZ;
}


struct score_neighbours: public Worker {

    // original tree
    const arma::Col<int> E;
    const std::vector<std::vector<double>> logP;
    const std::vector<double> logA;
    RVector<double> scores;

    // initialize with source and destination
    score_neighbours(const arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<double> logA, NumericVector scores): 
        E(E), logP(logP), logA(logA), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::vector<arma::Col<int>> Ep = nnin_cpp(E, i+1);
            scores[2*i] = score_tree_bp_wrapper2(Ep[0], logP, logA);
            scores[2*i+1] = score_tree_bp_wrapper2(Ep[1], logP, logA);
        }
    }
};

// Cached variant: efficiently scores NNI neighbours using NNICache deltas.
struct score_neighbours_cached: public Worker {
	const NNICache* cache; // read-only cache (shared across threads)
	double base_ll;        // current total log-likelihood
	RVector<double> scores;

	score_neighbours_cached(const NNICache* cache, NumericVector scores)
		: cache(cache), base_ll(cache->total_loglik()), scores(scores) {}

	void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; ++i) {
			const int edge_n = static_cast<int>(i) + 1; // 1-indexed internal-edge id
            
            double ll0 = cache->compute_new_loglik(edge_n, 0);
            double ll1 = cache->compute_new_loglik(edge_n, 1);
            
			scores[2*i]   = ll0;
			scores[2*i+1] = ll1;
		}
	}
};

// [[Rcpp::export]]
NumericVector nni_cpp_parallel_cached(arma::Col<int> E,
	const std::vector<std::vector<double>>& logP,
	const std::vector<double>& logA) {
	// Build a fresh cache for this tree
	NNICache cache(E, logP, logA);

	// Number of internal edges (same formula as the original function)
	int n = E.n_elem / 4 - 1;
	NumericVector scores(2 * n);

	// Parallel evaluation of all 2*n neighbors using cached deltas
	score_neighbours_cached worker(&cache, scores);
	parallelFor(0, n, worker);
	return scores;
}

struct score_neighbours_multi: public Worker {

    // original tree
    const arma::Col<int> E;
    const std::vector<std::vector<double>> logP;
    const std::vector<std::vector<double>> logA;
    RVector<double> scores;

    // initialize with source and destination
    score_neighbours_multi(const arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<std::vector<double>> logA, NumericVector scores): 
        E(E), logP(logP), logA(logA), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::vector<arma::Col<int>> Ep = nnin_cpp(E, i+1);
            scores[2*i] = score_tree_bp_wrapper_multi(Ep[0], logP, logA);
            scores[2*i+1] = score_tree_bp_wrapper_multi(Ep[1], logP, logA);
        }
    }
};

// [[Rcpp::export]]
NumericVector nni_cpp_parallel_multi(arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<std::vector<double>> logA) {

    E = reorderRcpp(E);

    int n = E.n_elem / 4 - 1;

    NumericVector scores(2*n);

    score_neighbours_multi score_neighbours_multi(E, logP, logA, scores);

    parallelFor(0, n, score_neighbours_multi);

    return scores;

}

/////////////////////////////////////// MCMC ////////////////////////////////////////

// Thread-safe variant: no R objects, no XPtr, safe inside RcppParallel workers.
// [[Rcpp::export]]
std::vector<arma::Col<int>> tree_mcmc_cpp_cached_threadsafe(
	arma::Col<int> E,
	const std::vector< std::vector<double> >& logP,
	const std::vector<double>& logA,
    int max_iter = 100, int seed = -1, bool reorder = true) {

	// Number of internal edges (unchanged by reordering)
	const int n = static_cast<int>(E.n_elem / 4) - 1;

	// Build cache locally (constructor reorders and 0-indexes internally)
	NNICache cache(E, logP, logA, reorder);

    std::vector<arma::Col<int>> tree_list(static_cast<size_t>(max_iter) + 1);

    // Starting log-likelihood and tree
    double l_0 = cache.total_loglik();
    tree_list[0] = cache.E + 1; // back to 1-indexed

	// RNG
	std::mt19937 gen;
	if (seed == -1) {
		std::random_device rd;
		gen.seed(rd());
	} else {
		gen.seed(seed);
	}
	std::uniform_int_distribution<> dis1(1, n);
	std::uniform_int_distribution<> dis2(0, 1);
	std::uniform_real_distribution<> dis3(0.0, 1.0);

    for (int i = 0; i < max_iter; ++i) {
		// propose: pick internal edge and which swap (0/1)
		const int r1 = dis1(gen);
		const int r2 = dis2(gen);

        // local log-likelihood delta using cached messages
        double new_ll = cache.compute_new_loglik(r1, r2, true);
        double dl;
        if (l_0 == -std::numeric_limits<double>::infinity()) {
             if (new_ll > -std::numeric_limits<double>::infinity()) dl = std::numeric_limits<double>::infinity();
             else dl = 0.0;
        } else {
             dl = new_ll - l_0;
        }

		// accept using log form (stable)
        if (std::log(dis3(gen)) < dl) {
            cache.commit_staged_nni();
            l_0 = new_ll;
		} else {
			cache.discard_staged_nni();
		}

        // store current tree (1-indexed)
        tree_list[static_cast<size_t>(i) + 1] = cache.E + 1;
	}

	return tree_list;
}

struct SeededTreeChainWorker : public Worker {
    const std::vector< arma::Col<int> >& start_edges;
    const std::vector< std::vector<double> >& logP;
    const std::vector<double>& logA;
    const std::vector<int>& max_iter_vec;
    const std::vector<int>& seeds;
    std::vector< std::vector<arma::Col<int>> >& chain_results;

    SeededTreeChainWorker(const std::vector< arma::Col<int> >& start_edges,
                      const std::vector< std::vector<double> >& logP,
                      const std::vector<double>& logA,
                      const std::vector<int>& max_iter_vec,
                      const std::vector<int>& seeds,
                      std::vector< std::vector<arma::Col<int>> >& chain_results)
        : start_edges(start_edges), logP(logP), logA(logA),
          max_iter_vec(max_iter_vec), seeds(seeds), chain_results(chain_results) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            chain_results[i] = tree_mcmc_cpp_cached_threadsafe(start_edges[i], logP, logA, max_iter_vec[i], seeds[i], false);
        }
    }
};

// [[Rcpp::export]]
std::vector< std::vector<arma::Col<int>> > tree_mcmc_parallel_seeded(std::vector< arma::Col<int> > start_edges,
    const std::vector< std::vector<double> >& logP,
    const std::vector<double>& logA,
    const std::vector<int>& max_iter_vec,
    const std::vector<int>& seeds) {

    const std::size_t nchains = start_edges.size();
    if (max_iter_vec.size() != nchains || seeds.size() != nchains) {
        stop("start_edges, max_iter_vec, and seeds must have the same length");
    }

    for (std::size_t i = 0; i < nchains; ++i) {
        start_edges[i] = reorderRcpp(start_edges[i]);
    }

    std::vector< std::vector<arma::Col<int>> > chain_results(nchains);

    SeededTreeChainWorker worker(start_edges, logP, logA, max_iter_vec, seeds, chain_results);
    parallelFor(0, static_cast<std::size_t>(nchains), worker);

    return chain_results;
}

struct MC3ChainResult {
    std::vector<arma::Col<int>> cold_trace;
    std::vector<arma::Col<int>> final_states;
    std::vector<int> swap_attempts;
    std::vector<int> swap_accepts;
};

MC3ChainResult tree_mc3_cpp_cached_threadsafe(
    const std::vector<arma::Col<int>>& start_edges,
    const std::vector<std::vector<double>>& logP,
    const std::vector<double>& logA,
    int max_iter,
    int seed,
    const std::vector<double>& temperatures,
    int swap_interval) {

    const std::size_t ntemps = temperatures.size();
    if (ntemps == 0 || start_edges.size() != ntemps) {
        stop("start_edges and temperatures must have the same positive length");
    }
    if (max_iter < 0) stop("max_iter must be non-negative");
    if (swap_interval < 1) stop("swap_interval must be positive");
    if (std::abs(temperatures[0] - 1.0) > 1e-12) {
        stop("The first MC3 temperature must be 1");
    }
    for (std::size_t t = 0; t < ntemps; ++t) {
        if (!std::isfinite(temperatures[t]) || temperatures[t] < 1.0 ||
            (t > 0 && temperatures[t] <= temperatures[t - 1])) {
            stop("MC3 temperatures must be finite, start at 1, and strictly increase");
        }
    }

    const int n = static_cast<int>(start_edges[0].n_elem / 4) - 1;
    if (n < 1) stop("MC3 requires a tree with at least one internal edge");

    std::vector<std::unique_ptr<NNICache>> caches;
    caches.reserve(ntemps);
    std::vector<double> loglik(ntemps);
    for (std::size_t t = 0; t < ntemps; ++t) {
        caches.emplace_back(new NNICache(start_edges[t], logP, logA, false));
        loglik[t] = caches[t]->total_loglik();
    }

    MC3ChainResult result;
    result.cold_trace.resize(static_cast<std::size_t>(max_iter) + 1);
    result.cold_trace[0] = caches[0]->E + 1;
    result.swap_attempts.assign(ntemps > 0 ? ntemps - 1 : 0, 0);
    result.swap_accepts.assign(ntemps > 0 ? ntemps - 1 : 0, 0);

    std::mt19937 gen;
    if (seed == -1) {
        std::random_device rd;
        gen.seed(rd());
    } else {
        gen.seed(seed);
    }
    std::uniform_int_distribution<> edge_dist(1, n);
    std::uniform_int_distribution<> move_dist(0, 1);
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    std::uniform_int_distribution<> pair_dist(0, ntemps > 1 ? static_cast<int>(ntemps) - 2 : 0);

    for (int iter = 0; iter < max_iter; ++iter) {
        for (std::size_t t = 0; t < ntemps; ++t) {
            const int edge_n = edge_dist(gen);
            const int which = move_dist(gen);
            const double new_ll = caches[t]->compute_new_loglik(edge_n, which, true);
            double delta;
            if (loglik[t] == -std::numeric_limits<double>::infinity()) {
                delta = new_ll > -std::numeric_limits<double>::infinity()
                    ? std::numeric_limits<double>::infinity() : 0.0;
            } else {
                delta = (new_ll - loglik[t]) / temperatures[t];
            }
            if (std::log(uniform(gen)) < delta) {
                caches[t]->commit_staged_nni();
                loglik[t] = new_ll;
            } else {
                caches[t]->discard_staged_nni();
            }
        }

        if (ntemps > 1 && (iter + 1) % swap_interval == 0) {
            const int pair = pair_dist(gen);
            const int next = pair + 1;
            const double beta_i = 1.0 / temperatures[static_cast<std::size_t>(pair)];
            const double beta_j = 1.0 / temperatures[static_cast<std::size_t>(next)];
            const double log_alpha = (beta_i - beta_j) *
                (loglik[static_cast<std::size_t>(next)] - loglik[static_cast<std::size_t>(pair)]);
            result.swap_attempts[static_cast<std::size_t>(pair)]++;
            if (std::log(uniform(gen)) < log_alpha) {
                std::swap(caches[static_cast<std::size_t>(pair)], caches[static_cast<std::size_t>(next)]);
                std::swap(loglik[static_cast<std::size_t>(pair)], loglik[static_cast<std::size_t>(next)]);
                result.swap_accepts[static_cast<std::size_t>(pair)]++;
            }
        }
        result.cold_trace[static_cast<std::size_t>(iter) + 1] = caches[0]->E + 1;
    }

    result.final_states.resize(ntemps);
    for (std::size_t t = 0; t < ntemps; ++t) {
        result.final_states[t] = caches[t]->E + 1;
    }
    return result;
}

struct SeededMC3Worker : public Worker {
    const std::vector<std::vector<arma::Col<int>>>& start_edges;
    const std::vector<std::vector<double>>& logP;
    const std::vector<double>& logA;
    const std::vector<int>& max_iter_vec;
    const std::vector<int>& seeds;
    const std::vector<double>& temperatures;
    int swap_interval;
    std::vector<MC3ChainResult>& results;

    SeededMC3Worker(
        const std::vector<std::vector<arma::Col<int>>>& start_edges,
        const std::vector<std::vector<double>>& logP,
        const std::vector<double>& logA,
        const std::vector<int>& max_iter_vec,
        const std::vector<int>& seeds,
        const std::vector<double>& temperatures,
        int swap_interval,
        std::vector<MC3ChainResult>& results)
        : start_edges(start_edges), logP(logP), logA(logA), max_iter_vec(max_iter_vec),
          seeds(seeds), temperatures(temperatures), swap_interval(swap_interval), results(results) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            results[i] = tree_mc3_cpp_cached_threadsafe(
                start_edges[i], logP, logA, max_iter_vec[i], seeds[i], temperatures, swap_interval);
        }
    }
};

// [[Rcpp::export]]
List tree_mc3_parallel_seeded_serial(
    std::vector<std::vector<arma::Col<int>>> start_edges,
    const std::vector<std::vector<double>>& logP,
    const std::vector<double>& logA,
    const std::vector<int>& max_iter_vec,
    const std::vector<int>& seeds,
    const std::vector<double>& temperatures,
    int swap_interval = 10) {

    const std::size_t nchains = start_edges.size();
    if (nchains == 0 || max_iter_vec.size() != nchains || seeds.size() != nchains) {
        stop("start_edges, max_iter_vec, and seeds must have the same positive length");
    }
    if (temperatures.empty() || std::abs(temperatures[0] - 1.0) > 1e-12) {
        stop("MC3 temperatures must start at 1");
    }
    if (swap_interval < 1) stop("swap_interval must be positive");
    for (std::size_t t = 0; t < temperatures.size(); ++t) {
        if (!std::isfinite(temperatures[t]) || temperatures[t] < 1.0 ||
            (t > 0 && temperatures[t] <= temperatures[t - 1])) {
            stop("MC3 temperatures must be finite, start at 1, and strictly increase");
        }
    }
    for (std::size_t i = 0; i < nchains; ++i) {
        if (max_iter_vec[i] < 0) stop("max_iter_vec values must be non-negative");
        if (start_edges[i].size() != temperatures.size()) {
            stop("Each MC3 ensemble must have one start state per temperature");
        }
        for (std::size_t t = 0; t < start_edges[i].size(); ++t) {
            start_edges[i][t] = reorderRcpp(start_edges[i][t]);
        }
    }

    std::vector<MC3ChainResult> results(nchains);
    SeededMC3Worker worker(
        start_edges, logP, logA, max_iter_vec, seeds, temperatures, swap_interval, results);
    parallelFor(0, nchains, worker);

    List traces(nchains);
    List final_states(nchains);
    IntegerMatrix attempts(nchains, temperatures.size() > 0 ? temperatures.size() - 1 : 0);
    IntegerMatrix accepts(nchains, temperatures.size() > 0 ? temperatures.size() - 1 : 0);
    for (std::size_t i = 0; i < nchains; ++i) {
        traces[i] = wrap(results[i].cold_trace);
        final_states[i] = wrap(results[i].final_states);
        for (std::size_t p = 0; p < results[i].swap_attempts.size(); ++p) {
            attempts(i, p) = results[i].swap_attempts[p];
            accepts(i, p) = results[i].swap_accepts[p];
        }
    }
    return List::create(
        _["traces"] = traces,
        _["final_states"] = final_states,
        _["swap_attempts"] = attempts,
        _["swap_accepts"] = accepts,
        _["temperatures"] = temperatures);
}

struct MC3CacheInitWorker : public Worker {
    const std::vector<std::vector<arma::Col<int>>>& start_edges;
    const std::vector<std::vector<double>>& logP;
    const std::vector<double>& logA;
    std::size_t ntemps;
    std::vector<std::unique_ptr<NNICache>>& caches;
    std::vector<double>& loglik;

    MC3CacheInitWorker(
        const std::vector<std::vector<arma::Col<int>>>& start_edges,
        const std::vector<std::vector<double>>& logP,
        const std::vector<double>& logA,
        std::size_t ntemps,
        std::vector<std::unique_ptr<NNICache>>& caches,
        std::vector<double>& loglik)
        : start_edges(start_edges), logP(logP), logA(logA), ntemps(ntemps),
          caches(caches), loglik(loglik) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t task = begin; task < end; ++task) {
            const std::size_t ensemble = task / ntemps;
            const std::size_t temperature = task % ntemps;
            caches[task].reset(new NNICache(
                start_edges[ensemble][temperature], logP, logA, false));
            loglik[task] = caches[task]->total_loglik();
        }
    }
};

struct MC3UpdateWorker : public Worker {
    std::size_t ntemps;
    int n_internal_edges;
    int block_start;
    int block_end;
    const std::vector<int>& max_iter_vec;
    const std::vector<double>& temperatures;
    std::vector<std::unique_ptr<NNICache>>& caches;
    std::vector<double>& loglik;
    std::vector<std::mt19937>& rngs;
    std::vector<std::vector<arma::Col<int>>>& cold_traces;

    MC3UpdateWorker(
        std::size_t ntemps,
        int n_internal_edges,
        int block_start,
        int block_end,
        const std::vector<int>& max_iter_vec,
        const std::vector<double>& temperatures,
        std::vector<std::unique_ptr<NNICache>>& caches,
        std::vector<double>& loglik,
        std::vector<std::mt19937>& rngs,
        std::vector<std::vector<arma::Col<int>>>& cold_traces)
        : ntemps(ntemps), n_internal_edges(n_internal_edges),
          block_start(block_start), block_end(block_end),
          max_iter_vec(max_iter_vec), temperatures(temperatures), caches(caches),
          loglik(loglik), rngs(rngs), cold_traces(cold_traces) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t task = begin; task < end; ++task) {
            const std::size_t ensemble = task / ntemps;
            const std::size_t temperature = task % ntemps;
            const int final_iter = std::min(block_end, max_iter_vec[ensemble]);
            if (block_start >= final_iter) continue;

            std::mt19937& gen = rngs[task];
            std::uniform_int_distribution<> edge_dist(1, n_internal_edges);
            std::uniform_int_distribution<> move_dist(0, 1);
            std::uniform_real_distribution<> uniform(0.0, 1.0);
            for (int iter = block_start; iter < final_iter; ++iter) {
                const int edge_n = edge_dist(gen);
                const int which = move_dist(gen);
                const double new_ll = caches[task]->compute_new_loglik(edge_n, which, true);
                double delta;
                if (loglik[task] == -std::numeric_limits<double>::infinity()) {
                    delta = new_ll > -std::numeric_limits<double>::infinity()
                        ? std::numeric_limits<double>::infinity() : 0.0;
                } else {
                    delta = (new_ll - loglik[task]) / temperatures[temperature];
                }
                if (std::log(uniform(gen)) < delta) {
                    caches[task]->commit_staged_nni();
                    loglik[task] = new_ll;
                } else {
                    caches[task]->discard_staged_nni();
                }
                if (temperature == 0) {
                    cold_traces[ensemble][static_cast<std::size_t>(iter) + 1] =
                        caches[task]->E + 1;
                }
            }
        }
    }
};

struct MC3SwapWorker : public Worker {
    std::size_t ntemps;
    int boundary_iter;
    const std::vector<int>& max_iter_vec;
    const std::vector<double>& temperatures;
    std::vector<std::unique_ptr<NNICache>>& caches;
    std::vector<double>& loglik;
    std::vector<std::mt19937>& swap_rngs;
    std::vector<std::vector<arma::Col<int>>>& cold_traces;
    std::vector<int>& swap_attempts;
    std::vector<int>& swap_accepts;

    MC3SwapWorker(
        std::size_t ntemps,
        int boundary_iter,
        const std::vector<int>& max_iter_vec,
        const std::vector<double>& temperatures,
        std::vector<std::unique_ptr<NNICache>>& caches,
        std::vector<double>& loglik,
        std::vector<std::mt19937>& swap_rngs,
        std::vector<std::vector<arma::Col<int>>>& cold_traces,
        std::vector<int>& swap_attempts,
        std::vector<int>& swap_accepts)
        : ntemps(ntemps), boundary_iter(boundary_iter), max_iter_vec(max_iter_vec),
          temperatures(temperatures), caches(caches), loglik(loglik),
          swap_rngs(swap_rngs), cold_traces(cold_traces),
          swap_attempts(swap_attempts), swap_accepts(swap_accepts) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t ensemble = begin; ensemble < end; ++ensemble) {
            if (max_iter_vec[ensemble] < boundary_iter) continue;

            std::mt19937& gen = swap_rngs[ensemble];
            std::uniform_int_distribution<> pair_dist(
                0, static_cast<int>(ntemps) - 2);
            std::uniform_real_distribution<> uniform(0.0, 1.0);
            const int pair = pair_dist(gen);
            const int next = pair + 1;
            const std::size_t first_task = ensemble * ntemps + static_cast<std::size_t>(pair);
            const std::size_t second_task = first_task + 1;
            const double beta_i = 1.0 / temperatures[static_cast<std::size_t>(pair)];
            const double beta_j = 1.0 / temperatures[static_cast<std::size_t>(next)];
            const double log_alpha = (beta_i - beta_j) *
                (loglik[second_task] - loglik[first_task]);
            const std::size_t stat_index = ensemble * (ntemps - 1) +
                static_cast<std::size_t>(pair);
            swap_attempts[stat_index]++;
            if (std::log(uniform(gen)) < log_alpha) {
                std::swap(caches[first_task], caches[second_task]);
                std::swap(loglik[first_task], loglik[second_task]);
                swap_accepts[stat_index]++;
            }
            // The serial implementation records the cold state after a swap.
            cold_traces[ensemble][static_cast<std::size_t>(boundary_iter)] =
                caches[ensemble * ntemps]->E + 1;
        }
    }
};

// [[Rcpp::export]]
List tree_mc3_parallel_seeded(
    std::vector<std::vector<arma::Col<int>>> start_edges,
    const std::vector<std::vector<double>>& logP,
    const std::vector<double>& logA,
    const std::vector<int>& max_iter_vec,
    const std::vector<int>& seeds,
    const std::vector<double>& temperatures,
    int swap_interval = 10) {

    const std::size_t nchains = start_edges.size();
    const std::size_t ntemps = temperatures.size();
    if (nchains == 0 || max_iter_vec.size() != nchains || seeds.size() != nchains) {
        stop("start_edges, max_iter_vec, and seeds must have the same positive length");
    }
    if (ntemps == 0 || std::abs(temperatures[0] - 1.0) > 1e-12) {
        stop("MC3 temperatures must start at 1");
    }
    if (swap_interval < 1) stop("swap_interval must be positive");
    for (std::size_t temperature = 0; temperature < ntemps; ++temperature) {
        if (!std::isfinite(temperatures[temperature]) || temperatures[temperature] < 1.0 ||
            (temperature > 0 && temperatures[temperature] <= temperatures[temperature - 1])) {
            stop("MC3 temperatures must be finite, start at 1, and strictly increase");
        }
    }
    for (std::size_t ensemble = 0; ensemble < nchains; ++ensemble) {
        if (max_iter_vec[ensemble] < 0) stop("max_iter_vec values must be non-negative");
        if (start_edges[ensemble].size() != ntemps) {
            stop("Each MC3 ensemble must have one start state per temperature");
        }
    }

    // A one-temperature ensemble is ordinary MCMC, including its exact seed
    // stream and output, rather than a special coupled-sampler approximation.
    if (ntemps == 1) {
        std::vector<arma::Col<int>> ordinary_starts(nchains);
        for (std::size_t ensemble = 0; ensemble < nchains; ++ensemble) {
            ordinary_starts[ensemble] = start_edges[ensemble][0];
        }
        const std::vector<std::vector<arma::Col<int>>> ordinary_traces =
            tree_mcmc_parallel_seeded(ordinary_starts, logP, logA, max_iter_vec, seeds);
        List traces(nchains);
        List final_states(nchains);
        for (std::size_t ensemble = 0; ensemble < nchains; ++ensemble) {
            traces[ensemble] = wrap(ordinary_traces[ensemble]);
            final_states[ensemble] = List::create(wrap(ordinary_traces[ensemble].back()));
        }
        return List::create(
            _["traces"] = traces,
            _["final_states"] = final_states,
            _["swap_attempts"] = IntegerMatrix(nchains, 0),
            _["swap_accepts"] = IntegerMatrix(nchains, 0),
            _["temperatures"] = temperatures);
    }

    for (std::size_t ensemble = 0; ensemble < nchains; ++ensemble) {
        for (std::size_t temperature = 0; temperature < ntemps; ++temperature) {
            start_edges[ensemble][temperature] = reorderRcpp(start_edges[ensemble][temperature]);
        }
    }
    const int n_internal_edges = static_cast<int>(start_edges[0][0].n_elem / 4) - 1;
    if (n_internal_edges < 1) stop("MC3 requires a tree with at least one internal edge");

    const std::size_t ntasks = nchains * ntemps;
    std::vector<std::unique_ptr<NNICache>> caches(ntasks);
    std::vector<double> loglik(ntasks);
    MC3CacheInitWorker init_worker(start_edges, logP, logA, ntemps, caches, loglik);
    parallelFor(0, ntasks, init_worker);

    std::vector<std::mt19937> rngs(ntasks);
    std::vector<std::mt19937> swap_rngs(nchains);
    std::random_device random_device;
    for (std::size_t ensemble = 0; ensemble < nchains; ++ensemble) {
        const std::uint32_t base_seed = seeds[ensemble] == -1
            ? random_device() : static_cast<std::uint32_t>(seeds[ensemble]);
        for (std::size_t temperature = 0; temperature < ntemps; ++temperature) {
            std::seed_seq stream_seed{
                base_seed,
                static_cast<std::uint32_t>(ensemble),
                static_cast<std::uint32_t>(temperature),
                UINT32_C(0x4d43334d)};
            rngs[ensemble * ntemps + temperature].seed(stream_seed);
        }
        std::seed_seq swap_seed{
            base_seed,
            static_cast<std::uint32_t>(ensemble),
            static_cast<std::uint32_t>(ntemps),
            UINT32_C(0x53574150)};
        swap_rngs[ensemble].seed(swap_seed);
    }

    std::vector<std::vector<arma::Col<int>>> cold_traces(nchains);
    int longest_chain = 0;
    for (std::size_t ensemble = 0; ensemble < nchains; ++ensemble) {
        cold_traces[ensemble].resize(static_cast<std::size_t>(max_iter_vec[ensemble]) + 1);
        cold_traces[ensemble][0] = caches[ensemble * ntemps]->E + 1;
        longest_chain = std::max(longest_chain, max_iter_vec[ensemble]);
    }
    std::vector<int> swap_attempts(nchains * (ntemps - 1), 0);
    std::vector<int> swap_accepts(nchains * (ntemps - 1), 0);

    for (int block_start = 0; block_start < longest_chain; block_start += swap_interval) {
        const int block_end = std::min(block_start + swap_interval, longest_chain);
        MC3UpdateWorker update_worker(
            ntemps, n_internal_edges, block_start, block_end, max_iter_vec,
            temperatures, caches, loglik, rngs, cold_traces);
        // Synchronous return from parallelFor is the pre-swap barrier.
        parallelFor(0, ntasks, update_worker);

        const int swap_boundary = block_start + swap_interval;
        if (swap_boundary <= longest_chain) {
            MC3SwapWorker swap_worker(
                ntemps, swap_boundary, max_iter_vec, temperatures, caches, loglik,
                swap_rngs, cold_traces, swap_attempts, swap_accepts);
            parallelFor(0, nchains, swap_worker);
        }
    }

    List traces(nchains);
    List final_states(nchains);
    IntegerMatrix attempts(nchains, ntemps - 1);
    IntegerMatrix accepts(nchains, ntemps - 1);
    for (std::size_t ensemble = 0; ensemble < nchains; ++ensemble) {
        traces[ensemble] = wrap(cold_traces[ensemble]);
        List ensemble_states(ntemps);
        for (std::size_t temperature = 0; temperature < ntemps; ++temperature) {
            ensemble_states[temperature] = wrap(caches[ensemble * ntemps + temperature]->E + 1);
        }
        final_states[ensemble] = ensemble_states;
        for (std::size_t pair = 0; pair < ntemps - 1; ++pair) {
            const std::size_t stat_index = ensemble * (ntemps - 1) + pair;
            attempts(ensemble, pair) = swap_attempts[stat_index];
            accepts(ensemble, pair) = swap_accepts[stat_index];
        }
    }
    return List::create(
        _["traces"] = traces,
        _["final_states"] = final_states,
        _["swap_attempts"] = attempts,
        _["swap_accepts"] = accepts,
        _["temperatures"] = temperatures);
}
