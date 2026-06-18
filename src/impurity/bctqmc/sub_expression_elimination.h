#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <string>
//#include <chrono>

using namespace std;
//using namespace std::chrono;


// ===================== Problem model =====================

enum class Kind : uint8_t { C=0, R=1, I=2, H=3 };

struct Token {
    Kind kind;
    int  idx;   // for C/R/I: global Gi index; for H: temp id
};

static inline string tok2str(const Token& t) {
    if (t.kind == Kind::H) return "H" + to_string(t.idx);
    char c = (t.kind==Kind::C ? 'C' : (t.kind==Kind::R ? 'R' : 'I'));
    return string(1, c) + to_string(t.idx);
}

struct Occ {
    int prod, start, end;  // inclusive [start..end]
};

struct TempInfo {
    int id;
    vector<Token> sig;     // the block (may include H tokens)
    vector<Occ>   occs;    // accepted occurrences
};

struct Plan {
    vector<TempInfo> temps;
    vector<vector<Token>> reduced;
};

struct CSEPlanner {
    vector<vector<Token>> products;  // input products
    Plan run();                      // planner
};

// ===================== Small signature (fast hashing) =====================

// Pack token into 32 bits: [2 bits kind | 30 bits idx]
static inline uint32_t packTok(Kind k, uint32_t idx) {
    return (uint32_t(idx) << 2) | (uint32_t(k) & 0x3u);
}

// Small fixed-capacity signature for short sequences (L <= 16 by default)
template<int MAXL=16>
struct SmallSig {
    uint8_t len;
    uint32_t data[MAXL];
    SmallSig() : len(0) {}
    inline void clear() { len = 0; }
    inline void push(uint32_t x) { data[len++] = x; }
    inline bool operator==(const SmallSig& o) const noexcept {
        if (len != o.len) return false;
        for (uint8_t i=0;i<len;++i) if (data[i]!=o.data[i]) return false;
        return true;
    }
};

struct SmallSigHash {
    size_t operator()(SmallSig<16> const& s) const noexcept {
        size_t h = 1469598103934665603ULL;
        for (uint8_t i=0;i<s.len;++i) { h ^= s.data[i]; h *= 1099511628211ULL; }
        h ^= size_t(s.len) * 0x9e3779b97f4a7c15ULL; h *= 1099511628211ULL;
        return h;
    }
};

// Seen-key (signature + occurrence set) to avoid recreating same TEMP across rounds
struct OccKey {
    int prod, start, end;
    bool operator==(const OccKey& o) const noexcept {
        return prod==o.prod && start==o.start && end==o.end;
    }
};
struct OccKeyHash {
    size_t operator()(const OccKey& x) const noexcept {
        return (size_t(uint32_t(x.prod))<<32) ^ (size_t(uint32_t(x.start))<<16) ^ size_t(uint32_t(x.end));
    }
};

struct SeenKey {
    SmallSig<16> sig;
    vector<OccKey> occs; // sorted
    bool operator==(const SeenKey& o) const noexcept {
        if (!(sig==o.sig) || occs.size()!=o.occs.size()) return false;
        for (size_t i=0;i<occs.size();++i) if (!(occs[i]==o.occs[i])) return false;
        return true;
    }
};
struct SeenKeyHash {
    size_t operator()(SeenKey const& k) const noexcept {
        size_t h = SmallSigHash{}(k.sig);
        for (auto& o : k.occs) {
            size_t x = (size_t(uint32_t(o.prod))<<32) ^ (size_t(uint32_t(o.start))<<16) ^ size_t(uint32_t(o.end));
            h ^= x; h *= 1099511628211ULL;
        }
        return h;
    }
};

// hash & eq for Token so it can be a key
struct TokenHash {
    size_t operator()(const Token& t) const noexcept {
        // pack (kind, idx) into a 64-bit hash
        return (size_t(uint32_t(t.idx)) << 3) ^ size_t(uint32_t(t.kind));
    }
};
struct TokenEq {
    bool operator()(const Token& a, const Token& b) const noexcept {
        return a.kind == b.kind && a.idx == b.idx;
    }
};

// ===================== Planner implementation =====================

Plan CSEPlanner::run() {
    Plan plan;
    if (products.empty()) return plan;

    // Working copy we mutate (insert H tokens)
    vector<vector<Token>> cur = products;

    // ---- Config knobs ----
    const int MIN_LEN = 2;    // ignore blocks shorter than this (try 3~4 for speed)
    const int MAX_ROUNDS = 10;
    
    // Hash of current products (structure only) — to detect no-op rounds
    auto hash_cur = [&](const vector<vector<Token>>& v)->size_t {
        size_t h = 1469598103934665603ULL;
        for (const auto& seq : v) {
            for (const auto& t : seq) {
                size_t x = (size_t(uint32_t(t.idx))<<2) ^ size_t(uint32_t(t.kind));
                h ^= x; h *= 1099511628211ULL;
            }
            h ^= 0x9e3779b97f4a7c15ULL; h *= 1099511628211ULL;
        }
        return h;
    };

    // Maximal extension (left & right) — must match across ALL occurrences
    auto maximal_extend = [&](vector<Token>& sig, vector<Occ>& occs) {
        // left
        for (;;) {
            bool ok = true, set = false; Token t0;
            for (auto& oc : occs) {
                if (oc.start == 0) { ok = false; break; }
                const Token& t = cur[oc.prod][oc.start - 1];
                if (!set) { t0=t; set=true; }
                else if (t.kind != t0.kind || t.idx != t0.idx) { ok=false; break; }
            }
            if (!ok) break;
            sig.insert(sig.begin(), t0);
            for (auto& oc : occs) --oc.start;
        }
        // right
        for (;;) {
            bool ok = true, set = false; Token t0;
            for (auto& oc : occs) {
                if (oc.end + 1 >= (int)cur[oc.prod].size()) { ok=false; break; }
                const Token& t = cur[oc.prod][oc.end + 1];
                if (!set) { t0=t; set=true; }
                else if (t.kind != t0.kind || t.idx != t0.idx) { ok=false; break; }
            }
            if (!ok) break;
            sig.push_back(t0);
            for (auto& oc : occs) ++oc.end;
        }
    };

    // Apply replacements: each accepted (prod,start..end) -> H(temp_id)
    struct Span { int start, end, tid; };
    auto apply_replacements = [&](const vector<TempInfo>& newTemps) {
        const int N = (int)cur.size();
        vector<vector<Span>> spans(N);
        for (const auto& t : newTemps) {
            for (const auto& oc : t.occs) spans[oc.prod].push_back({oc.start, oc.end, t.id});
        }
        for (int p=0;p<N;++p) {
            auto& v = spans[p];
            if (v.empty()) continue;
            sort(v.begin(), v.end(), [](const Span& a, const Span& b){
                return a.start < b.start;
            });
            const auto& seq = cur[p];
            vector<Token> out; out.reserve(seq.size());
            int i = 0; size_t k = 0;
            while (i < (int)seq.size()) {
                if (k < v.size() && i == v[k].start) {
                    out.push_back(Token{Kind::H, v[k].tid});
                    i = v[k].end + 1;
                    ++k;
                } else {
                    out.push_back(seq[i]);
                    ++i;
                }
            }
            cur[p].swap(out);
        }
    };

    // Global seen set to avoid re-creating the same TEMP across rounds
    unordered_set<SeenKey, SeenKeyHash> seen; seen.reserve(1024);

    int next_tid = 0;
    bool changed = true;
    int rounds = 0;

    while (changed && rounds < MAX_ROUNDS) {
        ++rounds;
        changed = false;
        size_t h_before = hash_cur(cur);

        // --------- 1) Mine substrings (len >= 2), TEMP included, with SmallSig ----------
        using Sig = SmallSig<16>;
        struct SigHash { size_t operator()(const Sig& s) const noexcept { return SmallSigHash{}(s); } };
        struct SigEq   { bool operator()(const Sig& a, const Sig& b) const noexcept { return a==b; } };

        unordered_map<Sig, vector<Occ>, SigHash, SigEq> bucket;
        size_t approx = 0;
        for (auto& s : cur) { int Lp=(int)s.size(); approx += (size_t)Lp*(Lp-1)/2; }
        bucket.reserve(approx);

        for (int p=0; p<(int)cur.size(); ++p){ // over all products that we need to simplify
	  const auto& s = cur[p];              // one product is s
	  int Lp = (int)s.size();              // its current size is Lp
	  for (int i=0;i<Lp;++i) {             // sub-product starts at i and ends at j (inclusive i and j)
	    Sig sig; // empty
	    for (int j=i;j<Lp;++j) {
	      sig.push(packTok(s[j].kind, uint32_t(s[j].idx))); // sub-product packed into integer
	      if (sig.len >= 2) bucket[sig].push_back({p,i,j}); // p--in which product it occurs, i--where is starts, j--where it ends
	      if (sig.len == 16) break; // safety for MAXL      // with this integer packing we can not go beyond 16, which corresponds to order 8.
	    }
	  }
        }

	//cout<<rounds<<") bucket size="<<bucket.size()<<endl;
        // --------- 2) Build maximal candidates; dedup per-round & skip already seen ----------
        struct Cand {
	  vector<Token> sig;
	  vector<Occ>   occs;   // sorted
	  long long     savings; // (occ-1)*(len-1)
        };

        struct PairKey {
	  Sig sig;
	  vector<Occ> occs; // sorted
        };
        struct PairHash {
	  size_t operator()(PairKey const& k) const noexcept {
	    size_t h = SmallSigHash{}(k.sig);
	    for (auto& oc : k.occs) {
	      size_t x = (size_t(uint32_t(oc.prod))<<32) ^ (size_t(uint32_t(oc.start))<<16) ^ size_t(uint32_t(oc.end));
	      h ^= x; h *= 1099511628211ULL;
	    }
	    return h;
	  }
        };
        struct PairEq {
	  bool operator()(PairKey const& a, PairKey const& b) const noexcept {
	    if (!(a.sig==b.sig) || a.occs.size()!=b.occs.size()) return false;
	    for (size_t i=0;i<a.occs.size();++i) {
	      if (a.occs[i].prod!=b.occs[i].prod || a.occs[i].start!=b.occs[i].start || a.occs[i].end!=b.occs[i].end) return false;
	    }
	    return true;
	  }
        };

        unordered_map<PairKey, Cand, PairHash, PairEq> uniq;
        uniq.reserve(bucket.size());
	
        for (auto& kv : bucket) {
            auto& occs0 = kv.second;
            if (occs0.size() < 2) continue;    // at least product of two, otherwise its not a product

            // Rebuild vector<Token> from Sig once (cheap)
            vector<Token> sig;
            sig.reserve(kv.first.len);
            for (uint8_t a=0;a<kv.first.len;++a) {
                uint32_t packed = kv.first.data[a];
                Kind k = Kind(packed & 0x3u);
                int  idx = int(packed >> 2);
                sig.push_back(Token{k, idx});
            }
	    
            // Sort occs (normalize)
            vector<Occ> occs = occs0;
            sort(occs.begin(), occs.end(), [](const Occ& A, const Occ& B){
                if (A.prod!=B.prod) return A.prod<B.prod;
                if (A.start!=B.start) return A.start<B.start;
                return A.end<B.end;
            });

            // Maximal left/right extension
            maximal_extend(sig, occs);

            // Savings with uniform matrix sizes: (occ-1)*(len-1)
            long long occ = (long long)occs.size();
            long long len = (long long)sig.size();
            if (len < MIN_LEN) continue; // filter small ones aggressively
            long long sav = (occ - 1) * (len - 1);
            if (sav <= 0) continue;

            // Build per-round uniq key (repack sig into SmallSig)
            Sig ss; ss.len = 0;
            for (auto& t : sig) {
                if (ss.len == 16) break;
                ss.push(packTok(t.kind, uint32_t(t.idx)));
            }
            PairKey key{ss, occs};

            // Global seen (signature + occ set) — skip if already materialized in past rounds
            SeenKey skey; skey.sig = ss;
            skey.occs.reserve(occs.size());
            for (auto& oc : occs) skey.occs.push_back(OccKey{oc.prod, oc.start, oc.end});
            if (seen.find(skey) != seen.end()) continue;

            auto it = uniq.find(key);
            if (it == uniq.end() || sav > it->second.savings) {
                uniq[key] = Cand{ std::move(sig), std::move(occs), sav };
            }
        }

        if (uniq.empty()) break;
	
        vector<Cand> cands; cands.reserve(uniq.size());
        for (auto& kv2 : uniq) cands.push_back(std::move(kv2.second));

        // --------- 3) Greedy select (non-overlapping); sort by savings desc, then length desc ----------
        sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b){
            if (a.savings != b.savings) return a.savings > b.savings;
            return a.sig.size() > b.sig.size();
        });

        vector<vector<char>> used(cur.size());
        for (int p=0;p<(int)cur.size();++p) used[p].assign(cur[p].size(), 0);

        vector<TempInfo> newTemps;
        newTemps.reserve(cands.size());

        for (auto& c : cands) {
            // accept non-overlapping occurrences in order
            const auto& occs = c.occs; // already sorted
            vector<Occ> accepted;
            for (auto& oc : occs) {
                bool ok = true;
                for (int k=oc.start; k<=oc.end; ++k)
                    if (used[oc.prod][k]) { ok=false; break; }
                if (!ok) continue;
                for (int k=oc.start; k<=oc.end; ++k) used[oc.prod][k] = 1;
                accepted.push_back(oc);
            }
            if ((int)accepted.size() >= 2) {
                TempInfo ti;
                ti.id   = next_tid++;
                ti.sig  = c.sig;
                ti.occs = std::move(accepted);

                // Record in global seen with the accepted occs (normalized)
                SeenKey skey;
                // pack sig
                for (auto& t : ti.sig) {
                    if (skey.sig.len == 16) break;
                    skey.sig.push(packTok(t.kind, uint32_t(t.idx)));
                }
                skey.occs.reserve(ti.occs.size());
                for (auto& oc : ti.occs)
                    skey.occs.push_back(OccKey{oc.prod, oc.start, oc.end});
                // sort occs in skey (should already be sorted)
                sort(skey.occs.begin(), skey.occs.end(),
                     [](const OccKey& A, const OccKey& B){
                         if (A.prod!=B.prod) return A.prod<B.prod;
                         if (A.start!=B.start) return A.start<B.start;
                         return A.end<B.end;
                     });
                seen.insert(std::move(skey));

                newTemps.push_back(std::move(ti));
            }
        }

        if (!newTemps.empty()) {
            size_t h0 = hash_cur(cur);
            apply_replacements(newTemps);
            size_t h1 = hash_cur(cur);
            if (h1 != h0) {
                for (auto& t : newTemps) plan.temps.push_back(std::move(t));
                changed = true;
            }
        }

        if (!changed) break;
    }

    plan.reduced = std::move(cur);
    return plan;
}

static void print_products(ostream& cout, const vector<vector<Token>>& prods, const string& title) {
  cout << "\n" << title << ":\n";
  for (size_t p=0;p<prods.size();++p) {
    cout << "P" << p << " = ";
    for (size_t i=0;i<prods[p].size();++i) {
      if (i) cout << " · ";
      cout << tok2str(prods[p][i]);
    }
    cout << "\n";
  }
}

/*
// ===================== Demo: build a test and run =====================
typedef std::chrono::duration<float> fsec;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    auto C = [](int g){ return Token{Kind::C, g}; };
    auto R = [](int g){ return Token{Kind::R, g}; };
    auto I = [](int g){ return Token{Kind::I, g}; };
    vector<vector<Token>> products = {
            { C(4) , C(5) , C(6) , C(7) , C(8) , C(9) , I(0) , R(1) },
            { C(4) , C(5) , C(6) , C(7) , C(8) , C(9) , R(0) , I(1) },
            { C(4) , C(5) , C(6) , C(7) , C(8) , C(9) , R(0) , R(1) },
            { R(4) , R(5) , R(6) , R(7) , R(8) , R(9) , C(0) , C(1) },
            { I(4) , R(5) , R(6) , R(7) , R(8) , R(9) , C(0) , C(1) },
            { R(4) , I(5) , R(6) , R(7) , R(8) , R(9) , C(0) , C(1) },
            { R(4) , R(5) , I(6) , R(7) , R(8) , R(9) , C(0) , C(1) },
            { R(4) , R(5) , R(6) , I(7) , R(8) , R(9) , C(0) , C(1) },
            { R(4) , R(5) , R(6) , R(7) , I(8) , R(9) , C(0) , C(1) },
            { R(4) , R(5) , R(6) , R(7) , R(8) , I(9) , C(0) , C(1) },
            { R(8) , R(9) , I(0) , R(1) , R(2) , R(3) , C(4) , C(5) },
            { R(8) , R(9) , R(0) , I(1) , R(2) , R(3) , C(4) , C(5) },
            { R(8) , R(9) , R(0) , R(1) , I(2) , R(3) , C(4) , C(5) },
            { R(8) , R(9) , R(0) , R(1) , R(2) , I(3) , C(4) , C(5) },
            { C(8) , C(9) , C(0) , C(1) , C(2) , C(3) , I(4) , R(5) },
            { C(8) , C(9) , C(0) , C(1) , C(2) , C(3) , R(4) , I(5) },
            { C(8) , C(9) , C(0) , C(1) , C(2) , C(3) , R(4) , R(5) },
            { R(8) , R(9) , R(0) , R(1) , R(2) , R(3) , C(4) , C(5) },
            { I(8) , R(9) , R(0) , R(1) , R(2) , R(3) , C(4) , C(5) },
            { R(8) , I(9) , R(0) , R(1) , R(2) , R(3) , C(4) , C(5) },
            { C(2) , C(3) , C(4) , R(5) , R(6) , R(7) , R(8) , R(9) },
            { R(2) , R(3) , R(4) , C(5) , C(6) , C(7) , C(8) , C(9) },
            { I(2) , R(3) , R(4) , C(5) , C(6) , C(7) , C(8) , C(9) },
            { R(2) , I(3) , R(4) , C(5) , C(6) , C(7) , C(8) , C(9) },
            { R(2) , R(3) , I(4) , C(5) , C(6) , C(7) , C(8) , C(9) },
            { C(2) , C(3) , C(4) , I(5) , R(6) , R(7) , R(8) , R(9) },
            { C(2) , C(3) , C(4) , R(5) , I(6) , R(7) , R(8) , R(9) },
            { C(2) , C(3) , C(4) , R(5) , R(6) , I(7) , R(8) , R(9) },
            { C(2) , C(3) , C(4) , R(5) , R(6) , R(7) , I(8) , R(9) },
            { C(2) , C(3) , C(4) , R(5) , R(6) , R(7) , R(8) , I(9) },
            { I(0) , R(1) , R(2) , R(3) , R(4) , R(5) , C(6) , C(7) },
            { R(0) , I(1) , R(2) , R(3) , R(4) , R(5) , C(6) , C(7) },
            { R(0) , R(1) , I(2) , R(3) , R(4) , R(5) , C(6) , C(7) },
            { R(0) , R(1) , R(2) , I(3) , R(4) , R(5) , C(6) , C(7) },
            { R(0) , R(1) , R(2) , R(3) , I(4) , R(5) , C(6) , C(7) },
            { R(0) , R(1) , R(2) , R(3) , R(4) , I(5) , C(6) , C(7) },
            { C(0) , C(1) , C(2) , C(3) , C(4) , C(5) , I(6) , R(7) },
            { C(0) , C(1) , C(2) , C(3) , C(4) , C(5) , R(6) , I(7) },
            { C(0) , C(1) , C(2) , C(3) , C(4) , C(5) , R(6) , R(7) },
            { R(0) , R(1) , R(2) , R(3) , R(4) , R(5) , C(6) , C(7) },
            { R(9) , I(0) , R(1) , C(2) , C(3) , C(4) , C(5) , C(6) },
            { R(9) , R(0) , I(1) , C(2) , C(3) , C(4) , C(5) , C(6) },
            { C(9) , C(0) , C(1) , I(2) , R(3) , R(4) , R(5) , R(6) },
            { C(9) , C(0) , C(1) , R(2) , I(3) , R(4) , R(5) , R(6) },
            { C(9) , C(0) , C(1) , R(2) , R(3) , I(4) , R(5) , R(6) },
            { C(9) , C(0) , C(1) , R(2) , R(3) , R(4) , I(5) , R(6) },
            { C(9) , C(0) , C(1) , R(2) , R(3) , R(4) , R(5) , I(6) },
            { C(9) , C(0) , C(1) , R(2) , R(3) , R(4) , R(5) , R(6) },
            { R(9) , R(0) , R(1) , C(2) , C(3) , C(4) , C(5) , C(6) },
            { I(9) , R(0) , R(1) , C(2) , C(3) , C(4) , C(5) , C(6) }
    };
    
    print_products(cout, products, "Original products");
    cout<<"# matrix multiplications="<< products.size()*(products[0].size()-1)<<endl;
    
    CSEPlanner planner{products};
    auto t0 = high_resolution_clock::now();
    Plan plan = planner.run();
    auto t1 = high_resolution_clock::now();


  // Print selected temporaries
  int matrix_multiply=0;
  cout << "\nTemporary sub-products found:\n";
  for (auto& t : plan.temps){
    cout << "H" << t.id << " = ";
    for (size_t i=0;i<t.sig.size();++i){
      if (i) cout << " · ";
      cout << tok2str(t.sig[i]);
    }
    //auto d = t.outDim;
    //cout << "   (dims " << d.m << "x" << d.n << ")\n   used at: ";
    cout << "     used at: ";
    for (auto& oc : t.occs) cout << "[P" << oc.prod << ":" << oc.start << "-" << oc.end << "] ";
    cout << "\n";
    matrix_multiply += t.sig.size()-1;
  }
  cout<<"matrix_multiplications for temps="<< matrix_multiply << endl;
  // Print reduced products
  cout << "\nReduced products:\n";
  for (size_t p=0;p<plan.reduced.size();++p){
    cout << "P" << p << " = ";
    const auto& seq = plan.reduced[p];
    for (size_t i=0;i<seq.size();++i){
      if (i) cout << " · ";
      cout << tok2str(seq[i]);
    }
    matrix_multiply += seq.size()-1;
    cout << "\n";
  }
  cout<<"total number of matrix multiplications ="<< matrix_multiply << endl;
  fsec dt = t1-t0;
  cout<<"Time : "<< dt.count() <<"s "<<endl;

  return 0;
}
*/
