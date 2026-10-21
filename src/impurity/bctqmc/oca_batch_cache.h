#ifndef EDMFT_OCA_BATCH_CACHE_H
#define EDMFT_OCA_BATCH_CACHE_H

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

namespace edmft::oca_dlr {

// Persistent topology only: never store propagator values, AAA poles or DLR
// transforms. A generation is immutable after publication; a tiny CURRENT
// pointer selects it atomically. Old generations are deliberately not deleted
// here: another solver may still be reading them. Cleanup is an explicit user
// operation, only when no jobs use this cache.
//
// Pair files contain COMPLETE, nonzero, unique key sums, concatenated from the
// spool's leaves. No key appears in two leaves. Thus the reader may split them
// at ANY record boundary and obey today's batch-size limit, independently of
// the batch size and MPI layout that originally built the cache.
template<class Record>
class OCABatchCache {
public:
  using Signature=std::array<uint64_t,8>;
  static constexpr uint64_t offset_basis=UINT64_C(1469598103934665603);
  struct Entry {
    uint64_t pair=0, count=0, hash=offset_basis;
  };
  static uint64_t hash_bytes(const void* data,size_t size,uint64_t hash=offset_basis) {
    const auto* bytes=static_cast<const unsigned char*>(data);
    for (size_t i=0;i<size;++i) { hash^=bytes[i]; hash*=UINT64_C(1099511628211); }
    return hash;
  }
private:
  std::filesystem::path base_,generation_;
  Signature signature_{};
  std::vector<Entry> entries_;
  std::map<uint64_t,Entry> written_;
  bool owns_unpublished_=false;
  static bool safe_generation(const std::string& s) {
    return s.rfind("build-",0)==0 && s.size()<=80 &&
      s.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-")==std::string::npos;
  }
  static uint64_t endian_tag() {
    const uint64_t word=1;
    return *reinterpret_cast<const unsigned char*>(&word)==1 ? 1 : 2;
  }
  static void durable(const std::filesystem::path& p) {
    const int fd=::open(p.c_str(),O_RDONLY);
    if (fd<0) throw std::runtime_error("Cannot open OCA cache for sync: "+p.string());
    const int result=::fsync(fd);
    ::close(fd);
    if (result!=0) throw std::runtime_error("Cannot sync OCA cache: "+p.string());
  }
  static std::string read_text(const std::filesystem::path& p) {
    // Metadata is at most a few MiB even at the packed-key bath limit. Bound
    // the read before parsing, so corrupt counts cannot cause huge allocations.
    if (std::filesystem::file_size(p)>16*1024*1024)
      throw std::runtime_error("Oversized OCA cache metadata");
    std::ifstream in(p,std::ios::binary);
    if (!in) throw std::runtime_error("Cannot read OCA cache metadata");
    std::ostringstream text; text<<in.rdbuf();
    if (in.bad()) throw std::runtime_error("Cannot read OCA cache metadata");
    return text.str();
  }
  static void write_text(const std::filesystem::path& p,const std::string& text) {
    std::ofstream out(p,std::ios::binary|std::ios::trunc);
    if (!out) throw std::runtime_error("Cannot create OCA cache metadata: "+p.string());
    out<<text; out.close();
    if (!out) throw std::runtime_error("Cannot write OCA cache metadata: "+p.string());
    durable(p);
  }
  std::filesystem::path pair_path(uint64_t pair) const {
    return generation_/("pair-"+std::to_string(pair)+".bin");
  }
  std::string serialize() const {
    std::ostringstream out;
    // Version must change whenever the packed key or any weight convention
    // changes. Native record layout is checked, not assumed portable.
    out<<"OCA_BATCH_CACHE 1 "<<sizeof(Record)<<" "<<endian_tag()<<"\n";
    for (auto value:signature_) out<<value<<" ";
    out<<"\n"<<generation_.filename().string()<<"\n"<<entries_.size()<<"\n";
    for (const auto& e:entries_) out<<e.pair<<" "<<e.count<<" "<<e.hash<<"\n";
    return out.str();
  }
public:
  OCABatchCache(std::string directory,Signature signature)
    : base_(std::move(directory)),signature_(signature) {
    static_assert(std::is_trivially_copyable_v<Record>);
  }
  OCABatchCache(const OCABatchCache&)=delete;
  OCABatchCache& operator=(const OCABatchCache&)=delete;
  ~OCABatchCache() {
    // Only the master owns its freshly created private build directory.
    // Never delete CURRENT, published data, or a user-supplied directory.
    if (owns_unpublished_) {
      std::error_code error; std::filesystem::remove_all(generation_,error);
    }
  }
  const auto& entries() const { return entries_; }
  const auto& directory() const { return base_; }
  uint64_t total_records() const {
    uint64_t total=0; for (const auto& e:entries_) total+=e.count; return total;
  }
  // The master does the small metadata/size checks once, then broadcasts this
  // manifest text. Payload checksums are verified by the assigned reader while
  // streaming, avoiding a second full traversal of a many-GB cache.
  std::string load(std::string& reason) {
    try {
      const std::string current=read_text(base_/"CURRENT");
      std::istringstream pointer(current);
      std::string name,extra; pointer>>name;
      if (!safe_generation(name) || (pointer>>extra))
        throw std::runtime_error("Invalid OCA cache generation pointer");
      const auto manifest=read_text(base_/name/"manifest");
      accept_manifest(manifest);
      if (generation_.filename()!=name) throw std::runtime_error("OCA cache generation mismatch");
      for (const auto& e:entries_)
        if (std::filesystem::file_size(pair_path(e.pair))!=e.count*sizeof(Record))
          throw std::runtime_error("OCA cache pair file size mismatch");
      reason="valid"; return manifest;
    } catch (const std::exception& e) {
      entries_.clear(); generation_.clear(); reason=e.what(); return {};
    }
  }
  void accept_manifest(const std::string& text) {
    std::istringstream in(text);
    std::string magic,name,extra;
    uint64_t version=0,record_size=0,endian=0,n=0;
    Signature sig{};
    in>>magic>>version>>record_size>>endian;
    for (auto& x:sig) in>>x;
    in>>name>>n;
    if (!in || magic!="OCA_BATCH_CACHE" || version!=1 || record_size!=sizeof(Record) ||
        endian!=endian_tag() || sig!=signature_ || !safe_generation(name) || n>65536)
      throw std::runtime_error("OCA cache format or input signature mismatch");
    std::vector<Entry> entries;
    uint64_t total=0;
    for (uint64_t i=0;i<n;++i) {
      Entry e; in>>e.pair>>e.count>>e.hash;
      if (!in || !e.count || e.count>uint64_t(std::numeric_limits<std::streamoff>::max())/sizeof(Record) ||
          e.count>std::numeric_limits<uint64_t>::max()-total ||
          (!entries.empty() && e.pair<=entries.back().pair))
        throw std::runtime_error("Invalid OCA cache pair manifest");
      entries.push_back(e); total+=e.count;
    }
    if (in>>extra) throw std::runtime_error("Trailing OCA cache metadata");
    generation_=base_/name; entries_=std::move(entries);
  }
  std::string begin() {
    std::filesystem::create_directories(base_);
    std::string pattern=(base_/"build-XXXXXX").string();
    std::vector<char> path(pattern.begin(),pattern.end()); path.push_back('\0');
    if (!::mkdtemp(path.data())) throw std::runtime_error("Cannot create persistent OCA cache generation");
    generation_=path.data(); owns_unpublished_=true;
    return generation_.filename().string();
  }
  void join(const std::string& name) {
    if (!safe_generation(name)) throw std::runtime_error("Invalid OCA cache build name");
    generation_=base_/name;
  }
  void append(uint64_t pair,const std::vector<Record>& records) {
    if (records.empty()) return;
    const auto bytes=records.size()*sizeof(Record);
    if (std::filesystem::space(generation_).available<UINT64_C(1073741824)+bytes)
      throw std::runtime_error("Insufficient persistent OCA cache space (1 GiB reserve)");
    auto [it,added]=written_.try_emplace(pair,Entry{pair,0,offset_basis});
    // Unique pair ownership is guaranteed by enumeration. Replicated mode
    // has exactly one writer (master); all other ranks still evaluate normally.
    std::ofstream out(pair_path(pair),std::ios::binary|(added?std::ios::trunc:std::ios::app));
    out.write(reinterpret_cast<const char*>(records.data()),bytes); out.close();
    if (!out) throw std::runtime_error("OCA persistent cache write failed");
    it->second.count+=records.size();
    it->second.hash=hash_bytes(records.data(),bytes,it->second.hash);
  }
  void finish_writer(int rank) {
    std::ostringstream out; out<<written_.size()<<"\n";
    for (const auto& [pair,e]:written_) {
      durable(pair_path(pair));
      out<<e.pair<<" "<<e.count<<" "<<e.hash<<"\n";
    }
    write_text(generation_/("rank-"+std::to_string(rank)+".index"),out.str());
  }
  // Called by master ONLY after all writers have closed/synced their data.
  // Publishing a new generation never changes data already opened by readers.
  void publish(int writer_ranks) {
    std::map<uint64_t,Entry> merged;
    for (int rank=0;rank<writer_ranks;++rank) {
      std::istringstream in(read_text(generation_/("rank-"+std::to_string(rank)+".index")));
      uint64_t n=0; in>>n;
      if (!in || n>65536) throw std::runtime_error("Invalid OCA writer index");
      for (uint64_t i=0;i<n;++i) {
        Entry e; in>>e.pair>>e.count>>e.hash;
        if (!in || !e.count || !merged.emplace(e.pair,e).second ||
            e.count>uint64_t(std::numeric_limits<std::streamoff>::max())/sizeof(Record) ||
            std::filesystem::file_size(pair_path(e.pair))!=e.count*sizeof(Record))
          throw std::runtime_error("Incomplete or multiply owned OCA cache pair");
      }
      std::string extra; if (in>>extra) throw std::runtime_error("Invalid OCA writer index tail");
    }
    entries_.clear(); for (const auto& [pair,e]:merged) entries_.push_back(e);
    write_text(generation_/"manifest",serialize()); durable(generation_);
    // Keep the pointer temp file INSIDE this private generation so exception
    // cleanup also removes it. Rename is atomic on this one filesystem.
    write_text(generation_/"pointer",generation_.filename().string()+"\n");
    std::filesystem::rename(generation_/"pointer",base_/"CURRENT");
    owns_unpublished_=false; // CURRENT now refers here, even if directory sync fails
    durable(base_);
  }
  std::vector<int> owners(int ranks) const {
    if (ranks<1) throw std::invalid_argument("OCA cache rank count must be positive");
    std::vector<size_t> order(entries_.size());
    for (size_t i=0;i<order.size();++i) order[i]=i;
    std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return entries_[a].count>entries_[b].count;});
    std::vector<uint64_t> load(ranks,0); std::vector<int> result(entries_.size());
    for (auto i:order) {
      int owner=std::min_element(load.begin(),load.end())-load.begin();
      result[i]=owner; load[owner]+=entries_[i].count;
    }
    return result;
  }
  void read_pair(size_t index,size_t cap,const std::function<void(std::vector<Record>&)>& consume) const {
    if (!cap) throw std::invalid_argument("OCA cache batch capacity must be positive");
    const auto& e=entries_.at(index);
    std::ifstream in(pair_path(e.pair),std::ios::binary);
    if (!in) throw std::runtime_error("Cannot read OCA cache pair (directory must be visible to every rank)");
    uint64_t remaining=e.count,hash=offset_basis;
    while (remaining) {
      const auto n=static_cast<size_t>(std::min<uint64_t>(remaining,cap));
      std::vector<Record> records(n);
      const auto bytes=n*sizeof(Record);
      in.read(reinterpret_cast<char*>(records.data()),bytes);
      if (!in || in.gcount()!=static_cast<std::streamsize>(bytes))
        throw std::runtime_error("Truncated OCA persistent cache payload");
      hash=hash_bytes(records.data(),bytes,hash);
      remaining-=n;
      // Validate the final checksum before consuming the last chunk. Earlier
      // chunks have only updated local accumulators; no result is published
      // until ALL pairs pass validation and the final MPI reduction succeeds.
      if (!remaining && (in.peek()!=std::char_traits<char>::eof() || hash!=e.hash))
        throw std::runtime_error("Corrupt OCA persistent cache payload (checksum/length mismatch)");
      consume(records);
    }
  }
};
} // namespace edmft::oca_dlr
#endif
