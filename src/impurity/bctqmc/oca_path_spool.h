#ifndef EDMFT_OCA_PATH_SPOOL_H
#define EDMFT_OCA_PATH_SPOOL_H

#include <array>
#include <cstdint>
#include <cstdlib>
#include <unistd.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace edmft::oca_dlr {

// Exact, external-memory reduction of partial diagram weights. Record contains
// a complete key and five weights. No record is evaluated until every partial
// weight for its key has met in the same leaf. Files are private temporary
// native-format files, never portable restart/cache files.
//
// First separate reduced bath pairs. If a pair exceeds max_keys, partition its
// file by successive bits of a mixed full-key hash. This rereads a
// disk file, NOT the expensive atomic diagram walk. Peak hash/plan size is
// bounded even when a single bath pair is enormous. Equal keys always follow
// exactly the same route; exact cancellations are removed only after reduction.
template<class Record, class Hash>
class OCAPathSpool {
  using Key=decltype(Record{}.key);
  using Weights=decltype(Record{}.weights);
  using Table=std::unordered_map<Key,Weights,Hash>;
  std::filesystem::path directory_;
  size_t max_keys_;
  std::function<uint64_t(const Key&)> group_;
  std::map<uint64_t,std::filesystem::path> roots_;
  uint64_t next_file_=0;
  uint64_t bytes_written_=0;
  uint64_t partial_records_=0;
  uint64_t partitions_=0;
  size_t peak_keys_=0;
  uint64_t next_space_check_=0;

  void check_space(size_t upcoming_bytes) {
    if (bytes_written_<next_space_check_) return;
    // Leave 1 GiB headroom, rather than filling the filesystem and harming
    // the user's other jobs. Quota/network errors are also caught by checked
    // writes. This is a safety margin, not a reservation across MPI ranks.
    const auto space=std::filesystem::space(directory_);
    if (space.available<UINT64_C(1073741824)+upcoming_bytes)
      throw std::runtime_error("Insufficient OCA scratch disk space (1 GiB reserve): "+
                               directory_.string());
    next_space_check_=bytes_written_+UINT64_C(8388608);
  }

  static bool nonzero(const Weights& weights) {
    for (double w:weights) if (w!=0.0) return true;
    return false;
  }
  static uint64_t mixed_hash(const Key& key) {
    // Mix once, then use independent digits at each level. Do not hash only
    // the pseudo indices: bath indices are also part of duplicate identity.
    uint64_t h=static_cast<uint64_t>(Hash{}(key));
    h=(h^(h>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    h=(h^(h>>27))*UINT64_C(0x94d049bb133111eb);
    return h^(h>>31);
  }
  std::filesystem::path fresh_file() {
    return directory_/std::to_string(next_file_++);
  }
  static bool read_record(std::ifstream& in,Record& record) {
    in.read(reinterpret_cast<char*>(&record),sizeof(record));
    if (in.gcount()==sizeof(record)) return true;
    if (in.gcount()!=0 || !in.eof())
      throw std::runtime_error("OCA scratch file is truncated or unreadable");
    return false;
  }
  static void remove_file(const std::filesystem::path& path) {
    if (!std::filesystem::remove(path))
      throw std::runtime_error("OCA scratch file disappeared: "+path.string());
  }
  void reduce(const std::filesystem::path& path,unsigned shift,
              const std::function<void(std::vector<Record>&)>& consume) {
    Table table;
    table.max_load_factor(0.7f);
    table.reserve(max_keys_);
    std::ifstream in(path,std::ios::binary);
    if (!in) throw std::runtime_error("Cannot read OCA scratch: "+path.string());
    Record record{};
    bool overflow=false;
    while (read_record(in,record)) {
      auto it=table.find(record.key);
      if (it==table.end()) {
        if (table.size()==max_keys_) { overflow=true; break; }
        table.emplace(record.key,record.weights);
      } else {
        for (size_t j=0;j<record.weights.size();++j)
          it->second[j]+=record.weights[j];
      }
    }
    in.close();
    peak_keys_=std::max(peak_keys_,table.size());
    if (!overflow) {
      std::vector<Record> records;
      records.reserve(table.size());
      for (const auto& [key,weights]:table) {
        if (!nonzero(weights)) continue;
        Record out{};
        out.key=key;
        out.weights=weights;
        records.push_back(out);
      }
      Table().swap(table); // release hash nodes BEFORE building sparse plans
      std::sort(records.begin(),records.end(),
                [](const Record& a,const Record& b){return a.key<b.key;});
      if (!records.empty()) consume(records);
      remove_file(path);
      return;
    }
    Table().swap(table);
    if (shift>=64)
      throw std::runtime_error("OCA scratch hash collision exceeds batch capacity");

    // All streams close before recursion, so file descriptors and buffers do
    // not grow with recursion depth. Children preserve per-key record order.
    // Binary splits keep accepted leaves near capacity. A wide radix can
    // create hundreds of tiny plans and needlessly repeat pole preparation.
    std::array<std::filesystem::path,2> children;
    {
      std::array<std::ofstream,2> outputs;
      in.open(path,std::ios::binary);
      if (!in) throw std::runtime_error("Cannot repartition OCA scratch");
      while (read_record(in,record)) {
        const auto digit=(mixed_hash(record.key)>>shift)&1;
        if (children[digit].empty()) {
          children[digit]=fresh_file();
          outputs[digit].open(children[digit],std::ios::binary);
          if (!outputs[digit]) throw std::runtime_error("Cannot create OCA scratch partition");
        }
        check_space(sizeof(record));
        outputs[digit].write(reinterpret_cast<const char*>(&record),sizeof(record));
        if (!outputs[digit]) throw std::runtime_error("OCA scratch write failed (disk full?)");
        bytes_written_+=sizeof(record);
      }
      in.close();
      for (size_t i=0;i<children.size();++i) if (!children[i].empty()) {
        outputs[i].close();
        if (!outputs[i]) throw std::runtime_error("OCA scratch close failed (disk full?)");
      }
    }
    ++partitions_;
    remove_file(path); // do not retain parent and all descendants on disk
    for (const auto& child:children)
      if (!child.empty()) reduce(child,shift+1,consume);
  }

public:
  OCAPathSpool(size_t max_keys,const std::string& scratch_directory,
               std::function<uint64_t(const Key&)> group)
    : max_keys_(max_keys),group_(std::move(group)) {
    static_assert(std::is_trivially_copyable_v<Record>);
    if (max_keys==0) throw std::invalid_argument("OCA batch capacity must be positive");
    const auto parent=scratch_directory.empty()
      ? std::filesystem::temp_directory_path() : std::filesystem::path(scratch_directory);
    // mkdtemp prevents overwriting another rank/job's files. Only this newly
    // created directory is eligible for destructor cleanup.
    std::string pattern=(parent/"oca-physical-XXXXXX").string();
    std::vector<char> name(pattern.begin(),pattern.end());
    name.push_back('\0');
    if (!::mkdtemp(name.data()))
      throw std::runtime_error("Cannot create OCA scratch directory under "+parent.string());
    directory_=name.data();
  }
  OCAPathSpool(const OCAPathSpool&)=delete;
  OCAPathSpool& operator=(const OCAPathSpool&)=delete;
  ~OCAPathSpool() {
    // Never remove a user-supplied directory: only our unique mkdtemp child.
    std::error_code error;
    if (!directory_.empty()) std::filesystem::remove_all(directory_,error);
  }
  const auto& directory() const { return directory_; }
  uint64_t bytes_written() const { return bytes_written_; }
  uint64_t partial_records() const { return partial_records_; }
  uint64_t partitions() const { return partitions_; }
  size_t peak_keys() const { return peak_keys_; }

  void append(const std::vector<Record>& records) {
    // Group a fixed-size flush, so opening/closing files costs once per bath
    // pair per flush, not once per contribution. At most max_keys records live
    // in this grouping copy; no file stays open between flushes.
    std::map<uint64_t,std::vector<Record>> groups;
    for (const auto& record:records) groups[group_(record.key)].push_back(record);
    for (const auto& [group,values]:groups) {
      auto [it,inserted]=roots_.try_emplace(group);
      if (inserted) it->second=fresh_file();
      std::ofstream out(it->second,std::ios::binary|std::ios::app);
      if (!out) throw std::runtime_error("Cannot open OCA scratch output");
      const size_t bytes=values.size()*sizeof(Record);
      check_space(bytes);
      out.write(reinterpret_cast<const char*>(values.data()),bytes);
      out.close();
      if (!out) throw std::runtime_error("OCA scratch write failed (disk full?)");
      bytes_written_+=bytes;
      partial_records_+=values.size();
    }
  }
  void consume(const std::function<void(std::vector<Record>&)>& callback) {
    for (const auto& [group,path]:roots_) reduce(path,0,callback);
    roots_.clear();
  }
};
} // namespace edmft::oca_dlr
#endif
