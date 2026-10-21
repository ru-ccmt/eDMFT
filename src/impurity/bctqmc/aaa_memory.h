// Read-only, optional process-memory diagnostics for AAA initialization.
#ifndef EDMFT_AAA_MEMORY_H
#define EDMFT_AAA_MEMORY_H
#include <ostream>
#if defined(__APPLE__)
#include <unistd.h>
#include <sys/resource.h>
// libproc.h pulls Mach IPC headers which are incompatible with some Homebrew
// GCC/macOS SDK combinations. The buffer ABI comes from sys/resource.h; this
// is the same C declaration as libproc.h, without its unrelated IPC headers.
extern "C" int proc_pid_rusage(int, int, rusage_info_t*);
#elif defined(__linux__)
#include <fstream>
#include <unistd.h>
#include <sys/resource.h>
#endif

namespace edmft::aaa::detail {
inline void print_memory(std::ostream* out, const char* stage)
{
  if (!out) return;
  constexpr double mib=1024.0*1024.0;
  *out << "AAA memory: stage=" << stage;
#if defined(__APPLE__)
  // RSS alone excludes compressed/swapped pages. Apple's physical footprint
  // is therefore reported separately, matching the relevant Activity Monitor
  // accounting more closely. Neither field is a sum over MPI ranks.
  struct rusage_info_v2 info{};
  if (proc_pid_rusage(getpid(),RUSAGE_INFO_V2,
                      reinterpret_cast<rusage_info_t*>(&info))==0)
    *out << " rss_mib=" << info.ri_resident_size/mib
         << " footprint_mib=" << info.ri_phys_footprint/mib;
  struct rusage usage{};
  if (getrusage(RUSAGE_SELF,&usage)==0)
    *out << " peak_rss_mib=" << usage.ru_maxrss/mib;
#elif defined(__linux__)
  // /proc is optional (e.g. a restricted container); omit unavailable fields.
  std::ifstream statm("/proc/self/statm");
  unsigned long long pages=0,resident=0;
  const long page_size=sysconf(_SC_PAGESIZE);
  if (page_size>0 && (statm>>pages>>resident))
    *out << " rss_mib=" << resident*static_cast<double>(page_size)/mib;
  struct rusage usage{};
  if (getrusage(RUSAGE_SELF,&usage)==0)
    *out << " peak_rss_mib=" << usage.ru_maxrss/1024.0;
#else
  *out << " unavailable_on_this_platform";
#endif
  // Flush at phase boundaries so a killed run still identifies its last
  // allocation stage. Peak RSS is a process-lifetime high water mark.
  *out << std::endl;
}
} // namespace edmft::aaa::detail
#endif
