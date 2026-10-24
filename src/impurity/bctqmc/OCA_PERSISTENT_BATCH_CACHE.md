Persistent bounded physical OCA enumeration
=========================================

The batched physical OCA evaluator can reuse its completed enumeration across
separate impurity-solver/DMFT invocations. This does NOT cache Green functions,
hybridization values, AAA fits, DLR transforms, or numerical diagram integrals.
Every invocation evaluates the saved topology with its current propagators.

Activation and location
-----------------------
Existing parameters control the feature; no new PARAMS entries are required:

  "oca_dlr_disk_cache": true
  "oca_dlr_cache_file": ".oca_dlr_paths.cache"

These defaults create a directory named
  .oca_dlr_paths.cache.physical.batched/
in the solver working directory. Preserve the WHOLE directory for reuse,
including CURRENT and its referenced build-* generation. It is a hidden
directory: ordinary shell '*' copying may omit it.

To reuse across different DMFT working directories, set oca_dlr_cache_file to
the same absolute persistent base path, or copy the entire cache directory.
The path is a BASE: the solver appends .physical.batched for this cache.
All MPI ranks must see that directory (shared filesystem or identically staged
files). Building a new cache requires shared writable storage. Existing .self
and full-plan .physical cache files are separate and are not changed here.

This feature applies when oca_dlr_physical_mode="batched", or when "auto"
selects the batched evaluator (currently Npseudo_unique >= 4096). The default
SC cluster full-plan path is unchanged. Setting oca_dlr_disk_cache=false
restores temporary-only spooling. oca_dlr_batch_directory still selects the
temporary raw-enumeration scratch directory, NOT the persistent cache directory.

What is saved and validated
---------------------------
Each record contains a packed topology key and all five merged diagram weights,
including Bulla contributions. Duplicate partial records are merged first;
only complete nonzero unique-key sums are saved. The cache validates:

  * CIX file contents;
  * actual loaded signed Sindex mapping, not just filename or dimensions;
  * bath/operator mappings and actual Bulla operator data;
  * actual pseudo-particle importance labels and OCA_cut_sum;
  * pseudo/bath dimensions and the binary format/convention revision.

Changing only hybridization, pseudo-propagator values, AAA/DLR parameters, MPI
rank count, or batch capacity does not invalidate this topology. If changing
propagators changes the importance labels, however, that DOES invalidate it.
CIX changes are treated conservatively: even a comment change rebuilds it.
The payload uses native binary layout; record size and endianness are checked.
Future changes to packed keys or weight conventions must bump the cache version.

Memory, communication and disk behavior
--------------------------------------
Reads obey the CURRENT oca_dlr_physical_batch_paths capacity, even when the
cache was created with a different capacity. Whole bath pairs are reassigned
greedily by saved record counts to the CURRENT rank count. There is no large
allgather and no per-batch MPI communication. Only small manifest/startup/end
messages and the existing final output-node reductions are collective.
Reading requires sequential disk I/O; it removes enumeration and raw-spool
merging, not the subsequent DLR evaluation or small per-batch plan construction.

On a cold build, temporary raw scratch and the growing merged persistent cache
coexist. The latter costs sizeof(OCAPhysicalDiskRecord) per unique path (56 bytes
with the current supported layout), plus small metadata. The solver checks for
a 1 GiB free-space reserve on writes; this is not a global disk-space guarantee.
Cache-hit runs need no raw spool directory.

Safety and diagnostics
----------------------
Look for:
  OCA bounded cache miss: ... reason=...
  OCA bounded cache published: ... pairs=... unique_paths=... data_GiB=...
or, on a subsequent invocation:
  OCA bounded cache hit: ... batch_paths=... ranks=...

Publication occurs after every rank has finished writing/evaluating its batches
and syncing the files. An atomic CURRENT replacement selects an immutable
generation. An interrupted build is not considered a valid new cache. A failed
build does not overwrite the previous published generation. Abrupt termination
can leave an unreferenced build-* directory, which can be removed after checking
that no running job owns/reads it.

Published older generations are deliberately retained: another invocation may
still be reading them. This means invalidations can accumulate disk usage.
When NO solver uses this cache, obsolete build-* generations may be removed;
keep the generation named by CURRENT, or remove the entire cache to rebuild.
Do not delete or edit cache files while jobs are using them.

Metadata/size/signature mismatches trigger rebuilding. Same-size payload
corruption is detected by a checksum during the streaming read and causes a
clear solver failure, not a silently incorrect result. Remove/rename the broken
cache after stopping its users, then rerun to rebuild.

A solver process already running an older executable cannot acquire this
feature mid-run. Its next invocation with the updated executable builds the
cache; following compatible invocations can reuse it.

Standalone infrastructure test: make test_oca_batch_cache && ./test_oca_batch_cache
