/*---------------------------------------------------------------------------*\
  fastSnappyHexMesh
  Copyright (C) 2026 Seyit Temiz
-------------------------------------------------------------------------------
License
    This file is part of fastSnappyHexMesh.

    fastSnappyHexMesh is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    fastSnappyHexMesh is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fastSnappyHexMesh.  If not, see
    <http://www.gnu.org/licenses/>.

  fastSnappyHexMesh profiling

Description
    Hierarchical wall-clock + memory profiler. A scope adds its wall time to
    the bucket named by the current scope stack ("castellate/surfaceRefine/
    refine/changeMesh"), so repeated iterations aggregate into one row.

    Per bucket:
      calls, wall [s] (total, max single call),
      heap  [MB] live malloc'ed bytes at scope exit (max over calls),
      dHeap [MB] largest heap growth during one call,
      hwm+  [MB] how much this scope raised the process RSS high-water mark.

    No collective communication: every rank keeps its own table. Master
    prints its table via Info; every rank writes <case>/fsProfile/rank<N>.

    Enable/disable with env FS_PROFILE=0 (default on).

\*---------------------------------------------------------------------------*/

#ifndef fsProfiler_H
#define fsProfiler_H

#include <string>

namespace Foam
{
namespace fsProfiler
{

//- RAII timing scope
class scope
{
    double t0_;
    double heap0_;
    double hwm0_;
    bool active_;

public:

    explicit scope(const char* name);
    explicit scope(const std::string& name);
    ~scope();

    scope(const scope&) = delete;
    void operator=(const scope&) = delete;
};

//- Current live heap [MB] (malloc in-use + mmapped)
double heapMB();

//- Current RSS and RSS high-water mark [MB]
double rssMB();
double hwmMB();

//- Return free memory to the system (tcmalloc only; no-op otherwise)
void releaseFreeMemory();

//- Low-memory mode: tcmalloc aggressive decommit (freed memory goes back to
//  the system at once). Measured at 3.5M cells: peak RSS -25..-33 %,
//  time +25..+30 %. Returns false if not running on tcmalloc.
bool setLowMemoryMode();

//- Name of the memory allocator in use ("tcmalloc (built in)" or "glibc")
const char* allocatorName();

//- Print master table to Info and write per-rank file into dir
void report(const std::string& dir);

} // End namespace fsProfiler
} // End namespace Foam

#define FS_PROF_CAT2(a, b) a##b
#define FS_PROF_CAT(a, b) FS_PROF_CAT2(a, b)
#define FS_PROF(name) \
    ::Foam::fsProfiler::scope FS_PROF_CAT(fsProfScope_, __LINE__)(name)

#endif

// ************************************************************************* //
