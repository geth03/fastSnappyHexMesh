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

  fastSnappyHexMesh profiling - see fsProfiler.H
\*---------------------------------------------------------------------------*/

#include "fsProfiler.H"
#include "Pstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "IOstreams.H"

#include <malloc.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>
#include <algorithm>

namespace
{

struct bucket
{
    long calls = 0;
    double wall = 0;
    double wallMax = 0;
    double heapEnd = 0;
    double heapDelta = 0;
    double hwmRise = 0;
    int depth = 0;
};

struct state
{
    std::vector<std::string> stack;
    std::vector<std::string> order;
    std::map<std::string, bucket> buckets;
    bool enabled = true;
    double tStart = 0;

    // Memory is sampled only for scopes up to this depth: reading the
    // tcmalloc heap counter walks its free lists and VmHWM parses
    // /proc/self/status - too expensive for small, frequent scopes.
    int memDepth = 4;

    state()
    {
        const char* e = std::getenv("FS_PROFILE");
        enabled = !(e && std::strcmp(e, "0") == 0);
        if (const char* m = std::getenv("FS_PROFILE_MEMDEPTH"))
        {
            memDepth = std::atoi(m);
        }
    }
};

state& S()
{
    static state s;
    return s;
}

double now()
{
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

double procStatusMB(const char* key)
{
    FILE* f = std::fopen("/proc/self/status", "r");
    if (!f) return 0;
    char line[256];
    const size_t n = std::strlen(key);
    double kb = 0;
    while (std::fgets(line, sizeof(line), f))
    {
        if (std::strncmp(line, key, n) == 0)
        {
            kb = std::atof(line + n);
            break;
        }
    }
    std::fclose(f);
    return kb/1024.0;
}

std::string path()
{
    std::string p;
    for (const std::string& s : S().stack)
    {
        if (!p.empty()) p += '/';
        p += s;
    }
    return p;
}

} // End anonymous namespace


// tcmalloc's C API, present only when tcmalloc is linked in (weak symbol)
extern "C" int MallocExtension_GetNumericProperty
(
    const char* property,
    size_t* value
) __attribute__((weak));


extern "C" void MallocExtension_ReleaseFreeMemory() __attribute__((weak));

extern "C" int MallocExtension_SetNumericProperty
(
    const char* property,
    size_t value
) __attribute__((weak));

bool Foam::fsProfiler::setLowMemoryMode()
{
    if (MallocExtension_SetNumericProperty)
    {
        return MallocExtension_SetNumericProperty
        (
            "tcmalloc.aggressive_memory_decommit",
            1
        );
    }
    return false;
}


void Foam::fsProfiler::releaseFreeMemory()
{
    if (MallocExtension_ReleaseFreeMemory)
    {
        MallocExtension_ReleaseFreeMemory();
    }
}


const char* Foam::fsProfiler::allocatorName()
{
    return MallocExtension_GetNumericProperty
        ? "tcmalloc (built in)"
        : "glibc malloc";
}


double Foam::fsProfiler::heapMB()
{
    if (MallocExtension_GetNumericProperty)
    {
        size_t bytes = 0;
        if
        (
            MallocExtension_GetNumericProperty
            (
                "generic.current_allocated_bytes",
                &bytes
            )
        )
        {
            return double(bytes)/(1024.0*1024.0);
        }
    }

    const struct mallinfo2 mi = mallinfo2();
    return double(mi.uordblks + mi.hblkhd)/(1024.0*1024.0);
}


double Foam::fsProfiler::rssMB()
{
    return procStatusMB("VmRSS:");
}


double Foam::fsProfiler::hwmMB()
{
    return procStatusMB("VmHWM:");
}


Foam::fsProfiler::scope::scope(const char* name)
:
    scope(std::string(name))
{}


Foam::fsProfiler::scope::scope(const std::string& name)
:
    t0_(0),
    heap0_(0),
    hwm0_(0),
    active_(S().enabled)
{
    if (!active_) return;

    if (S().tStart == 0) S().tStart = now();

    S().stack.push_back(name);
    const std::string p(path());
    if (!S().buckets.count(p))
    {
        S().order.push_back(p);
        S().buckets[p].depth = int(S().stack.size()) - 1;
    }

    if (int(S().stack.size()) - 1 <= S().memDepth)
    {
        heap0_ = heapMB();
        hwm0_ = hwmMB();
    }
    t0_ = now();
}


Foam::fsProfiler::scope::~scope()
{
    if (!active_) return;

    const double dt = now() - t0_;

    bucket& b = S().buckets[path()];
    b.calls++;
    b.wall += dt;
    b.wallMax = std::max(b.wallMax, dt);

    double heap = 0;
    double hwm = 0;
    if (b.depth <= S().memDepth)
    {
        heap = heapMB();
        hwm = hwmMB();
        b.heapEnd = std::max(b.heapEnd, heap);
        b.heapDelta = std::max(b.heapDelta, heap - heap0_);
        b.hwmRise += hwm - hwm0_;
    }

    // Live line for the outer two levels
    if (b.depth <= 1)
    {
        char buf[512];
        std::snprintf
        (
            buf, sizeof(buf),
            "[fsProf] %-50s %9.3f s  heap %8.1f MB  rss %8.1f MB  hwm %8.1f MB",
            path().c_str(), dt, heap, rssMB(), hwm
        );
        Info<< buf << endl;
    }

    S().stack.pop_back();
}


void Foam::fsProfiler::report(const std::string& dir)
{
    if (!S().enabled) return;

    const double total = now() - S().tStart;

    std::vector<std::string> lines;
    {
        char buf[512];
        std::snprintf
        (
            buf, sizeof(buf),
            "%-64s %6s %10s %6s %9s %9s %9s %9s",
            "# step", "calls", "wall_s", "%", "max1_s",
            "heap_MB", "dHeap_MB", "hwm+_MB"
        );
        lines.push_back(buf);
    }
    for (const std::string& p : S().order)
    {
        const bucket& b = S().buckets[p];

        // Indented leaf name for readability
        std::string leaf = p.substr(p.find_last_of('/') == std::string::npos
            ? 0 : p.find_last_of('/') + 1);
        std::string label(2*b.depth, ' ');
        label += leaf;

        char buf[512];
        std::snprintf
        (
            buf, sizeof(buf),
            "%-64s %6ld %10.3f %6.1f %9.3f %9.1f %9.1f %9.1f",
            label.c_str(), b.calls, b.wall, 100.0*b.wall/std::max(total, 1e-9),
            b.wallMax, b.heapEnd, b.heapDelta, b.hwmRise
        );
        lines.push_back(buf);
    }
    {
        char buf[256];
        std::snprintf
        (
            buf, sizeof(buf),
            "# total_wall_s %.3f  rss_MB %.1f  hwm_MB %.1f  heap_MB %.1f",
            total, rssMB(), hwmMB(), heapMB()
        );
        lines.push_back(buf);
    }

    Info<< nl << "fastSnappy profile (master rank)" << nl;
    for (const std::string& l : lines)
    {
        Info<< l.c_str() << nl;
    }
    Info<< endl;

    // Per-rank files, full paths so they can be merged
    Foam::mkDir(dir);
    const std::string fn =
        dir + "/rank" + std::to_string(Pstream::parRun() ? Pstream::myProcNo() : 0);
    FILE* f = std::fopen(fn.c_str(), "w");
    if (f)
    {
        std::fprintf
        (
            f, "# path calls wall_s max1_s heap_MB dHeap_MB hwm+_MB\n"
        );
        for (const std::string& p : S().order)
        {
            const bucket& b = S().buckets[p];
            std::fprintf
            (
                f, "%s %ld %.4f %.4f %.1f %.1f %.1f\n",
                p.c_str(), b.calls, b.wall, b.wallMax,
                b.heapEnd, b.heapDelta, b.hwmRise
            );
        }
        std::fprintf
        (
            f, "#total %.4f rss %.1f hwm %.1f\n", total, rssMB(), hwmMB()
        );
        std::fclose(f);
    }
}


// ************************************************************************* //
