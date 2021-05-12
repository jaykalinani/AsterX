#include "driver.hxx"
#include "io.hxx"
#include "io_silo.hxx"
#include "io_tsv.hxx"
#include "reduction.hxx"
#include "schedule.hxx"
#include "timer.hxx"

#include <CactusBase/IOUtil/src/ioGH.h>
#include <CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <AMReX.H>
#include <AMReX_Orientation.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <yaml-cpp/yaml.h>

#include <future>
#include <regex>

namespace CarpetX {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// Recovering

int recover_iteration = -1;

extern "C" int CarpetX_RecoverParameters() {
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("RecoverParameters");
  Interval interval(timer);

#ifdef HAVE_CAPABILITY_Silo
  // TODO: Stop at paramcheck time when Silo input parameters are
  // set, but Silo is not available
  recover_iteration = InputSiloParameters(recover_dir, recover_file);
#else
  CCTK_ERROR("No parameter recovery method available");
#endif

  return recover_iteration >= 0;
}

////////////////////////////////////////////////////////////////////////////////

void RecoverGridStructure(cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

#ifdef HAVE_CAPABILITY_Silo
  // TODO: Stop at paramcheck time when Silo input parameters are
  // set, but Silo is not available
  InputSiloGridStructure(cctkGH, recover_dir, recover_file, recover_iteration);
#else
  CCTK_ERROR("No grid structure recovery method available");
#endif
}

////////////////////////////////////////////////////////////////////////////////

void RecoverGH(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("RecoverGH");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root)
    CCTK_VINFO("RecoverGH: iteration %d, time %f", cctk_iteration,
               double(cctk_time));

  // Find input groups
  const vector<bool> group_enabled = [&] {
    vector<bool> enabled(CCTK_NumGroups(), false);
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
      const GHExt::GlobalData &globaldata = ghext->globaldata;
      if (globaldata.arraygroupdata.at(gi)) {
        // grid array
        enabled.at(gi) = globaldata.arraygroupdata.at(gi)->do_checkpoint;
      } else {
        // grid function
        const int level = 0;
        const GHExt::LevelData &leveldata = ghext->leveldata.at(level);
        assert(leveldata.groupdata.at(gi));
        enabled.at(gi) = leveldata.groupdata.at(gi)->do_checkpoint;
      }
    }
    return enabled;
  }();

#ifdef HAVE_CAPABILITY_Silo
  // TODO: Stop at paramcheck time when Silo input parameters are
  // set, but Silo is not available
  InputSilo(cctkGH, group_enabled, recover_dir, recover_file);
#else
  CCTK_ERROR("No grid hierarchy recovery method available");
#endif
}

void InputGH(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("InputGH");
  Interval interval(timer);

  if (filereader_ID_files[0] == '\0')
    return;

  // Find input groups
  const vector<bool> input_group = [&]() {
    vector<bool> enabled(CCTK_NumGroups(), false);
    const auto callback{
        [](const int index, const char *const optstring, void *const arg) {
          vector<bool> &enabled = *static_cast<vector<bool> *>(arg);
          enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
        }};
    CCTK_TraverseString(filereader_ID_vars, callback, &enabled,
                        CCTK_GROUP_OR_VAR);
    return enabled;
  }();

  const auto num_in_groups =
      count(input_group.begin(), input_group.end(), true);
  if (num_in_groups == 0)
    return;

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root) {
    CCTK_VINFO("InputGH: iteration %d, time %f", cctk_iteration,
               double(cctk_time));
    CCTK_VINFO("Input for groups:");
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi)
      if (input_group.at(gi))
        CCTK_VINFO("  %s", CCTK_FullGroupName(gi));
  }

#ifdef HAVE_CAPABILITY_Silo
  // TODO: Stop at paramcheck time when Silo input parameters are
  // set, but Silo is not available
  // TODO: handle multiple file names in `filereader_ID_files`
  InputSilo(cctkGH, input_group, filereader_ID_dir, filereader_ID_files);
#else
  CCTK_ERROR("No file reader method available");
#endif
}

////////////////////////////////////////////////////////////////////////////////

void OutputPlotfile(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputPlotfile");
  Interval interval(timer);

  const int numgroups = CCTK_NumGroups();
  vector<bool> group_enabled(numgroups, false);
  auto enable_group{[](int index, const char *optstring, void *callback) {
    vector<bool> &group_enabled = *static_cast<vector<bool> *>(callback);
    // CCTK_TraverseString expands the given groups into their variables
    // so I condense it down again to groups only
    group_enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
  }};
  CCTK_TraverseString(out_plotfile_groups, enable_group, &group_enabled,
                      CCTK_GROUP_OR_VAR);

  // TODO: Set the number of output files (the number of files written
  // in parallel) to reduce the total number of files
  // VisMF::SetNOutFiles(saveNFiles);

  for (int gi = 0; gi < numgroups; ++gi) {
    if (!group_enabled.at(gi))
      continue;
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;

    auto &restrict groupdata0 = *ghext->leveldata.at(0).groupdata.at(gi);
    if (groupdata0.mfab.size() > 0) {
      const int tl = 0;

      string groupname = CCTK_FullGroupName(gi);
      groupname = regex_replace(groupname, regex("::"), "-");
      for (auto &c : groupname)
        c = tolower(c);
      const string filename = [&]() {
        ostringstream buf;
        buf << out_dir << "/" << groupname << ".it" << setw(6) << setfill('0')
            << cctk_iteration;
        return buf.str();
      }();

      amrex::Vector<string> varnames(groupdata0.numvars);
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        ostringstream buf;
        buf << CCTK_VarName(groupdata0.firstvarindex + vi);
        for (int i = 0; i < tl; ++i)
          buf << "_p";
        varnames.at(vi) = buf.str();
      }

      amrex::Vector<const amrex::MultiFab *> mfabs(ghext->leveldata.size());
      amrex::Vector<amrex::Geometry> geoms(ghext->leveldata.size());
      amrex::Vector<int> iters(ghext->leveldata.size());
      amrex::Vector<amrex::IntVect> reffacts(ghext->leveldata.size());
      for (const auto &restrict leveldata : ghext->leveldata) {
        mfabs.at(leveldata.level) = &*leveldata.groupdata.at(gi)->mfab.at(tl);
        geoms.at(leveldata.level) = ghext->amrcore->Geom(leveldata.level);
        iters.at(leveldata.level) = cctk_iteration;
        reffacts.at(leveldata.level) = amrex::IntVect{2, 2, 2};
      }

      // TODO: Output all groups into a single file
      WriteMultiLevelPlotfile(filename, mfabs.size(), mfabs, varnames, geoms,
                              cctk_time, iters, reffacts);

      const bool is_root = CCTK_MyProc(nullptr) == 0;
      if (is_root) {
        const string visitname = [&]() {
          ostringstream buf;
          buf << out_dir << "/" << groupname << ".visit";
          return buf.str();
        }();
        ofstream visit(visitname, ios::app);
        assert(visit.good());
        // visit << filename << "/Header\n";
        visit << groupname << ".it" << setw(6) << setfill('0') << cctk_iteration
              << "/Header\n";
        visit.close();
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void OutputNorms(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (out_norm_vars[0] == '\0')
    return;

  static Timer timer("OutputNorms");
  Interval interval(timer);

  // Find output groups
  const vector<bool> group_enabled = [&] {
    vector<bool> enabled(CCTK_NumGroups(), false);
    const auto callback{
        [](const int index, const char *const optstring, void *const arg) {
          vector<bool> &enabled = *static_cast<vector<bool> *>(arg);
          enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
        }};
    CCTK_TraverseString(out_norm_vars, callback, &enabled, CCTK_GROUP_OR_VAR);
    if (verbose) {
      CCTK_VINFO("TSV output for groups:");
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi)
        if (group_enabled.at(gi))
          CCTK_VINFO("  %s", CCTK_FullGroupName(gi));
    }
    return enabled;
  }();
  const auto num_out_groups =
      count(group_enabled.begin(), group_enabled.end(), true);
  if (num_out_groups == 0)
    return;

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  const string sep = "\t";

  ofstream file;
  if (is_root) {
    ostringstream buf;
    buf << out_dir << "/norms.it" << setw(6) << setfill('0') << cctk_iteration
        << ".tsv";
    const string filename = buf.str();
    file = ofstream(filename);

    // get more precision for floats, could also use
    // https://stackoverflow.com/a/30968371
    file << setprecision(numeric_limits<CCTK_REAL>::digits10 + 1) << scientific;

    int col = 0;
    file << "# " << ++col << ":iteration";
    file << sep << ++col << ":time";
    file << sep << ++col << ":varname";
    file << sep << ++col << ":min";
    file << sep << ++col << ":max";
    file << sep << ++col << ":sum";
    file << sep << ++col << ":avg";
    file << sep << ++col << ":stddev";
    file << sep << ++col << ":volume";
    file << sep << ++col << ":L1norm";
    file << sep << ++col << ":L2norm";
    file << sep << ++col << ":maxabs";
    if (!out_norm_omit_unstable) {
      for (int d = 0; d < dim; ++d)
        file << sep << ++col << ":minloc[" << d << "]";
      for (int d = 0; d < dim; ++d)
        file << sep << ++col << ":maxloc[" << d << "]";
    }
    file << "\n";
  }

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    if (!group_enabled.at(gi))
      continue;
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;

    const int level = 0;
    const GHExt::LevelData &restrict leveldata = ghext->leveldata.at(level);
    const GHExt::LevelData::GroupData &restrict groupdata =
        *leveldata.groupdata.at(gi);

    const int tl = 0;
    for (int vi = 0; vi < groupdata.numvars; ++vi) {

      const reduction<CCTK_REAL, dim> red = reduce(gi, vi, tl);

      if (is_root) {
        file << cctk_iteration << sep << cctk_time << sep
             << CCTK_FullVarName(groupdata.firstvarindex + vi) << sep << red.min
             << sep << red.max << sep << red.sum << sep << red.avg() << sep
             << red.sdv() << sep << red.norm0() << sep << red.norm1() << sep
             << red.norm2() << sep << red.norm_inf();
        if (!out_norm_omit_unstable) {
          for (int d = 0; d < dim; ++d)
            file << sep << red.minloc[d];
          for (int d = 0; d < dim; ++d)
            file << sep << red.maxloc[d];
        }
        file << "\n";
      }
    }
  }

  if (is_root)
    file.close();
}

////////////////////////////////////////////////////////////////////////////////

void OutputScalars(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputScalars");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (!is_root)
    return;
  const string sep = "\t";

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_SCALAR)
      continue;

    string groupname = CCTK_FullGroupName(gi);
    groupname = regex_replace(groupname, regex("::"), "-");
    for (auto &c : groupname)
      c = tolower(c);
    ostringstream buf;
    buf << out_dir << "/" << groupname << ".tsv";
    const string filename = buf.str();

    ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::app);

    // TODO: write proper header (once)
#if 0
    // Output header
    file << "#"
         << sep
         << "iteration"
         << sep
         << "time";
    for (int vi = 0; vi < groupdata.numvars; ++vi) {
      file << sep << CCTK_VarName(group.firstvarindex + vi);
    file << "\n";
#endif

    // Output data
    const GHExt::GlobalData &restrict globaldata = ghext->globaldata;
    const GHExt::GlobalData::ArrayGroupData &restrict arraygroupdata =
        *globaldata.arraygroupdata.at(gi);
    const int tl = 0;
    file << cctkGH->cctk_iteration << sep << cctkGH->cctk_time;
    for (int vi = 0; vi < arraygroupdata.numvars; ++vi) {
      file << sep << arraygroupdata.data.at(tl).at(vi);
    }
    file << "\n";
  }
}

////////////////////////////////////////////////////////////////////////////////

struct parameters {};
YAML::Emitter &operator<<(YAML::Emitter &yaml, parameters) {
  yaml << YAML::LocalTag("parameters-1.0.0");
  yaml << YAML::BeginMap;
  int first = 1;
  for (;;) {
    const cParamData *data;
    // This call is most likely not thread safe
    int ierr = CCTK_ParameterWalk(first, nullptr, nullptr, &data);
    assert(ierr >= 0);
    if (ierr)
      break;
    const string fullname = string(data->thorn) + "::" + data->name;
    yaml << YAML::Key << fullname << YAML::Value << YAML::Flow;
    int type;
    const void *pvalue = CCTK_ParameterGet(data->name, data->thorn, &type);
    switch (type) {
    case PARAMETER_KEYWORD:
    case PARAMETER_STRING:
      yaml << *static_cast<const char *const *>(pvalue);
      break;
    case PARAMETER_INT:
      yaml << *static_cast<const CCTK_INT *>(pvalue);
      break;
    case PARAMETER_REAL:
      yaml << *static_cast<const CCTK_REAL *>(pvalue);
      break;
    case PARAMETER_BOOLEAN:
      yaml << static_cast<bool>(*static_cast<const CCTK_INT *>(pvalue));
      break;
    default:
      assert(0);
    }
    first = 0;
  }
  yaml << YAML::EndMap;
  return yaml;
}

void OutputMetadata(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!out_metadata)
    return;

  static Timer timer("OutputMetadata");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (!is_root)
    return;

  YAML::Emitter yaml;
  yaml << YAML::Comment("CarpetX");
  yaml << YAML::BeginDoc;
  yaml << YAML::LocalTag("carpetx-metadata-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "nghostzones";
  yaml << YAML::Value << YAML::Flow << YAML::BeginSeq;
  for (int d = 0; d < dim; ++d)
    yaml << cctk_nghostzones[d];
  yaml << YAML::EndSeq;
  // Note: origin_space and delta_space are nan at this point
  // yaml << YAML::Key << "origin_space";
  // yaml << YAML::Value << YAML::Flow << YAML::BeginSeq;
  // for (int d = 0; d < dim; ++d)
  //   yaml << cctk_origin_space[d];
  // yaml << YAML::EndSeq;
  // yaml << YAML::Key << "delta_space";
  // yaml << YAML::Value << YAML::Flow << YAML::BeginSeq;
  // for (int d = 0; d < dim; ++d)
  //   yaml << cctk_delta_space[d];
  // yaml << YAML::EndSeq;
  yaml << YAML::Key << "iteration";
  yaml << YAML::Value << cctk_iteration;
  yaml << YAML::Key << "time";
  yaml << YAML::Value << cctk_time;
  yaml << YAML::Key << "delta_time";
  yaml << YAML::Value << cctk_delta_time;
  yaml << YAML::Key << "ghext";
  yaml << YAML::Value << *ghext;
  yaml << YAML::Key << "parameters";
  yaml << YAML::Value << parameters();
  yaml << YAML::EndMap;
  yaml << YAML::EndDoc;

  ostringstream buf;
  buf << out_dir << "/metadata.yaml";
  const string filename = buf.str();

  ofstream file(filename.c_str(), std::ofstream::out);
  file << yaml.c_str();
}

////////////////////////////////////////////////////////////////////////////////

int OutputGH(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputGH");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root)
    CCTK_VINFO("OutputGH: iteration %d, time %f, run time %d s", cctk_iteration,
               double(cctk_time), CCTK_RunTime());

  if (out_every > 0 && cctk_iteration % out_every == 0) {
    vector<future<void> > tasks;
    const auto launch = [&](const auto &fun) {
      tasks.push_back(async(std::launch::async, fun));
    };

    launch([&]() { OutputMetadata(cctkGH); });

    OutputNorms(cctkGH);

    OutputPlotfile(cctkGH);

    OutputScalars(cctkGH);

#ifdef HAVE_CAPABILITY_Silo
    // TODO: Stop at paramcheck time when Silo output parameters are
    // set, but Silo is not available

    // Find output groups
    const vector<bool> group_enabled = [&] {
      vector<bool> enabled(CCTK_NumGroups(), false);
      const auto callback{
          [](const int index, const char *const optstring, void *const arg) {
            vector<bool> &enabled = *static_cast<vector<bool> *>(arg);
            enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
          }};
      CCTK_TraverseString(out_silo_vars, callback, &enabled, CCTK_GROUP_OR_VAR);
      if (verbose) {
        CCTK_VINFO("Silo output for groups:");
        for (int gi = 0; gi < CCTK_NumGroups(); ++gi)
          if (enabled.at(gi))
            CCTK_VINFO("  %s", CCTK_FullGroupName(gi));
      }
      return enabled;
    }();

    // Obtain the parameter file name
    const string parfilename = [&] {
      vector<char> buf(10000);
      int ilen = CCTK_ParameterFilename(buf.size(), buf.data());
      assert(ilen < int(buf.size() - 1));
      string parfilename(buf.data());
      // Remove directory prefix, if any
      auto slash = parfilename.rfind('/');
      if (slash != string::npos)
        parfilename = parfilename.substr(slash + 1);
      // Remove suffix, if it is there
      auto suffix = parfilename.rfind('.');
      if (suffix != string::npos && parfilename.substr(suffix) == ".par")
        parfilename = parfilename.substr(0, suffix);
      return parfilename;
    }();

    const string simulation_name = [&]() {
      string name = parfilename;
      const size_t last_slash = name.rfind('/');
      if (last_slash != string::npos && last_slash < name.length())
        name = name.substr(last_slash + 1);
      const size_t last_dot = name.rfind('.');
      if (last_dot != string::npos && last_dot > 0)
        name = name.substr(0, last_dot);
      return name;
    }();

    OutputSilo(cctkGH, group_enabled, out_dir, simulation_name);
#endif

    OutputTSVold(cctkGH);

    OutputTSV(cctkGH);

    for (auto &task : tasks)
      task.wait();
  }

  // TODO: This should be the number of variables output
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void Checkpoint(const cGH *const restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static int last_checkpoint_iteration = -1;

  if (cctkGH->cctk_iteration <= last_checkpoint_iteration) {
    CCTK_VINFO(
        "Already wrote checkpoint at iteration %d; skipping checkpointing",
        cctkGH->cctk_iteration);
    return;
  }
  last_checkpoint_iteration = cctkGH->cctk_iteration;

  static Timer timer("Checkpoint");
  Interval interval(timer);

#ifdef HAVE_CAPABILITY_Silo
  // TODO: Stop at paramcheck time when Silo output parameters are
  // set, but Silo is not available

  const vector<bool> checkpoint_group = [&] {
    vector<bool> enabled(CCTK_NumGroups(), false);
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
      const GHExt::GlobalData &globaldata = ghext->globaldata;
      if (globaldata.arraygroupdata.at(gi)) {
        // grid array
        enabled.at(gi) = globaldata.arraygroupdata.at(gi)->do_checkpoint;
      } else {
        // grid function
        const int level = 0;
        const GHExt::LevelData &leveldata = ghext->leveldata.at(level);
        assert(leveldata.groupdata.at(gi));
        enabled.at(gi) = leveldata.groupdata.at(gi)->do_checkpoint;
      }
    }
    return enabled;
  }();

  OutputSilo(cctkGH, checkpoint_group, checkpoint_dir, checkpoint_file);
#else
  CCTK_ERROR("No checkpointing method available");
#endif
}

extern "C" void CarpetX_CheckpointInitial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (checkpoint_ID) {
    CCTK_VINFO("Checkpointing initial conditions at iteration %d, time %f",
               cctk_iteration, double(cctk_time));
    Checkpoint(cctkGH);
  }
}

extern "C" void CarpetX_Checkpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static int last_checkpoint_runtime = -1; // seconds

  const int runtime = CCTK_RunTime(); // seconds

  const bool checkpoint_by_iteration =
      checkpoint_every > 0 && cctk_iteration % checkpoint_every == 0;
  const bool checkpoint_by_walltime =
      checkpoint_every_walltime_hours > 0 &&
      runtime >= last_checkpoint_runtime +
                     lrint(checkpoint_every_walltime_hours * 3600);
  last_checkpoint_runtime = runtime;

  if (checkpoint_by_iteration || checkpoint_by_walltime) {
    CCTK_VINFO("Checkpointing at iteration %d, time %f, run time %.2f h",
               cctk_iteration, double(cctk_time), double(runtime / 3600));
    Checkpoint(cctkGH);
  }
}

extern "C" void CarpetX_CheckpointTerminate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (checkpoint_on_terminate) {
    CCTK_VINFO("Checkpointing before terminating at iteration %d, time %f",
               cctk_iteration, double(cctk_time));
    Checkpoint(cctkGH);
  }
}

} // namespace CarpetX
