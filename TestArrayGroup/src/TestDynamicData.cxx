#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <cassert>
#include <cmath>
#include <memory>
#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <tuple>
#include <type_traits>
#include <vector>
#include <map>

extern "C" void TestArrayGroup_DynamicData(cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  std::map<std::string, int> error = {
      {"lsh", 0},  {"ash", 0},  {"gsh", 0},        {"lbnd", 0},
      {"ubnd", 0}, {"bbox", 0}, {"nghostzones", 0}};

  CCTK_INT gf_gi, array_gi, scalar_gi;
  cGroupDynamicData gf_data, array_data, scalar_data;

  gf_gi = CCTK_GroupIndex("TestArrayGroup::test_gf");
  array_gi = CCTK_GroupIndex("TestArrayGroup::test_array");
  scalar_gi = CCTK_GroupIndex("TestArrayGroup::test_scalar");

  CCTK_GroupDynamicData(cctkGH, gf_gi, &gf_data);
  CCTK_GroupDynamicData(cctkGH, array_gi, &array_data);
  CCTK_GroupDynamicData(cctkGH, scalar_gi, &scalar_data);

  // Validate grid function dynamic data
  cGroup gf_group;
  int ierr = CCTK_GroupData(gf_gi, &gf_group);
  if (ierr)
    CCTK_ERROR("TestArrayGroup: error in GroupData for grid functions");

  if (gf_data.dim != gf_group.dim || gf_data.dim <= 0)
    CCTK_ERROR(
        "TestArrayGroup: incorrect dimension in grid function dynamic data");
  if (gf_data.activetimelevels != CCTK_ActiveTimeLevelsGI(cctkGH, gf_gi))
    CCTK_ERROR("TestArrayGroup: incorrect activetimelevels in grid function "
               "dynamic data");
  for (int i = 0; i < gf_data.dim; i++) {
    if (gf_data.lsh[i] != cctkGH->cctk_lsh[i])
      error.at("lsh") += 1;
    if (gf_data.ash[i] != cctkGH->cctk_ash[i])
      error.at("ash") += 1;
    if (gf_data.gsh[i] != cctkGH->cctk_gsh[i])
      error.at("gsh") += 1;
    if (gf_data.lbnd[i] != cctkGH->cctk_lbnd[i])
      error.at("lbnd") += 1;
    if (gf_data.ubnd[i] != cctkGH->cctk_ubnd[i])
      error.at("lbnd") += 1;
    if (gf_data.bbox[2 * i] != cctkGH->cctk_bbox[2 * i])
      error.at("bbox") += 1;
    if (gf_data.bbox[2 * i + 1] != cctkGH->cctk_bbox[2 * i + 1])
      error.at("bbox") += 1;
    if (gf_data.nghostzones[i] != cctkGH->cctk_nghostzones[i])
      error.at("nghostzones") += 1;
  }
  if (error.at("lsh"))
    CCTK_ERROR(
        "TestArrayGroup: incorrect lsh data in grid function dynamic data");
  if (error.at("ash"))
    CCTK_ERROR(
        "TestArrayGroup: incorrect ash data in grid function dynamic data");
  if (error.at("gsh"))
    CCTK_ERROR(
        "TestArrayGroup: incorrect gsh data in grid function dynamic data");
  if (error.at("lbnd"))
    CCTK_ERROR(
        "TestArrayGroup: incorrect lbnd data in grid function dynamic data");
  if (error.at("ubnd"))
    CCTK_ERROR(
        "TestArrayGroup: incorrect ubnd data in grid function dynamic data");
  if (error.at("bbox"))
    CCTK_ERROR(
        "TestArrayGroup: incorrect bbox data in grid function dynamic data");
  if (error.at("nghostzones"))
    CCTK_ERROR("TestArrayGroup: incorrect nghostzones data in grid function "
               "dynamic data");

  // Validate grid scalar dynamic data
  for (auto &[key, val] : error) {
    val = 0;
  }

  if (scalar_data.dim != 0)
    CCTK_ERROR("TestArrayGroup: incorrect dimension in scalar dynamic data");
  if (scalar_data.activetimelevels != 1)
    CCTK_ERROR(
        "TestArrayGroup: incorrect activetimelevels in scalar dynamic data");
  for (int i = 0; i < gf_data.dim; i++) {
    if (scalar_data.lsh[i] != -1)
      error.at("lsh") += 1;
    if (scalar_data.ash[i] != -1)
      error.at("ash") += 1;
    if (scalar_data.gsh[i] != -1)
      error.at("gsh") += 1;
    if (scalar_data.lbnd[i] != -1)
      error.at("lbnd") += 1;
    if (scalar_data.ubnd[i] != -1)
      error.at("lbnd") += 1;
    if (scalar_data.bbox[2 * i] != -1)
      error.at("bbox") += 1;
    if (scalar_data.bbox[2 * i + 1] != -1)
      error.at("bbox") += 1;
    if (scalar_data.nghostzones[i] != -1)
      error.at("nghostzones") += 1;
  }
  if (error.at("lsh"))
    CCTK_ERROR("TestArrayGroup: incorrect lsh data in scalar dynamic data");
  if (error.at("ash"))
    CCTK_ERROR("TestArrayGroup: incorrect ash data in scalar dynamic data");
  if (error.at("gsh"))
    CCTK_ERROR("TestArrayGroup: incorrect gsh data in scalar dynamic data");
  if (error.at("lbnd"))
    CCTK_ERROR("TestArrayGroup: incorrect lbnd data in scalar dynamic data");
  if (error.at("ubnd"))
    CCTK_ERROR("TestArrayGroup: incorrect ubnd data in scalar dynamic data");
  if (error.at("bbox"))
    CCTK_ERROR("TestArrayGroup: incorrect bbox data in scalar dynamic data");
  if (error.at("nghostzones"))
    CCTK_ERROR(
        "TestArrayGroup: incorrect nghostzones data in scalar dynamic data");

  // Validate grid array dynamic data
  for (auto &[key, val] : error) {
    val = 0;
  }

  cGroup array_group;
  ierr = CCTK_GroupData(array_gi, &array_group);
  if (ierr)
    CCTK_ERROR("TestArrayGroup: error in GroupData for arrays");

  int sz[2] = {5, 6};

  if (array_data.dim != array_group.dim || array_data.dim <= 0)
    CCTK_ERROR("TestArrayGroup: incorrect dimension in array dynamic data");
  if (scalar_data.activetimelevels != 1)
    CCTK_ERROR(
        "TestArrayGroup: incorrect activetimelevels in array dynamic data");
  for (int i = 0; i < array_data.dim; i++) {
    if (array_data.lsh[i] != sz[i])
      error.at("lsh") += 1;
    if (array_data.ash[i] != sz[i])
      error.at("ash") += 1;
    if (array_data.gsh[i] != sz[i])
      error.at("gsh") += 1;
    if (array_data.lbnd[i] != 0)
      error.at("lbnd") += 1;
    if (array_data.ubnd[i] != sz[i] - 1)
      error.at("lbnd") += 1;
    if (array_data.bbox[2 * i] != 1)
      error.at("bbox") += 1;
    if (array_data.bbox[2 * i + 1] != 1)
      error.at("bbox") += 1;
    if (array_data.nghostzones[i] != 0)
      error.at("nghostzones") += 1;
  }
  if (error.at("lsh"))
    CCTK_ERROR("TestArrayGroup: incorrect lsh data in array dynamic data");
  if (error.at("ash"))
    CCTK_ERROR("TestArrayGroup: incorrect ash data in array dynamic data");
  if (error.at("gsh"))
    CCTK_ERROR("TestArrayGroup: incorrect gsh data in array dynamic data");
  if (error.at("lbnd"))
    CCTK_ERROR("TestArrayGroup: incorrect lbnd data in array dynamic data");
  if (error.at("ubnd"))
    CCTK_ERROR("TestArrayGroup: incorrect ubnd data in array dynamic data");
  if (error.at("bbox"))
    CCTK_ERROR("TestArrayGroup: incorrect bbox data in array dynamic data");
  if (error.at("nghostzones"))
    CCTK_ERROR(
        "TestArrayGroup: incorrect nghostzones data in array dynamic data");

  return;
}
