#include <string.h>
const char *ApplyRules(float mA, float mQ, float vDt, float vA, float vQ, float sDt, float sA, float sQ, float kDt, float kA, float kQ) {
  if ((sA > 1.900437) && (sQ <= 1.620706)) return "no tunnel";
  if ((vQ <= 6128.384798) && (sA > 1.680212 && sA <= 1.914077) && (kA <= 3.496445)) return "no tunnel";
  if ((sA > 1.658121 && sA <= 1.899986) && (sQ <= 1.380559)) return "no tunnel";
  if ((sA > 1.395776) && (sQ <= 0.366934) && (kA <= 3.962453)) return "no tunnel";
  if ((mA <= 317.144000) && (sDt <= 13.726500) && (sA > 1.615381 && sA <= 3.154511) && (sQ > 15.512480) && (kQ > -2.999931)) return "no tunnel";
  if ((vQ > 2092.184342) && (sA > 1.878432 && sA <= 3.154511)) return "no tunnel";
  if ((sA <= 1.402062)) return "tunnel";
  if ((sA <= 1.536331) && (kA > 1.095606)) return "tunnel";
  if ((vA <= 42240.751998) && (sA <= 1.632753) && (kDt > 46.624990) && (kA > 0.950843) && (kQ > -2.796590 && kQ <= 10670.354038)) return "tunnel";
  if ((mA > 263.456000) && (vDt > 4.418289) && (vA <= 38950.853523) && (sA <= 1.888413) && (kQ <= -2.358466)) return "tunnel";
  if ((sQ > 1.574043 && sQ <= 17.846329)) return "tunnel";
  if ((vQ > 56.362098 && vQ <= 60.987803) && (sQ > 0.298086)) return "tunnel";
  if ((mA > 237.683500) && (vDt <= 14.958344) && (vA <= 37338.200630) && (vQ > 46.533357 && vQ <= 58.356108) && (sQ > -0.008318)) return "tunnel";

  return "no tunnel";
}
