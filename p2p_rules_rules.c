#include <string.h>
const char *ApplyRules(float mQ, float vQ, float sQ) {
  if ((mQ > 89.002500) && (sQ > 1.433302 && sQ <= 21.677047)) return "tunnel";
  if ((vQ > 298.788268 && vQ <= 4038.600630) && (sQ <= 30.038342)) return "tunnel";
  if ((mQ <= 89.899500)) return "no tunnel";
  if ((vQ <= 139.285840)) return "no tunnel";

  return "no tunnel";
}
