#!/usr/bin/env python3
"""The rules a settings file has to obey to be worth exchanging.

WHY THIS EXISTS. A settings file is meant to travel: to another bracket, and to
another person. Every rule below is about the RECIPIENT, and every one of them
is the kind that fails silently -- the file loads, the picture is wrong, and
nothing says so. So each is pinned here rather than trusted.

  1. A key the file omits goes back to its DEFAULT, not to whatever the
     recipient's slider happened to be on. Otherwise the same file lands
     differently on two machines and neither of them is wrong.
  2. A key this build does not know is ignored, and the file still loads. A
     preset library has to survive the next release.
  3. Settings tied to one image's geometry -- the disc, the black point, the
     diamond ring -- are applied only when the file is loaded back into the
     folder it came from. Carried across, they do not carry a look, they carry
     a mistake.
  4. A file written before 0.22.84 says so, because Amplification was
     renormalised there and its number means something else now.
  5. The .params.json the export has always written beside an image loads too.
     People have those already, next to every render they liked.

    python tools/smoke_settings.py
"""
import json, os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from eclipseforgehdr import server
from eclipseforgehdr.render import DEFAULTS, classify_settings

app = server.app.test_client()
tmp = tempfile.mkdtemp()
server.STATE["folder"] = tmp

# a look: some portable, some local
look = dict(DEFAULTS)
look.update({"clarity": 0.80, "hillMix": 0.5, "hillGain": 0.12,
             "discLevel": 0.111, "ringBlend": 0.7, "bgBlack": 0.099})

r = app.post("/api/settings/save", json={"params": look, "name": "look1"})
j = r.get_json(); assert j["ok"], j
print("saved: %d portable, %d local -> %s" % (j["n_portable"], j["n_local"],
                                              os.path.basename(j["path"])))
path = j["path"]

# 1. same folder -> everything comes back
r = app.post("/api/settings/load", json={"path": path}); j = r.get_json()
assert j["ok"] and j["same_stack"], j
assert j["params"]["discLevel"] == 0.111 and j["params"]["clarity"] == 0.80
print("same folder : local settings applied (discLevel %.3f)" % j["params"]["discLevel"])

# 2. different folder -> local skipped, portable kept
server.STATE["folder"] = tempfile.mkdtemp()
r = app.post("/api/settings/load", json={"path": path}); j = r.get_json()
assert j["ok"] and not j["same_stack"], j
assert j["params"]["clarity"] == 0.80, "portable lost"
assert j["params"]["discLevel"] == DEFAULTS["discLevel"], "local leaked across stacks"
print("other folder: portable kept, local reset to default (%.3f)" % j["params"]["discLevel"])
print("              note: %s" % [n for n in j["notes"] if "skipped" in n][0][:70])

# 3. a key the file omits must go back to DEFAULT, not stay where it was
doc = json.load(open(path)); del doc["portable"]["clarity"]
json.dump(doc, open(path, "w"))
r = app.post("/api/settings/load", json={"path": path}); j = r.get_json()
assert j["params"]["clarity"] == DEFAULTS["clarity"], "omitted key kept a stale value"
print("omitted key : reset to default (%.2f), not left at 0.80" % j["params"]["clarity"])

# 4. unknown keys ignored, old version warns about hillGain
doc["portable"]["somethingFromTheFuture"] = 1.0
doc["version"] = "0.22.83"
json.dump(doc, open(path, "w"))
r = app.post("/api/settings/load", json={"path": path}); j = r.get_json()
assert j["ok"], j
assert any("unknown" in n for n in j["notes"]), j["notes"]
assert any("renormalised" in n for n in j["notes"]), j["notes"]
print("future key  : ignored, file still loads")
print("old version : %s" % [n for n in j["notes"] if "renormalised" in n][0][:88] + "...")

# 5. the plain .params.json sidecar the export has always written
side = os.path.join(tmp, "x.tif.params.json")
json.dump(look, open(side, "w"))
r = app.post("/api/settings/load", json={"path": side}); j = r.get_json()
assert j["ok"] and j["params"]["clarity"] == 0.80, j
assert any("no version" in n for n in j["notes"]), j["notes"]
print("export sidecar: loads too, flagged as carrying no version")


# 6. THE PROCESSING RECIPE IS A SEPARATE FILE, and the separation is the point.
#
# A look can be dropped on any bracket and seen at once. A recipe decides how
# the data is STACKED and means nothing until the next Start. Mixing them would
# make loading a look a decision about re-stacking, so:
#
#   6a. the look file carries no processing keys at all;
#   6b. the recipe is SET, not applied, and reports which keys differ from what
#       this folder was actually stacked with;
#   6c. a key the recipe omits is LEFT ALONE -- the look file's "unmentioned
#       goes back to default" rule is right there and wrong here;
#   6d. neither file loads through the other's endpoint.
doc = json.load(open(path))
assert "processing" not in doc, "the look file must not carry processing keys"
print("look file     : carries no processing keys")

r = app.post("/api/processing/save", json={
    "name": "recipe1",
    "processing": {"align_corr": "cross", "stack_combine": "clip",
                   "fnrgf_preset": "published"}})
j = r.get_json(); assert j["ok"], j
ppath = j["path"]
print("recipe saved  : %d setting(s) -> %s" % (j["n"], os.path.basename(ppath)))

pdoc = json.load(open(ppath))
pdoc["processing"]["notAKeyWeKnow"] = "x"
json.dump(pdoc, open(ppath, "w"))
r = app.post("/api/processing/load", json={"path": ppath}); j = r.get_json()
assert j["ok"], j
pr = j["processing"]
assert "notAKeyWeKnow" not in pr, "unknown key was not dropped"
assert pr["align_corr"] == "cross" and pr["stack_combine"] == "clip", pr
assert "tier_mode" not in pr, "a key the file omits must not be invented"
assert (j.get("ids") or {}).get("align_corr") == "alignCorr", j.get("ids")
assert set(j["changed"]) >= {"align_corr", "stack_combine", "fnrgf_preset"}, j["changed"]
assert any("does not know" in n for n in j["notes"]), j["notes"]
print("recipe loaded : %d set, %d flagged for a re-stack, unknown key dropped"
      % (j["n"], len(j["changed"])))
print("omitted key   : left alone, not reset (tier_mode absent from the reply)")

r = app.post("/api/processing/load", json={"path": path}); j = r.get_json()
assert not j["ok"] and "not a processing file" in j["error"], j
r = app.post("/api/settings/load", json={"path": ppath}); j = r.get_json()
assert j["ok"] and not j.get("processing"), \
    "a recipe loaded as a look must not smuggle processing keys through"
print("cross-loading : a look is refused as a recipe; a recipe carries no look")

print("\nOK — all six rules hold, and the two file types stay apart")
