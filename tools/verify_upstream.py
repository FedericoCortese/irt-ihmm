"""Verify the byte-preserved upstream source snapshot (standard library only)."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / "vendor/manifest.json").read_text())
for name, expected in manifest["sha256"].items():
    actual = hashlib.sha256((root / "vendor/psbpHMM" / name).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit(f"Upstream source hash mismatch: {name}")
print(f"Verified {len(manifest['sha256'])} source files at {manifest['commit']}")


def cpp_function(text, name):
    start = text.index("List " + name + "(")
    brace = text.index("{", start)
    end, depth = brace + 1, 1
    while depth:
        depth += (text[end] == "{") - (text[end] == "}")
        end += 1
    return text[start:end]


original = (root / "vendor/psbpHMM/src/tester.cpp").read_text()
derived = (root / "src/hoskovec.cpp").read_text()
expected_pi = cpp_function(original, "updatePi").replace("List updatePi(", "List hoskovec_pi(")
expected_alpha = cpp_function(original, "up_ajk").replace("List up_ajk(", "List hoskovec_alpha(")
expected_alpha = (expected_alpha.replace("int K, int n,", "int K, int n_sticks, int n,")
                  .replace("NumericVector ajk(K);", "NumericVector ajk(n_sticks);")
                  .replace("k < K;", "k < n_sticks;")
                  .replace("zi[t]>=(k+1) & zi[t-1]", "zi[t]>=(k+1) && zi[t-1]"))
if cpp_function(derived, "hoskovec_pi") != expected_pi:
    raise SystemExit("Unexpected change to imported probability kernel")
if cpp_function(derived, "hoskovec_alpha") != expected_alpha:
    raise SystemExit("Unexpected change to adapted intercept kernel")
print("Verified the two reused kernels against their documented transformations")
