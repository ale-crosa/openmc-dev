import re
from pathlib import Path
import matplotlib.pyplot as plt

surfaces = ["Schwarz_p", "Gyroid", "Diamond"]
pitches = [2, 5, 10, 20]
files = ["res_0.125.out", "res_0.062.out", "res_0.031.out", "res_0.015.out"]

pattern = re.compile(r"Combined k-effective\s*=\s*([0-9.]+)")

fig, axs = plt.subplots(1, 3, figsize=[12.8,4.8])

for surface, ax in zip(surfaces, axs):
    for fname in files:
        keffs = []

        for pitch in pitches:
            text = Path(f"{surface}/{pitch}", fname).read_text()

            match = pattern.search(text)
            if not match:
                raise ValueError(f"No k-effective found in {pitch}/{fname}")

            keffs.append(float(match.group(1)))

        ax.plot(pitches, keffs, marker="o", label=fname)
    ax.set_title(surface)
    ax.set_xlabel("pitch")
    ax.set_ylabel("keff")
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()