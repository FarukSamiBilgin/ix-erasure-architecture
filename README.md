# An erasure-biased architecture for interlayer-exciton position qubits in a gate-defined moiré trap lattice

> **Status: in revision.** An internal adversarial audit found issues that are being corrected; the summary below reports only results that survived it.

Companion repository: [heterotmd-purcell-emitters](https://github.com/FarukSamiBilgin/heterotmd-purcell-emitters).
System-level feasibility study: [moire-300-logical-qubit-feasibility](https://github.com/FarukSamiBilgin/moire-300-logical-qubit-feasibility).

## Question

Can the *position* of a single dark interlayer exciton (IX) in a gate-defined double trap on an aligned TMD heterobilayer serve as a qubit, and can its dominant failure (loss) be converted into a flagged erasure so that an erasure-tailored code does the heavy lifting?

## Architecture concept

- **Qubit:** dual-rail position of one dark H-stacked IX in a double trap (|L⟩/|R⟩); all logic is electrical.
- **Two-qubit gate:** a "polarize-and-phase" CZ in which a shared vertical-field pulse switches the pair from a quadrupole-parked sweet spot to a dipolar ZZ interaction, with a mid-gate polarity-flip echo.
- **Erasure conversion:** IX loss is detected by a photon-free electrical occupancy probe (a cavity click is an optional extra flag); the fast loop is photon-free, not readout-free.

## Results that survived the audit

1. **Double-trap tunnelling map:** an OFF-state tunnelling window t ≤ 0.05 µeV (bit flip < 10⁻⁵ per round) exists below the 200 meV trap-depth line, with a barrier-gate X gate of 10–213 ps and < 10⁻³ leakage (`figures/fig2_tunel.*`).
2. **Screening by image charges:** explicit image-charge sums show that a single bottom gate does **not** screen the vertical-dipole coupling (it stays at 4.1–6.8 µeV for gate distances 5–20 nm, the far field even enhanced up to 2×), whereas a dual-gate sandwich cuts beyond-neighbour crosstalk exponentially (`figures/fig3_ekran.*`).
3. **Echo against quasi-static disorder:** Monte Carlo gate simulations show that the mid-gate polarity-flip echo cancels quasi-static site-energy disorder exactly; the 10⁻³ CZ target then reduces to a single specification, an ns-correlated noise floor below ≈ 0.4 µeV. A Jaksch-blockade variant is percent-class at any realistic trim residual.

## Open problems (being corrected)

- **Trap geometry vs. gate resolution:** the trap geometry is below the resolution of the specified gate stack and needs gate amplitudes of several volts; the CZ timing is therefore not quoted here.
- **Electrode count:** one trim line per trap is required, and there is no analog memory at 20 nm pitch; crossbar addressing is not claimed. The system-level study in [moire-300-logical-qubit-feasibility](https://github.com/FarukSamiBilgin/moire-300-logical-qubit-feasibility) finds per-site/per-bond control layers fatally limited by disorder for this platform family; this architecture does not claim otherwise.
- **Threshold margin:** the comparison with reported erasure-code thresholds, the idle Z error per round, two-phonon channels, and the low-barrier tunnel-contact reset are being reworked.
- The permanent dipole used is a value measured in a different (2H MoSe₂/WSe₂) stack; its transfer to the proposed stack is an assumption.

## Repository layout

- `scripts/mimari_hesap.py` (architecture budgets), `scripts/fig_uret.py`, `scripts/requirements.txt`. Some comments are in Turkish.
- `results/mimari_hesap.json`: raw output of the budget script (includes quantities still under revision; the README above is authoritative on what survived).
- `figures/`: tunnelling map and screening figures.
- `paper/ABSTRACT.md`: extended abstract (manuscript unpublished).

## AI assistance

AI assistance: AI coding assistants were used for code development and for adversarial review. All load-bearing results were re-derived by independent methods and checked by automated tests; the author is responsible for every claim.

## Citation

Faruk Sami Bilgin (Independent Researcher, ORCID [0009-0001-9064-4410](https://orcid.org/0009-0001-9064-4410)), unpublished manuscript in revision, 2026. See `CITATION.cff`.

License: code MIT (`LICENSE`); text, figures and data CC BY 4.0 (`LICENSE-CONTENT`).
