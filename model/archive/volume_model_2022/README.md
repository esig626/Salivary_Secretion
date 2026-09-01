# Archived 2022 standalone volume model

This directory preserves a later standalone implementation of the salivary fluid-secretion and cell-volume equations that was previously stored in the separate `Dynamical_systems` repository.

The code was written in 2022 at the Institute of Metabolism and Systems Research, University of Birmingham. It retains the same salivary-physiology components seen in the earlier publication sequence, including Na⁺/K⁺ ATPase, NKCC1, Ae4, Nhe1, Ca²⁺-activated Cl⁻ and K⁺ channels, tight-junction ion transport, osmotic water fluxes, intracellular ion concentrations, and cell volume.

## Contents

- `Volume.m` implements a single-cell volume and secretion model with a prescribed change in intracellular Ca²⁺.
- `Volume2.m` is a geometry-aware variant that distributes luminal variables across cell-neighbour interfaces.
- `Volume_ex.m` is an executable `ode15s` example for `Volume.m`.
- `Par_V.mat` contains the archived parameters used by `Volume_ex.m`.
- `P2.mat` is an additional small archived MATLAB parameter/data file associated with the geometry-aware implementation.

## Large historical parameter archive

The original `Volume_Model` directory also contained `P.mat`, a roughly 579 kB MATLAB binary created in 2018. It is retained byte-for-byte in the archived `Dynamical_systems` repository rather than being regenerated or converted during this consolidation. `Volume_ex.m`, the runnable standalone example documented here, does not load `P.mat`.

## Status

These files are retained as a historical research implementation rather than presented as the canonical code for one of the peer-reviewed salivary-gland papers. Every source file and MATLAB binary transferred here has the same Git blob hash as the original. Their source was `esig626/Dynamical_systems`, commit `d95cacfad8890a35e9a533e0e847921b128f8693`.
