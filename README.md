# Mathematical models of primary saliva secretion

This repository contains MATLAB research code, parameter data, and reconstructed three-dimensional geometry used in mechanistic modelling of calcium-regulated primary saliva secretion in parotid acinar cells. The work spans single-cell models, anatomically reconstructed cells, and a seven-cell acinus developed across a sequence of peer-reviewed studies.

## Scientific question

How much spatial and multicellular biological detail is required to predict fluid secretion from parotid acinar cells?

The models couple intracellular Ca²⁺ signalling to ion transport, membrane fluxes, cell-volume regulation, luminal transport, and water flow. Increasingly detailed anatomical representations were used to test which spatial features materially alter secretory output.

## Main findings represented by this archive

- In the anatomically accurate single-cell model, spatial heterogeneity of intracellular Ca²⁺ in the apical or basal regions had little effect on the predicted primary saliva secretion rate.
- The secretion rate was governed mainly by mean intracellular Ca²⁺ rather than the frequency of Ca²⁺ oscillations.
- For flow-rate prediction alone, much of the underlying spatial and temporal complexity could be reduced to a substantially simpler representation.
- The multicellular extension coupled seven reconstructed acinar cells through a common lumen and intercellular signalling.
- In that model, acinus topology had little effect on total fluid secretion, suggesting that a detailed multicellular spatial representation is not required when total flow is the quantity of interest.

## Start here

The repository has been organised so that the original research code remains intact while the scientific workflow is easier to inspect.

| Path | Contents |
| --- | --- |
| `model/` | MATLAB source, example scripts, initial conditions, and the archived `PSec.mat` parameter data |
| `geometry/` | Reconstructed cell and lumen geometry, including STL files |
| `results/` | Archived figures and fitted-model outputs |
| `docs/` | Poster and supporting research material retained with the project |
| `CITATION.cff` | Citation metadata for the repository |
| `LICENSE` | MIT licence |

Useful entry points in `model/` include

- `ND_ex.m` — single-cell simulation example using `ND.m`
- `ND_ex_closed.m` — comparison script using the closed and open single-cell formulations
- `Open_ex.m` — parameter sweep for the open formulation
- `Cluster_Secretion_ex.m` — seven-cell acinus simulation using `Cell.m`, `Cluster.m`, `Cluster_Secretion.m`, and `PSec.mat`

For the multicellular example, set the MATLAB working directory to `model/` before running `Cluster_Secretion_ex.m` so that the archived parameter file is found in the expected location.

## Reproducibility status

This is an archive of original research code rather than a packaged software release. The MATLAB version and toolbox environment were not pinned when the code was deposited, and the repository does not contain an automated test suite. The present cleanup changes repository organisation and documentation only. The numerical model source has not been rewritten.

## Associated publications

The repository sits within the following sequence of work on mathematical modelling of salivary secretion.

1. Vera-Sigüenza E, Catalán MA, Peña-Münzenmayer G, Melvin JE, Sneyd J. **A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion.** *Bulletin of Mathematical Biology* 80, 255–282, 2018. https://doi.org/10.1007/s11538-017-0370-6

2. Vera-Sigüenza E, Pages N, Rugis J, Yule DI, Sneyd J. **A Mathematical Model of Fluid Transport in an Accurate Reconstruction of Parotid Acinar Cells.** *Bulletin of Mathematical Biology* 81, 699–721, 2019. https://doi.org/10.1007/s11538-018-0534-z

3. Pages N, Vera-Sigüenza E, Rugis J, Kirk V, Yule DI, Sneyd J. **A Model of Ca²⁺ Dynamics in an Accurate Reconstruction of Parotid Acinar Cells.** *Bulletin of Mathematical Biology* 81, 1394–1426, 2019. https://doi.org/10.1007/s11538-018-00563-z

4. Vera-Sigüenza E, Pages N, Rugis J, Yule DI, Sneyd J. **A Multicellular Model of Primary Saliva Secretion in the Parotid Gland.** *Bulletin of Mathematical Biology* 82, article 38, 2020. https://doi.org/10.1007/s11538-020-00712-3

5. Sneyd J, Vera-Sigüenza E, Rugis J, Pages N, Yule DI. **Calcium Dynamics and Water Transport in Salivary Acinar Cells.** *Bulletin of Mathematical Biology* 83, article 31, 2021. https://doi.org/10.1007/s11538-020-00841-9

## Citation

If you use a specific model, please cite the corresponding paper above. Repository-level citation metadata are also provided in `CITATION.cff`.

## Licence

The code in this repository is released under the MIT License.