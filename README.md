# Network Diffusion Drives Multistable Cooperation Patterns in Eco-Evolutionary Games

MATLAB source code accompanying the research manuscript **“Network Diffusion Drives Multistable Cooperation Patterns in Eco-Evolutionary Games”** (prepared for *Proceedings of the Royal Society B: Biological Sciences*).

This repository implements a networked eco-evolutionary public-goods game, evaluates the homogeneous equilibrium and its linear stability, calculates network-specific dispersion relations, and simulates the evolution of cooperation and public goods on different graph topologies. It includes examples for both diffusion-driven pattern formation and diffusion-mediated changes in the long-term dynamical state.

> **Scope and reproducibility.** The source code is provided for inspection and numerical reproduction. Pre-generated network `.mat` files are **not included** in this package; generate them first as described below. Network realizations and random initial perturbations can alter the exact trajectories. The included examples illustrate specific numerical configurations but do not automatically reconstruct every figure in the manuscript. MATLAB execution has not been independently verified in the environment used to prepare this README.

## 1. Requirements

- **MATLAB:** R2018b or newer is recommended for compatibility with `isfile`, `isfolder`, `uniquetol`, and the other functions used here. The exact minimum release and numerical equivalence across MATLAB releases have not been tested.
- **Toolboxes:** No third-party packages or specialized MATLAB toolboxes were identified in the supplied source files. The code uses standard MATLAB matrix operations, plotting, linear algebra, and numerical routines.
- **Operating system:** Windows is supported by the path conventions used in the supplied network-saving functions. For macOS/Linux, use the platform-independent network-data preparation example in Section 4 or change the helper functions to build paths with `fullfile`.
- **Memory and time:** Requirements depend mainly on the number of network nodes, the network spectrum calculation (`eig(full(L))`), and the number of explicit Euler steps (`T/deltaT`). Large networks and small time steps may be computationally demanding.

No installation or compilation is needed beyond adding the repository to MATLAB's search path.

## 2. Repository structure

```text
.
├── README.md
├── Main.m
├── Fucn_Get_Equilibrium.m
├── Fucn_Get_Jacobian.m
├── Func_Check_Turing_Preconditions.m
├── Func_Calculate_Critical_Dm.m
├── Func_Plot_Network_Turing_Dispersion.m
├── Func_ComputePattern_PGG.m
├── Func_generate_and_save.m
├── Func_Get_Initail_Condition.m
└── generateNetwork/
    ├── BAnetwork.m
    ├── ERnetwork.m
    ├── LAnetwork_Customizable.m
    ├── Generating_BA.m
    ├── Generating_BA_Repeat.m
    ├── Generating_ER.m
    ├── Generating_ER_Repeat.m
    ├── Generating_LA.m
    ├── func_SaveSortNetData.m
    ├── func_SaveSortNetData_Repeat.m
    └── LaplacianMatrixSort_<N_nodes>/       # Generated, not supplied
        └── <network_type><degree>.mat     # Generated, not supplied
```

The additional `Main_two_stage_perturbation.m` and `Func_Get_Initail_Condition.m` are the two-stage perturbation example and initialization function discussed alongside the original archive; the original archive contained the other 19 MATLAB files. **Function names are preserved verbatim**, including the existing spellings `Fucn_...` and `Initail`, to avoid breaking calls.

### Main scripts

| File     | Discription                                                                                                                                                                                                                                                                                           | **Reproducibility:**                                                                                                                                                                                                                                                                                                                                                             |
| -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Main.m` | This script simulates the spatiotemporal evolution of cooperation in a networked eco-evolutionary public goods game. It performs stability analysis and two-stage numerical simulations to investigate the emergence and robustness of heterogeneous cooperation patterns under random perturbations. | This script reproduces the main numerical results presented in **Figure 4** of the manuscript. Results for other figures can be obtained by modifying the relevant model parameters, network configurations, diffusion coefficients, and simulation settings in this script. Exact numerical reproduction requires the corresponding network data and random initial conditions. |

### Model and numerical functions

| File                                    | Inputs and outputs                                                              | Description                                                                                                                                                                        |
| --------------------------------------- | ------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Fucn_Get_Equilibrium.m`                | `P` → `[ue, me]`                                                                | Calculates the nontrivial homogeneous equilibrium of the local reaction system.                                                                                                    |
| `Fucn_Get_Jacobian.m`                   | `P, ue, me` → `[Fu, Fm, Gu, Gm]`                                                | Evaluates the four partial derivatives of the two local reaction functions at a specified equilibrium.                                                                             |
| `Func_Check_Turing_Preconditions.m`     | `P` → `is_stable`                                                               | Tests local diffusion-free equilibrium stability: `trace(J) < 0` and `det(J) > 0`; prints diagnostic values.                                                                       |
| `Func_Calculate_Critical_Dm.m`          | `P` → `dm_low`                                                                  | Estimates a public-goods diffusion threshold from the **continuous-spectrum** linear instability condition; does not alone establish an unstable mode on a finite graph.           |
| `Func_Plot_Network_Turing_Dispersion.m` | `P, L` → `[lambda_fine, sigma_max_fine, lambda_abs, sigma_net, lambda_cross]`   | Computes/plots the continuous dispersion curve and growth rates for discrete graph Laplacian eigenvalues; estimates zero crossings within the sampled spectral range.              |
| `Func_ComputePattern_PGG.m`             | `P, N_nodes, L, T, deltaT, u0, m0` → `[U_final, M_final, U_history, M_history]` | Integrates the coupled node-level reaction-diffusion equations using explicit forward Euler. Final-state vectors are `N_nodes × 1`; history matrices are nodes × recorded samples. |
| `Func_generate_and_save.m`              | `N_nodes, filename` → `[deltau, deltam]`                                        | Generates independent standard-normal perturbations and stores them in a `.mat` file for reuse.                                                                                    |
| `Func_Get_Initail_Condition.m`          | `P, N_nodes` → `[u0, m0]`                                                       | Initializes node states near the homogeneous equilibrium with independent Gaussian perturbations of standard deviation `0.0005`; it does not save these random draws.              |

### Network generation files

| File | Description |
| --- | --- |
| `BAnetwork.m` | Generates an undirected preferential-attachment (BA-like) adjacency matrix from `N`, attachment count `m`, and initial seed size `m0`. Also produces a degree-distribution plot. |
| `ERnetwork.m` | Generates an undirected Erdős–Rényi graph `G(N,p)` and returns its realized average degree. |
| `LAnetwork_Customizable.m` | Generates a periodic `M × M` lattice whose neighborhoods are defined by a Manhattan-distance radius `MD`. |
| `Generating_BA.m` / `Generating_BA_Repeat.m` | Batch scripts for BA graph realizations with editable size and attachment parameters. |
| `Generating_ER.m` / `Generating_ER_Repeat.m` | Batch scripts for ER graph realizations with editable graph size and link probability. |
| `Generating_LA.m` | Batch script for periodic-lattice graphs at several Manhattan radii. |
| `func_SaveSortNetData.m` | Sorts nodes by degree and saves `A` and `L` in `.mat` files, plus a text representation of `L`. |
| `func_SaveSortNetData_Repeat.m` | Similar to the preceding function, with a realization/repetition suffix in output filenames. |

## 3. Mathematical model and conventions

Each network node represents a subpopulation with cooperator density `u_i` and public-goods level `m_i`. The total local population density is fixed at `N = P.NA`; variations in cooperator density represent changes in **strategy composition**, not population migration.

The implemented reaction-diffusion dynamics are

$$

\frac{du_i}{dt}=F(u_i,m_i)+d_u\sum_j L_{ij}u_j,
\qquad
\frac{dm_i}{dt}=G(u_i,m_i)+d_m\sum_j L_{ij}m_j,

$$
where
$$

F(u,m)=\frac{u(N-u)}{N^2}\Big[(\eta_1u+N\eta_2)m
-\delta\big(\eta_3u+N(d_0-b_0)\big)\Big],

$$
$$
G(u,m)=\varepsilon\left(\frac{\beta u}{\alpha+u}-kNm\right),
$$

and

\[
\eta_1=a_1-a_0+b_0-b_1+c_0-c_1+d_1-d_0,\quad
\eta_2=d_0-d_1+b_1-b_0,\quad
\eta_3=b_0+c_0-a_0-d_0.
\]

**Important graph convention:** `L = A - D`, where `A` is adjacency and `D` is the diagonal degree matrix. Thus `L` has nonpositive eigenvalues for an undirected graph. The dispersion function plots their nonnegative magnitudes `|lambda|`. Do **not** substitute a positive-semidefinite Laplacian `D - A` without changing the diffusion convention accordingly.

A positive growth rate on the continuous spectral envelope is not sufficient to establish classical Turing instability on a given finite network. The local equilibrium must be stable without diffusion, and at least one nonzero eigenmode of that actual network must have a positive linear growth rate. Diffusion-mediated selection of a heterogeneous steady state from an oscillatory regime is a distinct dynamical mechanism and need not meet the classical Turing precondition.

## 4. Quick start

### 4.1 Set the MATLAB working folder

Open MATLAB, navigate to the **repository root** (the directory containing this README and the main scripts), and execute:

```matlab
projectDir = pwd;
addpath(genpath(projectDir));
```

All main scripts use relative paths. Run them with the repository root as the working directory.

### 4.2 Prepare a matching network `.mat` file

The main scripts expect a file such as

```text
generateNetwork/LaplacianMatrixSort_100/BA10.mat
```

containing at least `L` (`A` is also recommended). **This file is not supplied with the source archive.** The filename's degree label is based on the **rounded, realized mean node degree**, not necessarily the requested or nominal degree parameter.

The following MATLAB commands generate a BA-like network with 100 nodes and save it using portable `fullfile` paths. The actual rounded degree is printed; set the corresponding main script's `average_degree` variable to that value.

```matlab
% Run these commands from the repository root.
projectDir = pwd;
addpath(genpath(projectDir));
rng(42, 'twister');  % Reproduce this network realization.

N_nodes = 100;
A = full(BAnetwork(N_nodes, 5, 6));

% Sort nodes by descending realized degree, as in the supplied savers.
[~, order] = sort(sum(A, 2), 'descend');
A = A(order, order);
L = A - diag(sum(A, 2));
average_degree = round(mean(sum(A, 2)));

outputFolder = fullfile(projectDir, 'generateNetwork', ...
    sprintf('LaplacianMatrixSort_%d', N_nodes));
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
filename = fullfile(outputFolder, sprintf('BA%d.mat', average_degree));
save(filename, 'A', 'L');

fprintf('Saved network: %s\n', filename);
fprintf('Set average_degree = %d in the main script.\n', average_degree);
```

For ER or lattice graphs, replace the adjacency-generation line with one of:

```matlab
[A, ~] = ERnetwork(100, 0.10);            % ER graph, N = 100.
[A, ~] = LAnetwork_Customizable(10, 1);   % 10-by-10 periodic lattice, N = 100.
```

Then perform the same degree sorting and Laplacian saving steps, replacing the filename prefix `BA` by `ER` or `LA`. For the lattice, `N_nodes = M^2`. The batch scripts in `generateNetwork/` offer another option: edit their built-in `N`, `m`, `k`, `M`, or `MD` values and run them with `generateNetwork` as the current folder. The supplied saving functions use Windows-style relative paths, so the portable code above is recommended on macOS and Linux.

### 4.3 Run an example

Ensure the chosen script's three network settings match the saved network:

```matlab
network_type = 'BA';
average_degree = 10;  % Replace by the value printed during generation.
N_nodes = 100;
```

Then, from the repository root, run **one** of:

```matlab
Main_two_stage_perturbation
% OR
Main_oscillationToHeterogeneousSteady
% OR
Main_test_different_network
```

The three main scripts have different default graph sizes and model settings; edit the selected script to point to the network you generated. In particular, `Main_test_different_network.m` initially expects a **400-node BA graph labeled BA74**, whereas the other two examples initially expect a **100-node BA graph labeled BA10**.

### 4.4 Expected results

Depending on the selected main script:

1. MATLAB prints the homogeneous equilibrium, Jacobian entries, local stability diagnostics, and an estimate of the continuous-spectrum diffusion threshold.
2. A plot displays the continuous dispersion relation and sampled graph eigenmode growth rates.
3. The forward Euler solver returns the cooperator/public-goods state histories. The two scripts named `Main_oscillationToHeterogeneousSteady` and `Main_two_stage_perturbation` also plot cooperator-density trajectories.
4. Results remain in the MATLAB workspace. They are **not automatically exported** as figure files or numerical datasets; `Func_generate_and_save` can independently save the random perturbations it creates.

## 5. Changing parameters in the main scripts

### 5.0 Default example configurations

The main scripts have **different parameter presets**. The following values are illustrative starting points; do not mix them when trying to reproduce a particular numerical experiment.

| Main script | Default graph | Selected local payoff parameters | Diffusion protocol | Time integration |
| --- | --- | --- | --- | --- |
| `Main_test_different_network.m` | BA, 400 nodes, nominal file `BA74.mat` | `P.a0 = 10`, `P.c0 = 0.1` | Off in Stage 1; `du = 0.08`, `dm = 3` in Stage 2 | Stage 1: `T=399`, `deltaT=0.005`; Stage 2: `T=600`, `deltaT=0.002` |
| `Main_oscillationToHeterogeneousSteady.m` | BA, 100 nodes, nominal file `BA10.mat` | `P.a0 = 10`, `P.c0 = 1` | Off in Stage 1; `du = 0.08`, `dm = 0.85` in Stage 2 | Stage 1: `T=200`; Stage 2: `T=899`; both `deltaT=0.001` |
| `Main_two_stage_perturbation.m` | BA, 100 nodes, nominal file `BA10.mat` | `P.a0 = 2.9`, `P.c0 = 0.1` | `du = 0.03`, `dm = 0.85` in **both** stages | Stage 1: `T=399`; Stage 2: `T=299`; both `deltaT=0.01` |

For these example presets, the remaining shared values are `a1=3`, `b0=2.3`, `b1=1`, `c1=2`, `d0=1`, `d1=2`, `alpha=0.9`, `beta=2.5`, `delta=2`, `epsilon=1.5`, `k=0.2`, and `NA=3`.

### 5.1 Network configuration

| Setting | Meaning | Effect of changing it |
| --- | --- | --- |
| `network_type` | Topology label: `BA`, `ER`, or `LA`. | Selects a different adjacency/Laplacian matrix. A correspondingly named `.mat` file must exist. |
| `N_nodes` | Number of graph nodes. | Changes the number of coupled subpopulations, the matrix size, and vector sizes. Must match `size(L,1)`. |
| `average_degree` | Rounded realized mean degree used in the filename. | Selects a saved graph with a different level of connectivity; it is **not** a continuous parameter changed internally by the solver. |
| Network random seed | Seed applied *before graph generation*. | Allows the same random graph realization to be regenerated when all generation settings and algorithms are unchanged. |

Changing only `average_degree` in a main script **does not generate a new graph**; it changes the file being loaded. Compare several graphs by preparing each network `.mat` file first and running the main script with matching settings.

### 5.2 Local model parameters (`P`)

| Parameter(s) | Interpretation | Typical use |
| --- | --- | --- |
| `P.a0`, `P.a1`, `P.b0`, `P.b1`, `P.c0`, `P.c1`, `P.d0`, `P.d1` | Payoff coefficients defining the eco-evolutionary game. | Study the effect of game incentives on equilibria, stability, and long-term cooperation. |
| `P.alpha` | Saturation parameter in public-goods production. | Change the nonlinear response of resource production to cooperation. |
| `P.beta` | Public-goods production coefficient. | Modify the production strength in `G(u,m)`. |
| `P.delta` | Public-goods reference/threshold parameter entering the payoff-dependent cooperator reaction term. | Examine the effect of environmental feedback on strategy dynamics. |
| `P.epsilon` | Relative timescale of public-goods dynamics. | Modify how rapidly public goods respond to the cooperator state. |
| `P.k` | Public-goods consumption/decay coefficient in `G(u,m)`. | Modify public-goods removal strength. |
| `P.NA` | Fixed total population density **per node**, denoted `N` in the equations. | Change the local population normalization; **not** the number of graph nodes. |
| `P.du` | Diffusion coefficient for the cooperator strategy density. | Adjust network-mediated strategy diffusion. Set to zero to disable this coupling. |
| `P.dm` | Diffusion coefficient for public goods. | Adjust public-goods diffusion. Set to zero to disable this coupling. |

**After changing the payoff or environmental parameters**, recompute the equilibrium, the Jacobian, the stability diagnostics, and the dispersion relation. The chosen root of the equilibrium equation may become non-real, outside the admissible density interval, or numerically singular for some parameter sets; this implementation does not automatically select an alternative physical root.

### 5.3 Numerical simulation controls

| Setting | Meaning | Practical guidance |
| --- | --- | --- |
| `T` or `T_total` | Duration of each simulation stage. | Increase to study long transients; each stage starts its own time axis at zero. |
| `deltaT` | Explicit Euler integration time step. | Smaller steps generally improve numerical resolution but increase cost. Check sensitivity to `deltaT`; overly large steps can cause numerical instability. |
| `PS` | Amplitude multiplier for stored standard-normal initial perturbations in the two original main scripts. | Change initial-condition amplitude without regenerating the saved random vector. |
| `perturbation_scale` in `Func_Get_Initail_Condition.m` | Standard deviation of initial Gaussian perturbations (default `0.0005`). | Controls the initial displacement from the homogeneous equilibrium in the additional two-stage example. |
| `0.1 * rand(N_nodes, 1)` | Independent **nonnegative uniform** perturbations applied at the start of Stage 2 in `Main_two_stage_perturbation.m`. | Change `0.1` to control additional perturbation magnitude; this is not a zero-mean perturbation. |

In `Func_ComputePattern_PGG`, time integration occurs at every `deltaT` but storage is intended at integer model-time points. The implementation tests the exact condition `mod(timeGrid(n), 1) == 0`, which can be sensitive to floating-point time-grid representation. Use a step size dividing one time unit (as in the examples), and verify the number and contents of recorded snapshots before quantitative analysis. The solver does not enforce bounds such as `0 <= u_i <= P.NA` or positivity of `m_i`.

### 5.4 Three distinct two-stage experiments

**A. Switch diffusion on after an initial reaction-only stage**

Use `Main_oscillationToHeterogeneousSteady.m`:

```matlab
% Stage 1: P.du = 0;    P.dm = 0;
% Stage 2: P.du = 0.08; P.dm = 0.85;
T_total = [200, 899];
deltaT = 0.001;
PS = 0.000001;
```

The final state of Stage 1 directly initializes Stage 2. Its existing stability/dispersion diagnostics are called **before** the diffusion switch, while `P.du` and `P.dm` are zero; to analyze the post-switch dispersion relation, invoke the dispersion function again *after* setting the nonzero diffusion coefficients.

**B. Compare a selected network under reaction-only and diffusion-coupled dynamics**

Use `Main_test_different_network.m`:

```matlab
% Default graph: BA, N_nodes = 400, average_degree = 74.
% Dispersion settings: P.du = 0.08; P.dm = 3;
% Stage 1: P.du = 0;    P.dm = 0; T = 399; deltaT = 0.005;
% Stage 2: P.du = 0.08; P.dm = 3; T = 600; deltaT = 0.002;
```

To compare BA/ER/LA graphs, generate matching network files and change `network_type`, `average_degree`, and `N_nodes`. This script indexes `lambda_cross(1)` and `lambda_cross(2)` without checking that two crossings exist. For parameter combinations with fewer than two detected crossings, skip or guard that indexing before running the full script.

**C. Apply additional perturbations while diffusion stays active**

Use `Main_two_stage_perturbation.m`:

```matlab
P.du = 0.03;
P.dm = 0.85;

% Stage 1: T = 399; deltaT = 0.01;
% Stage 2: T = 299; deltaT = 0.01;
% New initial state for Stage 2:
u0 = U_final + 0.1 * rand(N_nodes, 1);
m0 = M_final + 0.1 * rand(N_nodes, 1);
```

Both stages keep the **same** diffusion coefficients. Do not interpret this experiment as switching diffusion on at the stage boundary.

## 6. Reproducibility and data handling

1. **Fix random seeds.** Add `rng(42, 'twister');` near the beginning of a main script, immediately after `clear; clc;` (or before any random sampling), to reproduce its random initialization, provided the same network data are loaded. Set a separate seed before generating random networks.
2. **Save the networks.** Archive each generated `A`/`L` pair along with its topology, node count, realized mean degree, network-generation settings, and random seed. This is important for exact graph-by-graph comparisons.
3. **Reuse initial perturbations when appropriate.** The two original main scripts look for `initial_perturbation_<N_nodes>.mat` and reuse `deltau`/`deltam` if it already exists. Deleting it causes new perturbations to be generated. The additional two-stage example instead calls `Func_Get_Initail_Condition` and does not save its draws by default.
4. **Record the simulation configuration.** At minimum retain `P`, `N_nodes`, `network_type`, `average_degree`, `T`, `deltaT`, `PS` (where applicable), perturbation amplitude, and the two stage-specific diffusion settings.
5. **Save outputs explicitly.** For example, after a run:

   ```matlab
   save('simulation_results.mat', 'P', 'L', ...
       'data1', 'data2', 'U_final2', 'M_final2');
   ```

   The above example applies to main scripts defining those named variables; adapt the list for `Main_test_different_network.m`, which instead defines `data`.

**Recorded-time interpretation:** The history matrices have node indices as rows and intended integer-time snapshots as columns. When two stages are concatenated, both stages' initial snapshots are included, so the join may contain a repeated physical time (and, in the additional-perturbation example, an instantaneous state jump). Plots using `1:size(data1,1)` show **recorded-sample indices**, not Euler step counts or an automatically reconstructed physical-time axis.

## 7. Common problems

| Symptom | Likely cause / action |
| --- | --- |
| `Unable to find file ... LaplacianMatrixSort_.../...mat` | Generate the network and confirm the filename matches the topology, node count, and **realized** rounded mean degree. Run the main script from the repository root. |
| `Unrecognized function or variable 'L'` | The loaded `.mat` file lacks the signed Laplacian. Save both `A` and `L`, with `L = A - diag(sum(A,2))`. |
| Matrix dimension mismatch | Check `size(L) == [N_nodes, N_nodes]` and column-vector initial conditions of length `N_nodes`. |
| `Index exceeds the number of array elements` at `lambda_cross(...)` | The chosen configuration does not yield two detected dispersion zero crossings; check `numel(lambda_cross)` before indexing. |
| Non-real, non-finite, or implausible equilibrium | The specified payoff/environmental parameters may not admit the selected physical equilibrium; inspect `gamma1`, the quadratic discriminant, and `0 <= ue <= P.NA`. |
| Solver diverges or produces negative densities | Reduce `deltaT`, check the model/parameter regime and initial perturbations, and examine physical admissibility. The explicit Euler solver does not constrain state variables. |
| Missing or zero-filled recorded history columns | The `mod(timeGrid(n), 1) == 0` recording condition may not recognize every expected integer time due to floating-point representation; use integer step counters for a more robust recording implementation. |
| Network-saving paths differ across operating systems | The supplied `func_SaveSortNetData*` functions use Windows-style path separators. Use the `fullfile`-based saving example in Section 4. |

## 8. Research-code notes

- The **continuous instability threshold** returned by `Func_Calculate_Critical_Dm` and the **finite-graph discrete-spectrum test** in `Func_Plot_Network_Turing_Dispersion` answer different questions. Check both the diffusion-free stability criterion and actual network eigenmodes before labeling a result a classical Turing instability.
- The final-state outputs of `Func_ComputePattern_PGG` are called `U_stable` and `M_stable` internally, but numerical integration ending at time `T` does **not** by itself prove asymptotic stability.
- **Known graph-helper limitations in the supplied source:** `func_SaveSortNetData_Repeat.m` reconstructs its `A` output with `A = L - diag(L)` rather than `A = L - diag(diag(L))`; review/correct this line before treating the saved `A` as an adjacency matrix. `LAnetwork_Customizable.m` uses `unique(rpn, 'row')`; if MATLAB rejects this option, use the documented `'rows'` option instead. These two lines have not been silently changed in this package.
- No automated unit/regression tests, reference network matrices, or checksum-verified benchmark trajectories are supplied. For archival reproducibility, publishing the graph realizations and seeds used for manuscript figures is recommended.
- The included scripts do not automatically perform systematic parameter sweeps, statistical ensemble aggregation, or figure-by-figure reproduction for the full manuscript. Those analyses require additional driver scripts and/or recorded configurations when applicable.

## 9. Citation and reuse

If using this implementation in research, cite the associated article after publication and link this repository to the article's **Data Accessibility / Code Availability** statement. Add the final bibliographic details and DOI when available. Licensing terms have not been specified in this code package; add an explicit `LICENSE` file before publishing if public reuse permissions are intended.
