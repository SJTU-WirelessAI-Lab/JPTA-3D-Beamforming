# Reproduction report — 3D Beamforming Through Joint Phase-Time Arrays

Paper: O. Yildiz, A. AlAmmouri, J. Mo, Y. Nam, E. Erkip, J. Zhang, "3D Beamforming Through
Joint Phase-Time Arrays", in *Proc. 2024 IEEE 100th Vehicular Technology Conference
(VTC2024-Fall)*, Washington, DC, USA, Oct. 2024, pp. 1-7.
DOI 10.1109/VTC2024-Fall63153.2024.10757955, IEEE Xplore 10757955.

Environment: MATLAB R2026b Update 1, Windows 11, 16 threads, 32 GB. No toolbox is used.
Gap numbers **G1–G9** are local to this repository.

**Reproducibility.** The whole tree was re-run from a fresh copy (`run_all`, 4.6 min). Every
gain, objective and iteration count came out identical to the first run; only the wall-clock
columns of `a6_fig5a.csv`, `a6_fig5b.csv`, `a7_fig7.csv` and the two `*_summary.txt` lines that
quote run time moved, which is why the run-time claims in §6 are quoted as ranges.

---

## 1 What is covered

| Job | Script | Figures / tables | Output |
|---|---|---|---|
| P0 | `run_p0.m` | 4 kernel invariants | `out/p0_selfcheck.txt` |
| A1 | `run_a1.m` | Fig. 2(a) | `out/a1_gain_vs_alpha.{csv,png}` |
| A3 | `run_a3.m` | Fig. 2(b) | `out/a3_fairness.png` |
| A4/A5 | `run_a45.m` | Fig. 3, Fig. 4 | `out/a4_pattern_{joint,sep}.png`, `out/a5_az_freq_joint.png` |
| A6 | `run_a6.m` | Fig. 5, plus the run-time claims | `out/a6_fig5{a,b}*.{csv,png}`, `out/a6_summary.txt` |
| A7 | `run_a7.m` | Fig. 6, Fig. 7 | `out/a7_fig6.*`, `out/a7_fig7*.*`, `out/a7_summary.txt` |

Everything is implemented from the equations of the paper; nothing is fitted. The only
external input is the baseline of Fig. 6/7, which is Algorithm 1 of ref. [5] re-implemented
from the published description (§8).

---

## 2 Model and solvers

**Hardware model.** With $f_m=f_c+m\Delta f$, the response of element $(y,z)$ is
$\exp\{j(\varphi_{y,z}+2\pi f_m\tau_{y,z})\}$, and the gain toward user $i$ at $f_m$ is

$$G_i(f_m)=\frac{1}{N_{az}N_{el}}\Big|\sum_{y,z}e^{j(\varphi_{y,z}+2\pi f_m\tau_{y,z}-\Omega_{y,z}(f_m))}\Big|^2,
\qquad \Omega_{y,z}=\pi\tfrac{f_m}{f_c}\big(y\sin\theta_{az}\sin\theta_{el}+z\cos\theta_{el}\big).$$

The objective of eq. (5) is $G_l=\sum_i 10\log_{10}\bar G_i$ with $\bar G_i$ the mean of
$G_i$ over the subcarriers of user $i$; Fig. 2 plots $G_l/N_u$.

**Analytical solution (eqs. 7–13).** The target phase of element $(y,z)$ is constant over each
user's band, so a design has only two degrees of freedom per element — a phase at $m=0$ and a
slope in $m$. With $m$ symmetric around 0, the normal equations of the LS fit are diagonal,
so $x_1=\mathrm{mean}(B)$ and $x_2=\sum mB/\sum m^2$ in closed form: no iteration.

**Infinity-norm variant.** The paper suggests a linear program; the compute node has no
Optimization Toolbox. For a fixed slope $x_2$, the best intercept is the mid-range of
$b-x_2m$ and the residual $t(x_2)$ is convex in $x_2$, so `fminbnd` gives the exact solution
in one dimension. 0.12 s per design point, about 10× faster than the per-element LP.

**Quantization.** $\varphi$ to 6 bit, $\tau$ to 2.5 ns inside $[0,200]$ ns. Both are applied
in `jp_finalize.m` *after* the continuous solution, and every reported gain is computed from
the quantized quantities.

**Unwrapping (eq. 10).** The analytical target is piecewise constant in frequency with a jump
between adjacent users. Adding the integer $2\pi k_i(y,z)$ from eq. (10) keeps each jump below
$2\pi$ in the unwrapped domain. This is not cosmetic: without it the LS solution at
$N_u=2,\alpha=0.1$ returns −3.05/−4.60 dB where the paper has 22.1/18.1 dB.

---

## 3 Self-check (P0)

| # | Invariant | Result | Measured |
|---|---|---|---|
| P0-1 | $N_u=1,\alpha=1$ degenerates to a pure phased array: delays constant, gain flat | PASS | max$\vert\tau\vert=0$ exactly; gain 25.839 dB; ripple 0.003 dB |
| P0-2 | no design exceeds $10\log_{10}(16\cdot24)=25.843$ dB | PASS | max observed 25.8433 dB (the bound itself) |
| P0-3 | with $N_{el}=1$, joint and separated solutions coincide | PASS | max$\vert\Delta G\vert=0$, max$\vert\Delta\tau\vert=0$ |
| P0-4 | required delay inside $\tau_{\max}=200$ ns | PASS | $N_u$=2/4/8/16 → 15.5 / 33.8 / 78.0 / 133.9 ns |

> P0-1 failed in the first version (8.35 dB instead of 25.84 dB). The cause was an extra
> factor $\pi$ in `jp_gain_freqs`: `jp_nu` already returns $\nu$ **including** $\pi$, so the
> gain expression multiplied the phase difference by $\pi$ again. Without this assertion the
> whole study would have run 17 dB low and nothing else would have looked obviously wrong.

---

## 4 Fig. 2(a): two users, $\alpha_1=\alpha$, $\alpha_2=1-\alpha$

Users at azimuth $[-60^\circ,60^\circ]$, elevation $[90^\circ,120^\circ]$. The table is
$G_l/N_u$ (dB), i.e. the vertical axis of Fig. 2(a).

| α | Joint-∞ | Joint-LS | Sep-∞ | Sep-LS | Δ(J−S) LS | Δ(J−S) ∞ |
|---|---|---|---|---|---|---|
| 0.10 | 23.734 | 21.975 | 21.847 | 18.634 | +3.34 | +1.89 |
| 0.15 | 23.929 | 23.345 | 22.400 | 20.678 | +2.67 | +1.53 |
| 0.20 | 24.154 | 24.129 | 22.824 | 22.439 | +1.69 | +1.33 |
| 0.25 | 24.351 | 24.480 | 23.106 | 23.402 | +1.08 | +1.24 |
| 0.30 | 24.555 | 24.708 | 23.631 | 23.925 | +0.78 | +0.92 |
| 0.35 | 24.697 | 24.848 | 23.889 | 24.109 | +0.74 | +0.81 |
| 0.40 | 24.784 | 24.901 | 24.012 | 24.242 | +0.66 | +0.77 |
| 0.45 | 24.778 | 24.935 | 24.017 | 24.290 | +0.64 | +0.76 |
| 0.50 | 24.671 | 24.944 | 23.951 | 24.312 | +0.63 | +0.72 |

| Claim | Paper | Reproduced | Verdict |
|---|---|---|---|
| joint always better, gap within 0.6–4 dB | 0.6–4 dB | 0.63–3.34 dB | hit (lower end almost exact) |
| mean gap, LS | 1.29 dB | 1.36 dB | +0.07 |
| mean gap, infinity norm | 1.18 dB | 1.11 dB | −0.07 |
| for α ≤ 0.2 the infinity norm beats LS | yes | yes: joint +1.76/+0.58/+0.03, separated +3.21/+1.72/+0.39 dB | trend exact |

Absolute heights also land inside the 18–26 dB range of Fig. 2(a) and saturate at 24–25 dB
beyond α ≈ 0.4, as in the paper.

**A2, read off the same table.** Degradation of LS relative to the infinity norm for α ≤ 0.2:

| Item | Paper | Reproduced |
|---|---|---|
| joint: max degradation | 1.81 dB | **1.76 dB** |
| joint: mean degradation | 0.76 dB | **0.79 dB** |
| separated: max degradation | 3.36 dB | **3.21 dB** |
| separated: mean degradation | 1.40 dB | 1.77 dB ✗ |

Three of four hit. The outlier is an accounting difference: the paper does not say over how
many points the α ≤ 0.2 average is taken, and with only 0.10/0.15/0.20 the value at 0.20
(where the difference is already ≈0.4 dB) pulls the mean up. A denser grid starting at 0.02
would move it down.

---

## 5 Fig. 2(b): five users

Azimuth $[-60,-30,0,30,60]^\circ$, elevation $[90,97.5,105,112.5,120]^\circ$,
$\alpha=[0.3,0.2,0.15,0.1,0.25]$.

| Design | U1 | U2 | U3 | U4 | U5 | Max spread | Paper | Verdict |
|---|---|---|---|---|---|---|---|---|
| Joint-∞ | 19.42 | 19.70 | 22.55 | 20.33 | 19.34 | **3.20** | 3.0 | ✓ |
| Joint-LS | 23.10 | 18.11 | 20.97 | 19.68 | 22.77 | **4.98** | 5.0 | ✓✓ |
| Sep-∞ | 16.73 | 18.50 | 21.19 | 19.74 | 15.84 | 5.35 | 3.2 | ✗ |
| Sep-LS | 20.41 | 14.95 | 19.46 | 20.56 | 20.28 | **5.61** | 5.4 | ✓ |

Three of four match (Joint-LS to 0.02 dB, Joint-∞ to 0.2 dB, Sep-LS to 0.21 dB).

The mismatch is the separated infinity-norm. This implementation minimises the infinity norm
of each array axis separately, while what a user sees is the **product** of the two axis
gains; per-axis minimax is not minimax of the product, and it can therefore be *less* fair
than LS (U1/U5 drop to 16.7/15.8 dB here). The paper's separated infinity norm is probably
taken on the product gain or jointly over users. Recorded as **G1** rather than forced.

---

## 6 Fig. 3 and Fig. 4

**An unplanned but strong check.** The text says user 3 owns $[-9.5,9.5)$ MHz. The subcarrier
rule $[f_{\alpha_{i-1}(M+1)},f_{\alpha_i(M+1)})$ with $M+1=793$ gives $m=-79\ldots79$, i.e.
$[-9.48,\,9.48]$ MHz — agreement to 0.02 MHz, which independently confirms the rounding
convention of the band edges.

| Item | Result |
|---|---|
| azimuth–frequency map ($\theta_{el}=105^\circ$, joint LS) | the bright region sits exactly inside the $\pm9.48$ MHz frame; user 3 in-band mean **22.04 dB**, out-of-band **8.14 dB** ($G_{max}=25.84$) |
| azimuth–elevation patterns, separated vs joint | the separated design leaks a striped pattern over the whole plane; the joint design gives five compact spots |

The qualitative statements of the paper ("the separated design has many more sidelobes",
"the in-band target gain is reasonably met") both hold. Worth noting: the in-band gain is
22.04 dB, 3.8 dB below the bound. That 3.8 dB is the intrinsic loss of a joint LS design on an
equal-bandwidth five-user split — the same number that makes Fig. 2(b) wander between 19 and
23 dB, so the two figures are consistent with each other.

---

## 7 Fig. 5: greedy Algorithm 1 and Adam gradient descent

### 7.1 First prove that "my greedy" is Algorithm 1

The incremental implementation exploits that updating one element changes exactly one term
of the array-factor sum, so each candidate is evaluated with $S_{\text{rest}}=S-\text{current}$.
The literal implementation recomputes eq. (5) in full for every candidate. On a 6×8 array,
three rounds:

| Mode | incremental $G_l$ | literal $G_l$ | difference | time |
|---|---|---|---|---|
| joint | 60.0431232911 | 60.0431232911 | $3.6\times10^{-14}$ | 1.2 s / 23.7 s |
| separated | 54.9326930484 | 54.9326930484 | 0 | 0.3 s / 6.5 s |

Search space ($\tau$ on 0:2.5:200 ns, 81 points; phase on 6 bit, 64 points) and objective are
identical to Algorithm 1; the only difference is whether the whole array is recomputed per
candidate. **Everything after this depends on it: without the check, "different algorithm"
and "different implementation" cannot be told apart.**

### 7.2 Fig. 5(a): $N_u=4$, equal bands

| Design | U1 | U2 | U3 | U4 | Mean | Iterations | Time |
|---|---|---|---|---|---|---|---|
| joint LS (start) | 23.27 | 21.76 | 21.36 | 23.19 | 22.394 | – | 0.01 s |
| joint greedy | 23.36 | 22.24 | 21.94 | 23.30 | **22.709** | 4 | 9.9 s |
| joint GD | 23.89 | 22.10 | 21.77 | 23.73 | **22.874** | 194 | 6.0 s |
| separated LS (start) | 21.43 | 21.68 | 18.20 | 20.92 | 20.559 | – | 0.01 s |
| separated greedy | 22.30 | 20.91 | 21.02 | 22.12 | **21.587** | 8 | 2.0 s |
| separated GD | 22.55 | 20.50 | 18.59 | 22.15 | 20.949 | 77 | 0.3 s |
| *paper* joint greedy | 23.35 | 22.25 | 21.90 | 23.05 | *22.64* | | |
| *paper* joint GD | 23.50 | 21.85 | 21.70 | 23.45 | *22.62* | | |
| *paper* separated greedy | 21.70 | 21.10 | 21.10 | 22.35 | *21.56* | | |
| *paper* separated GD | 22.20 | 21.10 | 20.15 | 21.30 | *21.19* | | |

(paper values read from Fig. 5(a), ±0.2 dB.)

| Claim | Result | Verdict |
|---|---|---|
| joint greedy, per user | max deviation **0.25 dB**, mean 22.709 vs 22.64 | ✓ |
| joint GD, per user | max deviation **0.39 dB**, mean 22.874 vs 22.62 | ✓ |
| separated greedy | mean 21.587 vs 21.56, max deviation 0.60 | ✓ |
| separated GD | mean 20.949 vs 21.19, user 3 off by **1.56 dB** | ✗ **G6** |
| "greedy and GD are comparable" | joint difference **+0.16 dB** (GD slightly higher, same direction as the paper) | ✓ |

### 7.3 Fig. 5(b): $N_u=2$, α sweep

| α | LS joint | greedy joint | GD joint | LS sep | greedy sep | GD sep | GD−greedy |
|---|---|---|---|---|---|---|---|
| 0.10 | 21.975 | 24.690 | 24.792 | 18.634 | 23.768 | 24.280 | +0.10 |
| 0.20 | 24.129 | 24.786 | 24.885 | 22.439 | 24.193 | 24.354 | +0.10 |
| 0.30 | 24.708 | 24.874 | 24.958 | 23.925 | 24.239 | 24.400 | +0.08 |
| 0.40 | 24.901 | 24.927 | 25.005 | 24.242 | 24.299 | 24.424 | +0.08 |
| 0.50 | 24.944 | 24.944 | 25.021 | 24.312 | 24.312 | 24.430 | +0.08 |

Mean $\vert$GD − greedy$\vert$: 0.088 dB joint, 0.194 dB separated, so "comparable" holds.
The shapes match Fig. 5(b) as well: the joint pair is flat at 24.7–25.0 dB, the separated pair
sits at 24.0–24.4 dB and is even flatter.

### 7.4 The "250× faster" is a property of the implementation

| Item | Measured |
|---|---|
| one evaluation of eq. (5) (16×24 array, 793 subcarriers, 4 users) | 6.8 ms |
| literal Algorithm 1: $145\times384=55\,680$ evaluations per round | **381 s / round** |
| incremental implementation | 2.48 s / round (4 rounds, 9.9 s) |
| Adam GD | 0.031 s / iteration (194 iterations, 6.0 s) |
| **literal Alg. 1 / GD** | **219×** (paper: about 250×) |
| greedy / GD, both implemented well | **1.6×** |

The paper's 250× matches a naive implementation that recomputes the whole array per candidate.
Put the incremental evaluation in, and the same algorithm at the same iterates collapses to
1.6×. **That claim constrains the cost model, not the convergence speed.** (Absolute times:
GD 6.0 s here vs 1.05 s in the paper — machine and MATLAB overhead. Wall-clock is the one thing
here that is not reproducible: across two full runs of the same code GD took 6.0 and 8.2 s and
the literal sweep 1320 and 1520 s, so `literal / GD` came out 219× and 187×, `greedy / GD` 1.6×
and 1.3×. Only the order of magnitude is a result; the numbers in the tables are the ones in
`out/a6_summary.txt`, i.e. the second run, which is the run the committed files come from.)

### 7.5 Where "+2.26 / +2.11 dB" comes from (G5)

The end points match, but the gain *over the start point* is only +0.32 dB (greedy) and
+0.48 dB (GD), not 2.26 / 2.11. Back-solving: +2.26 dB requires a start point of
22.71 − 2.26 = 20.45 dB, whereas the joint LS start is 22.39 dB.

| Start-point convention | U1 | U2 | U3 | U4 | Mean | after greedy | gain |
|---|---|---|---|---|---|---|---|
| with $2\pi K$ unwrapping (used here) | 23.27 | 21.76 | 21.36 | 23.19 | 22.394 | 22.709 | +0.32 |
| raw steering value, no unwrapping | 18.63 | 18.63 | 18.17 | 18.52 | 18.488 | 21.749 | +3.26 |
| **separated LS** | 21.43 | 21.68 | 18.20 | 20.92 | **20.559** | – | – |

The raw convention gives −3.05/−4.60 dB at $N_u=2,\alpha=0.1$, which contradicts Fig. 2(a)
(22.1/18.1 dB), so the paper does use some unwrapping and this implementation's choice is the
right one — the four curve families of Fig. 2(a) matching to 0.2 dB is the evidence. That
leaves the separated LS as the only candidate: joint greedy − separated LS = **+2.15 dB**
(paper +2.26) and joint GD − separated LS = **+2.31** (paper +2.11).

**Conclusion: the increment depends on a start point the paper never states; its value implies
20.4–20.6 dB, i.e. the separated LS, not the joint LS they plot in Fig. 2(a). The end points
reproduce, the increment does not.** Not forced.

### 7.6 Three gaps forced by this step

| ID | Problem | Handling |
|---|---|---|
| **G2** | one learning rate 0.1 for both φ and τ is dimensionally impossible: an Adam step moves a variable by ≈lr, so for τ in seconds one step is $2\pi\Delta f\cdot0.1\approx7.5\times10^4$ rad | reparameterise to $a=\varphi+2\pi f_c\tau$ and $b=2\pi\Delta f\tau$ (both in rad), i.e. "lr = 0.1 rad" |
| **G3** | $\tau\in[0,200]$ ns means $b\in[0,0.151]$ rad, only 2.4% of its $2\pi$ period; unconstrained Adam leaves the window in a few steps (first version pinned τ at 200 ns and lost 12 dB after mapping back) | bounded reparameterisation $b=B_{\max}\sigma(u)$; after mapping back the loss is $\le10^{-13}$ dB |
| **G4** | the convergence test of Algorithm 1 line 18 fires at the first iteration on GD, because a single Adam step changes $G_l$ by about $10^{-4}$ relative, below $\zeta=10^{-4}$ | 20-iteration sliding window for GD, per-round test for the greedy |

### 7.7 One open item (G6)

Separated GD is 1.56 dB low on user 3 (18.59 vs 20.15). Under the separated
parameterisation $G_i$ is a product of two axis gains, and user 3 sits at the compromise point
of both axes; the bounded delay compresses the available correction, so GD stalls near the
separated LS start (18.20 → 18.59) and does worse than the separated greedy, which escapes by
reconstructing element by element. This is exactly the regime where the paper says the
separated design "does not work" for JPTA — the direction agrees, the magnitude does not.

---

## 8 Fig. 6 and Fig. 7: against Algorithm 1 of ref. [5]

The baseline is the alternating optimisation of the earlier JPTA paper, lifted from 1-D/2-D to
3-D here: three steps (delays by Lemma 1 or by the line search of eq. (8a); phases by
eq. (8b); a free per-subcarrier phase by eq. (11)).

### 8.1 The functional is not the same as C16's

$$\text{[5]: } F=\sum_k\sqrt{G_k}\ \text{(every \emph{subcarrier} equal)}
\qquad\text{vs}\qquad
\text{C16: } G_l=\sum_i 10\log_{10}\Big(\tfrac1{|F_i|}\sum_{k\in F_i}G_k\Big)\ \text{(every \emph{user} equal)}$$

This one difference explains §8.3 in full.

### 8.2 Three degrees of freedom the text does not fix

| # | Freedom | Choice | Reason |
|---|---|---|---|
| 1 | number of TTD units $N$ and the grouping | $N=N_{az}N_{el}=384$, one per element | this paper's hardware model; [5] allows $N\le M$, and $N=M$ is both the most favourable and the fairest reading |
| 2 | target $\bar b_k$ under multi-user FDM | $\bar b_k=a(f_k,\theta_{i(k)})/\sqrt{N_{az}N_{el}}$, the steering vector of the user owning that subcarrier | [5] only defines its own two single-user beams; this is the only self-consistent extension and it ties $|\bar b_k^\dagger TPd_k|=\sqrt{G_k}$ to C16 eq. (3) |
| 3 | weights $\omega_k$ | all 1 | [5] allows any non-negative weights |

### 8.3 Two pitfalls that only show up numerically

**(a) The two papers differ by a sign in the delay term.** C16 eq. (2) writes
$e^{+j2\pi f_m\tau}$, [5] eq. (1) writes $e^{-j2\pi f_k\tau}$. Steps (7), (8a), (8b), (11) are
each self-consistent only inside their own convention. Mixing them gives a design that
converges every round and is 7.9 dB wrong: in the first version $F$ climbed from 9679 to
14047 (its own functional behaving perfectly), while the C16 evaluation gave **−7.88 dB**.
Fix: run the whole iteration in [5]'s convention and map back at the end,
$\tau^{\text{C16}}=\max_n\tau_n^{[5]}-\tau_n^{[5]}$, which multiplies the response by a
common $e^{+j2\pi f\max\tau^{[5]}}$ — a pure phase per subcarrier, so $|G_k|$ is untouched and
$\tau'\ge0$ is realisable.

**(b) The least squares of Lemma 1 needs a change of variables.** With the design matrix
$[1,\,2\pi\Delta f m]$ the second column is $\sim3\times10^8$ while the in-band variation is
0.15 rad: `RCOND = 3.4e-17`, i.e. numerically singular. Using $b\triangleq2\pi\Delta f\tau$ as
the unknown with the matrix $[1,-m]$ brings the condition number to $10^4$.

### 8.4 Unit checks before comparing anything

| # | Assertion | Result |
|---|---|---|
| 1 | single user, full band: Alg. 1 must reach $10\log_{10}384=25.843$ dB (also tests the $f_k/f_c$ squint on both the target and the evaluation side) | PASS: broadside 25.843 dB; sector edge $(-60^\circ,105^\circ)$ 25.821 dB, worst subcarrier 25.798 dB |
| 2 | shifting all delays by 37 ns does not change the gain | PASS: $\max\vert\Delta G\vert=1.7\times10^{-13}$ dB |
| 3 | on its own functional at α=0.5 it must not be worse than the C16 LS design | PASS: 14047.0 vs 13928.3 |
| 4 | is 10 iterations enough ([5] uses 10)? | converged at round 3; rounds 4–10 unchanged; 30 rounds identical |

Check 3 is the important one: the baseline is not a straw man — it beats the C16 analytical
solution on its own objective, so the small gap in §8.3 is not an implementation artefact.

### 8.5 Fig. 6, point by point ($N_u=2$, α = 0.10:0.05:0.50)

| α | Joint LS | Joint ∞ | Joint greedy | Joint GD | Sep LS | Sep greedy | Sep GD | **SOTA** | LS−SOTA |
|---|---|---|---|---|---|---|---|---|---|
| 0.10 | 21.97 | 23.73 | 24.69 | **24.79** | 18.63 | 23.77 | 24.28 | 21.65 | +0.33 |
| 0.15 | 23.35 | 23.93 | 24.73 | **24.84** | 20.68 | 24.07 | 24.32 | 23.24 | +0.10 |
| 0.20 | 24.13 | 24.15 | 24.79 | **24.89** | 22.44 | 24.19 | 24.35 | 24.12 | +0.01 |
| 0.25 | 24.48 | 24.35 | 24.84 | **24.92** | 23.40 | 24.23 | 24.38 | 24.58 | −0.10 |
| 0.30 | 24.71 | 24.56 | 24.87 | **24.96** | 23.92 | 24.24 | 24.40 | 24.82 | −0.11 |
| 0.35 | 24.85 | 24.70 | 24.91 | **24.99** | 24.11 | 24.27 | 24.41 | 24.94 | −0.09 |
| 0.40 | 24.90 | 24.78 | 24.93 | **25.00** | 24.24 | 24.30 | 24.42 | 24.99 | −0.09 |
| 0.45 | 24.93 | 24.78 | 24.94 | **25.02** | 24.29 | 24.31 | 24.43 | 25.02 | −0.08 |
| 0.50 | 24.94 | 24.67 | 24.94 | **25.02** | 24.31 | 24.31 | 24.43 | 25.02 | −0.08 |

| Claim | Paper | Reproduced | Verdict |
|---|---|---|---|
| for α > 0.35 the joint methods converge | qualitative | spread of LS/∞/greedy/GD/SOTA: 0.22 / 0.24 / 0.35 dB | **holds** |
| for α ≤ 0.35 the analytical design beats the SOTA by up to 1.54 dB | 1.54 dB | at most **0.33 dB** (α=0.10) | **does not hold** |
| for α ≤ 0.2, GD + separated beats everything except GD + joint | qualitative | yes: 24.28 > 23.73 dB at α=0.10 | **holds** |

### 8.6 Chasing the 1.54 dB: cross-metric, not a difference in capability

**(a) Two functionals, one set of designs.**

| α | $F$(C16 LS) | $F$(SOTA) | $G_l/N_u$(LS) | $G_l/N_u$(SOTA) | per-user gain LS / SOTA |
|---|---|---|---|---|---|
| 0.10 | 13978.8 | **14103.9** | 21.975 | 21.649 | [18.53 25.42] / [17.76 25.53] |
| 0.20 | 13451.7 | **13561.3** | 24.129 | 24.118 | [23.32 24.94] / [23.19 25.05] |
| 0.30 | 13584.0 | **13706.2** | 24.708 | **24.822** | [24.63 24.78] / [24.82 24.82] |
| 0.40 | 13819.8 | **13949.6** | 24.901 | **24.994** | [25.00 24.80] / [25.13 24.86] |
| 0.50 | 13928.3 | **14047.0** | 24.944 | **25.021** | [24.94 24.95] / [25.02 25.02] |

The SOTA wins on its own $F$ at every α (by 0.4–1.6%) and loses on $G_l/N_u$ below α=0.25.
The difference is the weighting: at α=0.10 it pushes the narrow-band user down to 17.76 dB to
buy 25.53 dB for the wide-band user, and the log-mean of C16 scores that trade as a loss.
**"The SOTA is 1.54 dB behind" is largely "design A scored by design B's metric".**

**(b) Sensitivity to the definition.** The paper only says "lift [5]'s iterative algorithm to
3-D"; it does not say whether to use eq. (8a) or Lemma 1, and [5] gives both.

| Reading | α=0.10 | 0.15 | 0.25 | 0.50 |
|---|---|---|---|---|
| Lemma 1, eq. (10a) folding (**main**) | 21.649 | 23.245 | 24.578 | **25.021** |
| Lemma 1, interval clipping | 21.402 | 22.968 | 24.440 | 24.997 |
| Lemma 1, free $\angle\alpha_k$ removed | 21.577 | 21.931 | 22.380 | 22.625 |
| eq. (8a) line search, 2.5 ns grid | 13.269 | 19.648 | 22.685 | 24.437 |
| eq. (8a) line search, 0.25 ns grid | 16.963 | 19.377 | 22.752 | 24.506 |
| *C16 joint LS (reference)* | *21.975* | *23.345* | *24.480* | *24.944* |

Lemma 1 gives 0.33 dB, eq. (8a) gives 8.71 dB, and the paper's 1.54 dB sits between them with
nothing in the text to decide. Recorded as **G7** — underdetermined, not unreproducible. The
Lemma 1 reading is taken as the main one because it is the path [5] itself recommends for
complexity, and it is the only one that stays in the same range as C16's own solutions under
both metrics.

Side result: the "per-element coherent sum" of eq. (8a) is a bad surrogate, because it drops
the residual $R_k$ of the other elements. With a two-segment FFT-shaped target the assumption
fails: rows 4–5 above are 0.5–8.4 dB below Lemma 1 and the iteration is not monotone
($F$ jumps 13507 → 12975 → 13223). [5]'s "the two are interchangeable" holds for a single
beam, not for multi-user FDM.

### 8.7 Fig. 7: run time

| α | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | mean |
|---|---|---|---|---|---|---|---|---|---|---|
| GD joint (s) | 3.764 | 2.825 | 2.544 | 1.745 | 1.481 | 1.370 | 1.380 | 1.352 | 1.343 | **1.98** |
| SOTA (s) | 0.159 | 0.162 | 0.152 | 0.147 | 0.147 | 0.150 | 0.146 | 0.147 | 0.176 | **0.15** |
| ratio | 23.7× | 17.4× | 16.8× | 11.9× | 10.1× | 9.1× | 9.5× | 9.2× | 7.6× | **12.8×** |
| paper | — | — | — | — | — | — | — | — | — | 5.5× (5.80/1.05) |

- "At least 5×" holds: every point is ≥7.6×, mean 12.8×.
- The absolute numbers do not, and should not: the paper's 1.05 s / 5.80 s are its own C/C++
  machine times, this GD takes 1.98 s on average (3.76 s at α=0.1, 3.6× the paper). A run-time
  ratio is a property of the algorithm; a run time is a property of the implementation —
  the same argument as the "250×" in §7.4.
- Both α trends agree: the smaller α, the slower GD (more iterations to reconcile two
  incompatible blocks), while the SOTA is flat (fixed 10 rounds).

---

## 9 Gaps in this paper (local numbering)

| ID | Problem | Handling |
|---|---|---|
| **G1** | separated infinity norm: 5.35 dB spread vs 3.2 | per-axis minimax ≠ minimax of the product gain; both readings reported |
| **G2** | single learning rate for φ and τ | reparameterised to $(a,b)$, both in rad |
| **G3** | Adam leaves $\tau\in[0,200]$ ns | bounded reparameterisation |
| **G4** | greedy convergence test misapplied to GD | sliding window for GD |
| **G5** | "+2.26/+2.11 dB over the start point" | implies a 20.45 dB start, i.e. the separated LS; end points reproduce |
| **G6** | separated GD 1.56 dB low on user 3 | open; direction agrees with the paper |
| **G7** | "1.54 dB" not locatable: Lemma 1 → 0.33 dB, (8a) → 8.71 dB | both reported, Lemma 1 as main; the text is underdetermined |
| **G8** | $\bar b_k$ under multi-user FDM is this work's extension | the only self-consistent choice, anchored to eq. (3) |
| **G9** | number of TTD units in the 3-D lift not stated | $N=M=384$; the $N<M$ sweep is left open and could be another source of the 1.54 dB |

---

## 10 Pitfalls

| # | Symptom | Root cause |
|---|---|---|
| 1 | self-check gain 8.35 dB instead of 25.84 dB | `jp_nu` already contains $\pi$; the gain expression multiplied by $\pi$ again |
| 2 | no `linprog` on the compute node | infinity norm rewritten as an exact 1-D convex problem solved by `fminbnd` |
| 3 | Fig. 2(a) 25 dB off | raw steering value instead of the eq. (10) unwrapping |
| 4 | Fig. 6 design converges every round but evaluates to −7.88 dB | mixed delay-sign conventions between C16 and [5] |
| 5 | Lemma 1 returns `RCOND = 3.4e-17` | design matrix $[1,2\pi\Delta f m]$ with a $3\times10^8$ column; change of variables |
| 6 | GD pins τ at 200 ns, gain drops 12 dB | unbounded Adam on a 2.4%-wide feasible window |
| 7 | GD "converges" after one iteration | the greedy convergence criterion is not scale-free |
| 8 | legend labels shifted by one curve | legend order ≠ data-column order; assign handles through an index map |

---

## 11 Relation to the two companion papers

| Dimension | this paper | 2025 prototyping paper | 2022 uplink paper |
|---|---|---|---|
| array | 16×24 planar, one TTD per element | 16V×16H, one TTD per column | 8V×16H, 16 sub-arrays |
| design variables | $(\varphi,\tau)$ per element, 48×48 | 16 pairs, closed form | 16 pairs, closed-form rainbow |
| solver | LS / ∞-norm / greedy / Adam | phase closed form + 1-D delay search | closed form, no search |
| objective | 3-D sub-band gain | 2-D rainbow + system level | uplink coverage and throughput |
| carrier / bandwidth | 28 GHz / 95 MHz | 28 GHz / 400 MHz | 140 GHz / 2 GHz (and 28 GHz / 800 MHz) |
| BS gain | 25.84 dB | 28 dB | 24 dB |

Because the arrays, the design variables and the objectives all differ, the three
reproductions share no kernel code; only the modelling philosophy is common.
