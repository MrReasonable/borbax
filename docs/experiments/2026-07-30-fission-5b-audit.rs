//! PRESERVED SCRATCH — Task 5b's fusion+fission abundance audit.
//!
//! **This is not build input and must not be moved into a crate.** It lives in
//! `docs/experiments/` deliberately, alongside the G2 locality writeup, for
//! three reasons:
//!
//!   1. It uses `.exp()` and `f64::max` — banned by §13.1 for simulation code
//!      and perfectly fine here, because this measures a model rather than
//!      being one. `cargo xtask` scans `crates` and `experiments` and would
//!      (correctly) refuse it in either. It refused it once already, which is
//!      how it ended up here.
//!   2. Under `crates/*/examples/` it is a cargo target, so `--all-targets`
//!      would compile it on every gate run.
//!   3. `docs/` contains no `Cargo.toml`, so this does not create a crate
//!      directory outside the scan roots — the hole `xtask`'s own
//!      escape-hatch check exists to catch.
//!
//! It existed untracked in a single worktree and in no commit on any branch,
//! so this file is a backup, not an endorsement. It has **not** been built or
//! verified against current `main`, from which its branch diverged by ~7,000
//! lines. Task 5b remains scaffolding: it gates nothing, the shipped one-line
//! abundance profile already satisfies §7.2 in 2000 of 2000 universes, and its
//! own steps need a beaker that arrives at Task 19.
//!
//! Its findings are summarised in project memory under Task 5b — chiefly that
//! the shipped profile correlates with its own universe's binding peak at
//! r = -0.024, i.e. carries no information about the chemistry it belongs to.

//! AUDIT SCRATCH — not part of the build. Reproduces and attacks Task 5b's
//! fusion+fission abundance model.
#![allow(clippy::all, clippy::pedantic, missing_docs, unused)]

use borbax_universe::element::generate_elements;

#[derive(Clone)]
pub struct Tab {
    pub n: usize,
    /// `e[i]` = energy_per_unit for cluster size `i+1`
    pub e: Vec<f64>,
    /// `inst[i]` = instability for cluster size `i+1`
    pub inst: Vec<f64>,
    pub valence: Vec<f64>,
    pub peak: usize,
}

pub fn load(seed: u64) -> Tab {
    let t = generate_elements(seed);
    let peak = t.pattern().peak;
    let mut e = Vec::new();
    let mut inst = Vec::new();
    let mut valence = Vec::new();
    for (_, el) in t.iter() {
        e.push(el.energy_per_unit.get());
        inst.push(el.instability);
        valence.push(f64::from(el.valence));
    }
    Tab {
        n: e.len(),
        e,
        inst,
        valence,
        peak,
    }
}

/// Recompute instability from an arbitrary energy series, exactly as
/// `generate_elements` does.
pub fn instability_from(e: &[f64]) -> (Vec<f64>, usize) {
    let mut peak_i = 0usize;
    let mut best = f64::MIN;
    for (i, &v) in e.iter().enumerate() {
        if v > best {
            best = v;
            peak_i = i;
        }
    }
    let peak_energy = e[peak_i];
    let min_scale = f64::EPSILON;
    let scale = {
        let s = peak_energy.abs();
        if s > min_scale { s } else { min_scale }
    };
    let inst = e
        .iter()
        .map(|&v| {
            let d = (peak_energy - v) / scale;
            if d > 1.0 { 1.0 } else { d }
        })
        .collect();
    (inst, peak_i + 1)
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Split {
    /// total outflow `d0*inst*a`, shared uniformly over the splits
    Normalised,
    /// each split fires at `d0*inst*a`; total outflow scales with the split count
    PerSplit,
    /// shed a monomer only
    Alpha,
    /// midpoint only
    Symmetric,
}

#[derive(Clone, Copy)]
pub struct Cfg {
    pub k0: f64,
    pub d0: f64,
    /// fusion gated on Q > 0
    pub fusion_gate: bool,
    /// fission only reverses fusion-admitted pairs
    pub fission_gate: bool,
    pub split: Split,
    /// `None` = use the table's own instability; `Some(c)` = constant
    pub const_inst: Option<f64>,
}

impl Default for Cfg {
    fn default() -> Self {
        Cfg {
            k0: 0.5,
            d0: 0.0,
            fusion_gate: true,
            fission_gate: false,
            split: Split::Normalised,
            const_inst: None,
        }
    }
}

#[inline]
fn q(e: &[f64], i: usize, j: usize) -> f64 {
    // sizes M = i+1, N = j+1, product K = i+j+2 at index i+j+1
    let m = (i + 1) as f64;
    let nn = (j + 1) as f64;
    let p = i + j + 1;
    (m + nn) * e[p] - m * e[i] - nn * e[j]
}

pub fn run(t: &Tab, cfg: Cfg, epochs: usize) -> Vec<f64> {
    let mut a = vec![0.0f64; t.n];
    a[0] = 1.0;
    advance(t, cfg, &mut a, epochs);
    a
}

pub fn advance(t: &Tab, cfg: Cfg, a: &mut [f64], epochs: usize) {
    let n = t.n;
    let inst: Vec<f64> = match cfg.const_inst {
        None => t.inst.clone(),
        Some(c) => vec![c; n],
    };
    let mut d = vec![0.0f64; n];
    for _ in 0..epochs {
        for v in d.iter_mut() {
            *v = 0.0;
        }
        // fusion: i outer ascending, j inner from i; deltas issued at the pair
        for i in 0..n {
            for j in i..n {
                let p = i + j + 1;
                if p >= n {
                    continue;
                }
                if cfg.fusion_gate && !(q(&t.e, i, j) > 0.0) {
                    continue;
                }
                let f = if i == j {
                    0.5 * cfg.k0 * a[i] * a[i]
                } else {
                    cfg.k0 * a[i] * a[j]
                };
                d[i] -= f;
                d[j] -= f;
                d[p] += f;
            }
        }
        // fission
        if cfg.d0 != 0.0 {
            for p in 1..n {
                let rate = cfg.d0 * inst[p] * a[p];
                if rate == 0.0 {
                    continue;
                }
                // unordered splits (i <= j) with i + j + 1 == p
                let mut splits: Vec<(usize, usize)> = Vec::new();
                let hi = p / 2; // i <= j  <=>  i <= (p-1)/2 ; p-1-i >= i
                for i in 0..=((p - 1) / 2) {
                    let j = p - 1 - i;
                    if cfg.fission_gate && !(q(&t.e, i, j) > 0.0) {
                        continue;
                    }
                    splits.push((i, j));
                }
                let chosen: Vec<(usize, usize)> = match cfg.split {
                    Split::Normalised | Split::PerSplit => splits,
                    Split::Alpha => splits.iter().copied().filter(|&(i, _)| i == 0).collect(),
                    Split::Symmetric => {
                        // the most nearly equal split
                        splits.iter().copied().max_by_key(|&(i, _)| i).into_iter().collect()
                    }
                };
                if chosen.is_empty() {
                    continue;
                }
                let s = chosen.len() as f64;
                for &(i, j) in &chosen {
                    let f = match cfg.split {
                        Split::PerSplit => rate,
                        _ => rate / s,
                    };
                    d[p] -= f;
                    d[i] += f;
                    d[j] += f;
                }
            }
        }
        for i in 0..n {
            a[i] += d[i];
        }
    }
}

pub fn mass(a: &[f64]) -> f64 {
    let mut s = 0.0;
    for (i, &v) in a.iter().enumerate() {
        s += (i + 1) as f64 * v;
    }
    s
}

/// argmax over N of the mass in bin N
pub fn mass_argmax(a: &[f64]) -> usize {
    let mut bi = 0usize;
    let mut b = f64::MIN;
    for (i, &v) in a.iter().enumerate() {
        let m = (i + 1) as f64 * v;
        if m > b {
            b = m;
            bi = i;
        }
    }
    bi + 1
}

/// mass-weighted mean cluster size
pub fn mass_mean(a: &[f64]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for (i, &v) in a.iter().enumerate() {
        let nn = (i + 1) as f64;
        num += nn * (nn * v);
        den += nn * v;
    }
    num / den
}

pub fn f_w(a: &[f64], val: &[f64]) -> f64 {
    let mut n1 = 0.0;
    let mut n2 = 0.0;
    for (i, &v) in a.iter().enumerate() {
        n1 += v * val[i];
        n2 += v * val[i] * val[i];
    }
    n2 / n1
}

pub fn top_bin_mass_share(a: &[f64]) -> f64 {
    let n = a.len();
    let lo = n - n / 10;
    let mut tot = 0.0;
    let mut top = 0.0;
    for (i, &v) in a.iter().enumerate() {
        let m = (i + 1) as f64 * v;
        tot += m;
        if i >= lo {
            top += m;
        }
    }
    top / tot
}

pub fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

pub fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let mx = mean(x);
    let my = mean(y);
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for i in 0..x.len() {
        sxy += (x[i] - mx) * (y[i] - my);
        sxx += (x[i] - mx) * (x[i] - mx);
        syy += (y[i] - my) * (y[i] - my);
    }
    sxy / (sxx.sqrt() * syy.sqrt())
}

/// partial correlation of x and y controlling for z
pub fn partial(x: &[f64], y: &[f64], z: &[f64]) -> f64 {
    let rxy = pearson(x, y);
    let rxz = pearson(x, z);
    let ryz = pearson(y, z);
    (rxy - rxz * ryz) / ((1.0 - rxz * rxz).sqrt() * (1.0 - ryz * ryz).sqrt())
}

fn main() {
    let which = std::env::args().nth(1).unwrap_or_else(|| "repro".into());
    match which.as_str() {
        "repro" => repro(),
        "matched" => matched(),
        "decompose" => decompose(),
        "profile" => profile(),
        "flux" => flux(),
        "bound" => bound(),
        "clamp" => clamp(),
        "control" => control(),
        "db" => db(),
        "ratio" => ratio(),
        "balance" => balance(),
        "track" => track(),
        "equil" => equil(),
        "deletion" => deletion(),
        "sweep" => sweep(),
        "intervene" => intervene(),
        "k0" => k0_range(),
        _ => eprintln!("unknown"),
    }
}

fn repro() {
    println!("# reproducing the plan's fission table, 40 seeds, k0=0.5, half kernel, Jacobi");
    println!("d0      epochs   argmax(mean)  massmean  peak(mean)  |mass-1|max   f_w(med)");
    for (d0, ep) in [
        (0.0, 400usize),
        (0.0, 6400),
        (0.2, 1600),
        (0.5, 1600),
        (0.5, 6400),
    ] {
        let mut am = Vec::new();
        let mut mm = Vec::new();
        let mut pk = Vec::new();
        let mut fw = Vec::new();
        let mut worst = 0.0f64;
        for s in 0..40u64 {
            let t = load(s);
            let cfg = Cfg {
                d0,
                ..Default::default()
            };
            let a = run(&t, cfg, ep);
            am.push(mass_argmax(&a) as f64);
            mm.push(mass_mean(&a));
            pk.push(t.peak as f64);
            fw.push(f_w(&a, &t.valence));
            worst = worst.max((mass(&a) - 1.0).abs());
        }
        fw.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "{:<7} {:<8} {:<13.2} {:<9.2} {:<11.2} {:<13.2e} {:.4}",
            d0,
            ep,
            mean(&am),
            mean(&mm),
            mean(&pk),
            worst,
            fw[fw.len() / 2]
        );
    }
}

fn equil() {
    let n_seeds: u64 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20);
    println!("# claim 1: is it a steady state or a slow crawl? {n_seeds} seeds, k0=0.5");
    let ckpt = [
        400usize, 1600, 6400, 25600, 102_400, 409_600, 1_638_400, 6_553_600,
    ];
    for d0 in [0.2f64, 0.5] {
        println!("--- d0 = {d0}");
        println!(
            "epochs      argmax(mean) massmean(mean) top-decile   L1/epoch     seeds w/ argmax moved since prev"
        );
        let tabs: Vec<Tab> = (0..n_seeds).map(load).collect();
        let cfg = Cfg {
            d0,
            ..Default::default()
        };
        let mut states: Vec<Vec<f64>> = tabs
            .iter()
            .map(|t| {
                let mut a = vec![0.0; t.n];
                a[0] = 1.0;
                a
            })
            .collect();
        let mut prev: Option<(Vec<Vec<f64>>, usize)> = None;
        let mut done = 0usize;
        for &ep in &ckpt {
            for (t, a) in tabs.iter().zip(states.iter_mut()) {
                advance(t, cfg, a, ep - done);
            }
            done = ep;
            let am: Vec<f64> = states.iter().map(|a| mass_argmax(a) as f64).collect();
            let mm: Vec<f64> = states.iter().map(|a| mass_mean(a)).collect();
            let td: Vec<f64> = states.iter().map(|a| top_bin_mass_share(a)).collect();
            let (drift, moved) = if let Some((ref p, pe)) = prev {
                let mut tot = 0.0;
                let mut mv = 0usize;
                for s in 0..states.len() {
                    let mut l1 = 0.0;
                    for i in 0..states[s].len() {
                        l1 += (states[s][i] - p[s][i]).abs();
                    }
                    tot += l1;
                    if mass_argmax(&states[s]) != mass_argmax(&p[s]) {
                        mv += 1;
                    }
                }
                (tot / states.len() as f64 / (ep - pe) as f64, mv)
            } else {
                (f64::NAN, 0)
            };
            println!(
                "{:<11} {:<12.3} {:<14.4} {:<12.6} {:<12.3e} {}/{}",
                ep,
                mean(&am),
                mean(&mm),
                mean(&td),
                drift,
                moved,
                states.len()
            );
            prev = Some((states.clone(), ep));
        }
    }
}

fn deletion() {
    println!("# claim 2 deletion test: replace instability(K) with a constant.");
    println!("# 40 seeds, k0=0.5, 1600 epochs. peak(mean) is the binding peak.");
    println!();
    let n_seeds = 40u64;
    let mut peaks = Vec::new();
    let mut nel = Vec::new();
    for s in 0..n_seeds {
        let t = load(s);
        peaks.push(t.peak as f64);
        nel.push(t.n as f64);
    }
    println!(
        "binding peak: mean {:.2}   n_elements mean {:.2}   r(peak, n) = {:.3}",
        mean(&peaks),
        mean(&nel),
        pearson(&peaks, &nel)
    );
    println!();
    println!("variant                         d0     argmax  massmean  r(am,peak)  r(am,peak | n)  r(am,n)");
    for (label, ci) in [
        ("real instability(K)", None),
        ("CONST 0.5", Some(0.5f64)),
        ("CONST = mean inst", Some(-1.0)),
    ] {
        for d0 in [0.2f64, 0.5, 1.0] {
            let mut am = Vec::new();
            let mut mm = Vec::new();
            for s in 0..n_seeds {
                let t = load(s);
                let c = match ci {
                    Some(x) if x < 0.0 => Some(mean(&t.inst)),
                    other => other,
                };
                let a = run(
                    &t,
                    Cfg {
                        d0,
                        const_inst: c,
                        ..Default::default()
                    },
                    1600,
                );
                am.push(mass_argmax(&a) as f64);
                mm.push(mass_mean(&a));
            }
            println!(
                "{:<31} {:<6} {:<7.2} {:<9.2} {:<11.3} {:<15.3} {:.3}",
                label,
                d0,
                mean(&am),
                mean(&mm),
                pearson(&am, &peaks),
                partial(&am, &peaks, &nel),
                pearson(&am, &nel)
            );
        }
    }
}

fn sweep() {
    let n_seeds: u64 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(120);
    println!("# claims 3 and 4: does d0 replace epochs? is the heavy end bounded across d0?");
    println!("# {n_seeds} seeds, k0=0.5, 6400 epochs");
    println!("d0      argmax  massmean  top-decile  maxbin-share  f_w(med)  f_w<=2  p_c(med)  minbin");
    let tabs: Vec<Tab> = (0..n_seeds).map(load).collect();
    for d0 in [
        0.0f64, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0,
    ] {
        let mut td = Vec::new();
        let mut am = Vec::new();
        let mut mm = Vec::new();
        let mut fw = Vec::new();
        let mut pc = Vec::new();
        let mut topshare = Vec::new();
        let mut minbin = f64::MAX;
        let mut bad = 0usize;
        for t in &tabs {
            let a = run(
                t,
                Cfg {
                    d0,
                    ..Default::default()
                },
                6400,
            );
            td.push(top_bin_mass_share(&a));
            am.push(mass_argmax(&a) as f64);
            mm.push(mass_mean(&a));
            let w = f_w(&a, &t.valence);
            fw.push(w);
            pc.push(1.0 / (w - 1.0));
            if !(w > 2.0) {
                bad += 1;
            }
            let m: f64 = a.iter().enumerate().map(|(i, &v)| (i + 1) as f64 * v).sum();
            let mx = a
                .iter()
                .enumerate()
                .map(|(i, &v)| (i + 1) as f64 * v / m)
                .fold(0.0f64, f64::max);
            topshare.push(mx);
            for &v in &a {
                if v < minbin {
                    minbin = v;
                }
            }
        }
        fw.sort_by(|a, b| a.partial_cmp(b).unwrap());
        pc.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "{:<7} {:<7.2} {:<9.2} {:<11.5} {:<13.4} {:<9.4} {:<7} {:<9.4} {:.2e}",
            d0,
            mean(&am),
            mean(&mm),
            mean(&td),
            mean(&topshare),
            fw[fw.len() / 2],
            bad,
            pc[pc.len() / 2],
            minbin
        );
    }
}

/// Rate-matched deletion test. For each seed, bisect `d0` under a CONSTANT
/// instability until the mass-weighted mean cluster size matches what the real
/// `instability(K)` produced at the reference `d0`. Then compare where the
/// argmax sits. Matching the mass-mean removes the confound that a constant
/// changes the total fission rate as well as its shape.
fn matched() {
    let n_seeds: u64 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(120);
    let epochs = 1600usize;
    println!("# rate-matched deletion test, {n_seeds} seeds, k0=0.5, {epochs} epochs");
    println!("# `d0*` bisected per seed so CONST reproduces the real model's mass-mean.");
    println!();
    for d0ref in [0.2f64, 0.5] {
        let mut real_am = Vec::new();
        let mut const_am = Vec::new();
        let mut peaks = Vec::new();
        let mut nel = Vec::new();
        let mut mm_real = Vec::new();
        let mut mm_const = Vec::new();
        let mut d0stars = Vec::new();
        for s in 0..n_seeds {
            let t = load(s);
            let a_real = run(
                &t,
                Cfg {
                    d0: d0ref,
                    ..Default::default()
                },
                epochs,
            );
            let target = mass_mean(&a_real);
            // bisect d0 for const_inst = 0.5; mass_mean is monotone decreasing in d0
            let (mut lo, mut hi) = (0.0f64, 4.0f64);
            let mut best = vec![];
            for _ in 0..40 {
                let mid = 0.5 * (lo + hi);
                let a = run(
                    &t,
                    Cfg {
                        d0: mid,
                        const_inst: Some(0.5),
                        ..Default::default()
                    },
                    epochs,
                );
                let m = mass_mean(&a);
                if m > target {
                    lo = mid;
                } else {
                    hi = mid;
                }
                best = a;
                if (m - target).abs() / target < 1e-4 {
                    break;
                }
            }
            d0stars.push(0.5 * (lo + hi));
            real_am.push(mass_argmax(&a_real) as f64);
            const_am.push(mass_argmax(&best) as f64);
            mm_real.push(target);
            mm_const.push(mass_mean(&best));
            peaks.push(t.peak as f64);
            nel.push(t.n as f64);
        }
        let mae = |x: &[f64], y: &[f64]| -> f64 {
            x.iter().zip(y).map(|(a, b)| (a - b).abs()).sum::<f64>() / x.len() as f64
        };
        let flat_real = vec![mean(&real_am); real_am.len()];
        let flat_const = vec![mean(&const_am); const_am.len()];
        println!("=== reference d0 = {d0ref}");
        println!(
            "  peak:  mean {:.2}  sd {:.2}   n_elements mean {:.1}  r(peak,n) {:.3}",
            mean(&peaks),
            sd(&peaks),
            mean(&nel),
            pearson(&peaks, &nel)
        );
        println!(
            "  REAL inst : argmax mean {:.2} sd {:.2} | massmean {:.2} | r(am,peak) {:.3} (p {:.3}) | r(am,peak|n) {:.3} | MAE-to-peak {:.2} vs MAE-to-own-mean {:.2}",
            mean(&real_am), sd(&real_am), mean(&mm_real),
            pearson(&real_am, &peaks), pval(pearson(&real_am, &peaks), real_am.len()),
            partial(&real_am, &peaks, &nel),
            mae(&real_am, &peaks), mae(&real_am, &flat_real)
        );
        println!(
            "  CONST 0.5 : argmax mean {:.2} sd {:.2} | massmean {:.2} | r(am,peak) {:.3} (p {:.3}) | r(am,peak|n) {:.3} | MAE-to-peak {:.2} vs MAE-to-own-mean {:.2}   [d0* mean {:.3}]",
            mean(&const_am), sd(&const_am), mean(&mm_const),
            pearson(&const_am, &peaks), pval(pearson(&const_am, &peaks), const_am.len()),
            partial(&const_am, &peaks, &nel),
            mae(&const_am, &peaks), mae(&const_am, &flat_const),
            mean(&d0stars)
        );
        println!();
    }
}

/// Claim 4 rests on "heavy clusters are the *least* stable, so fission removes
/// them fastest". Measure the instability profile either side of the peak.
fn profile() {
    let n_seeds: u64 = 500;
    println!("# `instability` either side of the binding peak, {n_seeds} seeds");
    println!("# instability = clamp((E_peak - E(N))/|E_peak|, <=1).  Light end vs heavy end.");
    println!();
    let mut below = Vec::new();
    let mut above = Vec::new();
    let mut inst_top = Vec::new();
    let mut inst_mono = Vec::new();
    let mut inst_dimer = Vec::new();
    let mut ratio = Vec::new();
    let mut clamped_below = Vec::new();
    for s in 0..n_seeds {
        let t = load(s);
        let p = t.peak; // 1-based
        let lo: Vec<f64> = t.inst[..p.saturating_sub(1)].to_vec();
        let hi: Vec<f64> = t.inst[p..].to_vec();
        if lo.is_empty() || hi.is_empty() {
            continue;
        }
        below.push(mean(&lo));
        above.push(mean(&hi));
        inst_top.push(t.inst[t.n - 1]);
        inst_mono.push(t.inst[0]);
        inst_dimer.push(t.inst[1]);
        ratio.push(mean(&lo) / mean(&hi).max(1e-300));
        clamped_below.push(lo.iter().filter(|&&x| x >= 1.0).count() as f64 / lo.len() as f64);
    }
    println!("mean instability BELOW the peak (light side): {:.4}", mean(&below));
    println!("mean instability ABOVE the peak (heavy side): {:.4}", mean(&above));
    println!("ratio light/heavy: mean {:.1}x, median {:.1}x", mean(&ratio), {
        let mut r = ratio.clone();
        r.sort_by(|a, b| a.partial_cmp(b).unwrap());
        r[r.len() / 2]
    });
    println!();
    println!("instability(monomer)      mean {:.4}", mean(&inst_mono));
    println!("instability(dimer)        mean {:.4}", mean(&inst_dimer));
    println!("instability(heaviest)     mean {:.4}", mean(&inst_top));
    println!(
        "fraction of light-side elements CLAMPED at 1.0: {:.4}",
        mean(&clamped_below)
    );
    println!();
    println!("seed  n   peak  inst(1)  inst(2)  inst(peak+1)  inst(n)   mean<peak  mean>peak");
    for s in 0..8u64 {
        let t = load(s);
        let p = t.peak;
        println!(
            "{:<5} {:<3} {:<5} {:<8.4} {:<8.4} {:<13.4} {:<9.4} {:<10.4} {:.4}",
            s,
            t.n,
            p,
            t.inst[0],
            t.inst[1],
            if p < t.n { t.inst[p] } else { f64::NAN },
            t.inst[t.n - 1],
            mean(&t.inst[..p - 1]),
            if p < t.n { mean(&t.inst[p..]) } else { f64::NAN }
        );
    }
}

/// The control the plan never ran: with NO fission at all, does the argmax
/// already track the binding peak? If it does, fission adds no peak-tracking.
fn control() {
    let n_seeds = 200u64;
    let tabs: Vec<Tab> = (0..n_seeds).map(load).collect();
    let peaks: Vec<f64> = tabs.iter().map(|t| t.peak as f64).collect();
    let nel: Vec<f64> = tabs.iter().map(|t| t.n as f64).collect();
    println!("# the missing control: fusion only (d0 = 0), {n_seeds} seeds, k0 = 0.5");
    println!(
        "# peak mean {:.2} sd {:.2}; n_elements mean {:.1} sd {:.1}; r(peak,n) {:.3}",
        mean(&peaks),
        sd(&peaks),
        mean(&nel),
        sd(&nel),
        pearson(&peaks, &nel)
    );
    println!();
    println!("gate   epochs   argmax(mean)  r(am,peak)  r(am,peak|n)  MAE->peak  r(am,n)");
    for gate in [true, false] {
        for ep in [100usize, 400, 1600, 6400] {
            let am: Vec<f64> = tabs
                .iter()
                .map(|t| {
                    mass_argmax(&run(
                        t,
                        Cfg {
                            d0: 0.0,
                            fusion_gate: gate,
                            ..Default::default()
                        },
                        ep,
                    )) as f64
                })
                .collect();
            println!(
                "{:<6} {:<8} {:<13.2} {:<11.3} {:<13.3} {:<10.2} {:.3}",
                gate,
                ep,
                mean(&am),
                pearson(&am, &peaks),
                partial(&am, &peaks, &nel),
                am.iter().zip(&peaks).map(|(a, b)| (a - b).abs()).sum::<f64>() / am.len() as f64,
                pearson(&am, &nel)
            );
        }
    }
}

/// Open item 4: the recorded `> 1.0` category error now has a second consumer.
/// How often does it bind on a bin that can actually fission, and what does it
/// change?
fn clamp() {
    let n_seeds = 200u64;
    println!("# open item 4: the `> 1.0` clamp's second consumer");
    let mut clamped_fissionable = Vec::new();
    let mut any = 0usize;
    for s in 0..n_seeds {
        let t = load(s);
        // recompute the raw deficit
        let pe = t.e[t.peak - 1];
        let scale = pe.abs().max(f64::EPSILON);
        let raw: Vec<f64> = t.e.iter().map(|&v| (pe - v) / scale).collect();
        let c = raw[1..].iter().filter(|&&x| x > 1.0).count();
        clamped_fissionable.push(c as f64);
        if c > 0 {
            any += 1;
        }
    }
    println!(
        "elements with units>=2 (i.e. fissionable) whose raw deficit exceeds 1.0:\n  mean {:.3} per universe, at least one in {}/{} universes",
        mean(&clamped_fissionable),
        any,
        n_seeds
    );
    println!();
    println!("# effect of removing the clamp on the fission channel");
    println!("clamp   d0     argmax  massmean  f_w(med)  f_w<=2");
    let tabs: Vec<Tab> = (0..n_seeds).map(load).collect();
    for clamped in [true, false] {
        for d0 in [0.2f64, 0.5] {
            let mut am = Vec::new();
            let mut mm = Vec::new();
            let mut fw = Vec::new();
            let mut bad = 0usize;
            for t in &tabs {
                let mut t2 = t.clone();
                if !clamped {
                    let pe = t.e[t.peak - 1];
                    let scale = pe.abs().max(f64::EPSILON);
                    t2.inst = t.e.iter().map(|&v| (pe - v) / scale).collect();
                }
                let a = run(
                    &t2,
                    Cfg {
                        d0,
                        ..Default::default()
                    },
                    1600,
                );
                am.push(mass_argmax(&a) as f64);
                mm.push(mass_mean(&a));
                let w = f_w(&a, &t.valence);
                fw.push(w);
                if !(w > 2.0) {
                    bad += 1;
                }
            }
            fw.sort_by(|a, b| a.partial_cmp(b).unwrap());
            println!(
                "{:<7} {:<6} {:<7.2} {:<9.2} {:<9.4} {}/{}",
                clamped,
                d0,
                mean(&am),
                mean(&mm),
                fw[fw.len() / 2],
                bad,
                n_seeds
            );
        }
    }
}

/// Is the plan right that a *properly reversible* fission gives the rejected
/// Boltzmann-in-total-binding distribution? Rate for channel (M,N) out of K is
/// `d0 * exp(-Q(M,N)/T)` on exactly the channels fusion admits — the genuine
/// reverse, obeying detailed balance.
fn db() {
    let n_seeds = 60u64;
    println!("# a detailed-balance fission: rate = d0*exp(-Q(M,N)/T) on the fusion channel set");
    println!("# (the genuine reverse; `instability(K)` is a different function of the table)");
    println!("T       d0     argmax  massmean  heaviest-is-argmax  f_w(med)  top-decile");
    let tabs: Vec<Tab> = (0..n_seeds).map(load).collect();
    for tmp in [1.0f64, 3.0, 10.0, 30.0, 100.0] {
        for d0 in [0.5f64] {
            let mut am = Vec::new();
            let mut mm = Vec::new();
            let mut fw = Vec::new();
            let mut td = Vec::new();
            let mut top = 0usize;
            for t in &tabs {
                let n = t.n;
                let mut a = vec![0.0f64; n];
                a[0] = 1.0;
                let mut d = vec![0.0f64; n];
                for _ in 0..4000 {
                    for v in d.iter_mut() {
                        *v = 0.0;
                    }
                    for i in 0..n {
                        for j in i..n {
                            let p = i + j + 1;
                            if p >= n {
                                continue;
                            }
                            let qq = q(&t.e, i, j);
                            if !(qq > 0.0) {
                                continue;
                            }
                            let f = if i == j {
                                0.5 * 0.5 * a[i] * a[i]
                            } else {
                                0.5 * a[i] * a[j]
                            };
                            d[i] -= f;
                            d[j] -= f;
                            d[p] += f;
                            // reverse on the same channel
                            let r = d0 * (-qq / tmp).exp() * a[p];
                            d[p] -= r;
                            d[i] += r;
                            d[j] += r;
                        }
                    }
                    for i in 0..n {
                        a[i] += d[i];
                    }
                    if a.iter().any(|&v| !v.is_finite()) {
                        break;
                    }
                }
                if a.iter().any(|&v| !v.is_finite()) {
                    continue;
                }
                let g = mass_argmax(&a);
                am.push(g as f64);
                mm.push(mass_mean(&a));
                fw.push(f_w(&a, &t.valence));
                td.push(top_bin_mass_share(&a));
                if g == n {
                    top += 1;
                }
            }
            if am.is_empty() {
                println!("{tmp:<7} {d0:<6} (all runs non-finite)");
                continue;
            }
            fw.sort_by(|a, b| a.partial_cmp(b).unwrap());
            println!(
                "{:<7} {:<6} {:<7.2} {:<9.2} {:<19} {:<9.4} {:.4}",
                tmp,
                d0,
                mean(&am),
                mean(&mm),
                format!("{}/{}", top, am.len()),
                fw[fw.len() / 2],
                mean(&td)
            );
        }
    }
}

/// Open item 1: `d0`'s stability bound, which the plan says has not been derived.
fn bound() {
    println!("# d0 stability bound: first epoch at which any bin goes negative / non-finite");
    println!("# 60 seeds, 400 epochs, k0 = 0.5 (half kernel)");
    println!("d0        seeds w/ negative bin   seeds w/ non-finite   worst bin value");
    let tabs: Vec<Tab> = (0..60u64).map(load).collect();
    for d0 in [
        0.5f64, 0.9, 0.99, 1.0, 1.0000001, 1.01, 1.05, 1.1, 1.25, 1.5, 2.0,
    ] {
        let mut neg = 0usize;
        let mut nf = 0usize;
        let mut worst = 0.0f64;
        for t in tabs.iter() {
            let mut a = vec![0.0; t.n];
            a[0] = 1.0;
            let cfg = Cfg {
                d0,
                ..Default::default()
            };
            let mut sawneg = false;
            let mut sawnf = false;
            for _ in 0..400 {
                advance(t, cfg, &mut a, 1);
                for &v in &a {
                    if v < 0.0 {
                        sawneg = true;
                        if v < worst {
                            worst = v;
                        }
                    }
                    if !v.is_finite() {
                        sawnf = true;
                    }
                }
                if sawnf {
                    break;
                }
            }
            if sawneg {
                neg += 1;
            }
            if sawnf {
                nf += 1;
            }
        }
        println!("{d0:<9} {neg:<23} {nf:<21} {worst:.3e}");
    }
    println!();
    println!("# combined: the fusion outflow from bin i is k0*a[i]*(a[i]+C) and the");
    println!("# fission outflow is d0*inst[i]*a[i], so positivity needs k0*(a+C) + d0*inst <= 1.");
    println!("# With C <= 1 and inst <= 1 that is k0 + d0 <= 1 at worst, NOT d0 <= 1 alone.");
    println!("k0 + d0   k0     d0     seeds w/ negative bin");
    for (k0, d0) in [
        (0.5f64, 0.5f64),
        (0.5, 0.51),
        (0.5, 0.6),
        (0.9, 0.1),
        (0.9, 0.11),
        (1.0, 0.01),
    ] {
        let mut neg = 0usize;
        for t in tabs.iter() {
            let mut a = vec![0.0; t.n];
            a[0] = 1.0;
            let cfg = Cfg {
                k0,
                d0,
                ..Default::default()
            };
            let mut sawneg = false;
            for _ in 0..400 {
                advance(t, cfg, &mut a, 1);
                if a.iter().any(|&v| v < 0.0 || !v.is_finite()) {
                    sawneg = true;
                    break;
                }
            }
            if sawneg {
                neg += 1;
            }
        }
        println!("{:<9.2} {:<6} {:<6} {}", k0 + d0, k0, d0, neg);
    }
}

/// Where does the fission MASS flux actually go? Claim 4 says the heavy end.
fn flux() {
    let n_seeds = 200u64;
    println!("# claim 4: 'heavy clusters are the least stable, so fission removes them fastest'");
    println!("# mass fission flux by cluster size, at the state the plan calls equilibrium");
    println!("# quintiles of the table by size (Q1 = lightest 20%).");
    println!();
    println!("d0     epochs  Q1flux   Q2flux   Q3flux   Q4flux   Q5flux  | Q1mass  Q5mass | flux above peak");
    for (d0, ep) in [(0.2f64, 1600usize), (0.5, 1600), (0.5, 6400)] {
        let mut qs = [0.0f64; 5];
        let mut q1m = 0.0;
        let mut q5m = 0.0;
        let mut above = 0.0;
        let mut nn = 0.0;
        for s in 0..n_seeds {
            let t = load(s);
            let a = run(
                &t,
                Cfg {
                    d0,
                    ..Default::default()
                },
                ep,
            );
            // instantaneous fission mass outflow per bin
            let mut f = vec![0.0f64; t.n];
            for p in 1..t.n {
                f[p] = d0 * t.inst[p] * a[p] * (p + 1) as f64;
            }
            let tot: f64 = f.iter().sum();
            if tot <= 0.0 {
                continue;
            }
            let mut local = [0.0f64; 5];
            let mut abv = 0.0;
            for (i, &v) in f.iter().enumerate() {
                let qi = (i * 5 / t.n).min(4);
                local[qi] += v / tot;
                if i + 1 > t.peak {
                    abv += v / tot;
                }
            }
            for k in 0..5 {
                qs[k] += local[k];
            }
            above += abv;
            let m: f64 = a
                .iter()
                .enumerate()
                .map(|(i, &v)| (i + 1) as f64 * v)
                .sum();
            let mut lm = [0.0f64; 5];
            for (i, &v) in a.iter().enumerate() {
                lm[(i * 5 / t.n).min(4)] += (i + 1) as f64 * v / m;
            }
            q1m += lm[0];
            q5m += lm[4];
            nn += 1.0;
        }
        println!(
            "{:<6} {:<7} {:<8.4} {:<8.4} {:<8.4} {:<8.4} {:<8.4}| {:<7.4} {:<7.4}| {:.4}",
            d0,
            ep,
            qs[0] / nn,
            qs[1] / nn,
            qs[2] / nn,
            qs[3] / nn,
            qs[4] / nn,
            q1m / nn,
            q5m / nn,
            above / nn
        );
    }
}

/// Claim 3: with a steady state the limit should be set by the fusion/fission
/// balance. If so it depends on `d0/k0` alone. Test it.
fn ratio() {
    let n_seeds = 40u64;
    println!("# claim 3: is the long-time state a function of d0/k0 alone?");
    println!("# If yes, `d0` is a physical ratio. If no, `k0` and `epochs` still set the answer.");
    println!("k0     d0     d0/k0  epochs  argmax(mean)  massmean  f_w(med)");
    for (k0, d0) in [
        (0.125f64, 0.05f64),
        (0.25, 0.1),
        (0.5, 0.2),
        (1.0, 0.4),
        (0.125, 0.125),
        (0.25, 0.25),
        (0.5, 0.5),
    ] {
        for ep in [1600usize, 25600] {
            let mut am = Vec::new();
            let mut mm = Vec::new();
            let mut fw = Vec::new();
            for s in 0..n_seeds {
                let t = load(s);
                let a = run(
                    &t,
                    Cfg {
                        k0,
                        d0,
                        ..Default::default()
                    },
                    ep,
                );
                am.push(mass_argmax(&a) as f64);
                mm.push(mass_mean(&a));
                fw.push(f_w(&a, &t.valence));
            }
            fw.sort_by(|a, b| a.partial_cmp(b).unwrap());
            println!(
                "{:<6} {:<6} {:<6.2} {:<7} {:<13.2} {:<9.3} {:.4}",
                k0,
                d0,
                d0 / k0,
                ep,
                mean(&am),
                mean(&mm),
                fw[fw.len() / 2]
            );
        }
    }
}

/// Open item 3: fission is ungated, fusion is gated on `Q > 0`. How much fission
/// flux runs through channels fusion cannot reverse? If it is large, the
/// "equilibrium" is a cyclic non-equilibrium steady state, not a balance.
fn balance() {
    let n_seeds = 200u64;
    println!("# open item 3: detailed balance");
    println!("# share of fission flux through pairs the fusion gate REFUSES");
    println!();
    let mut irrev_pairs = Vec::new();
    for s in 0..n_seeds {
        let t = load(s);
        let mut tot = 0usize;
        let mut bad = 0usize;
        for p in 1..t.n {
            for i in 0..=((p - 1) / 2) {
                let j = p - 1 - i;
                tot += 1;
                if !(q(&t.e, i, j) > 0.0) {
                    bad += 1;
                }
            }
        }
        irrev_pairs.push(bad as f64 / tot as f64);
    }
    println!(
        "share of (M,N) splits that fusion refuses: mean {:.4}, max {:.4}",
        mean(&irrev_pairs),
        irrev_pairs.iter().cloned().fold(0.0f64, f64::max)
    );
    println!();
    println!("# and what gating fission does to the result");
    println!("fission-gate  d0     argmax(mean)  massmean  f_w(med)  r(am,peak)");
    let tabs: Vec<Tab> = (0..n_seeds).map(load).collect();
    let peaks: Vec<f64> = tabs.iter().map(|t| t.peak as f64).collect();
    for fg in [false, true] {
        for d0 in [0.2f64, 0.5] {
            let cfg = Cfg {
                d0,
                fission_gate: fg,
                ..Default::default()
            };
            let mut am = Vec::new();
            let mut mm = Vec::new();
            let mut fw = Vec::new();
            for t in &tabs {
                let a = run(t, cfg, 1600);
                am.push(mass_argmax(&a) as f64);
                mm.push(mass_mean(&a));
                fw.push(f_w(&a, &t.valence));
            }
            fw.sort_by(|a, b| a.partial_cmp(b).unwrap());
            println!(
                "{:<13} {:<6} {:<13.2} {:<9.2} {:<9.4} {:.3}",
                fg,
                d0,
                mean(&am),
                mean(&mm),
                fw[fw.len() / 2],
                pearson(&am, &peaks)
            );
        }
    }
}

/// Which channel carries the peak-tracking: the fusion exothermicity gate, or
/// the fission energetics? Sweep `d0` for each variant and read `r(argmax, peak)`
/// at matched argmax level.
fn decompose() {
    let n_seeds: u64 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);
    let tabs: Vec<Tab> = (0..n_seeds).map(load).collect();
    let peaks: Vec<f64> = tabs.iter().map(|t| t.peak as f64).collect();
    let nel: Vec<f64> = tabs.iter().map(|t| t.n as f64).collect();
    println!("# channel decomposition, {n_seeds} seeds, k0=0.5, 1600 epochs");
    println!(
        "# peak: mean {:.2} sd {:.2}; n_elements mean {:.1}",
        mean(&peaks),
        sd(&peaks),
        mean(&nel)
    );
    println!();
    println!("fusion-gate  fission          d0      argmax(mean)  r(am,peak)  MAE->peak  MAE->median  massmean");
    for gate in [true, false] {
        for (flabel, ci) in [
            ("none", Some(f64::NAN)),
            ("real inst(K)", None),
            ("CONST 0.5", Some(0.5f64)),
        ] {
            let d0s: Vec<f64> = if ci == Some(f64::NAN) {
                vec![0.0]
            } else if ci.is_none() {
                vec![0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
            } else {
                vec![0.005, 0.01, 0.02, 0.03, 0.05, 0.1]
            };
            for d0 in d0s {
                let cfg = Cfg {
                    d0,
                    fusion_gate: gate,
                    const_inst: if ci == Some(f64::NAN) { None } else { ci },
                    ..Default::default()
                };
                let am: Vec<f64> = tabs.iter().map(|t| mass_argmax(&run(t, cfg, 1600)) as f64).collect();
                let mm: Vec<f64> = tabs.iter().map(|t| mass_mean(&run(t, cfg, 1600))).collect();
                let mut sorted = am.clone();
                sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let med = sorted[sorted.len() / 2];
                println!(
                    "{:<12} {:<16} {:<7} {:<13.2} {:<11.3} {:<10.2} {:<12.2} {:.2}",
                    gate,
                    flabel,
                    d0,
                    mean(&am),
                    pearson(&am, &peaks),
                    am.iter().zip(&peaks).map(|(a, b)| (a - b).abs()).sum::<f64>() / am.len() as f64,
                    am.iter().map(|a| (a - med).abs()).sum::<f64>() / am.len() as f64,
                    mean(&mm)
                );
            }
        }
    }
}

pub fn sd(v: &[f64]) -> f64 {
    let m = mean(v);
    (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (v.len() as f64 - 1.0)).sqrt()
}

/// two-tailed p for a Pearson r via the t approximation, normal tail
pub fn pval(r: f64, n: usize) -> f64 {
    if !r.is_finite() || n < 4 {
        return f64::NAN;
    }
    let t = r * ((n as f64 - 2.0).sqrt()) / (1.0 - r * r).sqrt();
    // normal approximation to the two-tailed p
    let z = t.abs();
    let p = (-0.5 * z * z).exp() / (z * (2.0 * std::f64::consts::PI).sqrt() + 1e-300);
    (2.0 * p).min(1.0)
}

/// Does the argmax track the *peak*, or just the cross-seed mean?
fn track() {
    let n_seeds: u64 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);
    println!("# claim 2, tracking test, {n_seeds} seeds, k0=0.5, 1600 epochs");
    println!("# 'tracks the peak' requires MAE(argmax, peak_s) < MAE(argmax, constant)");
    println!("d0     argmax(mean)  peak(mean)  MAE->peak_s  MAE->best const  r(am,peak)   p");
    for d0 in [0.1f64, 0.15, 0.2, 0.25, 0.3, 0.5] {
        let mut am = Vec::new();
        let mut pk = Vec::new();
        for s in 0..n_seeds {
            let t = load(s);
            let a = run(
                &t,
                Cfg {
                    d0,
                    ..Default::default()
                },
                1600,
            );
            am.push(mass_argmax(&a) as f64);
            pk.push(t.peak as f64);
        }
        // best constant predictor under MAE is the median of am
        let mut sorted = am.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = sorted[sorted.len() / 2];
        let mae_peak =
            am.iter().zip(&pk).map(|(a, b)| (a - b).abs()).sum::<f64>() / am.len() as f64;
        let mae_const = am.iter().map(|a| (a - med).abs()).sum::<f64>() / am.len() as f64;
        let r = pearson(&am, &pk);
        println!(
            "{:<6} {:<13.2} {:<11.2} {:<12.2} {:<16.2} {:<12.3} {:.4}",
            d0,
            mean(&am),
            mean(&pk),
            mae_peak,
            mae_const,
            r,
            pval(r, am.len())
        );
    }
}

fn intervene() {
    // Hold n fixed, move the binding peak by construction, ask whether the mass
    // argmax follows.
    println!("# claim 2, intervention: synthetic energy series, n = 100 fixed,");
    println!("# E(N) = a - b*(N - P)^2 / P^2 shaped like the real one, peak forced to P.");
    println!("P(forced)  d0=0.2 argmax  d0=0.5 argmax  d0=1.0 argmax   [n=100]");
    let n = 100usize;
    for p in [20usize, 35, 50, 65, 80] {
        let e: Vec<f64> = (1..=n)
            .map(|nn| {
                let x = (nn as f64 - p as f64) / p as f64;
                4.5 - 3.0 * x * x
            })
            .collect();
        let (inst, gotpeak) = instability_from(&e);
        let t = Tab {
            n,
            e,
            inst,
            valence: vec![2.0; n],
            peak: gotpeak,
        };
        let mut row = Vec::new();
        for d0 in [0.2f64, 0.5, 1.0] {
            let a = run(&t, Cfg { d0, ..Default::default() }, 1600);
            row.push(mass_argmax(&a) as f64);
        }
        println!(
            "{:<10} {:<15.1} {:<15.1} {:<15.1} peak read back = {}",
            p, row[0], row[1], row[2], gotpeak
        );
    }
}

fn k0_range() {
    println!("# k0: the plan's measurements are all at 0.5; the slot-reuse line gives 0.04..0.16");
    println!("k0      d0    ep       argmax  massmean  f_w(med)");
    for k0 in [0.04f64, 0.1, 0.16, 0.5, 1.0] {
        for (d0, ep) in [(0.5f64, 1600usize), (0.5, 6400)] {
            let mut am = Vec::new();
            let mut mm = Vec::new();
            let mut fw = Vec::new();
            for s in 0..40u64 {
                let t = load(s);
                let a = run(&t, Cfg { k0, d0, ..Default::default() }, ep);
                am.push(mass_argmax(&a) as f64);
                mm.push(mass_mean(&a));
                fw.push(f_w(&a, &t.valence));
            }
            fw.sort_by(|a, b| a.partial_cmp(b).unwrap());
            println!(
                "{:<7} {:<5} {:<8} {:<7.2} {:<9.2} {:.4}",
                k0,
                d0,
                ep,
                mean(&am),
                mean(&mm),
                fw[fw.len() / 2]
            );
        }
    }
}
