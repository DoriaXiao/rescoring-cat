"""
CAT Rescoring Explorer
Interactive companion to: "Recovering Accuracy From Miscalibrated Item Banks: Scoring Error Is Fixable, Selection Error Is Not"
Xiao, Domingue, Ram, & Frank (2026), Educational and Psychological Measurement, https://doi.org/10.1177/00131644261486854

Deploy on Hugging Face Spaces with Gradio SDK.
"""

import numpy as np
import gradio as gr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

# ── IRT Functions ─────────────────────────────────────────────

def p_2pl(theta, a, b):
    return 1.0 / (1.0 + np.exp(-a * (theta - b)))

def info_2pl(theta, a, b):
    p = p_2pl(theta, a, b)
    return a**2 * p * (1 - p)

def eap_score(responses, a, b, grid, prior_sd=1.5):
    ll = np.zeros(len(grid))
    for j in range(len(responses)):
        p = p_2pl(grid, a[j], b[j])
        p = np.clip(p, 1e-10, 1 - 1e-10)
        ll += responses[j] * np.log(p) + (1 - responses[j]) * np.log(1 - p)
    lp = ll + (-0.5 * (grid / prior_sd)**2)
    w = np.exp(lp - np.max(lp))
    w /= w.sum()
    return np.sum(grid * w)


# ── CAT Engine ────────────────────────────────────────────────

GRID = np.linspace(-6, 6, 241)

def run_cat(theta_true, pool_a, pool_b, select_a, select_b, score_a, score_b,
            cat_length, return_trajectory=False):
    n_items = len(pool_a)
    used = np.zeros(n_items, dtype=bool)
    responses = []
    items_given = []
    theta_hat = 0.0
    traj = [] if return_trajectory else None

    for k in range(cat_length):
        info_vals = info_2pl(theta_hat, select_a, select_b)
        info_vals[used] = -np.inf
        chosen = np.argmax(info_vals)
        used[chosen] = True
        items_given.append(chosen)

        p = p_2pl(theta_true, pool_a[chosen], pool_b[chosen])
        resp = np.random.binomial(1, p)
        responses.append(resp)

        theta_hat = eap_score(np.array(responses),
                              select_a[items_given], select_b[items_given], GRID)
        if return_trajectory:
            theta_score_k = eap_score(np.array(responses),
                                      score_a[items_given], score_b[items_given], GRID)
            traj.append((theta_hat, theta_score_k))

    theta_final = eap_score(np.array(responses),
                            score_a[items_given], score_b[items_given], GRID)
    return theta_final, items_given, responses, traj


def run_fixed(theta_true, pool_a, pool_b, items, score_a, score_b):
    responses = np.random.binomial(1, p_2pl(theta_true, pool_a[items], pool_b[items]))
    theta = eap_score(responses, score_a[items], score_b[items], GRID)
    return theta


# ── Misalignment Functions ────────────────────────────────────

def generate_pool(n=100):
    b = np.random.normal(0, 0.75, n)
    a = np.maximum(0.5, np.exp(np.random.normal(0.1, 0.3, n)))
    return a, b

def add_correlated(b, rho):
    mu, sd = b.mean(), b.std()
    residual = sd * np.sqrt(max(1 - rho**2, 0))
    return mu + rho * (b - mu) + np.random.normal(0, residual, len(b))

def add_noise(b, sigma):
    if sigma == 0:
        return b.copy()
    return b + np.random.normal(0, sigma, len(b))

def add_shift_noise(b, delta, sigma):
    return b + delta + np.random.normal(0, sigma, len(b))


# ── Simulation Runner ─────────────────────────────────────────

def simulate_cell(rho=None, sigma=None, delta=None, n_exam=50, n_reps=10,
                  cat_length=20, return_stepwise=False):
    """Run one simulation cell. Returns dict of metrics per condition."""
    np.random.seed(42)
    pool_a, pool_b = generate_pool(100)

    # Create misaligned bank
    if rho is not None:
        noisy_b = add_correlated(pool_b, rho)
    elif delta is not None:
        noisy_b = add_shift_noise(pool_b, delta, sigma or 0.3)
    elif sigma is not None:
        noisy_b = add_noise(pool_b, sigma)
    else:
        noisy_b = pool_b.copy()

    noisy_a = pool_a.copy()  # a held stable
    fixed_items = np.argsort(pool_a)[::-1][:cat_length]

    results = {c: [] for c in ["CAT-O", "CAT-R", "CAT-B", "Fixed-20", "Full-100"]}
    stepwise = {c: [[] for _ in range(cat_length)] for c in ["CAT-O", "CAT-R", "CAT-B"]}

    for rep in range(n_reps):
        thetas = np.random.normal(0, 1, n_exam)
        for th in thetas:
            # CAT-O
            th_o, items_o, resp_o, traj_o = run_cat(
                th, pool_a, pool_b, noisy_a, noisy_b, noisy_a, noisy_b,
                cat_length, return_trajectory=return_stepwise)
            results["CAT-O"].append((th, th_o))

            # CAT-R: rescore CAT-O items with true bank
            th_r = eap_score(np.array(resp_o),
                             pool_a[items_o], pool_b[items_o], GRID)
            results["CAT-R"].append((th, th_r))

            # CAT-B
            th_b, _, _, traj_b = run_cat(
                th, pool_a, pool_b, pool_a, pool_b, pool_a, pool_b,
                cat_length, return_trajectory=return_stepwise)
            results["CAT-B"].append((th, th_b))

            # Fixed-20
            th_fix = run_fixed(th, pool_a, pool_b, fixed_items, pool_a, pool_b)
            results["Fixed-20"].append((th, th_fix))

            # Full-100
            th_full = run_fixed(th, pool_a, pool_b, np.arange(100), pool_a, pool_b)
            results["Full-100"].append((th, th_full))

            if return_stepwise and traj_o:
                for k in range(cat_length):
                    # CAT-O: scored with noisy bank (traj_select)
                    stepwise["CAT-O"][k].append((th, traj_o[k][0]))
                    # CAT-R: rescore CAT-O's first k+1 items with TRUE bank
                    items_k = items_o[:k+1]
                    resp_k  = np.array(resp_o[:k+1])
                    theta_r_k = eap_score(resp_k, pool_a[items_k],
                                          pool_b[items_k], GRID)
                    stepwise["CAT-R"][k].append((th, theta_r_k))
                    # CAT-B: scored with true bank
                    stepwise["CAT-B"][k].append((th, traj_b[k][1]))

    # Compute metrics
    metrics = {}
    for cond, pairs in results.items():
        pairs = np.array(pairs)
        errors = pairs[:, 1] - pairs[:, 0]
        mse = np.mean(errors**2)
        bias = np.mean(errors)
        metrics[cond] = {
            "rmse": np.sqrt(mse),
            "bias": bias,
            "cor": np.corrcoef(pairs[:, 0], pairs[:, 1])[0, 1],
            "mse": mse,
            "bias_sq": bias**2,
            "variance": mse - bias**2,
            "pct_bias": (bias**2 / mse * 100) if mse > 1e-10 else 0
        }

    step_metrics = None
    if return_stepwise:
        step_metrics = {}
        for cond in ["CAT-O", "CAT-R", "CAT-B"]:
            step_metrics[cond] = []
            for k in range(cat_length):
                pairs = np.array(stepwise[cond][k])
                if len(pairs) > 0:
                    errs = pairs[:, 1] - pairs[:, 0]
                    step_metrics[cond].append(np.sqrt(np.mean(errs**2)))
                else:
                    step_metrics[cond].append(np.nan)

    return metrics, step_metrics


# ── Plot Styles ───────────────────────────────────────────────

COLORS = {
    "CAT-O": "#E8604C", "CAT-R": "#D4A843", "CAT-B": "#2D5FA0",
    "Fixed-20": "#555555", "Full-100": "#C48FB5"
}
LINES = {
    "CAT-O": "-", "CAT-R": "-", "CAT-B": "-",
    "Fixed-20": "--", "Full-100": ":"
}

def style_ax(ax, xlabel, ylabel, title=""):
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    if title:
        ax.set_title(title, fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.2)


# ══════════════════════════════════════════════════════════════
# TAB 1: Cross-Population Misalignment
# ══════════════════════════════════════════════════════════════

def tab1_sweep(n_exam, n_reps):
    n_exam, n_reps = int(n_exam), int(n_reps)
    rho_levels = np.linspace(0.05, 0.95, 10)
    all_metrics = {c: {"rmse": [], "bias": [], "pct_bias": []}
                   for c in COLORS}

    for rho in rho_levels:
        m, _ = simulate_cell(rho=rho, n_exam=n_exam, n_reps=n_reps)
        for c in COLORS:
            all_metrics[c]["rmse"].append(m[c]["rmse"])
            all_metrics[c]["bias"].append(m[c]["bias"])
            all_metrics[c]["pct_bias"].append(m[c]["pct_bias"])

    # RMSE plot
    fig1, ax1 = plt.subplots(figsize=(7, 4))
    for c in COLORS:
        ax1.plot(rho_levels, all_metrics[c]["rmse"],
                 color=COLORS[c], linestyle=LINES[c], linewidth=2, label=c,
                 marker="o", markersize=4)
    style_ax(ax1, "Cross-population alignment (ρ)", "RMSE",
             "Rescoring recovers most accuracy")
    ax1.legend(fontsize=9, ncol=5, loc="upper right")
    fig1.tight_layout()

    # Scoring vs Selection error
    fig2, ax2 = plt.subplots(figsize=(7, 4))
    scoring = np.array(all_metrics["CAT-O"]["rmse"]) - np.array(all_metrics["CAT-R"]["rmse"])
    selection = np.array(all_metrics["CAT-R"]["rmse"]) - np.array(all_metrics["CAT-B"]["rmse"])
    ax2.fill_between(rho_levels, 0, scoring, alpha=0.4, color="#E8604C", label="Scoring error (fixable)")
    ax2.fill_between(rho_levels, 0, selection, alpha=0.4, color="#D4A843", label="Selection error (not fixable)")
    style_ax(ax2, "Cross-population alignment (ρ)", "RMSE gap",
             "Scoring error >> Selection error")
    ax2.legend(fontsize=10)
    fig2.tight_layout()

    # Bias plot
    fig3, ax3 = plt.subplots(figsize=(7, 4))
    for c in ["CAT-O", "CAT-R", "CAT-B"]:
        ax3.plot(rho_levels, all_metrics[c]["bias"],
                 color=COLORS[c], linewidth=2, label=c, marker="o", markersize=4)
    ax3.axhline(0, color="grey", linestyle="--", alpha=0.5)
    style_ax(ax3, "Cross-population alignment (ρ)", "Mean Bias",
             "Rescoring removes most bias (compare CAT-O with CAT-R)")
    ax3.legend(fontsize=10)
    fig3.tight_layout()

    # Summary text
    avg_pct = np.mean(all_metrics["CAT-O"]["pct_bias"])
    avg_pct_r = np.mean(all_metrics["CAT-R"]["pct_bias"])
    avg_scoring_pct = np.mean(scoring / (scoring + selection + 1e-10) * 100)
    summary = (
        f"**Key findings:**\n"
        f"- Scoring error accounts for ~{avg_scoring_pct:.0f}% of miscalibration-attributable error (fixable by rescoring)\n"
        f"- Bias as % of MSE, averaged across ρ: {avg_pct:.2f}% for operational scores (CAT-O), {avg_pct_r:.2f}% after rescoring (CAT-R)\n"
        f"- CAT-R stays close to CAT-B across all ρ levels\n\n"
        f"*Simulation: {n_exam} examinees × {n_reps} reps per ρ level. Small runs are noisy; the paper's estimates "
        f"use 100 examinees × 50 replications per condition, so raise the sliders to approach them.*"
    )

    return fig1, fig2, fig3, summary


# ══════════════════════════════════════════════════════════════
# TAB 2: Random Noise & Structured Shifts
# ══════════════════════════════════════════════════════════════

def tab2_noise(delta, n_exam, n_reps):
    n_exam, n_reps = int(n_exam), int(n_reps)
    delta = float(delta)
    sigma_levels = np.linspace(0, 1.0, 10)
    all_metrics = {c: {"rmse": [], "bias": []} for c in COLORS}

    for sig in sigma_levels:
        if delta == 0:
            m, _ = simulate_cell(sigma=sig, n_exam=n_exam, n_reps=n_reps)
        else:
            m, _ = simulate_cell(delta=delta, sigma=sig, n_exam=n_exam, n_reps=n_reps)
        for c in COLORS:
            all_metrics[c]["rmse"].append(m[c]["rmse"])
            all_metrics[c]["bias"].append(m[c]["bias"])

    # RMSE
    fig1, ax1 = plt.subplots(figsize=(7, 4))
    for c in COLORS:
        ax1.plot(sigma_levels, all_metrics[c]["rmse"],
                 color=COLORS[c], linestyle=LINES[c], linewidth=2, label=c,
                 marker="o", markersize=4)
    title = f"RMSE by calibration noise (Δb = {delta:+.1f})" if delta != 0 else "RMSE by calibration noise (random only)"
    style_ax(ax1, "Calibration noise (σ)", "RMSE", title)
    ax1.legend(fontsize=9, ncol=5, loc="upper left")
    fig1.tight_layout()

    # Bias
    fig2, ax2 = plt.subplots(figsize=(7, 4))
    for c in ["CAT-O", "CAT-R", "CAT-B"]:
        ax2.plot(sigma_levels, all_metrics[c]["bias"],
                 color=COLORS[c], linewidth=2, label=c, marker="o", markersize=4)
    ax2.axhline(0, color="grey", linestyle="--", alpha=0.5)
    title2 = f"Bias by calibration noise (Δb = {delta:+.1f})" if delta != 0 else "Bias by calibration noise"
    style_ax(ax2, "Calibration noise (σ)", "Mean Bias", title2)
    ax2.legend(fontsize=10)
    fig2.tight_layout()

    gap_at_max = all_metrics["CAT-O"]["rmse"][-1] - all_metrics["CAT-R"]["rmse"][-1]
    summary = (
        f"**At σ = 1.0:** CAT-R improves over CAT-O by {gap_at_max:.3f} RMSE\n"
        f"{'Structured shift Δb = ' + str(delta) + ' adds systematic offset but rescoring still helps' if delta != 0 else 'Pure random noise: rescoring consistently reduces error'}\n\n"
        f"*Simulation: {n_exam} examinees × {n_reps} reps per σ level*"
    )

    return fig1, fig2, summary


# ══════════════════════════════════════════════════════════════
# TAB 3: Stepwise Trajectory
# ══════════════════════════════════════════════════════════════

def tab3_stepwise(rho, cat_length, n_exam, n_reps):
    n_exam, n_reps = int(n_exam), int(n_reps)
    rho = float(rho)
    cat_length = int(cat_length)

    m, step_m = simulate_cell(rho=rho, n_exam=n_exam, n_reps=n_reps,
                              cat_length=cat_length, return_stepwise=True)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    steps = np.arange(1, cat_length + 1)

    for c in ["CAT-O", "CAT-R", "CAT-B"]:
        ax.plot(steps, step_m[c], color=COLORS[c], linewidth=2.5,
                label=c, marker="o", markersize=4)

    # Full-100 reference line
    ax.axhline(m["Full-100"]["rmse"], color=COLORS["Full-100"],
               linestyle=":", linewidth=1.5, label="Full-100 ceiling")

    style_ax(ax, "Number of items administered", "RMSE",
             f"Stepwise accuracy (ρ = {rho:.2f}, {cat_length} items): short tests benefit most")
    ax.legend(fontsize=10)
    fig.tight_layout()

    mid = min(4, cat_length - 1)  # ~5 items or max available
    gap_mid = step_m["CAT-O"][mid] - step_m["CAT-R"][mid]
    gap_end = step_m["CAT-O"][-1] - step_m["CAT-R"][-1]
    summary = (
        f"**Rescoring advantage:**\n"
        f"- At {mid+1} items: CAT-R improves by {gap_mid:.3f} RMSE over CAT-O\n"
        f"- At {cat_length} items: CAT-R improves by {gap_end:.3f} RMSE over CAT-O\n"
        f"- Short tests benefit {gap_mid/max(gap_end, 0.001):.1f}× more from rescoring\n\n"
        f"This explains why younger children (who get shorter CATs) are most affected by miscalibration.\n\n"
        f"*Simulation: {n_exam} examinees × {n_reps} reps at ρ = {rho:.2f}, test length = {cat_length}*"
    )

    return fig, summary


# ══════════════════════════════════════════════════════════════
# TAB 4: Single-Point Explorer
# ══════════════════════════════════════════════════════════════

def tab4_explore(rho, cat_length, n_exam, n_reps):
    n_exam, n_reps = int(n_exam), int(n_reps)
    rho, cat_length = float(rho), int(cat_length)

    m, _ = simulate_cell(rho=rho, n_exam=n_exam, n_reps=n_reps,
                         cat_length=cat_length)

    # Bar chart
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    conds = ["CAT-O", "CAT-R", "CAT-B", "Fixed-20", "Full-100"]
    colors = [COLORS[c] for c in conds]

    # RMSE
    rmse_vals = [m[c]["rmse"] for c in conds]
    axes[0].bar(conds, rmse_vals, color=colors, alpha=0.85)
    style_ax(axes[0], "", "RMSE", "Accuracy")
    axes[0].tick_params(axis="x", rotation=30)

    # Bias
    bias_vals = [m[c]["bias"] for c in conds]
    axes[1].bar(conds, bias_vals, color=colors, alpha=0.85)
    axes[1].axhline(0, color="grey", linestyle="--", alpha=0.5)
    style_ax(axes[1], "", "Bias", "Systematic shift")
    axes[1].tick_params(axis="x", rotation=30)

    # Correlation
    cor_vals = [m[c]["cor"] for c in conds]
    axes[2].bar(conds, cor_vals, color=colors, alpha=0.85)
    style_ax(axes[2], "", "Correlation with θ_true", "Recovery")
    axes[2].tick_params(axis="x", rotation=30)
    axes[2].set_ylim(0, 1.05)

    fig.suptitle(f"ρ = {rho:.2f}, {cat_length} items", fontsize=13, fontweight="bold")
    fig.tight_layout()

    # Error decomposition
    total = m["CAT-O"]["rmse"] - m["CAT-B"]["rmse"]
    scoring = m["CAT-O"]["rmse"] - m["CAT-R"]["rmse"]
    selection = m["CAT-R"]["rmse"] - m["CAT-B"]["rmse"]
    pct_s = scoring / max(total, 0.001) * 100

    summary = (
        f"| Condition | RMSE | Bias | Correlation |\n"
        f"|-----------|------|------|-------------|\n"
    )
    for c in conds:
        summary += f"| {c} | {m[c]['rmse']:.3f} | {m[c]['bias']:+.4f} | {m[c]['cor']:.3f} |\n"

    summary += (
        f"\n**Error decomposition:**\n"
        f"- Miscalibration-attributable error: {total:.3f}\n"
        f"- Scoring error (fixable): {scoring:.3f} ({pct_s:.0f}%)\n"
        f"- Selection error (not fixable): {selection:.3f} ({100-pct_s:.0f}%)\n"
        f"- Bias as % of MSE: {m['CAT-O']['pct_bias']:.2f}% for operational scores (CAT-O), "
        f"{m['CAT-R']['pct_bias']:.2f}% after rescoring (CAT-R)\n\n"
        f"*Single-point runs are noisy at small sample sizes; raise the sliders for stable estimates.*\n"
    )

    return fig, summary


# ══════════════════════════════════════════════════════════════
# GRADIO APP
# ══════════════════════════════════════════════════════════════

DESCRIPTION = """
# 🎯 CAT Rescoring Explorer

**Interactive companion to:** *"Recovering Accuracy From Miscalibrated Item Banks: Scoring Error Is Fixable, Selection Error Is Not"*
Xiao, Domingue, Ram, & Frank (2026), *Educational and Psychological Measurement* · [https://doi.org/10.1177/00131644261486854](https://doi.org/10.1177/00131644261486854)

Explore how rescoring — recomputing ability estimates with improved item parameters — recovers accuracy from miscalibrated computerized adaptive tests (CATs).

**Five conditions:** CAT-O (operational, wrong ruler) → CAT-R (rescored, fixed ruler) → CAT-B (benchmark, never wrong) | Fixed-20 | Full-100

⚡ *Simulations run in real-time. Increase examinees/reps for smoother curves (slower).*
"""

with gr.Blocks(
    title="CAT Rescoring Explorer",
    theme=gr.themes.Soft(
        primary_hue="amber",
        secondary_hue="blue",
        font=["Source Sans Pro", "sans-serif"]
    )
) as app:

    gr.Markdown(DESCRIPTION)

    with gr.Tab("📊 Cross-Population Misalignment"):
        gr.Markdown(
            "**Paper claim:** Scoring error dominates selection error across all ρ levels. "
            "After rescoring, bias is below 0.5% of MSE. Rescoring recovers most miscalibration-attributable error."
        )
        with gr.Row():
            t1_exam = gr.Slider(10, 100, value=30, step=10, label="Examinees per cell")
            t1_reps = gr.Slider(2, 20, value=5, step=1, label="Replications")
        t1_btn = gr.Button("Run sweep across ρ", variant="primary")
        with gr.Row():
            t1_fig1 = gr.Plot(label="RMSE by ρ")
            t1_fig2 = gr.Plot(label="Scoring vs Selection Error")
        t1_fig3 = gr.Plot(label="Bias by ρ")
        t1_summary = gr.Markdown()
        t1_btn.click(tab1_sweep, [t1_exam, t1_reps],
                     [t1_fig1, t1_fig2, t1_fig3, t1_summary])

    with gr.Tab("🔊 Random Noise & Shifts"):
        gr.Markdown(
            "**Paper claim:** Rescoring is robust to both random calibration noise "
            "and structured difficulty shifts (DIF-like offsets)."
        )
        with gr.Row():
            t2_delta = gr.Radio([0, -0.3, 0.3, 1.0], value=0, label="Δb (difficulty shift)")
            t2_exam = gr.Slider(10, 100, value=30, step=10, label="Examinees")
            t2_reps = gr.Slider(2, 20, value=5, step=1, label="Replications")
        t2_btn = gr.Button("Run sweep across σ", variant="primary")
        with gr.Row():
            t2_fig1 = gr.Plot(label="RMSE by σ")
            t2_fig2 = gr.Plot(label="Bias by σ")
        t2_summary = gr.Markdown()
        t2_btn.click(tab2_noise, [t2_delta, t2_exam, t2_reps],
                     [t2_fig1, t2_fig2, t2_summary])

    with gr.Tab("📈 Stepwise: Short Tests Benefit Most"):
        gr.Markdown(
            "**Paper claim:** The rescoring advantage is largest for short tests (5–15 items), "
            "explaining why younger children are most affected by miscalibration."
        )
        with gr.Row():
            t3_rho = gr.Slider(0.1, 0.95, value=0.40, step=0.05, label="ρ")
            t3_len = gr.Slider(5, 50, value=20, step=5, label="Test length (items)")
            t3_exam = gr.Slider(10, 100, value=30, step=10, label="Examinees")
            t3_reps = gr.Slider(2, 20, value=5, step=1, label="Replications")
        t3_btn = gr.Button("Run stepwise simulation", variant="primary")
        t3_fig = gr.Plot(label="RMSE by number of items")
        t3_summary = gr.Markdown()
        t3_btn.click(tab3_stepwise, [t3_rho, t3_len, t3_exam, t3_reps],
                     [t3_fig, t3_summary])

    with gr.Tab("🔍 Single-Point Explorer"):
        gr.Markdown(
            "**Explore any scenario.** Set ρ and test length, see all metrics at once."
        )
        with gr.Row():
            t4_rho = gr.Slider(0.05, 0.99, value=0.5, step=0.05, label="ρ")
            t4_len = gr.Slider(5, 50, value=20, step=5, label="CAT length")
            t4_exam = gr.Slider(10, 100, value=30, step=10, label="Examinees")
            t4_reps = gr.Slider(2, 20, value=5, step=1, label="Replications")
        t4_btn = gr.Button("Simulate", variant="primary")
        t4_fig = gr.Plot(label="Metrics comparison")
        t4_summary = gr.Markdown()
        t4_btn.click(tab4_explore, [t4_rho, t4_len, t4_exam, t4_reps],
                     [t4_fig, t4_summary])

    gr.Markdown(
        "---\n"
        "**Paper:** Xiao, X., Domingue, B. W., Ram, N., & Frank, M. C. (2026). "
        "Recovering accuracy from miscalibrated item banks: Scoring error is fixable, selection error is not. "
        "*Educational and Psychological Measurement*. [https://doi.org/10.1177/00131644261486854](https://doi.org/10.1177/00131644261486854)\n\n"
        "**Code:** [github.com/DoriaXiao/rescoring-cat](https://github.com/DoriaXiao/rescoring-cat) | "
        "**LEVANTE:** [levante-network.org](https://levante-network.org)"
    )

if __name__ == "__main__":
    app.launch()
