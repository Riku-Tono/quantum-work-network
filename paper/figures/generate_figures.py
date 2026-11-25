"""Render the manuscript's vector figures from archived milestone CSVs."""

import csv
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


FIGURES = Path(__file__).resolve().parent
ROOT = FIGURES.parents[1]
INK = "#253244"
BLUE = "#17668A"
ORANGE = "#BC4A2F"
TEAL = "#007C70"
PURPLE = "#76519B"
GRAY = "#6E7884"
LIGHT = "#DCE3E9"
WORK_UNITS = r"($10^{-6}$ model units)"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "mathtext.fontset": "dejavusans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.65,
    "axes.edgecolor": GRAY,
    "axes.labelcolor": INK,
    "axes.titlesize": 9.5,
    "axes.titleweight": "semibold",
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "grid.color": LIGHT,
    "grid.linewidth": 0.55,
    "legend.fontsize": 7.4,
    "legend.frameon": False,
    "text.color": INK,
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})


def records(path, **filters):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as source:
        return [row for row in csv.DictReader(source)
                if all(row[key] == value for key, value in filters.items())]


def panel(ax, letter, title):
    ax.set_title(f"({letter})  {title}", loc="left", color=INK, pad=8)
    ax.set_axisbelow(True)


def save(fig, filename):
    fig.savefig(FIGURES / filename, metadata={"CreationDate": None})
    plt.close(fig)


def predictive_stopping():
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.8, 2.7))
    fig.subplots_adjust(left=.103, right=.978, bottom=.20, top=.89, wspace=.40)

    reference = records("Milestone_24/results/m24_precursors.csv",
                        signal="ergotropy_forecast")
    peak = [float(row["peak_elapsed_total"]) for row in reference]
    first_error = sorted(abs(float(row["first_target_elapsed_total"]) - time)
                         for row, time in zip(reference, peak, strict=True))
    matched_error = sorted(abs(float(row["target_error"])) for row in reference)
    panel(ax, "a", "Peak-timing error (M24)")
    for errors, label, color, linestyle in (
        (first_error, "First forecast", BLUE, "-"),
        (matched_error, "Hindsight-selected forecast", GRAY, "--"),
    ):
        ax.step([0, *errors, 18], [0, *range(1, len(errors) + 1), len(errors)],
                where="post", color=color, label=label, linestyle=linestyle, lw=1.5)
    ax.axvline(.10, color=LIGHT, lw=1.0)
    hits = sum(error <= .10 + 1e-12 for error in first_error)
    ax.scatter(.10, hits, color=BLUE, zorder=4, s=15)
    ax.annotate(f"{hits}/{len(reference)} within 0.10", (.10, hits), xytext=(12, 7),
                textcoords="offset points", color=BLUE, fontsize=7.6)
    ax.set(xlabel="Absolute error (model time)", ylabel="Trajectories within error",
           xlim=(0, 18), ylim=(0, 36), yticks=[0, 10, 20, 30])
    ax.set_xscale("symlog", linthresh=.10)
    ax.set_xticks([0, .10, 1, 10], ["0", "0.10", "1", "10"])
    ax.grid(axis="y")
    ax.legend(loc="lower right", borderaxespad=.35)

    panel(bx, "b", "Repeated work (M26)")
    history = records("Milestone_26/results/m26_time_comparison.csv")
    clocks = sorted({float(row["elapsed_total"]) for row in history})
    for policy, label, color, marker, linestyle in (
        ("A", "A: predictive", BLUE, "o", "-"),
        ("E", "E: hindsight", ORANGE, "s", "--"),
    ):
        sample = [[1e6 * float(row["cumulative_net"]) for row in history
                   if row["policy"] == policy and float(row["elapsed_total"]) == time]
                  for time in clocks]
        bx.fill_between(clocks, [min(values) for values in sample],
                        [max(values) for values in sample], color=color, alpha=.12, lw=0)
        bx.plot(clocks, [mean(values) for values in sample], color=color, marker=marker,
                ms=3.2, lw=1.5, linestyle=linestyle, label=label)
    bx.set(xlabel="Elapsed model time", ylabel=f"Completed net work {WORK_UNITS}",
           xlim=(39, 552), ylim=(10, 110), xticks=[50, 150, 300, 450, 540],
           yticks=[20, 40, 60, 80, 100])
    bx.grid(axis="y")
    bx.legend(loc="upper left", borderaxespad=.30)
    save(fig, "predictive_stopping.pdf")


def adaptive_reconnection():
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.8, 2.8))
    fig.subplots_adjust(left=.105, right=.977, bottom=.22, top=.89, wspace=.39)

    panel(ax, "a", "Matched-state phase rule (M28)")
    endpoints = records("Milestone_28/results/m28_endpoints.csv", group="F2-R",
                        precision="coarse")
    work = {(row["origin"], row["policy"]): 1e6 * float(row["cumulative_net"])
            for row in endpoints}
    origins = ("S0", "A04", "B04", "C04")
    ax.axhline(0, color=GRAY, lw=.75)
    for policy, label, offset, color, marker in (
        ("R1", "Fixed wait 1", -.18, BLUE, "o"),
        ("R3", "Fixed wait 3", 0, PURPLE, "s"),
        ("RH", "History-based", .18, TEAL, "D"),
    ):
        ax.plot([index + offset for index in range(len(origins))],
                [work[origin, "RP"] - work[origin, policy] for origin in origins],
                color=color, marker=marker, ms=4.4, linestyle="none", label=label)
    ax.set(ylabel=f"RP work difference {WORK_UNITS}", xlabel="F2-R anchor origin",
           xlim=(-.40, 3.40), ylim=(-12.5, 16), yticks=[-10, -5, 0, 5, 10, 15])
    ax.set_xticks(range(len(origins)), origins)
    ax.grid(axis="y")
    ax.legend(loc="lower right", borderaxespad=.28, handletextpad=.3, labelspacing=.32)

    panel(bx, "b", "Current-only policy P3 (M29)")
    retention = records("Milestone_29/results/m29_maintenance.csv", precision="coarse",
                        policy="P3")
    rows = records("Milestone_29/results/m29_endpoints.csv", precision="coarse")
    offsets = {"F1": (7, 5), "F2-R": (7, -11), "unseen": (7, 1)}
    for row in retention:
        group = row["group"]
        burden = {policy: mean(int(item["primitive_coherence_component_values"])
                               for item in rows
                               if item["group"] == group and item["policy"] == policy)
                  for policy in ("P3", "RP")}
        x = 100 * burden["P3"] / burden["RP"]
        y = 100 * float(row["ratio_of_means"])
        label = "Unused preparations" if group == "unseen" else f"{group} anchors"
        color = ORANGE if group == "unseen" else BLUE
        marker = "s" if group == "unseen" else "o"
        bx.plot(x, y, marker=marker, color=color, ms=5, zorder=3)
        bx.annotate(label, (x, y), xytext=offsets[group], textcoords="offset points",
                    fontsize=7.7, color=color)
    bx.axhline(95, color=GRAY, linestyle="--", lw=.9)
    bx.annotate("95% mean criterion", (104, 95), xytext=(0, -4), ha="right", va="top",
                textcoords="offset points", fontsize=7.5, color=GRAY)
    bx.plot(100, 100, marker="D", mfc="white", mec=INK, ms=4, zorder=3)
    bx.annotate("RP reference", (100, 100), xytext=(-6, 4), ha="right",
                textcoords="offset points", fontsize=7.7)
    bx.set(xlabel="Acquired real components (% of RP)",
           ylabel="Mean net work (% of RP)", xlim=(35, 106), ylim=(74, 110),
           xticks=[40, 60, 80, 100], yticks=[75, 85, 95, 105])
    bx.grid(axis="y")
    save(fig, "adaptive_reconnection.pdf")


def u05_reconnection():
    fig, (ax, bx) = plt.subplots(1, 2, sharey=True, figsize=(6.8, 3.25))
    fig.subplots_adjust(left=.14, right=.966, bottom=.16, top=.90, wspace=.19)
    trial = records("Milestone_30/u05/trial6_summary.csv", propagation="rk4",
                    precision="coarse")
    endpoint = records("Milestone_30/results/m30_endpoints.csv", propagation="rk4",
                       precision="coarse", origin="U05", elapsed="540")
    rows = ("P3", "P3+I", "P3+A", "P3+IA", "P1", "P1+I", "P1+A", "RP")
    colors = (GRAY, BLUE, PURPLE, TEAL, GRAY, BLUE, PURPLE, INK)
    markers = ("o", "o", "o", "o", "s", "s", "s", "D")

    for plot, data, letter, title in (
        (ax, trial, "a", "Common-state trial 6"),
        (bx, endpoint, "b", "Elapsed 540"),
    ):
        panel(plot, letter, title)
        work = {row["policy"]: float(row["net_work"]) * 1e6 for row in data}
        plot.axhline(3.5, color=LIGHT, lw=.9)
        plot.axvline(0, color=GRAY, lw=.8)
        for y, (policy, color, marker) in enumerate(zip(rows, colors, markers, strict=True)):
            plot.hlines(y, 0, work[policy], color=color, alpha=.66, lw=1.2)
            plot.plot(work[policy], y, marker=marker, color=color, ms=4.6)
        plot.set(xlabel=f"Net work {WORK_UNITS}", ylim=(7.6, -.6))
        plot.yaxis.set_tick_params(length=0)
    ax.set_yticks(range(len(rows)), rows)
    ax.set(xlim=(-375, 300), xticks=[-300, 0, 225])
    bx.set(xlim=(-1300, 2090), xticks=[-1000, 0, 1000, 2000])
    ax.annotate(r"$\tau=2.4$", (-298.68, 0), xytext=(5, 4), color=GRAY,
                textcoords="offset points", fontsize=7.5)
    ax.annotate(r"$\tau=6$ cap", (224.59, 1), xytext=(-6, 4), color=BLUE,
                textcoords="offset points", ha="right", fontsize=7.5)
    save(fig, "u05_reconnection.pdf")


if __name__ == "__main__":
    predictive_stopping()
    adaptive_reconnection()
    u05_reconnection()
