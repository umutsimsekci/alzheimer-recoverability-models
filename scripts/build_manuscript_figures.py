"""Render manuscript figures from completed study CSVs; no simulation or fitting."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
FAMILIES=(("M0",10.),("M1",10.),("M2",10.),("M3",1.),("M3",10.),("M3",50.))
COLORS={"none":"#777777","drive":"#0072B2","repair":"#D55E00","both":"#009E73"}
SPEED_COLORS={60:"#0072B2",240:"#D55E00",960:"#009E73"}
LABELS={"none":"No controls","drive":"Drive reduction","repair":"Repair pulse","both":"Both controls"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))


def family_label(model,tau):
    return model if model!="M3" else f"M3 (τ={tau:g})"


def save(fig,out,name):
    for suffix in ("png","svg","pdf"):
        fig.savefig(out/f"{name}.{suffix}",dpi=300,bbox_inches="tight",facecolor="white")
    plt.close(fig)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir",type=Path,default=ROOT/"results/study_full_v1")
    parser.add_argument("--output-dir",type=Path,default=ROOT/"results/manuscript_figures_v1")
    args=parser.parse_args(argv)
    src=args.source_dir.resolve();out=args.output_dir.resolve()
    summary=json.loads((src/"summary.json").read_text(encoding="utf-8"))
    if summary.get("evidence_class")!="synthetic_simulation" or summary.get("failures") or summary.get("status")!="generated_partial_protocol_unreviewed":
        raise ValueError("Require completed full synthetic A–C outputs without failures")
    if {j["stage"] for j in summary["jobs"] if j["status"]=="generated_unreviewed"}!={"A","B","C"}:
        raise ValueError("A/B/C completion records missing")
    declared=json.loads((src/"output_hashes.json").read_text(encoding="utf-8"))
    used=[src/"summary.json",src/"M1_tau10_nominal_trajectories.csv",src/"primary_outcomes.csv",src/"two_history_endpoints.csv"]
    used += [src/f"ramp_{model}_tau{tau:g}_leg{duration}.csv" for model,tau in (("M1",10.),("M3",50.)) for duration in (60,240,960)]
    for p in used:
        if declared.get(p.name)!=sha(p):
            raise ValueError(f"Source hash mismatch: {p.name}")
    primary=read_csv(src/"primary_outcomes.csv")
    if len(primary)!=9252:
        raise ValueError(f"Expected 9252 full primary records, found {len(primary)}")
    out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                         "svg.fonttype":"none","pdf.fonttype":42,"axes.labelsize":11,"axes.titlesize":12,"legend.fontsize":9})
    trajectories=read_csv(src/"M1_tau10_nominal_trajectories.csv")
    fig,axes=plt.subplots(1,2,figsize=(11,4),sharex=True,sharey=True)
    for ax,policy,title in zip(axes,("maintenance","full_withdrawal"),("A  Continued drive reduction","B  All controls off at t=45")):
        for arm in ("none","drive","repair","both"):
            rows=[r for r in trajectories if r["arm"]==arm and (r["policy"]==policy or r["policy"]=="identical_policies")]
            rows.sort(key=lambda r:float(r["t"]))
            if len(rows)!=481:
                raise ValueError(f"Nominal trajectory expected481: {arm}/{policy}")
            ax.plot([float(r["t"]) for r in rows],[float(r["x"]) for r in rows],color=COLORS[arm],lw=2 if arm=="both" else 1.5,
                    linestyle="--" if arm=="repair" else "-",label=LABELS[arm])
        ax.axhline(.15,color="#CC79A7",ls=":",lw=1.6,label="Arbitrary target x≤.15")
        ax.axvspan(30,45,color="#DDDDDD",alpha=.45,zorder=0)
        ax.axvline(30,color="#555555",ls=":",lw=.8)
        if policy=="full_withdrawal":ax.axvline(45,color="#555555",ls=":",lw=.8)
        ax.set(title=title,xlabel="Dimensionless model time",xlim=(0,120),ylim=(0,1))
    axes[0].set_ylabel("Dimensionless impairment x")
    fig.legend(*axes[0].get_legend_handles_labels(),loc="lower center",ncol=3,bbox_to_anchor=(.5,-.12),frameon=False)
    fig.suptitle("Nominal M1: maintenance and full withdrawal are different final systems",fontsize=13,y=1.04)
    fig.text(.5,-.18,"Synthetic trajectories; shading marks the repair-pulse interval t=30–45.",ha="center",fontsize=9)
    fig.tight_layout();save(fig,out,"figure_1_nominal_policies")

    counts=[]
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),sharey=True)
    x=np.arange(len(FAMILIES));width=.36
    for ax,policy,title in zip(axes,("maintenance","full_withdrawal"),("A  Continued drive reduction","B  Full withdrawal at t=45")):
        for shift,arm in ((-width/2,"drive"),(width/2,"both")):
            fractions=[]
            for model,tau in FAMILIES:
                rows=[r for r in primary if r["setting_id"]!="nominal" and r["model"]==model and float(r["tau"])==tau and r["arm"]==arm and r["policy"]==policy]
                if len(rows)!=256 or len({r["setting_id"] for r in rows})!=256:
                    raise ValueError(f"Invalid256 denominator: {model}/{tau}/{arm}/{policy}")
                successes=sum(r["primary_target"]=="True" for r in rows)
                fractions.append(successes/256)
                counts.append(dict(model=model,tau=tau,policy=policy,arm=arm,attainment_count=successes,denominator=256,fraction=successes/256))
            bars=ax.bar(x+shift,fractions,width,color=COLORS[arm],label=LABELS[arm])
            for bar,fraction in zip(bars,fractions):
                ax.text(bar.get_x()+bar.get_width()/2,max(.02,fraction+.025),f"{fraction*256:.0f}/256",ha="center",va="bottom",rotation=90,fontsize=7)
        ax.set_xticks(x,[family_label(m,t) for m,t in FAMILIES],rotation=25,ha="right")
        ax.set(title=title,ylim=(0,1.25));ax.grid(axis="y",alpha=.2);ax.set_axisbelow(True)
        ax.set_yticks(np.linspace(0,1,6))
    axes[0].set_ylabel("Fraction meeting sampled x≤.15 target\nthroughout t=90–120")
    fig.legend(*axes[1].get_legend_handles_labels(),loc="lower center",ncol=2,bbox_to_anchor=(.5,-.10),frameon=False)
    fig.suptitle("Target attainment in the arbitrary uniform parameter design",fontsize=13,y=1.02)
    fig.text(.5,-.16,"256 synthetic settings per bar; nominal excluded. M0/M3 ignore b. These are not patient probabilities.",ha="center",fontsize=9)
    fig.tight_layout();save(fig,out,"figure_2_design_target_fractions")

    history=read_csv(src/"two_history_endpoints.csv")
    if len(history)!=12:raise ValueError("Expected12 two-history records")
    fig=plt.figure(figsize=(13,4.8));grid=fig.add_gridspec(1,3,width_ratios=(1.1,1,1))
    left=fig.add_subplot(grid[0]);rampaxes=[fig.add_subplot(grid[1]),fig.add_subplot(grid[2])]
    for index,(model,tau) in enumerate(FAMILIES):
        for x0,color,offset,label in ((.04,"#0072B2",-.1,"Initial x=.04"),(.8,"#D55E00",.1,"Initial x=.8")):
            rows=[r for r in history if r["model"]==model and float(r["tau"])==tau and float(r["x0"])==x0]
            if len(rows)!=1:raise ValueError("Two-history key mismatch")
            r=rows[0];settled=r["convergence"]=="converged"
            left.scatter(index+offset,float(r["x_endpoint"]),color=color,facecolors=color if settled else "none",s=55,marker="o" if x0==.04 else "s",label=label if index==0 else None,zorder=3)
    left.set_xticks(np.arange(6),[family_label(m,t) for m,t in FAMILIES],rotation=30,ha="right")
    left.set(title="A  Same final drive a=.06",ylabel="Impairment at t=480",ylim=(0,1))
    left.legend(loc="upper left",frameon=False);left.grid(axis="y",alpha=.2)
    for ax,(model,tau),title in zip(rampaxes,(("M1",10.),("M3",50.)),("B  M1: feedback model","C  M3 (τ=50): slow-input model")):
        for duration in (60,240,960):
            rows=read_csv(src/f"ramp_{model}_tau{tau:g}_leg{duration}.csv")
            ax.plot([float(r["a_command"]) for r in rows],[float(r["x"]) for r in rows],color=SPEED_COLORS[duration],lw=1.5,label=f"Leg duration {duration}")
        ax.set(title=title,xlabel="Commanded drive a",xlim=(.015,.145),ylim=(0,1));ax.grid(alpha=.15)
        ax.legend(loc="upper left",frameon=False,fontsize=8)
    fig.suptitle("Settled histories and finite-rate loops are separate diagnostics",fontsize=13,y=1.02)
    fig.text(.5,-.08,"Synthetic dimensionless models. Hollow history markers denote unresolved convergence; loop width alone does not establish bistability.",ha="center",fontsize=9)
    fig.tight_layout();save(fig,out,"figure_3_history_and_rate_diagnostics")
    with (out/"figure_2_exact_counts.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(counts[0]));writer.writeheader();writer.writerows(counts)
    manifest=dict(evidence_class="figures_from_existing_synthetic_CSVs",new_simulation=False,new_fitting=False,
                  source_run=str(src),source_run_status=summary["status"],source_records_primary=len(primary),
                  per_fraction_denominator=256,nominal_excluded=True,
                  source_hashes={p.name:sha(p) for p in used},builder_sha256=sha(Path(__file__)),
                  matplotlib_version=matplotlib.__version__,created_unix=time.time(),
                  output_hashes={p.name:sha(p) for p in out.iterdir() if p.is_file()})
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print(json.dumps(dict(output_dir=str(out),figure_files=9,primary_records=len(primary),fractions=len(counts),source_hashes_verified=True)))


if __name__=="__main__":main()
