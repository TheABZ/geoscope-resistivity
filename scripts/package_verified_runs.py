"""Verify six dashboard exports, make comparison figures, and optionally a PDF.

python scripts/package_verified_runs.py --inputs /path/to/exports \
    --output results/verified_runs --pdf GeoScope_Project_Report.pdf

Survey labels below are supplied by the user and screenshots. Legacy exports
do not encode electrode spacings; labels are not independently inferred.
Requires numpy and matplotlib; optional PDF output also requires reportlab.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASES = [
    ("200m_4000", "200 m / 4,000 draws", "diagnostics.json", "posterior_samples.npz"),
    ("30m_4000", "30 m / 4,000 draws", "diagnostics (1).json", "posterior_samples (1).npz"),
    ("30m_8000", "30 m / 8,000 draws", "diagnostics (2).json", "posterior_samples (2).npz"),
]


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def collect(inputs, output):
    runs, hashes = [], []
    for key, label, diagnostic_name, sample_name in CASES:
        diagnostic = json.loads((inputs/diagnostic_name).read_text())
        with np.load(inputs/sample_name, allow_pickle=False) as archive:
            samples = archive["log_parameters"].copy()
            names = archive["parameter_names"].tolist()
        assert samples.shape == (diagnostic["chains"], diagnostic["draws_per_chain"], len(names))
        assert np.isfinite(samples).all()
        rows = {row["parameter"]: row for row in diagnostic["parameters"]}
        assert set(rows) == set(names)
        for i, name in enumerate(names):
            q = np.quantile(np.exp(samples[:, :, i]), [.025, .5, .975])
            np.testing.assert_allclose(q, [rows[name][k] for k in ["q025", "median", "q975"]], rtol=1e-10)
        passed = all(r["rhat"] <= 1.01 and r["ess_bulk"] >= 400 and r["ess_tail"] >= 400 for r in rows.values())
        assert passed == (diagnostic["status"] == "diagnostics_passed")
        runs.append({"key": key, "label": label, "diagnostic": diagnostic, "samples": samples, "names": names})
        (output/(key+"_diagnostics.json")).write_text(json.dumps(diagnostic, indent=2)+"\n")
        for name in (diagnostic_name, sample_name):
            hashes.append({"case": key, "file_name": name, "sha256": hashlib.sha256((inputs/name).read_bytes()).hexdigest()})
    prefix = np.array_equal(runs[1]["samples"], runs[2]["samples"][:, :4000, :])
    manifest = {"label_origin": "Survey extents supplied by user and screenshots; not encoded in legacy exports",
                "thirty_metre_4000_is_exact_prefix_of_8000": prefix,
                "comparison_design": "One baseline and one short-survey trajectory extended from 4000 to 8000 draws; not independent replications",
                "files": hashes}
    (output/"source_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    summaries = []
    for run in runs:
        d = run["diagnostic"]; rows=d["parameters"]
        summaries.append({"case": run["key"], "draws_per_chain": d["draws_per_chain"], "chains": d["chains"],
            "max_rhat": max(r["rhat"] for r in rows), "min_bulk_ess": min(r["ess_bulk"] for r in rows),
            "min_tail_ess": min(r["ess_tail"] for r in rows), "min_acceptance": min(d["acceptance_rates"]),
            "max_acceptance": max(d["acceptance_rates"]), "status": d["status"]})
    write_csv(output/"convergence_comparison.csv", summaries)
    write_csv(output/"baseline_parameters.csv", runs[0]["diagnostic"]["parameters"])
    return runs, summaries, manifest


def figures(runs, output):
    colors = ["#087e8b", "#cd8b21", "#c34f46"]
    with plt.rc_context({"font.family":"DejaVu Sans", "font.size":11, "svg.fonttype":"none", "svg.hashsalt":"geoscope-verified"}):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), layout="constrained")
        names = runs[0]["names"]
        x = np.arange(len(names)); width=.24
        for j, run in enumerate(runs):
            rows={r["parameter"]:r for r in run["diagnostic"]["parameters"]}
            axes[0].bar(x+(j-1)*width, [rows[n]["rhat"] for n in names], width, color=colors[j], label=run["label"])
            conservative=[min(rows[n]["ess_bulk"],rows[n]["ess_tail"]) for n in names]
            axes[1].bar(x+(j-1)*width, conservative, width, color=colors[j])
        axes[0].axhline(1.01,color="#303b45",ls="--",lw=1.5)
        axes[1].axhline(400,color="#303b45",ls="--",lw=1.5)
        axes[0].set(ylabel="Rank-normalized split R-hat", title="Between-chain agreement | target ≤1.01", ylim=(0,2.15))
        axes[1].set(yscale="log", ylabel="Minimum of bulk and tail ESS", title="Effective sample information | target ≥400", ylim=(3,2500))
        axes[0].legend(loc="upper left",fontsize=8)
        for ax in axes:
            ax.set_xticks(x,names); ax.grid(axis="y",alpha=.15); ax.set_axisbelow(True)
            ax.spines[["top","right"]].set_visible(False)
        fig.suptitle("GeoScope | Convergence checks change the conclusion",fontweight="bold",fontsize=17)
        for extension in ["png","svg"]:
            kwargs={"metadata":{"Date":None}} if extension=="svg" else {"dpi":180}
            fig.savefig(output/f"convergence_comparison.{extension}",**kwargs)
        plt.close(fig)
        fig, axes=plt.subplots(3,1,figsize=(12,7),layout="constrained")
        for ax,run in zip(axes,runs):
            samples=run["samples"]
            i=run["names"].index("rho_3")
            stride=max(1,samples.shape[1]//700)
            for chain in range(samples.shape[0]):
                ax.plot(np.arange(0,samples.shape[1],stride),np.exp(samples[chain,::stride,i]),lw=.8,alpha=.8,label=f"Chain {chain+1}")
            status="PASSED" if run["diagnostic"]["status"]=="diagnostics_passed" else "INCONCLUSIVE"
            ax.set(title=run["label"]+" | "+status, yscale="log", ylabel="rho_3 (ohm m)",ylim=(40,3500))
            ax.axhline(800,color="#48505a",ls="--",lw=1,alpha=.6)
            ax.grid(alpha=.15)
        axes[0].legend(ncol=4,fontsize=8,loc="lower right")
        axes[-1].set_xlabel("Retained iteration | plotting thinned; diagnostics use every retained draw")
        fig.suptitle("Deep-layer traces | dashed line = synthetic truth, 800 ohm m",fontsize=16,fontweight="bold")
        for extension in ["png","svg"]:
            kwargs={"metadata":{"Date":None}} if extension=="svg" else {"dpi":180}
            fig.savefig(output/f"deep_layer_traces.{extension}",**kwargs)
        plt.close(fig)


def pdf_report(path, runs, summaries):
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    fonts=Path("/usr/share/fonts/truetype/dejavu")
    if (fonts/"DejaVuSans.ttf").exists():
        pdfmetrics.registerFont(TTFont("GeoSans",str(fonts/"DejaVuSans.ttf")))
        pdfmetrics.registerFont(TTFont("GeoSans-Bold",str(fonts/"DejaVuSans-Bold.ttf")))
        pdfmetrics.registerFontFamily("GeoSans",normal="GeoSans",bold="GeoSans-Bold",italic="GeoSans",boldItalic="GeoSans-Bold")
        regular,bold="GeoSans","GeoSans-Bold"
    else:
        regular,bold="Helvetica","Helvetica-Bold"
    styles=getSampleStyleSheet()
    for key in ("BodyText","Normal"):
        styles[key].fontName=regular
    for key in ("Title","Heading1","Heading2"):
        styles[key].fontName=bold
    styles.add(ParagraphStyle(name="Deck",fontName=regular,fontSize=12,leading=17,textColor=colors.HexColor("#425564"),spaceAfter=14))
    styles["BodyText"].fontSize=10; styles["BodyText"].leading=14; styles["BodyText"].spaceAfter=9
    styles["Heading1"].textColor=colors.HexColor("#087e8b");styles["Heading1"].fontSize=19
    styles["Heading2"].textColor=colors.HexColor("#087e8b");styles["Heading2"].fontSize=13
    def p(text,style="BodyText"):return Paragraph(text,styles[style])
    def table(rows,widths):
        t=Table(rows,colWidths=widths,hAlign="LEFT")
        t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),colors.HexColor("#e9f3f4")),("TEXTCOLOR",(0,0),(-1,0),colors.HexColor("#234b52")),("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),("FONTNAME",(0,1),(-1,-1),"Helvetica"),("FONTSIZE",(0,0),(-1,-1),9),("TOPPADDING",(0,0),(-1,-1),8),("BOTTOMPADDING",(0,0),(-1,-1),8),("LINEBELOW",(0,0),(-1,0),.6,colors.HexColor("#b8ced1")),("VALIGN",(0,0),(-1,-1),"MIDDLE"),("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f6f8f9")])]))
        t.setStyle(TableStyle([("FONTNAME",(0,0),(-1,0),bold),("FONTNAME",(0,1),(-1,-1),regular),("FONTSIZE",(0,0),(-1,-1),8.5)]))
        return t
    path.parent.mkdir(parents=True,exist_ok=True)
    doc=SimpleDocTemplate(str(path),pagesize=A4,leftMargin=42,rightMargin=42,topMargin=42,bottomMargin=43,title="GeoScope: Verified Synthetic Results",author="Gabriel Chiazor Unomah")
    story=[p("GeoScope", "Title"),p("Verified synthetic results and sampling limitations", "Deck"),p("Gabriel Chiazor Unomah | 12 November 2026"),p("1. What this demonstration establishes", "Heading1"),p("GeoScope combines layered-earth electrical resistivity modelling, bounded inversion, and uncertainty diagnostics. Three dashboard exports were checked against their posterior samples. The 200 m baseline passed the declared sampling checks; the two 30 m exports did not. Failed-run quantiles are excluded from scientific interpretation."),p("The numerical experiment", "Heading2"),p("The intended synthetic Earth has resistivities of 100, 25, and 800 ohm m, with finite layer thicknesses of 5 and 20 m. The Wenner-alpha sounding uses a flat, horizontally layered, isotropic Earth. The reference setup uses 40 measurements, 3% log-noise, data seed 42, and sampler seed 43."),p("Audit result", "Heading2")]
    rows=[["Run", "Max R-hat", "Min bulk ESS", "Min tail ESS", "Decision"]]
    for run,s in zip(runs,summaries):
        rows.append([run["label"],f"{s['max_rhat']:.4f}",f"{s['min_bulk_ess']:.1f}",f"{s['min_tail_ess']:.1f}","Passed" if s["status"]=="diagnostics_passed" else "Inconclusive"])
    story += [table(rows,[145,70,90,90,116]),Spacer(1,12),p("Criteria: every parameter must have rank-normalized split R-hat at most 1.01 and both bulk and tail effective sample sizes (ESS) at least 400. All files contain four chains and 2,000 discarded warmup steps."),p("Important provenance detail", "Heading2"),p("The 4,000-draw 30 m sample array is an exact prefix of the 8,000-draw array for every chain and parameter. These are a trajectory and its extension, not two independent replications. Legacy exports do not encode electrode spacings; the 200 m and 30 m labels come from the user's stated settings and screenshots."),p("No field measurements or real-site groundwater claims are presented. This is a software demonstration, not a peer-reviewed validation study."),PageBreak()]
    story += [p("2. Passed baseline and failed sampling", "Heading1"),Image(str(OUTPUT/"convergence_comparison.png"),width=511,height=221),Spacer(1,9),p("Figure 1. The 200 m run satisfies the declared checks for all five parameters. Both 30 m runs fail. The right panel uses the smaller of bulk and tail ESS for each parameter; dashed lines mark the unchanged thresholds."),p("200 m baseline: conditional parameter uncertainty", "Heading2")]
    rows=[["Parameter", "Known truth", "Posterior median", "95% marginal interval"]]
    truth=[100,25,800,5,20]
    for row,value in zip(runs[0]["diagnostic"]["parameters"],truth):
        unit="ohm m" if row["parameter"].startswith("rho") else "m"
        rows.append([row["parameter"],f"{value:g} {unit}",f"{row['median']:.2f} {unit}",f"{row['q025']:.2f} - {row['q975']:.2f} {unit}"])
    story += [table(rows,[70,105,140,196]),Spacer(1,11),p("All five known values fall within their marginal intervals in this one baseline run. This is not an empirical coverage study. The baseline's deep resistivity is less tightly constrained than its near-surface resistivity."),p("The intervals are conditional on the fixed three-layer model, independent Gaussian errors in log apparent resistivity, known noise, and independent log-uniform priors: resistivity 10-3,000 ohm m and thickness 0.5-100 m. They exclude uncertainty in layer count, geology, geometry, and the noise model."),p("The 30 m chains accept only about 0.55-4.25% of proposals in the 4,000-draw export and 0.85-4.09% in the extension. These low rates, poor R-hat, and small ESS are consistent with poor mixing. They do not establish the true deep-layer resistivity."),PageBreak()]
    story += [p("3. What to share, and what remains open", "Heading1"),Image(str(OUTPUT/"deep_layer_traces.png"),width=511,height=298),Spacer(1,9),p("Figure 2. Log-scale traces of deep-layer resistivity. All panels use the same vertical limits; horizontal axes reflect each run's retained length. The failed short-survey chains explore different regions. The dashed line is the known synthetic truth, not a fitted constraint."),p("Supported conclusion", "Heading2"),p("For these supplied runs, extending the 30 m simulation from 4,000 to 8,000 draws did not resolve convergence failure. The 200 m posterior passed the selected diagnostics. Passing diagnostics does not prove exhaustive mode exploration; failure alone does not prove geological non-identifiability."),p("Practical next stage", "Heading2"),p("Present the 200 m case as a checked synthetic baseline and the 30 m case as an openly documented sampler limitation. Further short-survey interpretation needs a separately validated sampling strategy. Field-data validation, model-class comparison, prior sensitivity, and repeated-simulation coverage testing remain future work."),p("Software and audit trail", "Heading2"),p("GeoScope uses SimPEG 0.25.0 for forward physics, SciPy for bounded least squares, and ArviZ for diagnostics. The project has 31 passing tests in the tested Python 3.12/Linux environment. Version 0.2.1 withholds failed-run parameter estimates by default and records survey metadata in new dashboard exports."),p('Code, methods, figures, machine-readable summaries and SHA-256 source manifest: <link href="https://github.com/TheABZ/geoscope-resistivity" color="#087e8b">github.com/TheABZ/geoscope-resistivity</link>. The packaging script verifies shapes, finite samples, reported quantiles, status thresholds, and the 30 m prefix relationship. The underlying sampler and thresholds were not changed for this report.')]
    def footer(canvas,document):
        canvas.saveState();canvas.setStrokeColor(colors.HexColor("#dce6e8"));canvas.line(42,34,553,34);canvas.setFont("Helvetica",8);canvas.setFillColor(colors.HexColor("#596b77"));canvas.drawString(42,22,"GeoScope | Synthetic software demonstration | Gabriel Chiazor Unomah");canvas.drawRightString(553,22,str(document.page));canvas.restoreState()
    doc.build(story,onFirstPage=footer,onLaterPages=footer)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs",type=Path,required=True)
    parser.add_argument("--output",type=Path,default=Path("results/verified_runs"))
    parser.add_argument("--pdf",type=Path)
    args=parser.parse_args()
    OUTPUT=args.output; OUTPUT.mkdir(parents=True,exist_ok=True)
    runs,summaries,manifest=collect(args.inputs,OUTPUT)
    figures(runs,OUTPUT)
    if args.pdf:pdf_report(args.pdf,runs,summaries)
    print(json.dumps({"runs_checked":len(runs),"prefix_verified":manifest["thirty_metre_4000_is_exact_prefix_of_8000"],"output":str(OUTPUT),"pdf":str(args.pdf) if args.pdf else None}))
