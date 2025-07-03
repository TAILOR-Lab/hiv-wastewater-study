import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import os
import re
import matplotlib.dates as mdates
import matplotlib.lines as mlines

metadata_dir = "./metadata"
read_summary_file = "read_mapping_summary.csv"
masking_file = "HIV_paper_city_site_codes_FAKE.xlsx"

metadata_files = glob(os.path.join(metadata_dir, "p*_metadata.xlsx"))
frames = []
for f in metadata_files:
    pool_id = os.path.basename(f).split("_")[0]
    meta = pd.read_excel(f)
    meta["Pool_ID"] = pool_id
    frames.append(meta)

meta_df = pd.concat(frames, ignore_index=True)
meta_df.columns = meta_df.columns.str.strip()
meta_df["Site"] = meta_df["Site"].astype(str).str.strip().str.replace("\xa0", " ", regex=False)
meta_df["City"] = meta_df["City"].astype(str).str.strip().str.replace("\xa0", " ", regex=False)
meta_df["Date"] = pd.to_datetime(meta_df["Date"])
meta_df["Sample_ID"] = meta_df["Sample_ID"].astype(str).str.strip()

city_mask = pd.read_excel(masking_file, sheet_name=0)
site_mask = pd.read_excel(masking_file, sheet_name=1)
city_mask.columns = city_mask.columns.str.strip()
site_mask.columns = site_mask.columns.str.strip()
city_mask["City"] = city_mask["City"].astype(str).str.strip().str.replace("\xa0", " ", regex=False)
city_mask["Code"] = city_mask["Code"].astype(str).str.strip().str.replace("\xa0", " ", regex=False)
site_mask["Site"] = site_mask["Site"].astype(str).str.strip().str.replace("\xa0", " ", regex=False)
site_mask["Code"] = site_mask["Code"].astype(str).str.strip().str.replace("\xa0", " ", regex=False)

meta_df = meta_df.merge(city_mask.rename(columns={"Code": "City_Code"}), on="City", how="left")
meta_df = meta_df.merge(site_mask.rename(columns={"Code": "Site_Code"}), on="Site", how="left")
meta_df["Site_Label"] = "Site " + meta_df["Site_Code"].str.replace("Site ", "", regex=False)

read_df = pd.read_csv(read_summary_file)
read_df["Sample_ID"] = read_df["Sample_ID"].astype(str).str.strip()
sample_read_counts = read_df["Sample_ID"].value_counts().rename_axis("Sample_ID").reset_index(name="Read_Count")

meta_df = meta_df.merge(sample_read_counts, on="Sample_ID", how="left")
meta_df["Read_Count"] = meta_df["Read_Count"].fillna(0).astype(int)

meta_df["YearWeek_dt"] = meta_df["Date"] - pd.to_timedelta(meta_df["Date"].dt.weekday, unit='d')
meta_df["YearWeek_dt"] = meta_df["YearWeek_dt"].dt.to_period("W").apply(lambda r: r.start_time)

agg = meta_df.groupby(["Site_Label", "YearWeek_dt"], as_index=False)["Read_Count"].sum()
agg["Detected"] = agg["Read_Count"] > 0
agg["log_reads"] = agg["Read_Count"].apply(lambda x: np.log10(x) if x > 0 else 0)

def extract_site_parts(site_label):
    match = re.search(r"Site\s*([A-Z])(\d+)", site_label)
    if match:
        return (match.group(1), int(match.group(2)))
    return ("Z", 999)

site_order_df = (
    pd.DataFrame({"Site_Label": agg["Site_Label"].unique()})
    .assign(Sort_Key=lambda df: df["Site_Label"].apply(extract_site_parts))
    .sort_values(by=["Sort_Key"])
    .reset_index(drop=True)
)
site_order_df["Y_Position"] = site_order_df.index

agg = agg.merge(site_order_df[["Site_Label", "Y_Position"]], on="Site_Label", how="left")

def log_to_size_continuous(log_val, base=20, scale=200):
    return base + scale * log_val if log_val > 0 else 0

agg["Dot_Size"] = agg["log_reads"].apply(lambda x: log_to_size_continuous(x, base=20, scale=200))

sns.set(style="whitegrid", font_scale=1.25)
plt.figure(figsize=(15, len(site_order_df) * 0.38))

ndf = agg[agg["Detected"] == False]
plt.scatter(
    ndf["YearWeek_dt"], ndf["Y_Position"],
    marker="x", s=55, color="lightgrey", zorder=1, label="Non-detect"
)

df = agg[agg["Detected"] == True]
plt.scatter(
    df["YearWeek_dt"], df["Y_Position"],
    s=df["Dot_Size"],
    c="navy",
    marker="o", edgecolor="black", linewidths=0.4, alpha=0.97, zorder=2
)

plt.yticks(site_order_df["Y_Position"], site_order_df["Site_Label"], fontsize=10)
plt.gca().invert_yaxis()
plt.xlabel("Collection Date")
plt.ylabel("Site")
plt.title("Weekly Wastewater HIV Detection by Site", fontsize=16)

ax = plt.gca()
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
plt.xticks(rotation=45)

legend_sizes = [np.log10(v) for v in [10, 100, 1000]]
legend_dots = [log_to_size_continuous(val, base=20, scale=200) for val in legend_sizes]
legend_labels = ["~10", "~100", "~1,000"]
legend_handles = [
    mlines.Line2D([], [], color='navy', marker='o', linestyle='None', markersize=np.sqrt(s), label=l)
    for s, l in zip(legend_dots, legend_labels)
]
legend_handles.append(
    mlines.Line2D([], [], color='lightgrey', marker='x', linestyle='None', markersize=10, label='Non-detect')
)
plt.legend(
    handles=legend_handles,
    title="Read Count",
    bbox_to_anchor=(1.01, 0.5),
    loc="center left",
    borderaxespad=0.0,
    frameon=True,
    labelspacing=1.1
)

plt.tight_layout(rect=[0, 0, 0.85, 1])

plt.savefig("hiv_dot_plot_masked.png", dpi=300, bbox_inches="tight")
plt.savefig("hiv_dot_plot_masked.svg", bbox_inches="tight")
plt.savefig("hiv_dot_plot_masked.pdf", dpi=300, bbox_inches="tight")
plt.show()
