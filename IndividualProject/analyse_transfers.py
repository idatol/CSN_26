"""Descriptive analysis of one season of transfers (file made by explore.py)."""
from pathlib import Path

import pandas as pd

CSV = Path(__file__).parent / "season_transfers_18_19.csv"

# Entries in the club columns that are not real clubs
NOT_CLUBS = ["Without Club", "Unknown", "Career break", "Retired", "Ban", "Own Youth"]

df = pd.read_csv(CSV)


def section(title):
    print(f"\n=== {title} ===")


def share(label, n, total):
    print(f"{label:<42}{n:>7,}  ({n / total:.1%})")


# --- Overview ---------------------------------------------------------------
section("Overview")
# Count clubs by id: a few clubs appear under more than one name
club_ids = pd.concat([df.from_club_id, df.to_club_id])
print(f"Transfers:      {len(df):,}")
print(f"Unique players: {df.player_id.nunique():,}")
print(f"Unique clubs:   {club_ids.nunique():,}  (including the non-club entries)")

# --- Fees -------------------------------------------------------------------
paid = df[df.transfer_fee > 0]

section("Transfers by fee")
share("Fee above 0", len(paid), len(df))
share("Fee of 0 (may include loans)", (df.transfer_fee == 0).sum(), len(df))
share("Fee missing", df.transfer_fee.isna().sum(), len(df))

section("Players by fee")
share("Players with at least one fee above 0", paid.player_id.nunique(), df.player_id.nunique())

section("Size of the fees above 0 (million EUR)")
fees = paid.transfer_fee / 1e6
print(f"Total:   {fees.sum():,.1f}")
print(f"Mean:    {fees.mean():.2f}")
print(f"Median:  {fees.median():.2f}")
print(f"Min/max: {fees.min():.2f} / {fees.max():.2f}")
top10 = fees.nlargest(max(1, round(len(fees) * 0.1)))
print(f"The top 10% of fees make up {top10.sum() / fees.sum():.1%} of all money spent")

section("Fees above 0, by range (million EUR)")
bins = [0, 0.5, 1, 5, 10, 25, 50, float("inf")]
labels = ["up to 0.5", "0.5-1", "1-5", "5-10", "10-25", "25-50", "over 50"]
for label, n in pd.cut(fees, bins=bins, labels=labels).value_counts(sort=False).items():
    share(label, n, len(fees))

section("10 highest fees")
print(paid.nlargest(10, "transfer_fee")[
    ["player_name", "from_club_name", "to_club_name", "transfer_fee"]
].to_string(index=False, float_format="{:,.0f}".format))

# --- Entries that are not clubs ---------------------------------------------
from_none = df.from_club_name.isin(NOT_CLUBS)
to_none = df.to_club_name.isin(NOT_CLUBS)
no_club = from_none | to_none

section("Transfers without a real club")
share("Coming from a non-club", from_none.sum(), len(df))
share("Going to a non-club", to_none.sum(), len(df))
share("Either of these", no_club.sum(), len(df))
share("Players involved", df[no_club].player_id.nunique(), df.player_id.nunique())
print()
print(pd.DataFrame({
    "from": df.from_club_name[from_none].value_counts(),
    "to": df.to_club_name[to_none].value_counts(),
}).fillna(0).astype(int).to_string())

# --- Would they distort a network? ------------------------------------------
section("Most connected clubs (number of different clubs traded with)")
ends = df[["from_club_id", "to_club_id"]]
both = pd.concat([ends.set_axis(["club", "other"], axis=1),
                  ends[["to_club_id", "from_club_id"]].set_axis(["club", "other"], axis=1)])
names = pd.concat([df.set_index("from_club_id").from_club_name,
                   df.set_index("to_club_id").to_club_name])
names = names[~names.index.duplicated()]
degree = both.groupby("club").other.nunique().sort_values(ascending=False)
print(f"Mean: {degree.mean():.1f}   Median: {degree.median():.0f}")
for club, n in degree.head(10).items():
    flag = "  <-- not a club" if names[club] in NOT_CLUBS else ""
    print(f"{names[club]:<25}{n:>5}{flag}")

section("Effect of removing the non-clubs")
clean = df[~no_club]
share("Transfers kept", len(clean), len(df))
share("Players kept", clean.player_id.nunique(), df.player_id.nunique())
