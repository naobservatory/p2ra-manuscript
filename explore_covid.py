import pandas as pd

with open("prevalence-data/time_series_covid19_confirmed_US.csv") as inf:
    df = pd.read_csv(inf)

diff_df = df.iloc[:, 11:] = (
    df.iloc[:, 11:].diff(axis=1).fillna(df.iloc[:, 11:])
)
montgomery_county_ohio = df[
    (df["Admin2"] == "Montgomery") & (df["Province_State"] == "Ohio")
]
row_with_largest_values = (
    montgomery_county_ohio.iloc[:, 11:].sum(axis=1).idxmax()
)
largest_14_day_sum = (
    montgomery_county_ohio.iloc[:, 11:]
    .rolling(window=14, axis=1)
    .sum()
    .max(axis=1)
    .idxmax()
)

print(df.loc[largest_14_day_sum])
print(diff_df.loc[largest_14_day_sum])
