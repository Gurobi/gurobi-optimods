# Generates the workforce dataset

from pathlib import Path

import names
import numpy as np
import pandas as pd
import sqlalchemy

if (p := Path("customer-data.db")).exists():
    p.unlink()
if (p := Path("staff-availability.xlsx")).exists():
    p.unlink()

data = pd.DataFrame(
    dict(
        Date=pd.date_range(
            start=pd.Timestamp("2022-01-01"),
            end=pd.Timestamp("2022-12-31"),
            freq="D",
        ),
    )
).assign(
    Customers=lambda df: np.random.uniform(low=5, high=100, size=df.shape[0]),
    DayOfWeek=lambda df: df.Date.dt.dayofweek,
)
days = range(7)
muls = [1.1, 0.9, 1.5, 1.0, 1.2, 2.3, 0.5]
for day, mul in zip(days, muls):
    data.loc[data.DayOfWeek == day, "Customers"] *= mul
data = data.assign(Customers=lambda df: df["Customers"].round().astype(int)).drop(
    columns=["DayOfWeek"]
)

engine = sqlalchemy.create_engine("sqlite:///customer-data.db")
with engine.connect() as connection:
    data.to_sql("history", connection, index=False)

staff_names = list({names.get_first_name() for i in range(100)})[:10]
preferences = pd.DataFrame(
    [
        {
            "Date": date,
            "Staff": staff,
            "Response": "Yes" if np.random.random() < 0.7 else "No",
        }
        for date in pd.date_range(
            start=pd.Timestamp("2023-06-01"),
            periods=30,
            freq="D",
        )
        for staff in staff_names
    ]
)
(
    preferences.set_index(["Date", "Staff"])["Response"]
    .unstack()
    .reset_index()
    .to_excel("staff-availability.xlsx", index=False)
)
