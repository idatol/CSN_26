from pathlib import Path

import duckdb

DB = Path(__file__).parent / "transfermarkt-datasets.duckdb"

# read_only avoids accidentally modifying the database
#with duckdb.connect(str(DB), read_only=True) as con:

    # .df() returns the query result as a pandas DataFrame
    # Initial query to get the top 10 transfers by transfer fee
#    transfers = con.sql("""
#        SELECT player_name, from_club_name, to_club_name, transfer_date, transfer_fee
#        FROM transfers
#        ORDER BY transfer_fee DESC NULLS LAST
#        LIMIT 10
#    """).df()

print(transfers)

with duckdb.connect(str(DB), read_only=True) as con:
    #Query to get all transfers for season 18/19, and save it to a file
    season_transfers = con.sql("""
        SELECT player_id, player_name, from_club_id, from_club_name,
               to_club_id, to_club_name, transfer_date, transfer_fee
        FROM transfers
        WHERE transfer_season = '18/19'
        ORDER BY transfer_fee DESC NULLS LAST
    """).df()

    print(season_transfers)
    season_transfers.to_csv(DB.parent / "season_transfers_18_19.csv", index=False)

with duckdb.connect(str(DB), read_only=True) as con:
    # Count transfers by fee status for season 18/19
    fee_counts = con.sql("""
        SELECT
            COUNT(*) AS total,
            COUNT(*) FILTER (transfer_fee > 0) AS paid,
            COUNT(*) FILTER (transfer_fee = 0) AS zero_fee,
            COUNT(*) FILTER (transfer_fee IS NULL) AS missing_fee
        FROM transfers
        WHERE transfer_season = '18/19'
    """).df()

    print(fee_counts)
    fee_counts.to_csv(DB.parent / "fee_counts_18_19.csv", index=False)