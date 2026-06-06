import sqlite3
import json
import csv
from pathlib import Path
from datetime import datetime

DB_DIR = Path(__file__).parent
DB_PATH = DB_DIR / "experiments.db"


def get_connection(db_path: str | Path = DB_PATH) -> sqlite3.Connection:
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


def init_db(db_path: str | Path = DB_PATH) -> None:
    conn = get_connection(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS experiments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp TEXT NOT NULL DEFAULT (datetime('now')),
            params TEXT NOT NULL DEFAULT '{}',
            tissue_type TEXT NOT NULL DEFAULT 'multilayer',
            implant_depth REAL NOT NULL,
            frequency REAL NOT NULL,
            path_loss REAL NOT NULL,
            snr REAL NOT NULL,
            energy_consumption REAL DEFAULT 0.0,
            environment TEXT NOT NULL DEFAULT 'simulation'
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            exp_id INTEGER NOT NULL,
            metric_name TEXT NOT NULL,
            metric_value REAL NOT NULL,
            unit TEXT NOT NULL DEFAULT '',
            timestamp TEXT NOT NULL DEFAULT (datetime('now')),
            FOREIGN KEY (exp_id) REFERENCES experiments(id) ON DELETE CASCADE
        )
    """)

    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_experiments_timestamp
        ON experiments(timestamp)
    """)

    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_experiments_depth_freq
        ON experiments(implant_depth, frequency)
    """)

    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_results_exp_id
        ON results(exp_id)
    """)

    conn.commit()
    conn.close()


def insert_experiment(
    params: dict | None = None,
    tissue_type: str = "multilayer",
    implant_depth: float = 0.0,
    frequency: float = 0.0,
    path_loss: float = 0.0,
    snr: float = 0.0,
    energy_consumption: float = 0.0,
    environment: str = "simulation",
    db_path: str | Path = DB_PATH,
) -> int:
    conn = get_connection(db_path)
    cursor = conn.cursor()

    params_json = json.dumps(params) if params else "{}"
    timestamp = datetime.now().isoformat()

    cursor.execute(
        """
        INSERT INTO experiments
            (timestamp, params, tissue_type, implant_depth, frequency,
             path_loss, snr, energy_consumption, environment)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            timestamp,
            params_json,
            tissue_type,
            implant_depth,
            frequency,
            path_loss,
            snr,
            energy_consumption,
            environment,
        ),
    )

    exp_id = cursor.lastrowid
    conn.commit()
    conn.close()
    return exp_id


def insert_result(
    exp_id: int,
    metric_name: str,
    metric_value: float,
    unit: str = "",
    db_path: str | Path = DB_PATH,
) -> int:
    conn = get_connection(db_path)
    cursor = conn.cursor()

    timestamp = datetime.now().isoformat()
    cursor.execute(
        """
        INSERT INTO results (exp_id, metric_name, metric_value, unit, timestamp)
        VALUES (?, ?, ?, ?, ?)
        """,
        (exp_id, metric_name, metric_value, unit, timestamp),
    )

    res_id = cursor.lastrowid
    conn.commit()
    conn.close()
    return res_id


def query_experiments(
    tissue_type: str | None = None,
    depth_min: float | None = None,
    depth_max: float | None = None,
    freq_min: float | None = None,
    freq_max: float | None = None,
    limit: int = 100,
    db_path: str | Path = DB_PATH,
) -> list[sqlite3.Row]:
    conn = get_connection(db_path)
    cursor = conn.cursor()

    conditions = []
    params_list = []

    if tissue_type:
        conditions.append("tissue_type = ?")
        params_list.append(tissue_type)
    if depth_min is not None:
        conditions.append("implant_depth >= ?")
        params_list.append(depth_min)
    if depth_max is not None:
        conditions.append("implant_depth <= ?")
        params_list.append(depth_max)
    if freq_min is not None:
        conditions.append("frequency >= ?")
        params_list.append(freq_min)
    if freq_max is not None:
        conditions.append("frequency <= ?")
        params_list.append(freq_max)

    where_clause = " AND ".join(conditions) if conditions else "1=1"
    sql = f"SELECT * FROM experiments WHERE {where_clause} ORDER BY timestamp DESC LIMIT ?"
    params_list.append(limit)

    cursor.execute(sql, params_list)
    rows = cursor.fetchall()
    conn.close()
    return rows


def query_results_for_experiment(
    exp_id: int,
    db_path: str | Path = DB_PATH,
) -> list[sqlite3.Row]:
    conn = get_connection(db_path)
    cursor = conn.cursor()
    cursor.execute(
        "SELECT * FROM results WHERE exp_id = ? ORDER BY timestamp ASC", (exp_id,)
    )
    rows = cursor.fetchall()
    conn.close()
    return rows


def export_experiments_to_csv(
    csv_path: str | Path | None = None,
    db_path: str | Path = DB_PATH,
) -> Path:
    if csv_path is None:
        csv_path = DB_DIR / "exported_experiments.csv"
    csv_path = Path(csv_path)

    conn = get_connection(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM experiments ORDER BY timestamp DESC")
    rows = cursor.fetchall()

    if not rows:
        conn.close()
        return csv_path

    column_names = [desc[0] for desc in cursor.description]

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(column_names)
        for row in rows:
            writer.writerow(list(row))

    conn.close()
    return csv_path


def get_statistics(db_path: str | Path = DB_PATH) -> dict:
    conn = get_connection(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT COUNT(*) FROM experiments")
    total_exp = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM results")
    total_results = cursor.fetchone()[0]

    cursor.execute("""
        SELECT
            MIN(path_loss), MAX(path_loss), AVG(path_loss),
            MIN(snr), MAX(snr), AVG(snr)
        FROM experiments
    """)
    stats = cursor.fetchone()

    cursor.execute("""
        SELECT tissue_type, COUNT(*) as cnt
        FROM experiments GROUP BY tissue_type ORDER BY cnt DESC
    """)
    tissue_counts = {row["tissue_type"]: row["cnt"] for row in cursor.fetchall()}

    conn.close()

    return {
        "total_experiments": total_exp,
        "total_results": total_results,
        "path_loss_min": round(stats[0], 2) if stats[0] is not None else None,
        "path_loss_max": round(stats[1], 2) if stats[1] is not None else None,
        "path_loss_avg": round(stats[2], 2) if stats[2] is not None else None,
        "snr_min": round(stats[3], 2) if stats[3] is not None else None,
        "snr_max": round(stats[4], 2) if stats[4] is not None else None,
        "snr_avg": round(stats[5], 2) if stats[5] is not None else None,
        "tissue_type_distribution": tissue_counts,
    }


def main():
    print("初始化数据库...")
    init_db()
    print(f"数据库已创建: {DB_PATH}")

    for params, depth_m, freq_hz, pl, snr_val, env in [
        ({"n_freq": 50, "n_depth": 20}, 0.03, 402e6, 45.2, 12.5, "simulation"),
        ({"n_freq": 50, "n_depth": 20}, 0.03, 915e6, 55.8, 5.2, "simulation"),
        ({"n_freq": 50, "n_depth": 20}, 0.05, 402e6, 52.1, 8.3, "simulation"),
        ({"n_freq": 50, "n_depth": 20}, 0.05, 2.4e9, 68.4, -2.1, "simulation"),
        ({"n_freq": 50, "n_depth": 20}, 0.08, 433e6, 60.0, 4.7, "simulation"),
    ]:
        exp_id = insert_experiment(
            params=params,
            tissue_type="multilayer",
            implant_depth=depth_m,
            frequency=freq_hz,
            path_loss=pl,
            snr=snr_val,
            energy_consumption=1e-8,
            environment=env,
        )
        insert_result(exp_id, "path_loss_exponent", 4.5, "dB/cm")
        insert_result(exp_id, "channel_capacity", 500e3, "bps")

    print("已插入 5 条示例实验记录")

    stats = get_statistics()
    print("\n数据库统计:")
    for key, value in stats.items():
        print(f"  {key}: {value}")

    csv_path = export_experiments_to_csv()
    print(f"\n实验数据已导出到: {csv_path}")

    print("\n最近实验记录:")
    for row in query_experiments(limit=5):
        print(
            f"  ID={row['id']}, 深度={row['implant_depth']*1000:.0f}mm, "
            f"频率={row['frequency']/1e6:.0f}MHz, "
            f"PL={row['path_loss']:.1f}dB, SNR={row['snr']:.1f}dB"
        )


if __name__ == "__main__":
    main()
