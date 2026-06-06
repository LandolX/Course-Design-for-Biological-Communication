import subprocess, json, time, re, os, sys

arxiv_ids = [
    "1704.04206", "1808.04273", "1808.05147", "1811.00417", "1812.05492",
    "1911.08291", "2007.01799", "2112.12485", "2207.01875", "2303.08015",
    "2305.05527", "2311.16356", "2403.20029", "2405.14044", "2406.09875",
    "2409.18616", "2410.15943", "2503.13738", "2504.12123", "2505.22849",
    "2506.14360", "2506.17112", "2506.22137", "2508.19739", "2510.11743",
    "2511.02074", "2603.28412"
]

def fetch_abstract(aid):
    url = f"https://arxiv.org/abs/{aid}"
    cmd = f'curl -s --max-time 15 "{url}"'
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=20)
        html = r.stdout
        title_m = re.search(r'<title>\s*\[([^\]]+)\]\s*(.+?)\s*</title>', html)
        title = title_m.group(2).strip() if title_m else "N/A"
        abs_m = re.search(r'<span class="descriptor">Abstract:</span>\s*(.+?)(?:</blockquote>|<span|\s*$)', html, re.DOTALL)
        abstract = abs_m.group(1).strip() if abs_m else "N/A"
        abstract = re.sub(r'<[^>]+>', '', abstract).strip()
        abstract = abstract[:800]
        return {"id": aid, "title": title, "abstract": abstract}
    except Exception as e:
        return {"id": aid, "title": "N/A", "abstract": f"Error: {e}"}

results = []
for i, aid in enumerate(arxiv_ids):
    print(f"[{i+1}/{len(arxiv_ids)}] Fetching {aid}...", flush=True)
    r = fetch_abstract(aid)
    results.append(r)
    print(f"  Title: {r['title'][:80]}", flush=True)
    if i < len(arxiv_ids) - 1:
        time.sleep(3.5)

out_path = os.path.join(os.path.dirname(__file__), "paper_metadata.json")
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
print(f"\nSaved to {out_path}")
