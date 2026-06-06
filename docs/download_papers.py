#!/usr/bin/env python3
"""Download arXiv papers related to in-body communication systems for drug delivery/in-vivo detection."""

import subprocess
import re
import os

OUTPUT_DIR = "/Users/landolx/Desktop/生物通信工程/课程设计/docs"

PAPERS = [
    # === Latest papers (2024-2026) ===
    ("2603.28412", "Joint Detection and Identification for Scalable Control of Nanorobot Swarms under Harsh Communication Constraints (2026)"),
    ("2511.02074", "Neural Network based Distance Estimation for Branched Molecular Communication Systems (2025)"),
    ("2510.11743", "Mixture of Inverse Gaussians for Hemodynamic Transport (MIGHT) in Vascular Networks (2025)"),
    ("2508.19739", "Molecular Communication for Gastroretentive Drug Delivery (2025)"),
    ("2506.22137", "On Drug Delivery System Parameter Optimisation via Semantic Information Theory (2025)"),
    ("2506.17112", "Closed-Loop Molecular Communication with Local and Global Degradation (2025)"),
    ("2506.14360", "Identification for Molecular Communication Based on Diffusion Channel (2025)"),
    ("2505.22849", "Flexure-FET-Based Receiver for Interference Mitigation in Molecular Communication (2025)"),
    ("2504.12123", "The CAM Model: An in vivo Testbed for Molecular Communication Systems (2025)"),
    ("2503.13738", "General Molecular Communication Model in Multi-Layered Spherical Channels (2025)"),
    ("2410.15943", "Molecular Signal Reception in Complex Vessel Networks (2024)"),
    ("2409.18616", "Enhanced Drug Delivery via Localization-Enabled Relaying in MC Nanonetworks (2024)"),
    ("2406.09875", "The Chorioallantoic Membrane Model: A 3D in vivo Testbed for MC Systems (2024)"),
    ("2405.14044", "Single Input Multi Output Model of Molecular Communication via Diffusion (2024)"),
    ("2403.20029", "Analysis of Signal Distortion in Molecular Communication Channels (2024)"),

    # === Core papers (2022-2023) ===
    ("2311.16356", "What Really is Molecule in Molecular Communications? (2023)"),
    ("2305.05527", "Microparticle-based Controlled Drug Delivery Systems (2023)"),
    ("2303.08015", "Molecular Communication for Quorum Sensing Inspired Cooperative Drug Delivery (2023)"),
    ("2207.01875", "The End-to-End MC Model of Extracellular Vesicle-based Drug Delivery (2022)"),
    ("2112.12485", "On the Reception Process of Molecular Communication-Based Drug Delivery (2021)"),

    # === Foundational papers (2018-2020) ===
    ("2007.01799", "Transfer Function Models for Cylindrical MC Channels with Diffusion and Laminar Flow (2020)"),
    ("1911.08291", "Toward a Wired Ad Hoc Nanonetwork (2019)"),
    ("1812.05492", "Channel Modeling for Diffusive Molecular Communication: A Tutorial Review (2018)"),
    ("1811.00417", "Diffusive Mobile MC for Controlled-Release Drug Delivery with Absorbing Receiver (2018)"),
    ("1808.05147", "Magnetic Nanoparticle Based Molecular Communication in Microfluidic Environments (2019)"),
    ("1808.04273", "A Molecular Communications Model for Drug Delivery (2018)"),
    ("1704.04206", "Molecular Communication using Magnetic Nanoparticles (2017)"),
]

def download_paper(arxiv_id, title):
    pdf_url = f"https://arxiv.org/pdf/{arxiv_id}.pdf"
    safe_title = re.sub(r'[^\w\s-]', '', title)
    safe_title = re.sub(r'[-\s]+', '_', safe_title).strip('_')[:60]
    filename = f"{arxiv_id}_{safe_title}.pdf"
    filepath = os.path.join(OUTPUT_DIR, filename)

    if os.path.exists(filepath) and os.path.getsize(filepath) > 1000:
        size_kb = os.path.getsize(filepath) / 1024
        print(f"  EXISTS ({size_kb:.0f} KB): {filename}")
        return True

    result = subprocess.run(
        ["curl", "-sL", "-o", filepath, pdf_url, "--max-time", "120"],
        capture_output=True, text=True, timeout=130
    )
    if os.path.exists(filepath) and os.path.getsize(filepath) > 1000:
        size_kb = os.path.getsize(filepath) / 1024
        print(f"  OK ({size_kb:.0f} KB): {filename}")
        return True
    else:
        print(f"  FAILED: {arxiv_id} ({title[:50]})")
        if os.path.exists(filepath):
            os.remove(filepath)
        return False

def main():
    print(f"{'='*60}")
    print(f"Papers to download: {len(PAPERS)}")
    print(f"Target folder: {OUTPUT_DIR}")
    print(f"{'='*60}")

    success = 0
    for arxiv_id, title in PAPERS:
        print(f"\n[{arxiv_id}] {title}")
        if download_paper(arxiv_id, title):
            success += 1

    print(f"\n{'='*60}")
    print(f"Summary: {success}/{len(PAPERS)} papers downloaded successfully")
    print(f"{'='*60}")
    print(f"\nContents of docs folder:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        if f.endswith('.pdf'):
            size_kb = os.path.getsize(os.path.join(OUTPUT_DIR, f)) / 1024
            print(f"  {f} ({size_kb:.0f} KB)")

if __name__ == "__main__":
    main()
