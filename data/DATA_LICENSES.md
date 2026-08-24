# Third-Party Observational Data Manifest and Usage Licenses

This document outlines the official provenance, intellectual property, usage terms, and mandatory academic citations for all third-party observational astrophysics and cosmology datasets utilized by the **Reotransductor** framework.

The core source code of Reotransductor is distributed under the open-source **MIT License**. Third-party observational datasets remain subject to the respective terms and citation policies of their scientific collaborations.

---

## 1. ESA Planck Satellite (2018 Data Release 3)
* **Dataset:** Binned Temperature Power Spectrum $C_\ell^{TT}$ ($2 \le \ell \le 2508$).
* **File:** `data/planck_2018/planck_2018_tt_binned.json` (Source: `COM_PowerSpect_CMB-TT-binned_R3.01.txt`).
* **Source Organization:** European Space Agency (ESA) & Planck Collaboration.
* **Persistent Identifier / DOI:** [10.1051/0004-6361/201833910](https://doi.org/10.1051/0004-6361/201833910)
* **License & Terms:** Open Data Access (ESA Open Data Policy). Free for scientific and educational use with mandatory citation.
* **Required Citation:**
  > Planck Collaboration, Aghanim, N., Akrami, Y., et al. (2020). *Planck 2018 results. VI. Cosmological parameters*. Astronomy & Astrophysics, 641, A6.

---

## 2. DESI (Dark Energy Spectroscopic Instrument - 2024 DR1 / DR2)
* **Dataset:** Baryon Acoustic Oscillation (BAO) measurements across multiple redshift tracers ($0.1 < z < 4.2$).
* **Files:** `data/desi_2024/desi_2024_dr1_bao.json`, `data/desi_2024/desi_dr2_bao.json`.
* **Source Organization:** Dark Energy Spectroscopic Instrument (DESI) Collaboration / Lawrence Berkeley National Laboratory (LBNL).
* **Persistent Identifier / DOI:** [10.48550/arXiv.2404.03002](https://doi.org/10.48550/arXiv.2404.03002)
* **License & Terms:** Public Scientific Data Release.
* **Required Citation:**
  > DESI Collaboration, Adame, A. G., Aguilar, J., et al. (2024). *DESI 2024 VI: Cosmological Constraints from the Measurements of Baryon Acoustic Oscillations*. arXiv:2404.03002.

---

## 3. SDSS-III BOSS (Baryon Oscillation Spectroscopic Survey - DR12)
* **Dataset:** Consensus BAO distance scale measurements ($z = 0.38, 0.51, 0.61$).
* **File:** `data/desi_2024/sdss_boss_dr12_bao.json`.
* **Source Organization:** Sloan Digital Sky Survey (SDSS-III) Collaboration.
* **Persistent Identifier / DOI:** [10.1093/mnras/stx721](https://doi.org/10.1093/mnras/stx721)
* **License & Terms:** Public Domain / SDSS Open Access Policy.
* **Required Citation:**
  > Alam, S., Ata, M., Bailey, S., et al. (2017). *The clustering of galaxies in the completed SDSS-III Baryon Oscillation Spectroscopic Survey: cosmological analysis of the DR12 galaxy sample*. Monthly Notices of the Royal Astronomical Society, 470(3), 2617-2652.

---

## 4. SPARC (Spitzer Photometry & Accurate Rotation Curves - 2020)
* **Dataset:** High-precision rotation curves and mass models for 175 disc galaxies.
* **File:** `data/sparc_2020/sparc_rotation_curves.json` (Source: `sparc_table1.dat`, `sparc_table2.dat`).
* **Principal Investigators:** Federico Lelli, Stacy S. McGaugh, James M. Schombert.
* **Persistent Identifier / DOI:** [10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157)
* **License & Terms:** Academic Open Access (AAS / Astronomical Journal).
* **Required Citation:**
  > Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves*. The Astronomical Journal, 152(6), 157.

---

## 5. Pantheon+ Supernova Collaboration (2022)
* **Dataset:** Standardized Type Ia Supernovae distance moduli and redshift catalogue (1,701 light curves across 1,550 unique SNe Ia, $0.001 < z < 2.26$).
* **Files:** `data/pantheon_2022/pantheon_plus_supernovae.json`, `data/pantheon_2022/pantheon_plus_full_1701.json` (Source: `Pantheon+SH0ES.dat`).
* **Principal Investigators:** Dan Scolnic, Dillon Brout, Adam G. Riess, et al.
* **Persistent Identifier / DOI:** [10.3847/1538-4357/ac8b04](https://doi.org/10.3847/1538-4357/ac8b04)
* **License & Terms:** Open Data Release (AAS / Astrophysical Journal).
* **Required Citation:**
  > Brout, D., Scolnic, D., Popovic, B., et al. (2022). *The Pantheon+ Analysis: Cosmological Constraints*. The Astrophysical Journal, 938(2), 110.

---

## 6. NANOGrav 15-Year Data Release (2023)
* **Dataset:** Spatial cross-correlation and timing residuals across 67 millisecond pulsars.
* **File:** `data/nanograv_2023/nanograv_15yr_pulsars.json`.
* **Source Organization:** North American Nanohertz Observatory for Gravitational Waves (NANOGrav).
* **Persistent Identifier / DOI:** [10.3847/2041-8213/acdac6](https://doi.org/10.3847/2041-8213/acdac6)
* **License & Terms:** Open Access Data Release (Astrophysical Journal Letters).
* **Required Citation:**
  > Agazie, G., Anumarlapudi, A., Archibald, A. M., et al. (2023). *The NANOGrav 15 yr Data Set: Evidence for a Gravitational-wave Background*. The Astrophysical Journal Letters, 951(1), L8.
