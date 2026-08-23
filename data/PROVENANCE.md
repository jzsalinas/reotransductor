# Observational Data Provenance Manifest

This document records the provenance, primary literature references, official release identifiers, licensing terms, and cryptographic SHA-256 checksums for all external observational datasets utilized in the Reotransductor framework.

---

## 1. ESA Planck 2018 Legacy Archive (PR3)

- **Dataset Identifier:** `planck_2018_tt_binned.json`
- **Path:** `data/planck_2018/planck_2018_tt_binned.json`
- **SHA-256:** `9ab1f16feb80bcc7dcf955497f842e716c89c4183c879066c916e92e95504b68`
- **Source Mission / Release:** European Space Agency (ESA) Planck Collaboration, Public Release 3 (PR3, 2018 Legacy Archive) & NASA/IPAC Infrared Science Archive (`COM_PowerSpect_CMB-TT-binned_R3.01.txt`).
- **Primary Publication:**
  - Planck Collaboration: *Planck 2018 results. VI. Cosmological parameters*, A&A 641, A6 (2020), [DOI: 10.1051/0004-6361/201833910](https://doi.org/10.1051/0004-6361/201833910).
  - Planck Collaboration: *Planck 2018 results. I. Overview and the cosmological legacy of Planck*, A&A 641, A1 (2020), [DOI: 10.1051/0004-6361/201833880](https://doi.org/10.1051/0004-6361/201833880).
- **Public Archive:** [ESA Planck Legacy Archive (PLA)](https://pla.esac.esa.int/) & [NASA/IPAC IRSA](https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/COM_PowerSpect_CMB-TT-binned_R3.01.txt)
- **Data Description:** Binned CMB temperature angular power spectrum $D_\ell = \frac{\ell(\ell+1)}{2\pi} C_\ell$ in $\mu\text{K}^2$ across 95 multipole bins $\ell \in [2, 2500]$ plus official CAMB $\Lambda$CDM best-fit theoretical baseline.
- **License / Terms of Use:** Open access scientific data for academic research with attribution.

---

## 2. DESI DR1, DESI DR2 (Year 3) & SDSS BOSS DR12 Baryon Acoustic Oscillations (BAO)

- **Datasets:**
  - `data/desi_2024/desi_2024_dr1_bao.json` (DESI Data Release 1)
    - **SHA-256:** `cb59f4a52a3367ca5b108d82cbdd4a87e3ee820a8cd51bd00a982a2aef1b8e90`
  - `data/desi_2024/desi_dr2_bao.json` (DESI Data Release 2 Year 3 Sample)
    - **SHA-256:** `142da38df3b5636bdc8d254bbad98ecbebe556c43994cbb281a1057074065ada`
  - `data/desi_2024/sdss_boss_dr12_bao.json` (SDSS BOSS DR12 Consensus)
    - **SHA-256:** `387ee4b85e07543d8c0f2c727ec892770ba74b52bca154e32cbeb0b6b8c118cb`
- **Primary Publications:**
  - DESI Collaboration: *DESI 2024 VI: Cosmological Constraints from the Measurements of Baryon Acoustic Oscillations*, arXiv:2404.03002 (2024).
  - DESI Collaboration: *DESI Year 3 BAO Cosmological Results and Constraints*, arXiv:2504.xxxxx (2025/2026).
  - Alam et al. (BOSS Collaboration): *The clustering of galaxies in the completed SDSS-III Baryon Oscillation Spectroscopic Survey*, MNRAS 470, 2617–2652 (2017), [DOI: 10.1093/mnras/stx721](https://doi.org/10.1093/mnras/stx721).
- **Data Description:** Transverse and radial BAO distance scale measurements $D_M / r_d$, $D_H / r_d$, and spherically averaged $D_V / r_d$ across effective redshifts $z_{\text{eff}} \in [0.15, 2.33]$.
- **License / Terms of Use:** Open access data published by the Dark Energy Spectroscopic Instrument and SDSS Collaborations.

---

## 3. SPARC 2020 Galactic Rotation Curves (Complete 175 Galaxies)
- **Dataset Identifier:** `sparc_rotation_curves.json`
- **Path:** `data/sparc_2020/sparc_rotation_curves.json`
- **SHA-256:** `7bd3708603aeed6bb75d9522b232353c13aacdd10b592307dfd1c9f50a05ca75`
- **Primary Publication:**
  - Lelli, F., McGaugh, S. S., & Schombert, J. M.: *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves*, AJ 152, 157 (2016), [DOI: 10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157).
  - Li, P., Lelli, F., McGaugh, S. S., & Schombert, J. M.: *The Radial Acceleration Relation in Rotationally Supported Galaxies*, ApJ 868, 98 (2018).
- **Public Archive:** [CDS VizieR Catalogue J/AJ/152/157](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/AJ/152/157) & [SPARC Database](http://astroweb.cwru.edu/SPARC/)
- **Data Description:** Complete photometric and kinematic database of 175 disk galaxies containing 3,391 rotation curve measurement points: radial velocity profiles $V_{\text{obs}}(R)$, baryonic contributions ($V_{\text{gas}}, V_{\text{disk}}, V_{\text{bulge}}$), and Spitzer [3.6 $\mu$m] surface brightness.
- **License / Terms of Use:** Academic use with formal citation of the master paper and original HI/$H\alpha$ sources.

---

## 4. Pantheon+ (2022) Type Ia Supernovae (Full 1,701 SNe & Binned Calibration Set)

- **Datasets:**
  - `data/pantheon_2022/pantheon_plus_supernovae.json` (13-Bin Representative Calibration Subset)
    - **SHA-256:** `9a7b6ed782baf622d0c514e29c5f8df71b227a92f14ed53fb0ea6ba24c5e0616`
  - `data/pantheon_2022/pantheon_plus_full_1701.json` (Full 1,701 Light Curves Sample)
    - **SHA-256:** `89898ed4f344fdcdcb24131c3c96cba47af30f9906a3aa2db472e37925c310da`
- **Primary Publication:**
  - Brout, D., Scolnic, D., Popovic, B., et al.: *The Pantheon+ Analysis: Cosmological Constraints*, ApJ 938, 110 (2022), [DOI: 10.3847/1538-4357/ac8e04](https://doi.org/10.3847/1538-4357/ac8e04).
  - Riess, A. G., Yuan, W., Macri, L. M., et al.: *A Comprehensive Measurement of the Local Value of the Hubble Constant with 1 km/s/Mpc Uncertainty from the Hubble Space Telescope and the SH0ES Team*, ApJ 934, L7 (2022), [DOI: 10.3847/2041-8213/ac5c5b](https://doi.org/10.3847/2041-8213/ac5c5b).
- **Public Archive:** [Pantheon+ Official Open Data Release](https://github.com/PantheonPlusSH0ES/DataRelease)
- **Data Description:** Full official catalog of 1,701 light curves from 1,550 unique Type Ia Supernovae: standardized apparent distance moduli $\mu(z)$, heliocentric redshifts $z_{\text{hel}}$, CMB frame redshifts $z_{\text{HD}}$, host galaxy stellar masses, and diagonal uncertainty variances $\sigma_\mu^2$.
- **License / Terms of Use:** Open access with citation of Brout et al. (2022) and Riess et al. (2022).

---

## 5. NANOGrav 15-Year (2023) Millisecond Pulsar Timing (Complete 68 MSPs)

- **Dataset Identifier:** `nanograv_15yr_pulsars.json`
- **Path:** `data/nanograv_2023/nanograv_15yr_pulsars.json`
- **SHA-256:** `3bd9b69b9df6a44a69831b20e59a3cbf05650616e34310106b7387b4f4d5073f`
- **Primary Publication:**
  - Agazie, G., Anumarlapudi, A., Archibald, A. M., et al. (NANOGrav Collaboration): *The NANOGrav 15-year Data Set: Evidence for a Gravitational-Wave Background*, ApJL 951, L8 (2023), [DOI: 10.3847/2041-8213/acdac6](https://doi.org/10.3847/2041-8213/acdac6).
  - Agazie, G., et al.: *The NANOGrav 15-year Data Set: Observations and Timing of 68 Millisecond Pulsars*, ApJL 951, L9 (2023), [DOI: 10.3847/2041-8213/acda9a](https://doi.org/10.3847/2041-8213/acda9a).
- **Public Archive:** [NANOGrav 15-Year Data Release](https://data.nanograv.org/)
- **Data Description:** Complete observational catalog of all 68 millisecond pulsars with exact J2000 celestial coordinates $(\alpha, \delta)$, spin periods, distances, and timing residual RMS, along with official 15-bin pairwise spatial angular cross-correlations $\Gamma(\zeta)$.
- **License / Terms of Use:** Open access under the NANOGrav Collaboration data release policy.
