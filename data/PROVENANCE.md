# Observational Data Provenance Manifest

This document records the provenance, primary literature references, official release identifiers, licensing terms, and cryptographic SHA-256 checksums for all external observational datasets utilized in the Reotransductor framework.

---

## 1. ESA Planck 2018 Legacy Archive (PR3)

- **Dataset Identifier:** `planck_2018_tt_binned.json`
- **Path:** `data/planck_2018/planck_2018_tt_binned.json`
- **SHA-256:** `ab8b4296b2a119a1456505b8e77ec71b29e3d037cdf32bc368323a6bff858531`
- **Source Mission / Release:** European Space Agency (ESA) Planck Collaboration, Public Release 3 (PR3, 2018 Legacy Archive).
- **Primary Publication:**
  - Planck Collaboration: *Planck 2018 results. VI. Cosmological parameters*, A&A 641, A6 (2020), [DOI: 10.1051/0004-6361/201833910](https://doi.org/10.1051/0004-6361/201833910).
  - Planck Collaboration: *Planck 2018 results. I. Overview and the cosmological legacy of Planck*, A&A 641, A1 (2020), [DOI: 10.1051/0004-6361/201833880](https://doi.org/10.1051/0004-6361/201833880).
- **Public Archive:** [ESA Planck Legacy Archive (PLA)](https://pla.esac.esa.int/)
- **Data Description:** Binned CMB temperature angular power spectrum $D_\ell = \frac{\ell(\ell+1)}{2\pi} C_\ell$ in $\mu\text{K}^2$ across multipoles $\ell \in [2, 2500]$.
- **License / Terms of Use:** Open access scientific data for academic research with attribution.

---

## 2. DESI 2024 DR1 & SDSS BOSS DR12 Baryon Acoustic Oscillations (BAO)

- **Datasets:**
  - `data/desi_2024/desi_2024_dr1_bao.json`
    - **SHA-256:** `cb59f4a52a3367ca5b108d82cbdd4a87e3ee820a8cd51bd00a982a2aef1b8e90`
  - `data/desi_2024/sdss_boss_dr12_bao.json`
    - **SHA-256:** `387ee4b85e07543d8c0f2c727ec892770ba74b52bca154e32cbeb0b6b8c118cb`
- **Primary Publications:**
  - DESI Collaboration: *DESI 2024 VI: Cosmological Constraints from the Measurements of Baryon Acoustic Oscillations*, arXiv:2404.03002 (2024).
  - Alam et al. (BOSS Collaboration): *The clustering of galaxies in the completed SDSS-III Baryon Oscillation Spectroscopic Survey: cosmological analysis of the DR12 galaxy sample*, MNRAS 470, 2617–2652 (2017), [DOI: 10.1093/mnras/stx721](https://doi.org/10.1093/mnras/stx721).
- **Data Description:** Transverse and radial BAO distance scale measurements $D_M / r_d$, $D_H / r_d$, and spherically averaged $D_V / r_d$ across effective redshifts $z_{\text{eff}} \in [0.15, 2.33]$.
- **License / Terms of Use:** Open access data published by the Dark Energy Spectroscopic Instrument and SDSS Collaborations.

---

## 3. SPARC 2020 Galactic Rotation Curves
- **Dataset Identifier:** `sparc_rotation_curves.json`
- **Path:** `data/sparc_2020/sparc_rotation_curves.json`
- **SHA-256:** `91302d6be323811df4e2e3983863319eb3142bcf6c9e35c3877640243c5b339d`
- **Primary Publication:**
  - Lelli, F., McGaugh, S. S., & Schombert, J. M.: *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves*, AJ 152, 157 (2016), [DOI: 10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157).
  - Li, P., Lelli, F., McGaugh, S. S., & Schombert, J. M.: *The Radial Acceleration Relation in Rotationally Supported Galaxies*, ApJ 868, 98 (2018).
- **Public Archive:** [SPARC Database (Case Western Reserve University)](http://astroweb.cwru.edu/SPARC/)
- **Data Description:** Radial velocity profiles $V_{\text{obs}}(R)$, baryonic contributions ($V_{\text{gas}}, V_{\text{disk}}, V_{\text{bulge}}$), and surface brightness for rotationally supported galaxies.
- **License / Terms of Use:** Academic use with formal citation of the master paper and original HI/$H\alpha$ sources.

---

## 4. Pantheon+ (2022) Type Ia Supernovae

- **Dataset Identifier:** `pantheon_plus_supernovae.json`
- **Path:** `data/pantheon_2022/pantheon_plus_supernovae.json`
- **SHA-256:** `ff9387b0fa53a15a4c95134afe0c6925e3c3560fb728427165fc65583118e758`
- **Primary Publication:**
  - Brout, D., Scolnic, D., Popovic, B., et al.: *The Pantheon+ Analysis: Cosmological Constraints*, ApJ 938, 110 (2022), [DOI: 10.3847/1538-4357/ac8e04](https://doi.org/10.3847/1538-4357/ac8e04).
  - Riess, A. G., Yuan, W., Macri, L. M., et al.: *A Comprehensive Measurement of the Local Value of the Hubble Constant with 1 km/s/Mpc Uncertainty from the Hubble Space Telescope and the SH0ES Team*, ApJ 934, L7 (2022), [DOI: 10.3847/2041-8213/ac5c5b](https://doi.org/10.3847/2041-8213/ac5c5b).
- **Public Archive:** [Pantheon+ Open Data Archive](https://github.com/PantheonPlusSH0ES/DataRelease)
- **Data Description:** Standardized apparent distance moduli $\mu(z)$, heliocentric redshifts $z_{\text{hel}}$, CMB frame redshifts $z_{\text{HD}}$, host galaxy environmental density classifications, and photometric covariance matrices for 1,701 light curves of 1,550 distinct Type Ia Supernovae.
- **License / Terms of Use:** Open access with citation of Brout et al. (2022) and Riess et al. (2022).

---

## 5. NANOGrav 15-Year (2023) Millisecond Pulsar Timing

- **Dataset Identifier:** `nanograv_15yr_pulsars.json`
- **Path:** `data/nanograv_2023/nanograv_15yr_pulsars.json`
- **SHA-256:** `2e7ee39cef62b8b902bca1054e406f8d144f6a5b1b8d210de87ae25c7bee8cd2`
- **Primary Publication:**
  - Agazie, G., Anumarlapudi, A., Archibald, A. M., et al. (NANOGrav Collaboration): *The NANOGrav 15-year Data Set: Evidence for a Gravitational-Wave Background*, ApJL 951, L8 (2023), [DOI: 10.3847/2041-8213/acdac6](https://doi.org/10.3847/2041-8213/acdac6).
  - Agazie, G., et al.: *The NANOGrav 15-year Data Set: Observations and Timing of 68 Millisecond Pulsars*, ApJL 951, L9 (2023), [DOI: 10.3847/2041-8213/acda9a](https://doi.org/10.3847/2041-8213/acda9a).
- **Public Archive:** [NANOGrav 15-Year Data Release](https://data.nanograv.org/)
- **Data Description:** Binned pairwise spatial angular cross-correlations $\Gamma(\zeta)$, timing residual noise metrics, and characteristic strain power parameters $A_{\text{GWB}} = 2.4 \times 10^{-15}$ at reference frequency $f_{\text{ref}} = 1/\text{yr}$ for 68 millisecond pulsars.
- **License / Terms of Use:** Open access under the NANOGrav Collaboration data release policy.
