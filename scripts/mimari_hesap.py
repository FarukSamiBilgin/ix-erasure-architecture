# -*- coding: utf-8 -*-
"""
IX mimari makalesi — sayisal omurga hesaplari (acik-hesap listesi,
merkez kalemler). Cikti: ../data/mimari_hesap.json

H1' : cift-tuzak tunelleme/X haritasi (1D sonlu-fark; t(V0,a), P_flip, X-kapisi)
H2' : kapi-ekranlamali dipolar J_eff (4-yuk imaj toplamlari; J_eff(a,R,d_m))
H4' : (a) ZZ "kutupla-ve-fazla" CZ — OU-MC, eko var/yok, sigma_E in {3.7,10,47}
      (b) Jaksch blokaj-CZ varyanti — motor-v1 3-seviye sablonu IX
          parametreleriyle (V=59.2/473 ueV, tau=139 ns) + sigma_E detuning MC
ERASE: kavite-silinme/flas gecikme butcesi (Purcell T1'=2.2 ns, eta-endeksli)
EC   : silinme-butcesi aritmetigi p_L(d; p_e) (analitik; Stim MC = H8' gelecek)

Birimler: enerji meV (ueV belirtilen yerde), uzunluk nm, zaman ps.
hbar = 0.65821195 meV*ps;  e^2/(4 pi eps0) = 1439.964 meV*nm.
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBAR = 0.65821195          # meV ps
COUL = 1439.964            # e^2/(4 pi eps0), meV nm
EPS_R = 4.5                # hBN
P_DIP = 0.43               # e nm (OLCULMUS H-yigili tuzakli IX dipolu, Baek)
M0_INV = 76.1996           # hbar^2/m0 = 76.1996 meV nm^2  (2*38.0998)
TAU_BRIGHT_NS = 139.0      # kalibre parlak IX omru (Purcell makalesi)

OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
R = {}

# =====================================================================
# H1' — cift-tuzak tunelleme haritasi
# =====================================================================


def cift_tuzak_t(V0, a, m_rel=0.9, sigma_w=None, dV_merkez=0.0, n=6001):
    """1D FD: V(x) = -V0[G(x-a/2)+G(x+a/2)] + dV_merkez*G(x); Gauss sigma_w.
    Dondurur: t = (E1-E0)/2 [meV], E0, E10 orbital bosluk E2-E1."""
    if sigma_w is None:
        sigma_w = a / 4.0
    L = a / 2 + 6 * sigma_w + 8.0
    x = np.linspace(-L, L, n)
    dx = x[1] - x[0]
    g = lambda c: np.exp(-(x - c) ** 2 / (2 * sigma_w ** 2))
    V = -V0 * (g(a / 2) + g(-a / 2)) + dV_merkez * g(0.0)
    h2m = (M0_INV / 2) / m_rel          # hbar^2/(2m) meV nm^2
    diag = V + 2 * h2m / dx ** 2
    off = -h2m / dx ** 2 * np.ones(n - 1)
    w = eigh_tridiagonal(diag, off, select="i", select_range=(0, 3),
                         eigvals_only=True)
    return 0.5 * (w[1] - w[0]), w[0], w[2] - w[1]


def h1_harita():
    V0_list = [50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0]
    a_list = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
    out = {"V0_meV": V0_list, "a_nm": a_list, "t_ueV": [], "E0_meV": [],
           "sigma_w": "a/4", "m_rel": 0.9}
    for V0 in V0_list:
        satir_t, satir_E = [], []
        for a in a_list:
            t, E0, _ = cift_tuzak_t(V0, a)
            satir_t.append(t * 1e3)       # ueV
            satir_E.append(E0)
        out["t_ueV"].append(satir_t)
        out["E0_meV"].append(satir_E)
    # duyarlilik: genislik ve kutle (v2 noktasi V0=150, a=20)
    duy = {}
    for et, kw in [("sw_a3", dict(sigma_w=20 / 3)), ("sw_a5", dict(sigma_w=4.0)),
                   ("m08", dict(m_rel=0.8)), ("m10", dict(m_rel=1.0))]:
        t, _, _ = cift_tuzak_t(150.0, 20.0, **kw)
        duy[et] = t * 1e3
    t_ref, _, _ = cift_tuzak_t(150.0, 20.0)
    duy["ref"] = t_ref * 1e3
    out["duyarlilik_v2_ueV"] = duy

    # X-kapisi: merkez-kapiyla bariyer indirme (dV_merkez<0 cukur acar),
    # t_on hedef 5-52 ueV (SWAP 207-20 ps)
    xk = {"dV_meV": [], "t_on_ueV": [], "tX_ps": [], "E_gap_meV": []}
    for dV in (0.0, -20.0, -40.0, -60.0, -80.0, -100.0, -110.0, -120.0,
               -130.0, -140.0):
        t, _, gap = cift_tuzak_t(150.0, 20.0, dV_merkez=dV)
        xk["dV_meV"].append(dV)
        xk["t_on_ueV"].append(t * 1e3)
        xk["tX_ps"].append(math.pi * HBAR / (2 * t) if t > 0 else None)
        xk["E_gap_meV"].append(gap)
    out["X_kapisi"] = xk

    # P_flip/tur haritasi: Gamma_flip = (t/delta)^2 * gamma_phi ;
    # gamma_phi = sigma_eff/hbar (hizli-defaz orani tavani), tur = 5 ns
    sigma_eff = 0.8e-3   # meV (DFS sonrasi orta deger)
    gphi = sigma_eff / HBAR          # 1/ps
    tur_ps = 5e3
    pf = {"t_ueV": [0.05, 0.1, 0.2, 0.3, 0.5], "delta_ueV": [20.0, 30.0, 43.0],
          "P_flip_tur": []}
    for tt in pf["t_ueV"]:
        satir = []
        for dd in pf["delta_ueV"]:
            satir.append((tt / dd) ** 2 * gphi * tur_ps)
        pf["P_flip_tur"].append(satir)
    pf["not"] = "Gamma_flip=(t/delta)^2*gamma_phi; gamma_phi=sigma_eff/hbar, sigma_eff=0.8 ueV; tur=5 ns"
    out["P_flip"] = pf
    return out


# =====================================================================
# H2' — kapi-ekranlamali dipolar etkilesim (4-yuk imaj toplamlari)
# =====================================================================


def dipol_cift_U(r_lat, d_m, delta=P_DIP, eps=EPS_R):
    """Iki dusey dipol (p = e*delta), merkezleri metal kapinin d_m ustunde,
    yanal ayrim r_lat. Topraklanmis duzlem imaji: q -> -q, z -> -z.
    U = U(1,2) + U(1, im2)  [meV, yuk e birimli]."""
    z0 = d_m
    q1 = [(+1.0, 0.0, z0 + delta / 2), (-1.0, 0.0, z0 - delta / 2)]
    q2 = [(+1.0, r_lat, z0 + delta / 2), (-1.0, r_lat, z0 - delta / 2)]
    im2 = [(-s, xx, -zz) for (s, xx, zz) in q2]
    U = 0.0
    for s1, x1, zz1 in q1:
        for s2, x2, zz2 in q2 + im2:
            U += s1 * s2 / math.hypot(x2 - x1, zz2 - zz1)
    return COUL / eps * U


def dipol_V_noktasal(r, p=P_DIP, eps=EPS_R):
    return COUL * p * p / (eps * r ** 3)


def h2_ekran():
    out = {"kontrol": {}, "V_r": {}, "J_eff": {}}
    # kontrol: d_m -> sonsuz, nokta-dipol limiti
    for r in (10.0, 20.0, 40.0):
        U_inf = dipol_cift_U(r, 1e6)
        out["kontrol"][f"r={r:g}"] = {
            "U_4yuk_ueV": U_inf * 1e3,
            "nokta_dipol_ueV": dipol_V_noktasal(r) * 1e3}
    # V(r; d_m)
    r_list = [10, 15, 20, 25, 30, 35, 40, 50, 60, 80]
    for d_m in (5.0, 10.0, 20.0, 1e6):
        key = "bare" if d_m > 1e5 else f"dm={d_m:g}"
        out["V_r"][key] = {"r_nm": r_list,
                           "V_ueV": [dipol_cift_U(float(r), d_m) * 1e3
                                     for r in r_list]}
    # J_eff(a, R=2a; d_m) = |V(R-a) - 2V(R) + V(R+a)|
    out["J_eff"]["tanim"] = "J=|V(R-a)-2V(R)+V(R+a)|, R=2a"
    for d_m in (5.0, 10.0, 20.0, 1e6):
        key = "bare" if d_m > 1e5 else f"dm={d_m:g}"
        rec = {}
        for a in (10.0, 20.0, 35.0):
            Rr = 2 * a
            J = abs(dipol_cift_U(Rr - a, d_m) - 2 * dipol_cift_U(Rr, d_m)
                    + dipol_cift_U(Rr + a, d_m))
            rec[f"a={a:g}"] = {"J_ueV": J * 1e3,
                               "t_CZ_ps": math.pi * HBAR / J if J > 0 else None}
        out["J_eff"][key] = rec
    # NN-otesi artik (bare ve dm=10): V(2R)/V(R) ve (J ile kiyas)
    for d_m, key in ((10.0, "dm=10"), (1e6, "bare")):
        V1 = dipol_cift_U(40.0, d_m) * 1e3
        V2 = dipol_cift_U(80.0, d_m) * 1e3
        out["J_eff"].setdefault("NN_otesi", {})[key] = {
            "V(40)_ueV": V1, "V(80)_ueV": V2, "oran": V2 / V1}
    return out


# =====================================================================
# H4'a — ZZ kutupla-ve-fazla CZ: OU Monte Carlo, eko var/yok
# =====================================================================


def ou_yol(n, dt, sigma, tau_c, rng):
    """Ornstein-Uhlenbeck yolu (stasyoner baslangic)."""
    x = np.empty(n)
    x[0] = rng.normal(0, sigma)
    if tau_c <= 0:
        return rng.normal(0, sigma, n)
    a = math.exp(-dt / tau_c)
    b = sigma * math.sqrt(1 - a * a)
    for i in range(1, n):
        x[i] = a * x[i - 1] + b * rng.normal()
    return x


def zz_cz_mc(J_ueV=5.82, sigma_eff_ueV=0.5, tau_c_ps=1e6, eko=True,
             n_traj=400, n_step=200, seed=7):
    """CZ: H = eps1 Z1/2 + eps2 Z2/2 + (J/4) Z1Z2 (diyagonal; yol-basina
    analitik faz). Eko: tau/2'de E_z polarite cevirme -> eps isaret degistirir
    (p -> -p), ZZ (p^2) etkilenmez. Fidelite: |<CZ psi0|psi>|^2, psi0=|++>."""
    J = J_ueV * 1e-3
    sig = sigma_eff_ueV * 1e-3
    tau = math.pi * HBAR / J          # ps
    dt = tau / n_step
    rng = np.random.default_rng(seed)
    F = 0.0
    yarim = n_step // 2
    isaret = np.ones(n_step)
    if eko:
        isaret[yarim:] = -1.0
    for _ in range(n_traj):
        e1 = ou_yol(n_step, dt, sig, tau_c_ps, rng)
        e2 = ou_yol(n_step, dt, sig, tau_c_ps, rng)
        a1 = float(np.sum(isaret * e1) * dt / HBAR)   # tek-qubit faz hatasi
        a2 = float(np.sum(isaret * e2) * dt / HBAR)
        # psi = 1/2 sum_{ab} exp(i phi_ab) |ab>; phi = a1*s1/2 + a2*s2/2 + pi/4*s1s2 (J tam)
        # hedef CZ fazlari: diag(1,1,1,-1) konvansiyonunda |11> pi fazi J'den tam gelir.
        amp = 0.0 + 0.0j
        for s1 in (+1, -1):
            for s2 in (+1, -1):
                ph = 0.5 * a1 * s1 + 0.5 * a2 * s2
                amp += 0.25 * np.exp(1j * ph)
        F += abs(amp) ** 2
    return 1.0 - F / n_traj, tau


def h4a_zz():
    out = {"J_ueV": 5.82, "satirlar": []}
    DFS = {"3.7": (0.4, 1.2), "10": (1.0, 3.3), "47": (4.7, 15.7)}
    for sE, (lo, hi) in DFS.items():
        for sig in (lo, hi):
            for tau_c, tc_et in ((1e3, "1 ns"), (1e4, "10 ns"), (1e6, "1 us")):
                for eko in (False, True):
                    e, tau = zz_cz_mc(sigma_eff_ueV=sig, tau_c_ps=tau_c,
                                      eko=eko, n_traj=600)
                    out["satirlar"].append(dict(
                        sigma_E_ueV=float(sE), sigma_eff_ueV=sig,
                        tau_c=tc_et, eko=eko, hata=e))
    out["t_CZ_ps"] = math.pi * HBAR / (5.82e-3)
    out["not"] = ("sigma_eff = sigma_E*Dx/lambda_c (DFS 3-10x); eko = kapi-ortasi "
                  "E_z polarite cevirme (p->-p: eps isaret, ZZ degismez)")
    return out


# =====================================================================
# H4'b — Jaksch blokaj-CZ (3 seviye), IX parametreleri + sigma_E detuning
# =====================================================================


def lindblad_rk4(rho, H, Ls, T, dt):
    def drho(r):
        d = -1j / HBAR * (H @ r - r @ H)
        for L in Ls:
            LdL = L.conj().T @ L
            d += L @ r @ L.conj().T - 0.5 * (LdL @ r + r @ LdL)
        return d
    n = max(1, int(round(T / dt)))
    h = T / n
    for _ in range(n):
        k1 = drho(rho)
        k2 = drho(rho + 0.5 * h * k1)
        k3 = drho(rho + 0.5 * h * k2)
        k4 = drho(rho + h * k3)
        rho = rho + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return rho


def blokaj_cz(V_ueV, tau_ns, Omega_ueV, det1_ueV=0.0, det2_ueV=0.0):
    """Motor-v1 kolA_cz sablonu (faz-duyarli metrik) + site detuningleri."""
    V = V_ueV * 1e-3
    Om = Omega_ueV * 1e-3
    G = HBAR / (tau_ns * 1e3)
    d = 3
    sp = np.zeros((d, d), complex); sp[2, 1] = 1
    nE = np.zeros((d, d)); nE[2, 2] = 1
    I3 = np.eye(d, dtype=complex)
    Hd = 0.5 * Om * (sp + sp.T.conj())
    D1 = det1_ueV * 1e-3 * np.kron(nE, I3)
    D2 = det2_ueV * 1e-3 * np.kron(I3, nE)
    Hint = V * np.kron(nE, nE) + D1 + D2
    L1 = np.sqrt(G / HBAR) * np.kron(sp.T.conj(), I3)
    L2 = np.sqrt(G / HBAR) * np.kron(I3, sp.T.conj())
    dt_cap = 0.2 * HBAR / max(V, Om, 1e-6)      # motor v2 dt kurali

    def seg(rho, hc, ht, T):
        H = np.kron(hc, I3) + np.kron(I3, ht) + Hint
        return lindblad_rk4(rho, H, [L1, L2], T, min(T / 600, dt_cap))

    tpi = math.pi * HBAR / Om
    psi0 = np.zeros(d * d, complex)
    for (a, b) in ((0, 0), (0, 1), (1, 0), (1, 1)):
        psi0[a * d + b] = 0.5
    rho = np.outer(psi0, psi0.conj())
    rho = seg(rho, Hd, 0 * Hd, tpi)
    rho = seg(rho, 0 * Hd, Hd, 2 * tpi)
    rho = seg(rho, Hd, 0 * Hd, tpi)
    hedef = np.zeros(d * d, complex)
    for (a, b), s in {(0, 0): 1, (0, 1): -1, (1, 0): -1, (1, 1): -1}.items():
        hedef[a * d + b] = 0.5 * s
    F = float(np.real(hedef.conj() @ rho @ hedef))
    return 1 - F, 4 * tpi


def h4b_blokaj():
    out = {"tau_bright_ns": TAU_BRIGHT_NS, "satirlar": []}
    rng = np.random.default_rng(11)
    for V in (59.2, 473.0):
        # once sigma=0'da Omega optimizasyonu
        best = None
        for Om in np.geomspace(1.0, min(0.5 * V, 60.0), 10):
            e, tg = blokaj_cz(V, TAU_BRIGHT_NS, float(Om))
            if best is None or e < best[0]:
                best = (e, float(Om), tg)
        e0, Om0, tg0 = best
        out["satirlar"].append(dict(V_ueV=V, sigma_ueV=0.0, Omega_ueV=Om0,
                                    t_gate_ps=tg0, hata=e0, n_mc=1))
        # sigma_E detuning MC (trim sonrasi kalinti / spec / olculmus)
        for sig in (0.5, 3.7, 10.0, 47.0):
            es = []
            for _ in range(12):
                d1, d2 = rng.normal(0, sig, 2)
                e, _ = blokaj_cz(V, TAU_BRIGHT_NS, Om0, d1, d2)
                es.append(e)
            out["satirlar"].append(dict(V_ueV=V, sigma_ueV=sig, Omega_ueV=Om0,
                                        t_gate_ps=tg0, hata=float(np.mean(es)),
                                        hata_std=float(np.std(es)), n_mc=12))
    out["not"] = ("Jaksch pi-2pi-pi; surulen |e> = parlak IX (tau=139 ns); "
                  "detuning = site enerji kaymasi (Stark-trim oncesi sigma_E, "
                  "sonrasi kalinti). Motor-v1 faz-duyarli metrik + dt kurali.")
    return out


# =====================================================================
# Silinme / flas gecikme butcesi (analitik)
# =====================================================================


def silinme_butce():
    T1p_ns = 2.2            # Purcell'li flas omru (Q=7.9k tasarim)
    beta = {"1.5K": 0.985, "10K": 0.956}
    eta_snspd = 0.8
    out = {"T1_flas_ns": T1p_ns, "beta": beta, "eta_SNSPD": eta_snspd,
           "pencere": [], "amortize": [], "R_e": {}}
    for Tet, b in beta.items():
        for eta in (0.05, 0.1, 0.2):
            p = b * eta * eta_snspd
            # dolu/bos ayrimi icin hata<1e-2: N_gamma = ln(100)/p foton;
            # cevrim ~ 2*T1' (yeniden-uyarim dahil)
            Ng = math.log(100.0) / p
            tw = Ng * 2 * T1p_ns
            out["pencere"].append(dict(T=Tet, eta=eta, p_klik=p,
                                       N_foton=Ng, t_pencere_ns=tw))
    # amortize raf kaybi: pencere/K tur basina; kural raf>=30*pencere, K>=25
    for tw, raf_us, K in ((515.0, 3.0, 25), (515.0, 10.0, 25),
                          (130.0, 1.0, 25), (1040.0, 30.0, 25)):
        kayip = tw * 1e-3 / (raf_us * K)
        out["amortize"].append(dict(t_pencere_ns=tw, raf_us=raf_us, K=K,
                                    kayip_tur=kayip))
    # R_e: tek-flas klik olasiligi (kayip-bayragi): pencere icinde >=1 klik
    for Tet, b in beta.items():
        rec = {}
        for eta in (0.05, 0.1, 0.2):
            p = b * eta * eta_snspd
            # pencerede N_c = t_w/(2T1') cevrim; R_e = 1-(1-p)^(N_c) ama
            # kayip-bayragi tek IX'in TOPLAM yakalanma olasiligi ~ p (tek foton,
            # kacis halinde yeniden-uyarim yok) -> durust bant [p, 1-(1-p)^5]
            rec[f"eta={eta}"] = {"alt(tek foton)": p,
                                 "ust(5 cevrim)": 1 - (1 - p) ** 5}
        out["R_e"][Tet] = rec
    out["not"] = ("Birincil kayip-bayragi FOTONSUZ doluluk sayaci (ancilla "
                  "dipol-flasi, kayma 59 ueV >> sigma_E); kavite kliki bonus. "
                  "R_e=0.4-0.7 bandi eta>=0.5 toplama VE cevrimli yakalama ister "
                  "-> V/xi-iyilestirmeli kaviteye kosullu (Purcell makalesi yol B).")
    return out


# =====================================================================
# EC butcesi aritmetigi (analitik; Stim MC = H8' gelecek isi)
# =====================================================================


def ec_butce():
    out = {"satirlar": []}
    for tur_ns, raf_us in ((5.0, 1.0), (5.0, 10.0), (10.0, 1.0)):
        p_e = tur_ns * 1e-3 / raf_us
        for p_th, kod in ((0.0415, "kod-kapasite R_e=0.98"),
                          (0.082, "devre-duzeyi biased-erasure"),
                          (0.103, "XZZX-fusion")):
            rec = dict(tur_ns=tur_ns, raf_us=raf_us, p_e=p_e, esik=p_th,
                       kod=kod, marj=p_th / p_e)
            for dd in (3, 5):
                rec[f"p_L_d{dd}"] = (p_e / p_th) ** dd
            out["satirlar"].append(rec)
    out["kaynak_sayimi"] = {"d3_tuzak": 26, "d5_tuzak": 74,
                            "d3_alan_20nm": "~200x200 nm",
                            "giris_merdiveni_n3": 8}
    out["not"] = "p_L ~ (p_e/p_th)^d (silinme-baskin); Pauli-Z katkisi bias eta=10-100 ile XZZX'e devredilir"
    return out


# =====================================================================
if __name__ == "__main__":
    print("H1' cift-tuzak...")
    R["H1_cift_tuzak"] = h1_harita()
    print("H2' ekranlama...")
    R["H2_ekran"] = h2_ekran()
    print("H4'a ZZ-CZ OU-MC...")
    R["H4a_zz_cz"] = h4a_zz()
    print("H4'b blokaj-CZ...")
    R["H4b_blokaj_cz"] = h4b_blokaj()
    print("silinme butcesi...")
    R["silinme"] = silinme_butce()
    R["EC"] = ec_butce()
    (OUT / "mimari_hesap.json").write_text(
        json.dumps(R, indent=1, ensure_ascii=False), encoding="utf-8")
    print("yazildi:", OUT / "mimari_hesap.json")

    # ozet
    h1 = R["H1_cift_tuzak"]
    print("\n--- t(V0,a) [ueV], satir=V0, sutun=a:", h1["a_nm"])
    for V0, satir in zip(h1["V0_meV"], h1["t_ueV"]):
        print(f"V0={V0:5.0f}:", ["%.3g" % v for v in satir])
    print("duyarlilik v2:", {k: "%.3g" % v for k, v in
                             h1["duyarlilik_v2_ueV"].items()})
    print("X kapisi:", [(d, "%.3g" % t, "%.3g" % (x or 0)) for d, t, x in
                        zip(h1["X_kapisi"]["dV_meV"], h1["X_kapisi"]["t_on_ueV"],
                            h1["X_kapisi"]["tX_ps"])])
    print("\n--- H2 kontrol:", json.dumps(R["H2_ekran"]["kontrol"], indent=1))
    print("J_eff:", json.dumps(R["H2_ekran"]["J_eff"], indent=1)[:900])
    print("\n--- H4a (ilk 8):")
    for r in R["H4a_zz_cz"]["satirlar"][:8]:
        print(r)
    print("\n--- H4b:")
    for r in R["H4b_blokaj_cz"]["satirlar"]:
        print({k: (round(v, 6) if isinstance(v, float) else v)
               for k, v in r.items()})
