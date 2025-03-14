import pycatima as catima
import numpy as np
from simcc import GetMaterial, GetMFP, GetMCHistories, GetMCDeltaE, GetMCProb
import time


def GetCAtimaCompound(zts, m_fractions):
    return [[0, z, m] for m, z in zip(m_fractions, zts)]


def GetMCEloss(A, Z, Q, energy, material, length, N=10000, random_state=None, histories=None):

    assert random_state is not None
    # 電荷状態を外から与えられるようにする
    config = catima.Config()

    # 有効Zを有効にする
    config.z_effective = 1

    # 電荷履歴がある場合に履歴の最後の厚みを取得
    if histories is not None:
        l1 = histories[0][-1][1]
    else:
        l1 = 0

    # 物質情報を取得
    result = GetMaterial(material)
    zts, m_fractions, density, solid_gas = result["zts"], result["m_fractions"], result["density"], result["solid_gas"]
    compound = GetCAtimaCompound(zts, m_fractions)

    # 分割数を検討
    mat = catima.Material(compound, density=density)
    mat.thickness_cm(length)
    layer = catima.Layers()
    layer.add(mat)
    res = catima.calculate_layers(catima.Projectile(A=A, Z=Z, Q=0, T=energy), layer, config)
    Eloss = res.total_result.Ein - res.total_result.Eout
    stragglingColRate = res.total_result.sigma_E / Eloss
    nlayer = int(np.ceil(Eloss / 5.0))

    # 分割層ごとのエネルギーロスを計算
    layer = catima.Layers()
    mat = catima.Material(compound, density=density)
    mat.thickness_cm(length / nlayer)
    for ilayer in range(nlayer):
        layer.add(mat)
    res = catima.calculate_layers(catima.Projectile(A=A, Z=Z, Q=0, T=energy), layer, config)
    print(f"{energy:.3f} -> {res.total_result.Eout:.3f} MeV/u")
    assert res.total_result.Eout >= 50

    # 有効Zを無効にしてQは外から与える
    config.z_effective = 0

    # 分割層ごと電荷変化MFPを取得
    MFPs = []
    for ilayer in range(nlayer):
        assert res.results[ilayer].Eout >= 50
        ene = (res.results[ilayer].Ein + res.results[ilayer].Eout) * 0.5
        MFPs.append(GetMFP(Z, ene, material))

    # 分割層ごとに電荷履歴を計算
    print("GetMCHistories", end=" ")
    start = time.time()
    for ilayer in range(nlayer):
        print(".", end="")
        histories = GetMCHistories(MFPs[ilayer], Q, length / nlayer, random_state=random_state, histories=histories, N=N)
    print(f"Elapsed time {time.time()-start:.1f} s")

    # 分割層ごとにdEdxを計算
    Eout = energy

    dEcc = np.zeros(N)
    print("GetMCDeltaE", end=" ")
    start = time.time()
    for ilayer in range(nlayer):
        print(".", end="")
        dedx_list = {}
        for Q in range(Z - 6, Z + 1):
            layer = catima.Layers()
            mat = catima.Material(compound, density=density)
            mat.thickness_cm(length / nlayer)
            layer.add(mat)
            res = catima.calculate_layers(catima.Projectile(A=A, Z=Z, Q=Q, T=Eout), layer, config)
            dedx_list[Q] = (res.total_result.Ein - res.total_result.Eout) / (length / nlayer)

        dE = GetMCDeltaE(histories, l1 + ilayer * length / nlayer, l1 + (ilayer + 1) * length / nlayer, dedx_list)
        Eout = Eout - np.mean(dE)
        dEcc += dE
    dEcolRate = 1 + random_state.normal(loc=0, scale=stragglingColRate, size=N)
    dEcol = np.mean(dEcc) * dEcolRate
    dEtotal = dEcc * dEcolRate
    print(f"Elapsed time {time.time()-start:.1f} s")
    print(f"{energy:.3f} -> {energy-np.mean(dEcc):.3f} MeV/u")

    # 電荷分布をlog scaleで20層に分けて計算
    print("GetMCProb", end=" ")
    start = time.time()
    length_log = np.array([0] + list(np.logspace(np.log10(length / 1000), np.log10(length), num=20)))
    charges = {"length": length_log}
    for l2 in length_log:
        print(".", end="")
        charge_prob = np.array(GetMCProb(histories, l1 + l2))
        for dQ in range(0, 7):
            if Z - dQ not in charges:
                charges[Z - dQ] = []
            charges[Z - dQ].append(charge_prob[dQ])
    print(f"Elapsed time {time.time()-start:.1f} s")
    return dEtotal, dEcol, dEcc, charges, histories


def GetAnalyticalEloss(A, Z, energy, material, length, z_effective=1, density=0):

    # 物質情報を取得
    result = GetMaterial(material)
    zts, m_fractions, nominal_density, solid_gas = result["zts"], result["m_fractions"], result["density"], result["solid_gas"]
    compound = GetCAtimaCompound(zts, m_fractions)

    # 分割数を検討
    mat = catima.Material(compound, density=nominal_density if density == 0 else density)
    mat.thickness_cm(length)
    layer = catima.Layers()
    layer.add(mat)

    config = catima.Config()
    if z_effective == 1:
        config.z_effective = 1
        res = catima.calculate_layers(catima.Projectile(A=A, Z=Z, Q=0, T=energy), layer, config)
    else:
        config.z_effective = 0
        res = catima.calculate_layers(catima.Projectile(A=A, Z=Z, Q=z_effective, T=energy), layer, config)

    Eloss = res.total_result.Ein - res.total_result.Eout

    return Eloss, res.total_result.sigma_E


def GetAnalyticalRange(A, Z, energy, material, density=0):
    result = GetMaterial(material)
    zts, m_fractions, nominal_density, solid_gas = result["zts"], result["m_fractions"], result["density"], result["solid_gas"]
    compound = GetCAtimaCompound(zts, m_fractions)
    mat = catima.Material(compound)
    config = catima.Config()
    config.z_effective = 1
    res = catima.calculate(catima.Projectile(A=A, Z=Z, Q=0, T=energy), mat, config)
    density = nominal_density if density == 0 else density
    return res.range / density, res.sigma_r / density
