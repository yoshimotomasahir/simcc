import bisect
import numpy as np
import scipy

periodic_table = ["H", "He",\
"Li", "Be", "B", "C", "N", "O", "F", "Ne",\
"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",\
"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",\
"Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",\
"Cs", "Ba", \
"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",\
"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",\
"Fr", "Ra", \
"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",\
"Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

z2symbol = {}
for z, symbol in zip(range(1, len(periodic_table) + 1), periodic_table):
    z2symbol[z] = symbol

symbol2z = {}
for z, symbol in zip(range(1, len(periodic_table) + 1), periodic_table):
    symbol2z[symbol] = z

standard_weights = [1.008, 4.003,\
6.941, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00, 20.18,\
22.99, 24.31, 26.98, 28.09, 30.97, 32.07, 35.45, 39.95,\
39.10, 40.08, 44.96, 47.87, 50.94, 52.00, 54.94, 55.85, 58.93, 58.69, 63.55, 65.38, 69.72, 72.63, 74.92, 78.97, 79.90, 83.80,\
85.47, 87.62, 88.91, 91.22, 92.91, 95.95, 99, 101.1, 102.9, 106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6, 126.9, 131.3,\
132.9, 137.3,\
138.9, 140.1, 140.9, 144.2, 145, 150.4, 152.0, 157.3, 158.9, 162.5, 164.9, 167.3, 168.9, 173.0, 175.0,\
178.5, 180.9, 183.8, 186.2, 190.2, 192.2, 195.1, 197.0, 200.6, 204.4, 207.2, 209.0, 210, 210, 222,\
223, 226,\
227, 232.0, 231.0, 238.0, 237, 239, 243, 247, 247, 252, 252, 257, 258, 259, 262,\
267, 268, 271, 272, 277, 276, 281, 280, 285, 278, 289, 289, 293, 293, 294]

z2weight = {}
for z, standard_weight in zip(range(1, len(standard_weights) + 1), standard_weights):
    z2weight[z] = standard_weight


def ParseCS(zp, energy, zt, solid_gas, n_line):
    """
    断面積ファイルから断面積を得る
    """
    import os

    with open(f"{os.path.dirname(__file__)}/ChargeStates/CS_{zp}_{solid_gas}.txt") as f:
        line = f.readlines()[n_line]
    cs = {}
    strs = line.split()
    assert str(zp) == strs[0]
    assert str(energy) == strs[1]
    assert str(zt) == strs[2]
    cs[f"{zp-0}->{zp-1}"] = float(strs[3])
    cs[f"{zp-1}->{zp-0}"] = float(strs[4])
    cs[f"{zp-1}->{zp-2}"] = float(strs[5])
    cs[f"{zp-2}->{zp-1}"] = float(strs[6])
    cs[f"{zp-2}->{zp-3}"] = float(strs[7])
    cs[f"{zp-3}->{zp-2}"] = float(strs[8])
    cs[f"{zp-3}->{zp-4}"] = float(strs[9])
    cs[f"{zp-4}->{zp-3}"] = float(strs[10])
    cs[f"{zp-4}->{zp-5}"] = float(strs[11])
    cs[f"{zp-5}->{zp-4}"] = float(strs[12])
    cs[f"{zp-5}->{zp-6}"] = float(strs[13])
    cs[f"{zp-6}->{zp-5}"] = float(strs[14])
    return cs


def GetCS(zp, energy, zt, solid_gas="solid"):
    """
    Projectile (Z=zp, energy [MeV/u]) を Target (Z=zt) に入射したときの電荷変化断面積を得る
    zpは29以上94以下
    energy は 50 MeV/u以上、1000 MeV/u以下
    断面積の単位はbarn
    """
    assert energy >= 50 and energy <= 1000
    assert zp >= 29 and zp <= 94
    assert zt >= 1 and zt <= 94
    assert solid_gas in ["solid", "gas"]


    if len([a for a in range(50, 1001, 5) if np.isclose(a, energy)]) > 0:
        n_line = 1 + (zt - 1) * 191 + (int(np.round(energy)) - 50) // 5
        return ParseCS(zp, int(np.round(energy)), zt, solid_gas, n_line)
    else:
        n_line = 1 + (zt - 1) * 191 + (int(energy) - 50) // 5
        e1 = ((int(energy)) // 5 + 0) * 5
        e2 = ((int(energy)) // 5 + 1) * 5
        cs1 = ParseCS(zp, e1, zt, solid_gas, n_line + 0)
        cs2 = ParseCS(zp, e2, zt, solid_gas, n_line + 1)
        cs = {}
        for key, value1, value2 in zip(cs1.keys(), cs1.values(), cs2.values()):
            cs[key] = value1 * (e2 - energy) / (e2 - e1) + value2 * (energy - e1) / (e2 - e1)
        return cs


def GetPureMFP(zp, energy, zt, solid_gas="solid", density=1):
    """
    Projectile (Z=zp, energy [MeV/u]) を Target (Z=zt) に入射したときの電荷変化平均自由行程を得る
    densityの単位はg/cm3で、値を正しく入れたときの平均自由行程の単位はcm
    densityを指定しない場合は、平均自由行程の単位はg/cm2
    """
    cs = GetCS(zp, energy, zt, solid_gas)
    u = 1.66054e-24  # g
    mfp = {}
    for key, value in cs.items():
        mfp[key] = z2weight[zt] * u / (value * density) * 1e24
    return mfp


def GetMixedMFP(zp, energy, zts, m_fractions, solid_gas="solid", density=1):
    """
    Projectile (Z=zp, energy [MeV/u]) を Target 混合物 に入射したときの電荷変化平均自由行程を得る
    混合物は Zの組と 個数比で指定する。質量比ではない
    densityの単位はg/cm3で、値を正しく入れたときの平均自由行程の単位はcm
    densityを指定しない場合は、平均自由行程の単位はg/cm2
    """
    MFPs = {}
    MFP = {}
    assert len(zts) == len(m_fractions)
    fractions = []
    for zt, m_f in zip(zts, m_fractions):
        fractions.append(z2weight[zt] * m_f)
    sum_fractions = sum(fractions)
    for zt in zts:
        MFPBuf = GetPureMFP(zp, energy, zt, solid_gas, 1)
        for a in MFPBuf.items():
            if a[0] not in MFPs:
                MFPs[a[0]] = []
            MFPs[a[0]].append(a[1])
    for a in MFPs.items():
        sum_cs = 0
        for i, b in enumerate(a[1]):
            sum_cs += 1 / b * fractions[i] / sum_fractions
        MFP[a[0]] = 1 / sum_cs / density
    return MFP


def GetMaterial(material, density_factor=1):

    if type(material) == int:
        material = str(material)

    if len(material.split()) == 3 and material.split()[2] == "Torr":
        density_factor = float(material.split()[1]) / 760
    if len(material.split()) >= 2:
        material = material.split()[0]

    densities = {"CH4": 0.667, "PureAr": 1.662, "PureXe": 5.46}
    if material == "CH4" or material == "Methane":
        solid_gas = "gas"
        density = 0.001 * densities["CH4"]  # PDG 0.667 (20deg 1atm) #CATIMAと密度が異なる
        zts, m_fractions = [6, 1], [1, 4]

    elif material == "CF4":
        solid_gas = "gas"
        density = 0.001 * 3.78  # PDG 3.78 (20deg 1atm)
        zts, m_fractions = [6, 9], [1, 4]

    elif material == "iC4H10" or material == "Isobutane":
        solid_gas = "gas"
        density = 0.001 * 2.49
        zts, m_fractions = [6, 1], [4, 10]

    elif material == "P10" or material == "CH4Ar9":
        solid_gas = "gas"
        density = 0.001 * (densities["CH4"] * 0.1 + densities["PureAr"] * 0.9)
        zts, m_fractions = [6, 1, 18], [1, 1 * 4, 9]

    elif material == "(CH4)1Xe9" or material == "CH4Xe9" or material == "Xe9":
        solid_gas = "gas"
        density = 0.001 * (densities["CH4"] * 0.1 + densities["PureXe"] * 0.9)
        zts, m_fractions = [6, 1, 54], [1, 1 * 4, 9]

    elif material == "(CH4)2Xe8" or material == "CH4Xe4" or material == "Xe8":
        solid_gas = "gas"
        density = 0.001 * (densities["CH4"] * 0.2 + densities["PureXe"] * 0.8)
        zts, m_fractions = [6, 1, 54], [2, 2 * 4, 8]

    elif material == "(CH4)3Xe7" or material == "Xe7":
        solid_gas = "gas"
        density = 0.001 * (densities["CH4"] * 0.3 + densities["PureXe"] * 0.7)
        zts, m_fractions = [6, 1, 54], [3, 3 * 4, 7]

    elif material == "(CH4)4Xe6" or material == "Xe6":
        solid_gas = "gas"
        density = 0.001 * (densities["CH4"] * 0.4 + densities["PureXe"] * 0.6)
        zts, m_fractions = [6, 1, 54], [4, 4 * 4, 6]

    elif material == "(CH4)5Xe5" or material == "CH4Xe" or material == "Xe5":
        solid_gas = "gas"
        density = 0.001 * (densities["CH4"] * 0.5 + densities["PureXe"] * 0.5)
        zts, m_fractions = [6, 1, 54], [5, 5 * 4, 5]

    elif material == "Mylar":
        solid_gas = "solid"
        density = 1.380
        zts, m_fractions = [1, 6, 8], [8, 10, 4]  # C10H8O4

    elif material == "Kapton":
        solid_gas = "solid"
        density = 1.420
        zts, m_fractions = [1, 6, 7, 8], [10, 22, 2, 5]  # C22H10O5N2 CATIMAと組成が異なる

    elif material == "Pla" or material == "Plastics" or material == "Plastic":
        solid_gas = "solid"
        density = 1.032
        zts, m_fractions = [1, 6], [10, 9]  # C9H10

    elif material == "Lucite" or material == "Acrylic" or material == "PMMA":
        solid_gas = "solid"
        density = 1.18
        zts, m_fractions = [1, 6, 8], [8, 5, 2]  #

    else:
        ZTarget = None
        if material == "GasHe":
            solid_gas = "gas"
            zts, m_fractions = [2], [1]
            density = 0.001 * 0.179
        elif material == "GasNe":
            solid_gas = "gas"
            zts, m_fractions = [10], [1]
            density = 0.001 * 0.839
        elif material == "GasAr":
            solid_gas = "gas"
            zts, m_fractions = [18], [1]
            density = 0.001 * densities["PureAr"]  # PDG 1.662 (20deg 1atm)
        elif material == "GasKr":
            solid_gas = "gas"
            zts, m_fractions = [36], [1]
            density = 0.001 * 3.4
        elif material == "GasXe" or material == "(CH4)0Xe10":
            solid_gas = "gas"
            zts, m_fractions = [54], [1]
            density = 0.001 * densities["PureXe"]  # PDG 5.483 (20deg 1atm)
        elif material == "Be" or material == "Beryllium":
            solid_gas = "solid"
            zts, m_fractions = [4], [1]
            density = 1.848
        elif material == "Diamond":
            solid_gas = "solid"
            zts, m_fractions = [6], [1]
            density = 3.52
        elif material == "Carbon":
            solid_gas = "solid"
            zts, m_fractions = [6], [1]
            density = 1.8
        elif material == "Al" or material == "Aluminium" or material == "Aluminum":
            solid_gas = "solid"
            zts, m_fractions = [13], [1]
            density = 2.702
        elif material == "Cu" or material == "Copper":
            solid_gas = "solid"
            zts, m_fractions = [29], [1]
            density = 8.96
        elif material == "Ta" or material == "Tantalum":
            solid_gas = "solid"
            zts, m_fractions = [73], [1]
            density = 16.65
        elif material == "W" or material == "Tungsten":
            solid_gas = "solid"
            zts, m_fractions = [74], [1]
            density = 19.3
        elif material == "Pt" or material == "Platinum":
            solid_gas = "solid"
            zts, m_fractions = [78], [1]
            density = 21.45
        elif material == "Au" or material == "Gold":
            solid_gas = "solid"
            zts, m_fractions = [79], [1]
            density = 19.32
        else:
            solid_gas = None
            zts, m_fractions = [int(material)], [1]
            density = 1
    return {"zts": zts, "m_fractions": m_fractions, "density": density * density_factor, "solid_gas": solid_gas}


def GetMFP(zp, energy, material, density_factor=1, solid_gas=None):
    """
    物質の電荷変化平均自由行程、固体かガスか、密度を得る
    solid_gas: "gas" or "solid"
    density: g/cm3 もしmaterialを数値で与えた場合は1g/cm3になる
    """
    result = GetMaterial(material, density_factor)
    zts, m_fractions, density, solid_gas2 = result["zts"], result["m_fractions"], result["density"], result["solid_gas"]
    return GetMixedMFP(zp, energy, zts, m_fractions, solid_gas if solid_gas2 is None else solid_gas2, density)


def GetAnalyticalEqProbFromCS(cs):
    """
    断面積の比から平衡状態における電荷分布を得る
    """
    cs_list = list(cs.values())
    r = []
    for i in range(0, len(cs_list), 2):
        r.append(cs_list[i + 1] / cs_list[i])
    d = 1
    for i in range(len(r)):
        psum = 1
        for j in range(i, len(r)):
            psum *= r[j]
        d += psum
    f = [0 for _ in r]
    f.append(1 / d)
    for i in range(len(r) - 1, -1, -1):
        f[i] = r[i] * f[i + 1]
    return f


def GetAnalyticalEqProb(MFP):
    """
    平均自由行程比から平衡状態における電荷分布を得る
    """
    cs = {}
    for key, value in MFP.items():
        cs[key] = 1 / value
    return GetAnalyticalEqProbFromCS(cs)


def GetAnalyticalEqCharge(MFP):
    """
    平均自由行程比から平衡状態における平均電荷を得る
    """
    dist = GetAnalyticalEqProb(MFP)
    zp = int(list(MFP.keys())[0].split("-")[0])
    sum0 = 0.0
    sum1 = 0.0
    for i, f in enumerate(dist):
        sum0 = sum0 + f
        sum1 = sum1 + f * (zp - i)
    return sum1 / sum0


def GetAnalyticalEqNcc(MFP, length):
    """
    平衡状態における長さあたりの平均電荷変化回数を得る
    """
    dist = GetAnalyticalEqProb(MFP)
    zp = int(list(MFP.keys())[0].split("-")[0])
    ncc_sum = 0
    for i, f in enumerate(dist):
        if f"{zp-i}->{zp-i+1}" in MFP:
            ncc_sum += f * length / MFP[f"{zp-i}->{zp-i+1}"]
        if f"{zp-i}->{zp-i-1}" in MFP:
            ncc_sum += f * length / MFP[f"{zp-i}->{zp-i-1}"]
    return ncc_sum


def GetMCHistories(MFP, initialQ, length, N=10000, random_state=None, histories=None, ignored=False):
    """
    平均自由行程からモンテカルロ法による電荷変化履歴を得る
    lengthは、MFPと同じ単位を与える。MFPの単位がcmならlengthもcmで
    ignored が Trueのときは、電荷変化履歴は計算するが長さは変えない
    initialQ は 初期電荷か、電荷リストを与える
    random_state は以下のようにrsを生成して使い回す
    rs = np.random.RandomState(1)
    """
    assert random_state is not None  # 乱数発生器の指定を確認
    rs = random_state

    zp = int(list(MFP.keys())[0].split("-")[0])  # 初期電荷を取得

    if histories == None:
        histories = [[] for _ in range(N)]
        if type(initialQ) is list:
            initialQs = np.array(initialQ, dtype=np.int32)
        else:
            initialQs = np.full(N, initialQ, dtype=np.int32)
        offset_length = 0
    else:
        assert len(histories) == N
        initialQs = np.array([histories[n][-1][0] for n in range(N)], dtype=np.int32)
        offset_length = histories[0][-1][1]

    length = length + offset_length  # シミュレーションの終点を設定

    # MFP のキャッシュ (辞書参照回数を削減)
    mfp_cache = {int(key.split("->")[0])*10+(int(key.split("->")[1])-int(key.split("->")[0])) : MFP[key] for key in MFP}
    Qmax = max([int(key.split("->")[0]) for key in MFP])
    Qmin = min([int(key.split("->")[0]) for key in MFP])
    mfp_cache[Qmax*10+1] = float("inf")
    mfp_cache[Qmin*10-1] = float("inf")
    lookup = np.vectorize(mfp_cache.get)

    current_length = np.full(N, offset_length, dtype=np.float64)
    Q = initialQs.copy()

    current_length_list = [current_length.copy()]
    Q_list = [Q]
    dQ_list = []
    while True:
        lp = rs.exponential(lookup(Q * 10 + 1))
        lm = rs.exponential(lookup(Q * 10 - 1))
        current_length += np.where(lp < lm, lp, lm)
        current_length_list.append(current_length.copy())
        Q = np.where(lp < lm, Q + 1, Q - 1)
        Q_list.append(Q.copy())
        dQ_list.append(np.where(lp < lm, +1, -1).copy())
        if np.all(current_length > length):
            break

    current_length_array = np.array(current_length_list).T
    Q_array = np.array(Q_list).T
    dQ_array = np.array(dQ_list).T
    for n in range(N):
        history = [[int(Q_array[n][0]),float(current_length_array[n][0]),"pre",zp]]
        for Q, current_length, dQ in zip(Q_array[n][1:], current_length_array[n][1:],dQ_array[n]):
            if current_length > length:
                history.append([int(Q - dQ), float(length), "post"])
                break
            else:
                history.append([int(Q), float(current_length), "+" if dQ > 0 else "-"])

        if ignored:
            # 計算した最後のQを距離0で追加する
            if len(histories[n]) > 0:
                history[-1][1] = histories[n][-1][1]
            else:
                history[-1][1] = 0
            history[-1][2] = "ignored"
            histories[n] += [history[0], history[-1]]
        else:
            histories[n] += history
    return histories


def CheckLength(histories, l1, l2):
    import math

    assert l1 <= l2
    if histories[0][0][1] > l1 and math.isclose(l1, histories[0][0][1]):
        l1 = histories[0][0][1]
    if l2 > histories[0][-1][1] and math.isclose(l2, histories[0][-1][1]):
        l2 = histories[0][-1][1]
    assert histories[0][0][1] <= l1
    assert l2 <= histories[0][-1][1]
    return l1, l2


def GetMCProbImpl(histories, length):
    """
    モンテカルロ法による電荷変化履歴からlengthの位置における電荷を得る
    長さを無視して電荷状態だけ計算する場合があるため、履歴の下流から探索する
    """
    import math

    Qs = []
    _, length = CheckLength(histories, 0, length)
    for history in histories:
        positions = [h[1] for h in history]
        index = bisect.bisect_right(positions, length)

        h = history[index - 1]
        if index > 0:
            if math.isclose(h[1], length) or h[1] < length:
                Qs.append(h[0])
                continue
    assert len(Qs) == len(histories)
    return Qs


def GetMCProb(histories, length):
    Qs = GetMCProbImpl(histories, length)
    return [np.sum(np.array(Qs) == Q) / len(histories) for Q in range(histories[0][0][3], histories[0][0][3] - 7, -1)]


def GetMCDeltaE(histories, l1, l2, dedx_list={}):
    """
    電荷ごとのdE/dxから、電荷変化履歴におけるl1からl2までの実際のdEを得る
    dedx_list[zp],dedx_list[zp-1],...,dedx_list[zp-6] まで計算して与える
    """
    dEs = []
    l1, l2 = CheckLength(histories, l1, l2)
    for history in histories:
        dE = 0
        for i, h in enumerate(history):
            if i == 0:
                continue
            if h[1] < l1:
                continue
            if l2 < h[1]:
                if history[i - 1][1] < l1:
                    dE += dedx_list[history[i - 1][0]] * (l2 - l1)
                else:
                    dE += dedx_list[history[i - 1][0]] * (l2 - history[i - 1][1])
                break
            else:
                if history[i - 1][1] < l1:
                    dE += dedx_list[history[i - 1][0]] * (history[i][1] - l1)
                else:
                    dE += dedx_list[history[i - 1][0]] * (history[i][1] - history[i - 1][1])
        dEs.append(dE)
    return dEs


def GetMCMeanProb(histories, l1, l2):
    """
    電荷変化履歴におけるl1からl2までの電荷状態の存在確率を得る
    """
    Z = histories[0][0][3]
    Qs = [Z, Z - 1, Z - 2, Z - 3, Z - 4, Z - 5, Z - 6]
    Ps = dict(zip(Qs, [[] for _ in Qs]))
    l1, l2 = CheckLength(histories, l1, l2)
    length = l2 - l1
    for history in histories:
        P = dict(zip(Qs, [0 for _ in Qs]))
        for i, h in enumerate(history):
            if i == 0:
                continue
            if h[1] < l1:
                continue
            if l2 < h[1]:
                if history[i - 1][1] < l1:
                    P[history[i - 1][0]] += (l2 - l1) / length
                else:
                    P[history[i - 1][0]] += (l2 - history[i - 1][1]) / length
                break
            else:
                if history[i - 1][1] < l1:
                    P[history[i - 1][0]] += (history[i][1] - l1) / length
                else:
                    P[history[i - 1][0]] += (history[i][1] - history[i - 1][1]) / length
        for k, v in P.items():
            Ps[k].append(v)
    return Ps


def GetMCNcc(histories, l1, l2):
    """
    電荷変化履歴におけるl1からl2までの電荷状態の変化回数を得る
    """
    l1, l2 = CheckLength(histories, l1, l2)
    import math

    ChargeChangings = []
    Qs = GetMCProbImpl(histories, l1)
    for history, Charge in zip(histories, Qs):
        ChargeChanging = 0
        l = l1

        for h in history:
            if h[1] < l1:
                continue  # 範囲外
            if l2 < h[1]:
                break  # 範囲外
            if math.isclose(h[1], l):  # 同じ場所
                Charge = h[0]
            else:
                l = h[1]
                if Charge != h[0]:  # chargeが違う
                    Charge = h[0]
                    ChargeChanging += 1
        ChargeChangings.append(ChargeChanging)
    return ChargeChangings


def GetMCMeanCharge(histories, l1, l2):
    """
    電荷変化履歴におけるl1からl2までの平均電荷を得る
    """
    dedx_list = {}
    Z = histories[0][0][3]
    for dQ in range(7):
        dedx_list[Z - dQ] = (Z - dQ) / (l2 - l1)
    return GetMCDeltaE(histories, l1, l2, dedx_list=dedx_list)


def GetAnalyticalProbImpl(MFP, x, P0):
    """
    平均自由行程データ MFP から距離 x 進んだときの状態確率を計算する関数の実装部分
    """
    states = sorted(set(int(k.split("->")[0]) for k in MFP.keys()))
    n = len(states)
    Zmax = max(states)

    # 平均自由行程行列
    MFP_matrix = np.zeros((n, n))
    for k, v in MFP.items():
        i, j = Zmax - int(k.split("->")[0]), Zmax - int(k.split("->")[1])
        MFP_matrix[j, i] = v

    # 遷移行列 L の構築
    T_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if MFP_matrix[i, j] == 0:
                continue
            T_matrix[i, j] = -1 / MFP_matrix[i, j]
    for i in range(n):
        T_matrix[i, i] = -np.sum(T_matrix[:, i])

    #  行列指数関数の計算
    exp_Tx = scipy.linalg.expm(-T_matrix * x)
    P1 = np.dot(exp_Tx, P0)

    return P1


def GetAnalyticalProb(MFP, x, charge_state=0):
    """
    平均自由行程データ MFP から距離 x 進んだときの状態確率を計算する関数
    """
    states = sorted(set(int(k.split("->")[0]) for k in MFP.keys()))
    n = len(states)

    P0 = np.zeros(n)
    try:
        iter(charge_state)
        for i, v in enumerate(charge_state):
            P0[i] = v
    except:
        P0[charge_state] = 1.0  # 初期状態を設定

    return GetAnalyticalProbImpl(MFP, x, P0).tolist()


def GetAnalyticalEqThick(MFP, charge_state=0, threshold=1 / np.exp(6)):
    """
    平衡状態に到達するまでの厚さ x を求める
    """

    def condition(x, A, B_func, threshold):
        B_x = B_func(x)
        return np.sum(np.abs(A - B_x)) / len(A) - threshold

    # 平衡状態の電荷分布
    EqDist = np.array(GetAnalyticalEqProb(MFP))

    P0 = np.zeros(len(EqDist))
    P0[charge_state] = 1.0  # 初期状態を設定

    # 初期状態で既に平衡状態の場合
    if condition(0, EqDist, lambda x: GetAnalyticalProbImpl(MFP, x, P0), threshold) < 0:
        return 0

    solution = scipy.optimize.fsolve(condition, 0.000001, args=(EqDist, lambda x: GetAnalyticalProbImpl(MFP, x, P0), threshold))
    return solution[0]
