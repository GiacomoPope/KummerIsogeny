import time
import cProfile
import pstats
from tabulate import tabulate

from sage.all import *

from kummer_line import KummerLine
from kummer_isogeny import KummerLineIsogeny_Velu, KummerLineIsogeny_VeluSqrt

proof.all(False)


def print_comparison(ells):
    data = []
    for ell in ells:
        if ell < 17:
            continue
        ell = ZZ(ell)
        cofactor = (p + 1) // ell
        R = cofactor * P
        ker = K(R)
        img = K(Q)

        repeat = 50

        velu_codomain = 10**100
        for _ in range(repeat):
            t0 = time.process_time_ns()
            psi = KummerLineIsogeny_Velu(K, ker, ell)
            velu_codomain = min(time.process_time_ns() - t0, velu_codomain)
        velu_codomain = velu_codomain // 1000

        velu_image = 10**100
        for _ in range(repeat):
            t0 = time.process_time_ns()
            _ = psi(img)
            velu_image = min(time.process_time_ns() - t0, velu_image)
        velu_image = velu_image // 1000

        sqrtvelu_codomain = 10**100
        for _ in range(repeat):
            t0 = time.process_time_ns()
            phi = KummerLineIsogeny_VeluSqrt(K, ker, ell)
            sqrtvelu_codomain = min(time.process_time_ns() - t0, sqrtvelu_codomain)
        sqrtvelu_codomain = sqrtvelu_codomain // 1000

        sqrtvelu_image = 10**100
        for _ in range(repeat):
            t0 = time.process_time_ns()
            _ = phi(img)
            sqrtvelu_image = min(time.process_time_ns() - t0, sqrtvelu_image)
        sqrtvelu_image = sqrtvelu_image // 1000

        j1 = psi.codomain().j_invariant()
        j2 = phi.codomain().j_invariant()

        assert j1 == j2

        # Ratios of image vs codomain
        ratio_velu = velu_image / velu_codomain
        ratio_sqrt = sqrtvelu_image / sqrtvelu_codomain
        ratio_velu = f"{ratio_velu:0.2f}"
        ratio_sqrt = f"{ratio_sqrt:0.2f}"

        ratio_codomain = velu_codomain / sqrtvelu_codomain
        ratio_image = velu_image / sqrtvelu_image
        ratio_codomain = f"{ratio_codomain:0.2f}"
        ratio_image = f"{ratio_image:0.2f}"

        ell_data = [
            ell,
            velu_codomain,
            sqrtvelu_codomain,
            velu_image,
            sqrtvelu_image,
            ratio_velu,
            ratio_sqrt,
            ratio_codomain,
            ratio_image,
        ]
        data.append(ell_data)

    print(
        tabulate(
            data,
            headers=[
                "ell",
                "velu codomain (us)",
                "sqrt codomain (us)",
                "velu image (us)",
                "sqrt image (us)",
                "image / codomain velu",
                "image / codomain sqrt",
                "velu/sqrt codomain",
                "velu/sqrt image",
            ],
        )
    )


def profile_codomain(ell):
    p = K.base_ring().characteristic()

    cofactor = (p + 1) // ell
    R = cofactor * P
    ker = K(R)

    p_codomain = cProfile.Profile()
    p_codomain.enable()

    for _ in range(100):
        _ = KummerLineIsogeny_VeluSqrt(K, ker, ell)

    p_codomain.disable()
    p_codomain.dump_stats("p_codomain.cProfile")
    p = pstats.Stats("p_codomain.cProfile")
    p.strip_dirs().sort_stats("cumtime").print_stats(30)


def profile_image(ell):
    p = K.base_ring().characteristic()
    ell = ZZ(ell)

    cofactor = (p + 1) // ell
    R = cofactor * P
    ker = K(R)
    img = K(Q)

    phi = KummerLineIsogeny_VeluSqrt(K, ker, ell)

    p_image = cProfile.Profile()
    p_image.enable()

    for _ in range(100):
        phi(img)

    p_image.disable()
    p_image.dump_stats("p_image.cProfile")
    p = pstats.Stats("p_image.cProfile")
    p.strip_dirs().sort_stats("cumtime").print_stats(30)


if __name__ == "__main__":
    # FESTA-128 Prime, lots of smooth factors to check
    p = ZZ(
        0x176C11CF13E54B11406FCEC87BD4C1480F2BF6B3CF47C54370FEBD1C756E54F72C1501712922BAF5993402979D50DD13D09A841FED4773CFDB168F19A73E323F656921D7DCD797059B7B9AC3245C4D7BE6B343FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    )

    """
    d1 = 3**6 * 19**2 * 29**2 * 37**2 * 83**2 * 139**2 * 167**2 * 251**2 * 419**2 * 421**2 * 701**2 * 839**2 * 1009**2 * 1259**2 * 3061**2 * 3779**2)
    d2 = 5**4 * 7**3 * 11**2 * 13**2 * 17**2 * 41**2 * 43**2 * 71**2 * 89**2 * 127**2 * 211**2 * 281**2 * 503**2 * 631**2 * 2309**2 * 2521**2 * 2647**2 * 2729**2)
    dA1_sqrt = 59 * 6299 * 6719 * 9181
    dA2_sqrt = 3023 * 3359 * 4409 * 5039 * 19531 * 22679 * 41161
    """
    d1 = ZZ(0x16C1B23A83CA861E79769376C62F334BB8F6346C820A0E9BD6A652D42ECF246B9)
    d2 = ZZ(0xDB6B524664D2555466FA41A3469CB451C0B756A33B3FDD8F09EDC26E4E9035DDF)
    dA1_sqrt = ZZ(0x14D9C07F458B)
    dA2_sqrt = ZZ(0xD4A42112B63144D6E2F7023)

    # Compute curve
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    i = F.gens()[0]

    E = EllipticCurve(F, [0, 6, 0, 1, 0])
    E.set_order((p + 1) ** 2)

    # Define generators of E[p+1] to avoid recomputing them
    Py0 = 0x60AD711E082F16D9F75123E93D365C356C60E0A2F454FEF4F2720D95BB6DD1AF983E68757D6B422D9BE60F933E030546AFCF9D083F12C07FFC8AD7593F819CFB975CAF1282A2B653571E0C3F9359840C4A8D51DB0A434DBDD40D8530CB6F92836432D45146D266255C0C877512793091A39901D6D6D59DC989F750BA92F27AE45DE486BA29D841B9736DE33B0E5A9B0FE5EB907BA28EDD0A91E2DAC6C1AEB08D068
    Py1 = 0x14E5027864CC66B9C0BDC6EC56776A6A87442C11799A8D2B817E8B319943BD8F87C4435DF720AB3A1CF86B473E8953C258DB24AFAEA3137DA7923B85D3692D7537261F191E2BABD8161B1515C83E509454ABCBF93AE99980C8E4DBB814FADF41625109C80D82BDC043E7321A390895202DD46146C2EA31717E379BE69A8D7D30CEAE1E26ADD83B8871946BBD3D5AA46C8EE4900D0CF20F3D8EF644B337E88D2D3C4
    Px0 = 2
    Px1 = 1
    P = E(Px0 + i * Px1, Py0 + i * Py1)

    Qx0 = 20
    Qx1 = 2
    Qy0 = 0x37C3CAB67CE8AB5020CF6B1B8D85A35D7C123AD99E57200C4DC3907A45E62348721688FCDD463F954BDEB67092B2674F4C8E607D111A1BC9DBEF4DAA928E1E92FD721FC9EE8DB01E0FCAA70F1BEA582C1D06C26A55EF673ACF348D2A1849575F808D8C8D8884E2F7DFAEBDC41D27ED0D0CE639FFFBB8F0A1FF76A2474C890AD10B00506F8D4A4E5750978FC82EA8A5FBFC8811156268604614394D3EE1205E4686C
    Qy1 = 0xCF92282ED56AFBB8C29EB1A65D45D1C9D1FF571FDD4F36ED7287F580B11B663C02A310AB99A22A1F839734D6B3F973A23889F14D32D1F355ACA203EEDF5A2FEB7BADED8DE0303922FBC9F616D1C43FCE6C7EDA1B326DEBAD961F977D2AAA016930D06414CAC009EBC73104AC66CDFEA2FF914DCA19DB704F98219172B4F846CA54F5C0B7FBAC1D5CECED54674866F95D02E1D533EE50159CBE3CCD1F390F6913DB5
    Q = E(Qx0 + i * Qx1, Qy0 + i * Qy1)

    K = KummerLine(E)

    ells = [ell for ell, _ in factor(p + 1)]
    print_comparison(ells)

    # profile_image(839)
    # profile_image(41161)
