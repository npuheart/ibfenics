class IU:
    g = 0.001
    s = 1
    cm = 0.01 
    # g = 1; # 基本单位
    # s = 1; #
    # cm  = 1;  #
    m = 100 * cm; # 基本单位(米、千克)
    kg = 1000 * g; #
    dyn = g * cm / s / s;  # 力(达因、牛顿)
    N = kg * m / s / s;  #
    Pa = N / m / m; #/ 压强(帕斯卡、毫米汞柱)
    kPa = 1000*Pa
    mmHg = 133.3223684 * Pa; #
    J = N * m; # 能量(焦耳、卡)
    cal = 4.1868 * J; #
    W = J / s; # 功率(瓦)

class UserIU:
    def __init__(self, g=0.001, cm=0.01, s=1):
        m = 100 * cm; # 基本单位(米、千克)
        kg = 1000 * g; #
        dyn = g * cm / s / s;  # 力(达因、牛顿)
        N = kg * m / s / s;  #
        Pa = N / m / m; #/ 压强(帕斯卡、毫米汞柱)
        kPa = 1000*Pa
        MPa = 1000*kPa
        mmHg = 133.3223684 * Pa; #
        J = N * m; # 能量(焦耳、卡)
        cal = 4.1868 * J; #
        W = J / s; # 功率(瓦)
        # export:
        self.g      = g
        self.cm     = cm
        self.s      = s
        self.m      = m
        self.kg     = kg
        self.dyn    = dyn
        self.N      = N
        self.Pa     = Pa
        self.kPa    = kPa
        self.MPa    = MPa
        self.mmHg   = mmHg
        self.J      = J
        self.cal    = cal
        self.W      = W
        