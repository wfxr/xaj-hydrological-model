package HydrologicalModel;

import static HydrologicalModelHelper.ParameterValidation.*;
import static com.google.common.base.Preconditions.*;

/**
 * 新安江水文模型
 * <p>
 * Created by Wenxuan on 2016/1/11.
 */
@SuppressWarnings("SpellCheckingInspection")
public class XajModel {
    /*土壤蓄水量参数*/
    private double WUM;     // 上层土壤蓄水容量
    private double WLM;     // 下层土壤蓄水容量
    private double WDM;     // 底层土壤蓄水容量
    private double WM;      // 土壤总蓄水容量

    /*蒸散发参数*/
    private double K;       // 蒸发皿折算系数
    private double C;       // 深层蒸散发系数

    /*产流参数*/
    private double B;       // 蓄水容量曲线的指数
    private double Imp;     // 不透水面积与全流域面积的比值

    /*水源划分参数*/
    private double SM;      // 流域平均自由水蓄水容量
    private double EX;      // 自由水蓄水容量曲线指数
    private double KSS;     // 自由水蓄水库对壤中流的出流系数
    private double KG;      // 自由水蓄水库对地下径流出流系数

    /*汇流参数*/
    private double KKSS;    // 地下水库消退系数
    private double KKG;     // 壤中流水库消退系数
    private double Area;    // 水文单元面积


    /**
     * 产流计算
     *
     * @param P   降雨量序列
     * @param EI  蒸发皿蒸发量序列
     * @param wu0 初始时刻上层土壤含水量
     * @param wl0 初始时刻下层土壤含水量
     * @param wd0 初始时刻深层土壤含水量
     * @return 产流计算结果
     */
    public RunoffGenerationResult ComputeRunoffGeneration(double[] P, double[] EI, double wu0, double wl0, double wd0) {
        int n = CheckIdenticalLength(P, EI);

        RunoffGenerationResult result = new RunoffGenerationResult(n);
        LayeredSoilParam w0 = new LayeredSoilParam(wu0, wl0, wd0);  // 时段初土壤分层含水量
        double w0Sum = w0.getSum();                   // 时段初土壤总含水量
        double wmmax = WM * (1 + B) / (1 - Imp);      // 流域点最大蓄水容量
        for (int i = 0; i < n; ++i) {
            double p = P[i];                          // 时段内降雨量
            double ep = EI[i] * K;                    // 计算时段内实际蒸发量
            double e = ComputeE(ep, p, w0.U, w0.L);   // 计算时段内土壤蒸发量

            double pe = p - e;                        // 计算时段内净雨量
            double r = ComputeR(pe, wmmax, w0Sum);    // 计算时段内产流量
            w0 = ComputeLayeredW(w0, pe, r);          // 计算时段末土壤分层含水量
            w0Sum = w0.getSum();                      // 计算时段末土壤总含水量

            result.Set(i, pe, r, w0Sum);
        }
        return result;
    }

    /**
     * 水源划分计算
     *
     * @param runoffYieldResult 产流计算结果
     * @param s0                初始时刻流域平均自由含水量
     * @param dt                时段长度
     * @return 水源划分计算结果
     */
    public SourcePartitionResult ComputeSourcePartition(RunoffGenerationResult runoffYieldResult, double s0, double dt) {
        int n = runoffYieldResult.Length;
        SourcePartitionResult result = new SourcePartitionResult(n);

        double KSSD = (1 - Math.pow(1 - (KG + KSS), dt / 24)) / (1 + KG / KSS);     // 转化后的壤中流出流系数
        double KGD = KSSD * KG / KSS;     // 转化后的地下水出流系数
        double Smax = (1 + EX) * SM;      // 流域点自由蓄水量最大值
        for (int i = 0; i < n; ++i) {
            double pe = runoffYieldResult.PE[i];
            double w = runoffYieldResult.W[i];
            double r = runoffYieldResult.R[i];

            double au = Smax * (1 - Math.pow(1 - s0 / SM, 1 / (1 + EX)));
            double FR;              // 产流面积
            double coef;
            if (pe > 0) {
                FR = r / pe - Imp;
                coef = (pe + au < Smax) ? SM - SM * Math.pow(1 - (pe + au) / Smax, 1 + EX) : SM;
            } else {
                FR = 1 - Imp - Math.pow(1 - w / WM, B / (1 + B)) * (1 - Imp);
                coef = s0;
            }

            double rimp = pe * Imp;
            double rs = Math.max((pe + s0 - coef) * FR, 0.0);
            double rss = coef * KSSD * FR;
            double rg = coef * KGD * FR;
            s0 = coef * (1 - KSSD - KGD);

            result.Set(i, rimp, rs, rss, rg, s0);
        }

        return result;
    }

    /**
     * 汇流计算
     *
     * @param sourcePartitionResult 水源划分计算结果
     * @param qrss0                 初始时刻壤中流流量
     * @param qrg0                  初始时刻地下径流流量
     * @param dt                    时段长度
     * @return 汇流计算结果
     */
    public RunoffConcentrationResult ComputeRunoffConcentration(SourcePartitionResult sourcePartitionResult, double qrss0, double qrg0, double dt) {
        double[] QRS = UnitHydrograph(sourcePartitionResult.RS, sourcePartitionResult.RIMP, dt);    // 计算地表径流
        double[] QRSS = LinearReservoir(sourcePartitionResult.RSS, KKSS, qrss0, dt);    // 计算壤中流
        double[] QRG = LinearReservoir(sourcePartitionResult.RG, KKG, qrg0, dt);        // 计算地下径流
        double[] Q = ComputeQ(QRS, QRSS, QRG);                                          // 计算总径流

        return new RunoffConcentrationResult(QRS, QRSS, QRG, Q);
    }

    /**
     * 设置流域土壤含水量参数
     *
     * @param wum 上层土壤含水量
     * @param wlm 下层土壤含水量
     * @param wdm 底层土壤含水量
     */
    public void SetSoilWaterStorageParam(double wum, double wlm, double wdm) {
        WUM = wum;
        WLM = wlm;
        WDM = wdm;
        WM = WUM + WLM + WDM;
    }

    /**
     * 设置流域蒸散发参数
     *
     * @param k 蒸发皿折算系数
     * @param c 深层蒸散发系数
     */
    public void SetEvapotranspirationParam(double k, double c) {
        K = k;
        C = c;
    }

    /**
     * 设置流域产流计算参数
     *
     * @param b   蓄水容量曲线的指数
     * @param imp 不透水面积与全流域面积的比值
     */
    public void SetRunoffGenerationParam(double b, double imp) {
        B = b;
        Imp = imp;
    }

    /**
     * 设置流域水源划分参数
     *
     * @param sm  流域平均自由水蓄水容量
     * @param ex  自由水蓄水容量曲线指数
     * @param kss 自由水蓄水库对壤中流的出流系数
     * @param kg  自由水蓄水库对地下径流的出流系数
     */
    public void SetSourcePartitionParam(double sm, double ex, double kss, double kg) {
        checkArgument(kss + kg < 1, "Parameter validation failed: kss + kg < 1");
        SM = sm;
        EX = ex;
        KSS = kss;
        KG = kg;
    }

    /**
     * 设置流域汇流计算参数
     *
     * @param kkss 地下水库消退系数
     * @param kkg  壤中流水库消退系数
     * @param area 水文单元面积
     */
    public void SetRunoffConcentrationParam(double kkss, double kkg, double area) {
        KKSS = kkss;
        KKG = kkg;
        Area = area;
    }

    /**
     * 计算时段内的土壤蒸发量
     *
     * @param ep  时段内实际蒸发量
     * @param p   时段内降雨量
     * @param wu0 时段初上层土壤含水量
     * @param wl0 时段初下层土壤含水量
     * @return 时段内的土壤蒸发量
     */
    private double ComputeE(double ep, double p, double wu0, double wl0) {
        double eu, el, ed;
        if (p + wu0 >= ep) {
            eu = ep;
            el = ed = 0;
        } else {
            eu = p + wu0;
            if (wl0 >= C * WLM) {
                el = (ep - eu) * wl0 / WLM;
                ed = 0;
            } else if (wl0 >= C * (ep - eu)) {
                el = C * (ep - eu);
                ed = 0;
            } else {
                el = C * wl0;
                ed = C * (ep - eu) - el;
            }
        }
        return eu + el + ed;
    }

    /**
     * 计算时段内的产流量
     *
     * @param pe    时段内净雨量
     * @param wmmax 流域点最大蓄水容量
     * @param w0    时段初土壤总含水量
     * @return 时段内产流量
     */
    private double ComputeR(double pe, double wmmax, double w0) {
        double a = wmmax * (1 - Math.pow(1 - w0 / WM, 1 / (1 + B)));  // 时段初流域蓄水量w相应的纵坐标

        double r;  // 产流量
        if (pe <= 0)
            r = 0;
        else if (pe + a < wmmax)
            r = pe - WM + w0 + WM * Math.pow(1 - (pe + a) / wmmax, 1 + B);
        else
            r = pe - (WM - w0);

        return r;
    }

    /**
     * 计算时段末的土壤含水量
     *
     * @param w0 时段初土壤含水量
     * @param pe 时段内净雨量
     * @param r  时段内产流量
     * @return 时段末的土壤含水量
     */
    private LayeredSoilParam ComputeLayeredW(LayeredSoilParam w0, double pe, double r) {
        double wu0 = w0.U;
        double wl0 = w0.L;
        double wd0 = w0.D;
        double wu, wl, wd;
        double dw = pe - r;
        if (dw > 0) {
            if (wu0 + dw < WUM) {
                wu = wu0 + dw;
                wl = wl0;
                wd = wd0;
            } else {
                wu = WUM;
                wl = wl0 + dw - (WUM - wu0);
                if (wl < WLM)
                    wd = wd0;
                else {
                    wl = WLM;
                    wd = wd0 + dw - (WUM - wu0) - (WLM - wl0);
                }
            }
        } else {
            if (wu0 + dw > 0) {
                wu = wu0 + dw;
                wl = wl0;
                wd = wd0;
            } else {
                wu = 0;
                wl = wu0 + dw + wl0;
                if (wl > 0)
                    wd = wd0;
                else {
                    wl = 0;
                    wd = wu0 + wl0 + dw + wd0;
                }
            }
        }
        return new LayeredSoilParam(wu, wl, wd);
    }

    /**
     * 汇流阶段用以计算壤中流或地下径流的线性水库法
     *
     * @param R   产流序列
     * @param KK  消退系数
     * @param qr0 初始时刻产流量
     * @param dt  时段长度
     * @return 径流序列
     */
    private double[] LinearReservoir(double[] R, double KK, double qr0, double dt) {
        double[] QR = new double[R.length];
        double KKD = Math.pow(KK, dt / 24);
        for (int i = 0; i < R.length; ++i) {
            qr0 = qr0 * KKD + R[i] * Area / (3.6 * dt) * (1 - KKD);
            QR[i] = qr0;
        }
        return QR;
    }

    /**
     * 汇流阶段用以计算地表径流的单位线法
     *
     * @param RS   地表径流产流序列
     * @param RIMP 不透水面积产流序列
     * @param dt   时段长度
     * @return 径流序列
     */
    private double[] UnitHydrograph(double[] RS, double[] RIMP, double dt) {
        int n = CheckIdenticalLength(RS, RIMP);
        double[] result = new double[n];
        double[] UH = new double[]{0.3, 0.6, 0.1};
        for (int i = 0; i < n; ++i) {
            double rs = RS[i];
            double rimp = RIMP[i];
            for (int j = 0; j < UH.length && i + j < n; ++j)
                result[i + j] += (rs + rimp) * UH[j] * Area / (3.6 * dt);
        }

        return result;
    }

    /**
     * 计算径流总量序列
     *
     * @param QRS  地表径流序列
     * @param QRSS 壤中流序列
     * @param QRG  地下径流序列
     * @return 径流总量序列
     */
    private static double[] ComputeQ(double[] QRS, double[] QRSS, double[] QRG) {
        int n = CheckIdenticalLength(QRS, QRSS, QRG);
        double[] Q = new double[n];
        for (int i = 0; i < n; ++i)
            Q[i] = QRS[i] + QRSS[i] + QRG[i];

        return Q;
    }
}

/**
 * 土壤分层参数
 */
class LayeredSoilParam {
    /**
     * @param u 上层参数值
     * @param l 下层参数值
     * @param d 深层参数值
     */
    public LayeredSoilParam(double u, double l, double d) {
        Set(u, l, d);
    }

    /**
     * @param u 上层参数值
     * @param l 下层参数值
     * @param d 深层参数值
     */
    public void Set(double u, double l, double d) {
        U = u;
        L = l;
        D = d;
    }

    /**
     * @return 分层参数值总和
     */
    public double getSum() {
        return U + L + D;
    }

    /**
     * 上层参数值
     */
    public double U;

    /**
     * 下层参数值
     */
    public double L;

    /**
     * 深层参数值
     */
    public double D;
}